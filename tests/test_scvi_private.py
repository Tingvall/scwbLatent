"""Private scVI boundary rejects unauthenticated metadata shapes before fitting."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "inst" / "python"))
from scvi_private import (  # noqa: E402
    PrivateScVIError,
    _device_observer,
    _metadata,
    _request,
    _select_device,
    _training_evidence,
    compute,
)


def _fixture(tmp_path):
    cells = ["cell_2", "cell_10", "cell_1"]
    snapshot = {
        "schema": "scworkbench.ordered_cell_metadata.snapshot.v1",
        "matrix_input": "matrix",
        "matrix_id": "RNA_counts",
        "input_descriptor_identity": {"scheme": "test", "sha256": "a" * 64},
        "input_content_identity": {"scheme": "test", "sha256": "b" * 64},
        "cell_count": 3,
        "cell_id": cells,
        "cell_order_identity": {"scheme": "test", "sha256": "c" * 64},
        "fields": [
            {
                "parameter": "batch_key",
                "name": "batch",
                "field_id": "core:batch",
                "value_type": "categorical",
                "storage_type": "TEXT",
                "values": ["b", "a", "a"],
            },
            {
                "parameter": "categorical_covariates",
                "name": "condition",
                "field_id": "core:condition",
                "value_type": "categorical",
                "storage_type": "TEXT",
                "values": ["stim", "ctrl", "ctrl"],
            },
            {
                "parameter": "continuous_covariates",
                "name": "quality",
                "field_id": "core:quality",
                "value_type": "continuous",
                "storage_type": "REAL",
                "values": [0.1, 0.3, 0.5],
            },
        ],
    }
    request = {
        "schema": "scwblatent.scvi_private.v1",
        "operation": "fit",
        "features": "features.txt",
        "chunks": [],
        "metadata": "metadata.json",
        "n_latent": 2,
        "n_layers": 1,
        "n_hidden": 16,
        "dropout": 0.1,
        "max_epochs": 1,
        "lr": 0.001,
        "early_stopping": True,
        "seed": 1,
        "batch_key": "batch",
        "use_gpu": False,
        "categorical_covariates": ["condition"],
        "continuous_covariates": ["quality"],
    }
    (tmp_path / "request.json").write_text(json.dumps(request), encoding="utf-8")
    (tmp_path / "metadata.json").write_text(json.dumps(snapshot), encoding="utf-8")
    return cells, request, snapshot


def test_ordered_typed_metadata_is_preserved(tmp_path):
    cells, request, _ = _fixture(tmp_path)
    parsed = _request(tmp_path)
    assert parsed == request
    obs = _metadata(tmp_path, request, cells)
    assert obs.index.tolist() == cells
    assert obs["batch"].tolist() == ["b", "a", "a"]
    assert obs["condition"].tolist() == ["stim", "ctrl", "ctrl"]
    assert obs["quality"].tolist() == [0.1, 0.3, 0.5]


@pytest.mark.parametrize(
    "change",
    [
        lambda snapshot: snapshot.update(cell_id=["cell_10", "cell_2", "cell_1"]),
        lambda snapshot: snapshot.update(cell_id=["cell_2", "foreign", "cell_1"]),
        lambda snapshot: snapshot["fields"][0].update(values=["b", None, "a"]),
        lambda snapshot: snapshot["fields"][1].update(name="batch"),
        lambda snapshot: snapshot["fields"][2].update(value_type="categorical"),
        lambda snapshot: snapshot["fields"][2].update(values=[0.1, "bad", 0.5]),
    ],
)
def test_rejects_reordered_foreign_missing_or_mistyped_covariates(tmp_path, change):
    cells, request, snapshot = _fixture(tmp_path)
    change(snapshot)
    (tmp_path / "metadata.json").write_text(json.dumps(snapshot), encoding="utf-8")
    with pytest.raises(PrivateScVIError):
        _metadata(tmp_path, request, cells)


def test_request_rejects_bundle_authority_and_nonboolean_gpu(tmp_path):
    _, request, _ = _fixture(tmp_path)
    for field, value in (("bundle_path", "/private/bundle"), ("use_gpu", "yes")):
        malformed = dict(request, **{field: value})
        (tmp_path / "request.json").write_text(json.dumps(malformed), encoding="utf-8")
        with pytest.raises(ValueError):
            _request(tmp_path)


def test_request_accepts_gpu_and_legacy_cpu_default(tmp_path):
    _, request, _ = _fixture(tmp_path)
    request["use_gpu"] = True
    (tmp_path / "request.json").write_text(json.dumps(request), encoding="utf-8")
    assert _request(tmp_path)["use_gpu"] is True
    del request["use_gpu"]
    (tmp_path / "request.json").write_text(json.dumps(request), encoding="utf-8")
    assert _request(tmp_path)["use_gpu"] is False


def test_explicit_device_selection_never_falls_back():
    calls = []
    cuda = SimpleNamespace(
        is_available=lambda: True,
        device_count=lambda: 1,
        synchronize=lambda index: calls.append(("sync", index)),
    )
    torch = SimpleNamespace(
        cuda=cuda,
        empty=lambda count, device: (
            calls.append(("empty", device))
            or SimpleNamespace(add_=lambda value: calls.append(("add", value)))
        ),
    )
    assert _select_device(torch, False) == "cpu"
    assert calls == []
    assert _select_device(torch, True) == "gpu"
    assert calls == [("empty", "cuda:0"), ("add", 1), ("sync", 0)]
    cuda.is_available = lambda: False
    with pytest.raises(PrivateScVIError, match="GPU requested.*no CUDA"):
        _select_device(torch, True)
    cuda.is_available = lambda: True
    torch.empty = lambda count, device: (_ for _ in ()).throw(RuntimeError("driver"))
    with pytest.raises(PrivateScVIError, match="GPU requested.*unusable"):
        _select_device(torch, True)


def test_training_device_evidence_requires_observed_batch_and_one_device():
    import torch

    cuda = SimpleNamespace(
        get_device_capability=lambda index: (8, 6),
        get_device_name=lambda index: "Synthetic GPU",
    )
    fake_torch = SimpleNamespace(
        __version__="test-torch",
        version=SimpleNamespace(cuda="12.4"),
        cuda=cuda,
    )
    progress = SimpleNamespace(current=SimpleNamespace(completed=2))
    trainer = SimpleNamespace(
        strategy=SimpleNamespace(root_device=torch.device("cuda:0")),
        fit_loop=SimpleNamespace(epoch_progress=progress),
        num_devices=1,
    )
    model = SimpleNamespace(trainer=trainer)
    observer = SimpleNamespace(device=torch.device("cuda:0"), batches=4)
    evidence = _training_evidence(model, observer, fake_torch, True, 5)
    assert evidence["observed_training_device"] == "cuda:0"
    assert evidence["completed_epochs"] == 2
    assert evidence["gpu_name"] == "Synthetic GPU"
    trainer.num_devices = 2
    with pytest.raises(PrivateScVIError, match="one-device"):
        _training_evidence(model, observer, fake_torch, True, 5)
    trainer.num_devices = 1
    observer.device = torch.device("cpu")
    with pytest.raises(PrivateScVIError, match="one-device"):
        _training_evidence(model, observer, fake_torch, True, 5)


def test_batch_observer_records_actual_tensor_device():
    import torch

    observer = _device_observer(torch)
    module = SimpleNamespace(device=torch.device("cpu"))
    observer.on_train_batch_start(None, module, {"X": torch.ones(2)}, 0)
    assert observer.device == torch.device("cpu")
    assert observer.batches == 1
    module.device = torch.device("cuda:0")
    with pytest.raises(PrivateScVIError, match="batch and model devices"):
        observer.on_train_batch_start(None, module, {"X": torch.ones(2)}, 1)


def test_gpu_unavailable_rejects_before_matrix_loading(tmp_path, monkeypatch):
    import scvi_private

    _, request, _ = _fixture(tmp_path)
    request["use_gpu"] = True
    (tmp_path / "request.json").write_text(json.dumps(request), encoding="utf-8")
    monkeypatch.setattr(
        scvi_private,
        "_matrix",
        lambda *args: pytest.fail("matrix was loaded before GPU preflight"),
    )
    with pytest.raises(PrivateScVIError, match="GPU requested"):
        compute(tmp_path)
    assert not (tmp_path / "model_output").exists()
    assert not (tmp_path / "latent.tsv").exists()
