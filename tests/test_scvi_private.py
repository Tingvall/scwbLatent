"""Private scVI boundary rejects unauthenticated metadata shapes before fitting."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "inst" / "python"))
from scvi_private import PrivateScVIError, _metadata, _request  # noqa: E402


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
            {"parameter": "batch_key", "name": "batch", "field_id": "core:batch",
             "value_type": "categorical", "storage_type": "TEXT",
             "values": ["b", "a", "a"]},
            {"parameter": "categorical_covariates", "name": "condition",
             "field_id": "core:condition", "value_type": "categorical",
             "storage_type": "TEXT", "values": ["stim", "ctrl", "ctrl"]},
            {"parameter": "continuous_covariates", "name": "quality",
             "field_id": "core:quality", "value_type": "continuous",
             "storage_type": "REAL", "values": [0.1, 0.3, 0.5]},
        ],
    }
    request = {
        "schema": "scwblatent.scvi_private.v1", "operation": "fit",
        "features": "features.txt", "chunks": [], "metadata": "metadata.json",
        "n_latent": 2, "n_layers": 1, "n_hidden": 16, "dropout": 0.1,
        "max_epochs": 1, "lr": 0.001, "early_stopping": True,
        "seed": 1, "batch_key": "batch",
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


@pytest.mark.parametrize("change", [
    lambda snapshot: snapshot.update(cell_id=["cell_10", "cell_2", "cell_1"]),
    lambda snapshot: snapshot.update(cell_id=["cell_2", "foreign", "cell_1"]),
    lambda snapshot: snapshot["fields"][0].update(values=["b", None, "a"]),
    lambda snapshot: snapshot["fields"][1].update(name="batch"),
    lambda snapshot: snapshot["fields"][2].update(value_type="categorical"),
    lambda snapshot: snapshot["fields"][2].update(values=[0.1, "bad", 0.5]),
])
def test_rejects_reordered_foreign_missing_or_mistyped_covariates(
    tmp_path, change
):
    cells, request, snapshot = _fixture(tmp_path)
    change(snapshot)
    (tmp_path / "metadata.json").write_text(json.dumps(snapshot), encoding="utf-8")
    with pytest.raises(PrivateScVIError):
        _metadata(tmp_path, request, cells)


def test_request_rejects_bundle_authority_and_unsupported_gpu(tmp_path):
    _, request, _ = _fixture(tmp_path)
    for field, value in (("bundle_path", "/private/bundle"),
                         ("use_gpu", True)):
        malformed = dict(request, **{field: value})
        (tmp_path / "request.json").write_text(
            json.dumps(malformed), encoding="utf-8"
        )
        with pytest.raises(ValueError):
            _request(tmp_path)
