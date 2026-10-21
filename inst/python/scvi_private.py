"""Private scVI fit from core-authenticated matrix and metadata snapshots.

This helper has no bundle, registry, provenance, or durable publication handle.
The R pack verifies its response and core owns the atomic output closure.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np

from peakvi_private import (
    _closed,
    _float,
    _int,
    _latent,
    _matrix,
    _model_tree,
    _stage_path,
)

SCHEMA = "scwblatent.scvi_private.v1"


class PrivateScVIError(ValueError):
    """Invalid private scVI request, metadata, or trained response."""


def _request(stage: Path) -> dict[str, Any]:
    request = _closed(
        json.loads(_stage_path(stage, "request.json", kind="file").read_text()),
        {
            "schema",
            "operation",
            "features",
            "chunks",
            "metadata",
            "n_latent",
            "n_layers",
            "n_hidden",
            "dropout",
            "max_epochs",
            "lr",
            "early_stopping",
            "seed",
            "batch_key",
            "categorical_covariates",
            "continuous_covariates",
            "use_gpu",
        },
        "scVI request",
    )
    required = {
        "schema",
        "operation",
        "features",
        "chunks",
        "metadata",
        "n_latent",
        "n_layers",
        "n_hidden",
        "dropout",
        "max_epochs",
        "lr",
        "early_stopping",
        "seed",
        "batch_key",
        "categorical_covariates",
        "continuous_covariates",
        "use_gpu",
    }
    if set(request) == required - {"use_gpu"}:
        request["use_gpu"] = False
    if (
        set(request) != required
        or request["schema"] != SCHEMA
        or request["operation"] != "fit"
    ):
        raise PrivateScVIError("unsupported or incomplete private scVI request")
    request["n_latent"] = _int(request["n_latent"], "n_latent", 2, 1024)
    request["n_layers"] = _int(request["n_layers"], "n_layers", 1, 16)
    request["n_hidden"] = _int(request["n_hidden"], "n_hidden", 1, 8192)
    request["max_epochs"] = _int(request["max_epochs"], "max_epochs", 1, 10000)
    request["seed"] = _int(request["seed"], "seed", 0, 2**31 - 1)
    request["dropout"] = _float(request["dropout"], "dropout", 0, 1)
    request["lr"] = _float(request["lr"], "lr", 1e-12, 1)
    if not isinstance(request["early_stopping"], bool):
        raise PrivateScVIError("early_stopping must be boolean")
    if not isinstance(request["use_gpu"], bool):
        raise PrivateScVIError("use_gpu must be boolean")
    if request["batch_key"] is not None and (
        not isinstance(request["batch_key"], str) or not request["batch_key"]
    ):
        raise PrivateScVIError("batch_key must be null or nonempty text")
    for parameter in ("categorical_covariates", "continuous_covariates"):
        values = request[parameter]
        if (
            not isinstance(values, list)
            or any(not isinstance(value, str) or not value for value in values)
            or len(set(values)) != len(values)
        ):
            raise PrivateScVIError(f"{parameter} must be unique nonempty names")
    return request


def _metadata(stage: Path, request: dict[str, Any], cells: list[str]):
    import pandas as pd

    snapshot = _closed(
        json.loads(_stage_path(stage, request["metadata"], kind="file").read_text()),
        {
            "schema",
            "matrix_input",
            "matrix_id",
            "input_descriptor_identity",
            "input_content_identity",
            "cell_count",
            "cell_id",
            "cell_order_identity",
            "fields",
        },
        "ordered cell metadata",
    )
    if (
        set(snapshot)
        != {
            "schema",
            "matrix_input",
            "matrix_id",
            "input_descriptor_identity",
            "input_content_identity",
            "cell_count",
            "cell_id",
            "cell_order_identity",
            "fields",
        }
        or snapshot["schema"] != "scworkbench.ordered_cell_metadata.snapshot.v1"
        or snapshot["matrix_input"] != "matrix"
        or snapshot["cell_id"] != cells
        or snapshot["cell_count"] != len(cells)
    ):
        raise PrivateScVIError("metadata snapshot is not aligned to matrix cells")
    fields = snapshot["fields"]
    if not isinstance(fields, list) or len(fields) > 16:
        raise PrivateScVIError("metadata fields exceed the closed resource bound")
    obs = pd.DataFrame(index=pd.Index(cells, name="cell_id"))
    seen: set[str] = set()
    parameter_names: dict[str, list[str]] = {
        "batch_key": [],
        "categorical_covariates": [],
        "continuous_covariates": [],
    }
    for field in fields:
        item = _closed(
            field,
            {"parameter", "name", "field_id", "value_type", "storage_type", "values"},
            "metadata field",
        )
        if set(item) != {
            "parameter",
            "name",
            "field_id",
            "value_type",
            "storage_type",
            "values",
        }:
            raise PrivateScVIError("metadata field is incomplete")
        parameter = item["parameter"]
        name = item["name"]
        values = item["values"]
        if (
            parameter not in parameter_names
            or not isinstance(name, str)
            or not name
            or name in seen
            or not isinstance(item["field_id"], str)
            or not item["field_id"]
            or not isinstance(values, list)
            or len(values) != len(cells)
        ):
            raise PrivateScVIError("metadata field identity or values are invalid")
        seen.add(name)
        parameter_names[parameter].append(name)
        continuous = parameter == "continuous_covariates"
        if item["value_type"] != ("continuous" if continuous else "categorical"):
            raise PrivateScVIError("metadata semantic value type changed")
        if continuous:
            if item["storage_type"] not in {"INTEGER", "REAL"} or any(
                value is not None
                and (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not np.isfinite(value)
                )
                for value in values
            ):
                raise PrivateScVIError("continuous metadata is not finite numeric")
            obs[name] = pd.Series(values, index=obs.index, dtype="float64")
        else:
            if item["storage_type"] not in {"INTEGER", "REAL", "TEXT"} or any(
                value is not None and not isinstance(value, str) for value in values
            ):
                raise PrivateScVIError("categorical metadata is not nullable text")
            obs[name] = pd.Series(values, index=obs.index, dtype="string")
    expected_batch = [] if request["batch_key"] is None else [request["batch_key"]]
    if (
        parameter_names["batch_key"] != expected_batch
        or parameter_names["categorical_covariates"]
        != request["categorical_covariates"]
        or parameter_names["continuous_covariates"] != request["continuous_covariates"]
    ):
        raise PrivateScVIError("metadata fields differ from declared covariates")
    if obs.isna().any(axis=None):
        raise PrivateScVIError("scVI training requires observed covariate values")
    return obs


def _select_device(torch, use_gpu: bool) -> str:
    if not use_gpu:
        return "cpu"
    if not torch.cuda.is_available() or torch.cuda.device_count() < 1:
        raise PrivateScVIError(
            "GPU requested but no CUDA device is available to PyTorch"
        )
    try:
        probe = torch.empty(1, device="cuda:0")
        probe.add_(1)
        torch.cuda.synchronize(0)
        del probe
    except Exception as error:
        raise PrivateScVIError(
            f"GPU requested but CUDA device is unusable: {error}"
        ) from error
    return "gpu"


def _device_observer(torch):
    from lightning.pytorch import Callback

    class DeviceObserver(Callback):
        def __init__(self):
            self.device = None
            self.batches = 0

        def on_train_batch_start(self, trainer, pl_module, batch, batch_idx):
            if not isinstance(batch, Mapping):
                raise PrivateScVIError(
                    "scVI training batch has no observable tensor device"
                )
            devices = {
                value.device
                for value in batch.values()
                if isinstance(value, torch.Tensor)
            }
            if len(devices) != 1 or next(iter(devices)) != pl_module.device:
                raise PrivateScVIError("scVI training batch and model devices differ")
            device = next(iter(devices))
            if self.device is not None and device != self.device:
                raise PrivateScVIError("scVI training device changed between batches")
            self.device = device
            self.batches += 1

    return DeviceObserver()


def _training_evidence(
    model, observer, torch, use_gpu: bool, max_epochs: int
) -> dict[str, Any]:
    trainer = model.trainer
    root_device = trainer.strategy.root_device
    expected = "cuda" if use_gpu else "cpu"
    completed = trainer.fit_loop.epoch_progress.current.completed
    if (
        observer.batches < 1
        or observer.device is None
        or observer.device != root_device
        or root_device.type != expected
        or trainer.num_devices != 1
        or not isinstance(completed, int)
        or not 1 <= completed <= max_epochs
    ):
        raise PrivateScVIError(
            "scVI trainer did not prove the requested one-device fit"
        )
    index = root_device.index if use_gpu else None
    if use_gpu and index is None:
        raise PrivateScVIError("scVI trainer did not report its CUDA device index")
    capability = torch.cuda.get_device_capability(index) if use_gpu else None
    return {
        "observed_training_device": str(observer.device),
        "trainer_root_device": str(root_device),
        "training_device_type": root_device.type,
        "training_device_index": index,
        "training_batches_observed": observer.batches,
        "completed_epochs": completed,
        "torch_version": torch.__version__,
        "torch_cuda_version": torch.version.cuda,
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES", ""),
        "gpu_name": torch.cuda.get_device_name(index) if use_gpu else "",
        "gpu_compute_capability": f"{capability[0]}.{capability[1]}" if use_gpu else "",
    }


def compute(stage: Path) -> dict[str, Any]:
    if stage.is_symlink():
        raise PrivateScVIError("private stage must not be a symlink")
    stage = stage.resolve(strict=True)
    request = _request(stage)
    import torch

    accelerator = _select_device(torch, request["use_gpu"])
    adata, nnz, largest_chunk_cells = _matrix(stage, request)
    if adata.X.nnz and not np.allclose(adata.X.data, np.rint(adata.X.data), atol=1e-8):
        raise PrivateScVIError("matrix violates the declared integer-count policy")
    adata.obs = _metadata(stage, request, adata.obs_names.tolist())
    import scvi

    scvi.settings.seed = request["seed"]
    torch.manual_seed(request["seed"])
    np.random.seed(request["seed"])
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=request["batch_key"],
        categorical_covariate_keys=request["categorical_covariates"] or None,
        continuous_covariate_keys=request["continuous_covariates"] or None,
    )
    model = scvi.model.SCVI(
        adata,
        n_hidden=request["n_hidden"],
        n_latent=request["n_latent"],
        n_layers=request["n_layers"],
        dropout_rate=request["dropout"],
    )
    observer = _device_observer(torch)
    model.train(
        max_epochs=request["max_epochs"],
        plan_kwargs={"lr": request["lr"]},
        early_stopping=request["early_stopping"],
        accelerator=accelerator,
        devices=1,
        callbacks=[observer],
    )
    training = _training_evidence(
        model, observer, torch, request["use_gpu"], request["max_epochs"]
    )
    latent = np.asarray(model.get_latent_representation(), dtype=float)
    if latent.shape != (adata.n_obs, request["n_latent"]):
        raise PrivateScVIError("scVI latent dimensions differ from request")
    model_dir = stage / "model_output"
    if model_dir.exists() or model_dir.is_symlink():
        raise PrivateScVIError("model output already exists")
    model.save(str(model_dir), overwrite=False, save_anndata=False)
    _, files = _model_tree(stage, "model_output")
    if not any(file["path"] == "model.pt" for file in files):
        raise PrivateScVIError("contained scVI model lacks model.pt")
    latent_path, latent_sha256 = _latent(stage, adata.obs_names.tolist(), latent)
    return {
        "schema": SCHEMA,
        "status": "ok",
        "operation": "fit",
        "trained": True,
        "latent_path": latent_path,
        "latent_sha256": latent_sha256,
        "latent_rows": int(latent.shape[0]),
        "latent_columns": int(latent.shape[1]),
        "model_output_dir": "model_output",
        "model_files": files,
        "matrix_nnz": int(nnz),
        "largest_input_chunk_cells": int(largest_chunk_cells),
        "scvi_version": scvi.__version__,
        **training,
        "effective_params": {
            key: request[key]
            for key in (
                "n_latent",
                "n_layers",
                "n_hidden",
                "dropout",
                "max_epochs",
                "lr",
                "early_stopping",
                "seed",
                "batch_key",
                "categorical_covariates",
                "continuous_covariates",
                "use_gpu",
            )
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="private scVI fit")
    parser.add_argument("--stage-root", required=True, type=Path)
    args = parser.parse_args()
    stage = args.stage_root
    if stage.is_symlink():
        raise PrivateScVIError("private stage must not be a symlink")
    stage = stage.resolve(strict=True)
    response_path = stage / "response.json"
    if response_path.exists() or response_path.is_symlink():
        raise PrivateScVIError("response already exists")
    preexisting = (stage / "model_output").exists()
    try:
        response = compute(stage)
        temporary = stage / ".response.tmp"
        with temporary.open("x", encoding="utf-8") as stream:
            json.dump(response, stream, allow_nan=False, separators=(",", ":"))
        os.replace(temporary, response_path)
        return 0
    except BaseException:
        (stage / ".response.tmp").unlink(missing_ok=True)
        (stage / "latent.tsv").unlink(missing_ok=True)
        output = stage / "model_output"
        if not preexisting and output.is_dir() and not output.is_symlink():
            shutil.rmtree(output)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
