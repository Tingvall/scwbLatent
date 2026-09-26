"""Private, unregistered CPU PeakVI fit/project computation for scwbLatent.

All inputs are typed files below one caller-owned private stage. This process
has no bundle, registry, provenance, or publication authority. The R pack and
core authenticate inputs and the output closure before any durable write.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import Any

import numpy as np
from scipy import io as scipy_io
from scipy import sparse

SCHEMA = "scwblatent.peakvi_private.v1"
MAX_FEATURES = 100_000
MAX_CELLS = 100_000
MAX_CHUNKS = 10_000
MAX_CHUNK_CELLS = 4_096
MAX_NNZ = 100_000_000


class PrivatePeakVIError(ValueError):
    """Invalid private-stage input or scientific output."""


def _closed(value: Any, allowed: set[str], label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise PrivatePeakVIError(f"{label} must be an object")
    extra = set(value) - allowed
    if extra:
        raise PrivatePeakVIError(
            f"{label} contains unsupported fields: {sorted(extra)}"
        )
    return value


def _stage_path(stage: Path, name: Any, *, kind: str) -> Path:
    if not isinstance(name, str) or not name or "\\" in name:
        raise PrivatePeakVIError(
            "private path must be nonempty stage-relative POSIX text"
        )
    if any(part in ("", ".", "..") for part in name.split("/")):
        raise PrivatePeakVIError("private path escapes the stage")
    relative = Path(name)
    if relative.is_absolute():
        raise PrivatePeakVIError("private path escapes the stage")
    path = stage / relative
    for parent in (path, *path.parents):
        if parent == stage.parent:
            break
        if parent.is_symlink():
            raise PrivatePeakVIError("private stage path contains a symlink")
    try:
        resolved = path.resolve(strict=True)
        resolved.relative_to(stage)
    except (OSError, ValueError) as exc:
        raise PrivatePeakVIError(
            "private stage path is absent or escapes the stage"
        ) from exc
    if (kind == "file" and not resolved.is_file()) or (
        kind == "directory" and not resolved.is_dir()
    ):
        raise PrivatePeakVIError(f"private stage path is not a {kind}")
    return resolved


def _ids(path: Path, label: str, ceiling: int) -> list[str]:
    values = path.read_text(encoding="utf-8").splitlines()
    if (
        not values
        or len(values) > ceiling
        or any(not value or "\t" in value or "\r" in value for value in values)
        or len(set(values)) != len(values)
    ):
        raise PrivatePeakVIError(f"{label} must contain unique nonempty identifiers")
    return values


def _int(value: Any, label: str, minimum: int, maximum: int) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not minimum <= value <= maximum
    ):
        raise PrivatePeakVIError(
            f"{label} must be an integer in [{minimum}, {maximum}]"
        )
    return value


def _float(value: Any, label: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PrivatePeakVIError(f"{label} must be numeric")
    result = float(value)
    if not np.isfinite(result) or not minimum <= result <= maximum:
        raise PrivatePeakVIError(
            f"{label} must be finite and in [{minimum}, {maximum}]"
        )
    return result


def _request(stage: Path) -> dict[str, Any]:
    path = _stage_path(stage, "request.json", kind="file")
    request = _closed(
        json.loads(path.read_text(encoding="utf-8")),
        {
            "schema",
            "operation",
            "features",
            "chunks",
            "n_latent",
            "n_layers",
            "n_hidden",
            "dropout",
            "max_epochs",
            "lr",
            "early_stopping",
            "seed",
            "model_input_dir",
        },
        "PeakVI request",
    )
    if request.get("schema") != SCHEMA or request.get("operation") not in {
        "fit",
        "project",
    }:
        raise PrivatePeakVIError("unsupported private PeakVI schema or operation")
    operation = request["operation"]
    if operation == "fit" and "model_input_dir" in request:
        raise PrivatePeakVIError("fit must not receive a model input")
    if operation == "project" and "model_input_dir" not in request:
        raise PrivatePeakVIError("project requires a private model input")
    fit_only = {
        "n_layers",
        "n_hidden",
        "dropout",
        "max_epochs",
        "lr",
        "early_stopping",
        "seed",
    }
    if operation == "project" and fit_only & set(request):
        raise PrivatePeakVIError("project contains unsupported fit parameters")
    required = {"features", "chunks", "n_latent"}
    if operation == "fit":
        required |= fit_only
    if not required <= set(request):
        raise PrivatePeakVIError(
            f"missing PeakVI request fields: {sorted(required - set(request))}"
        )
    request["n_latent"] = _int(request["n_latent"], "n_latent", 2, 4096)
    if operation == "fit":
        request["n_layers"] = _int(request["n_layers"], "n_layers", 1, 16)
        request["n_hidden"] = _int(request["n_hidden"], "n_hidden", 1, 8192)
        request["dropout"] = _float(request["dropout"], "dropout", 0, 1)
        request["max_epochs"] = _int(request["max_epochs"], "max_epochs", 1, 10000)
        request["lr"] = _float(request["lr"], "lr", 1e-12, 1)
        request["seed"] = _int(request["seed"], "seed", 0, 2**31 - 1)
        if not isinstance(request["early_stopping"], bool):
            raise PrivatePeakVIError("early_stopping must be boolean")
    return request


def _matrix(stage: Path, request: dict[str, Any]):
    import anndata as ad  # Import only after closed request/path validation.
    import pandas as pd

    features = _ids(
        _stage_path(stage, request["features"], kind="file"), "features", MAX_FEATURES
    )
    chunks = request["chunks"]
    if not isinstance(chunks, list) or not 1 <= len(chunks) <= MAX_CHUNKS:
        raise PrivatePeakVIError("chunks must be a bounded nonempty array")
    parts = []
    cells: list[str] = []
    total_nnz = 0
    largest_chunk_cells = 0
    for ordinal, item in enumerate(chunks, 1):
        spec = _closed(item, {"matrix", "cells"}, f"chunk {ordinal}")
        if not {"matrix", "cells"} <= set(spec):
            raise PrivatePeakVIError("each chunk requires matrix and cells paths")
        chunk_cells = _ids(
            _stage_path(stage, spec["cells"], kind="file"),
            f"chunk {ordinal} cells",
            MAX_CHUNK_CELLS,
        )
        matrix_path = _stage_path(stage, spec["matrix"], kind="file")
        try:
            rows, columns, declared_nnz, matrix_format, _, symmetry = scipy_io.mminfo(
                str(matrix_path)
            )
        except (OSError, ValueError) as exc:
            raise PrivatePeakVIError(
                f"chunk {ordinal} has an invalid Matrix Market header"
            ) from exc
        if (
            matrix_format != "coordinate"
            or symmetry != "general"
            or (rows, columns) != (len(features), len(chunk_cells))
            or declared_nnz > min(MAX_NNZ - total_nnz, len(features) * len(chunk_cells))
        ):
            raise PrivatePeakVIError(
                f"chunk {ordinal} Matrix Market header exceeds typed bounds"
            )
        matrix = scipy_io.mmread(str(matrix_path))
        matrix = (
            matrix.tocsc() if sparse.issparse(matrix) else sparse.csc_matrix(matrix)
        )
        if matrix.shape != (len(features), len(chunk_cells)):
            raise PrivatePeakVIError(f"chunk {ordinal} dimensions disagree with axes")
        if matrix.nnz and (
            not np.isfinite(matrix.data).all() or (matrix.data < 0).any()
        ):
            raise PrivatePeakVIError(
                "PeakVI matrix entries must be finite nonnegative values"
            )
        total_nnz += matrix.nnz
        largest_chunk_cells = max(largest_chunk_cells, len(chunk_cells))
        if total_nnz > MAX_NNZ or len(cells) + len(chunk_cells) > MAX_CELLS:
            raise PrivatePeakVIError(
                "PeakVI private matrix exceeds its resource ceiling"
            )
        parts.append(matrix)
        cells.extend(chunk_cells)
    if len(set(cells)) != len(cells):
        raise PrivatePeakVIError("chunk cell axes overlap")
    feature_by_cell = sparse.hstack(parts, format="csc")
    data = feature_by_cell.transpose().tocsr()
    obs = pd.DataFrame(index=pd.Index(cells, name="cell_id"))
    var = pd.DataFrame(index=pd.Index(features, name="feature_id"))
    return ad.AnnData(X=data, obs=obs, var=var), total_nnz, largest_chunk_cells


def _model_tree(stage: Path, relative: str) -> tuple[Path, list[dict[str, Any]]]:
    root = _stage_path(stage, relative, kind="directory")
    entries = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise PrivatePeakVIError("model tree contains a symlink")
        if path.is_dir():
            continue
        if not path.is_file() or path.stat().st_size <= 0:
            raise PrivatePeakVIError("model tree contains a nonregular or empty file")
        entries.append(
            {
                "path": path.relative_to(root).as_posix(),
                "size": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    if not entries:
        raise PrivatePeakVIError("model tree is empty")
    return root, entries


def _latent(stage: Path, cells: list[str], matrix: np.ndarray) -> tuple[str, str]:
    if (
        matrix.ndim != 2
        or matrix.shape[0] != len(cells)
        or not np.isfinite(matrix).all()
    ):
        raise PrivatePeakVIError("PeakVI returned an invalid latent matrix")
    path = stage / "latent.tsv"
    if path.exists() or path.is_symlink():
        raise PrivatePeakVIError("latent output already exists")
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["cell_id", *[f"latent_{i}" for i in range(1, matrix.shape[1] + 1)]]
        )
        for cell, row in zip(cells, matrix, strict=True):
            writer.writerow([cell, *[repr(float(x)) for x in row]])
    return path.name, hashlib.sha256(path.read_bytes()).hexdigest()


def compute(stage: Path) -> dict[str, Any]:
    if stage.is_symlink():
        raise PrivatePeakVIError("private stage must not be a symlink")
    stage = stage.resolve(strict=True)
    if not stage.is_dir():
        raise PrivatePeakVIError("private stage must be a regular directory")
    request = _request(stage)
    adata, nnz, largest_chunk_cells = _matrix(stage, request)
    import scvi  # Import only after typed file validation and dependency preflight.

    if request["operation"] == "fit":
        import torch

        seed = request["seed"]
        scvi.settings.seed = seed
        torch.manual_seed(seed)
        np.random.seed(seed)
        scvi.model.PEAKVI.setup_anndata(adata)
        model = scvi.model.PEAKVI(
            adata,
            n_latent=request["n_latent"],
            n_layers_encoder=request["n_layers"],
            n_layers_decoder=request["n_layers"],
            n_hidden=request["n_hidden"],
            dropout_rate=request["dropout"],
        )
        model.train(
            max_epochs=request["max_epochs"],
            lr=request["lr"],
            early_stopping=request["early_stopping"],
            accelerator="cpu",
            devices=1,
        )
        latent = np.asarray(model.get_latent_representation(), dtype=float)
        model_dir = stage / "model_output"
        if model_dir.exists() or model_dir.is_symlink():
            raise PrivatePeakVIError("model output already exists")
        model.save(str(model_dir), overwrite=False, save_anndata=False)
        _, files = _model_tree(stage, "model_output")
        model_result = {"model_output_dir": "model_output", "model_files": files}
        trained = True
    else:
        model_dir, _ = _model_tree(stage, request["model_input_dir"])
        # `load` and `load_query_data` both reuse fitted weights. This helper
        # deliberately never calls train() on a projected query.
        try:
            model = scvi.model.PEAKVI.load(
                str(model_dir), adata=adata, accelerator="cpu", device="auto"
            )
            load_method = "load"
        except Exception as first_error:
            try:
                model = scvi.model.PEAKVI.load_query_data(
                    adata=adata,
                    reference_model=str(model_dir),
                    inplace_subset_query_vars=True,
                    accelerator="cpu",
                    device="auto",
                )
                load_method = "load_query_data"
            except Exception as second_error:
                raise PrivatePeakVIError(
                    f"PeakVI transfer load failed: {first_error}; {second_error}"
                ) from second_error
        latent = np.asarray(model.get_latent_representation(), dtype=float)
        model_result = {"model_load_method": load_method}
        trained = False
    if latent.shape != (adata.n_obs, request["n_latent"]):
        raise PrivatePeakVIError("latent shape differs from declared fitted dimension")
    latent_path, latent_sha256 = _latent(stage, adata.obs_names.tolist(), latent)
    return {
        "schema": SCHEMA,
        "status": "ok",
        "operation": request["operation"],
        "latent_path": latent_path,
        "latent_sha256": latent_sha256,
        "latent_rows": int(latent.shape[0]),
        "latent_columns": int(latent.shape[1]),
        "scvi_version": scvi.__version__,
        "trained": trained,
        "matrix_nnz": int(nnz),
        "largest_input_chunk_cells": int(largest_chunk_cells),
        "effective_params": {
            key: request[key]
            for key in (
                (
                    "n_latent",
                    "n_layers",
                    "n_hidden",
                    "dropout",
                    "max_epochs",
                    "lr",
                    "early_stopping",
                    "seed",
                )
                if trained
                else ("n_latent",)
            )
        },
        **model_result,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="private CPU PeakVI fit/project computation"
    )
    parser.add_argument("--stage-root", required=True, type=Path)
    args = parser.parse_args()
    if args.stage_root.is_symlink():
        raise PrivatePeakVIError("private stage must not be a symlink")
    stage = args.stage_root.resolve(strict=True)
    response_path = stage / "response.json"
    model_output_preexisted = (stage / "model_output").exists()
    latent_preexisted = (stage / "latent.tsv").exists() or (
        stage / "latent.tsv"
    ).is_symlink()
    temporary_preexisted = (stage / ".response.tmp").exists() or (
        stage / ".response.tmp"
    ).is_symlink()
    if response_path.exists() or response_path.is_symlink():
        raise PrivatePeakVIError("response output already exists")
    if latent_preexisted or temporary_preexisted:
        raise PrivatePeakVIError("private output already exists")
    try:
        response = compute(stage)
        temporary = stage / ".response.tmp"
        with temporary.open("x", encoding="utf-8") as stream:
            json.dump(response, stream, allow_nan=False, separators=(",", ":"))
        os.replace(temporary, response_path)
        return 0
    except BaseException:
        if not temporary_preexisted:
            (stage / ".response.tmp").unlink(missing_ok=True)
        if not latent_preexisted:
            (stage / "latent.tsv").unlink(missing_ok=True)
        model_output = stage / "model_output"
        if (
            not model_output_preexisted
            and model_output.is_dir()
            and not model_output.is_symlink()
        ):
            shutil.rmtree(model_output)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
