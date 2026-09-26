"""Focused, real CPU PeakVI private-stage fit/transfer proof."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import io as scipy_io
from scipy import sparse

HELPER = Path(__file__).parents[1] / "inst" / "python" / "peakvi_private.py"
SCHEMA = "scwblatent.peakvi_private.v1"


def _stage(root: Path, cells: list[str], features: list[str], *, seed: int) -> dict:
    root.mkdir()
    (root / "features.txt").write_text("\n".join(features) + "\n", encoding="utf-8")
    (root / "chunks").mkdir()
    random = np.random.default_rng(seed)
    chunks = []
    for ordinal, selected in enumerate(
        (cells[: len(cells) // 2], cells[len(cells) // 2 :]), 1
    ):
        matrix = (random.random((len(features), len(selected))) < 0.3).astype(float)
        matrix[0, :] = 1
        matrix[:, 0] = 1
        matrix_name = f"chunks/{ordinal:06d}.mtx"
        cells_name = f"chunks/{ordinal:06d}.cells.txt"
        scipy_io.mmwrite(str(root / matrix_name), sparse.csc_matrix(matrix))
        (root / cells_name).write_text("\n".join(selected) + "\n", encoding="utf-8")
        chunks.append({"matrix": matrix_name, "cells": cells_name})
    return {
        "schema": SCHEMA,
        "features": "features.txt",
        "chunks": chunks,
        "n_latent": 2,
    }


def _run(stage: Path, request: dict) -> dict:
    (stage / "request.json").write_text(json.dumps(request), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(HELPER), "--stage-root", str(stage)],
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )
    assert result.returncode == 0, result.stderr[-6000:]
    return json.loads((stage / "response.json").read_text(encoding="utf-8"))


def _latent(stage: Path) -> list[list[str]]:
    with (stage / "latent.tsv").open(newline="", encoding="utf-8") as stream:
        return list(csv.reader(stream, delimiter="\t"))


def test_real_cpu_fit_then_disjoint_query_projection_without_refit(tmp_path):
    features = [f"peak{i:03d}" for i in range(40)]
    reference_cells = [f"reference{i:03d}" for i in range(32)]
    query_cells = [f"query{i:03d}" for i in range(12)]
    reference_stage = tmp_path / "fit"
    fit_request = _stage(reference_stage, reference_cells, features, seed=7)
    fit_request.update(
        operation="fit",
        n_layers=1,
        n_hidden=16,
        dropout=0.1,
        max_epochs=1,
        lr=1e-3,
        early_stopping=False,
        seed=7,
    )
    fit = _run(reference_stage, fit_request)
    assert fit["status"] == "ok"
    assert fit["trained"] is True
    assert fit["latent_rows"] == len(reference_cells)
    assert fit["latent_columns"] == 2
    assert fit["largest_input_chunk_cells"] == 16
    assert fit["model_output_dir"] == "model_output"
    assert fit["model_files"]
    assert (
        fit["latent_sha256"]
        == hashlib.sha256((reference_stage / "latent.tsv").read_bytes()).hexdigest()
    )
    for entry in fit["model_files"]:
        model_file = reference_stage / "model_output" / entry["path"]
        assert entry["size"] == model_file.stat().st_size
        assert entry["sha256"] == hashlib.sha256(model_file.read_bytes()).hexdigest()
    fit_latent = _latent(reference_stage)
    assert fit_latent[0] == ["cell_id", "latent_1", "latent_2"]
    assert [row[0] for row in fit_latent[1:]] == reference_cells
    assert np.isfinite(
        [[float(value) for value in row[1:]] for row in fit_latent[1:]]
    ).all()

    query_stage = tmp_path / "project"
    project_request = _stage(query_stage, query_cells, features, seed=11)
    shutil.copytree(reference_stage / "model_output", query_stage / "model_input")
    project_request.update(operation="project", model_input_dir="model_input")
    project = _run(query_stage, project_request)
    assert project["status"] == "ok"
    assert project["trained"] is False
    assert project["model_load_method"] in {"load", "load_query_data"}
    assert project["latent_rows"] == len(query_cells)
    assert project["latent_columns"] == fit["latent_columns"]
    assert (
        project["latent_sha256"]
        == hashlib.sha256((query_stage / "latent.tsv").read_bytes()).hexdigest()
    )
    assert [row[0] for row in _latent(query_stage)[1:]] == query_cells
    assert not (query_stage / "model_output").exists()


def test_rejects_escape_and_bundle_authority_before_outputs(tmp_path):
    stage = tmp_path / "bad"
    request = _stage(
        stage, [f"cell{i}" for i in range(4)], [f"peak{i}" for i in range(4)], seed=2
    )
    request.update(
        operation="fit",
        n_layers=1,
        n_hidden=8,
        dropout=0.1,
        max_epochs=1,
        lr=1e-3,
        early_stopping=False,
        seed=1,
        bundle_path="/forbidden",
    )
    (stage / "request.json").write_text(json.dumps(request), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(HELPER), "--stage-root", str(stage)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "unsupported fields" in result.stderr
    assert not (stage / "response.json").exists()
    assert not (stage / "latent.tsv").exists()
    assert not (stage / "model_output").exists()
    request.pop("bundle_path")
    request["features"] = "../features.txt"
    (stage / "request.json").write_text(json.dumps(request), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(HELPER), "--stage-root", str(stage)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "escapes the stage" in result.stderr
    assert not (stage / "response.json").exists()
    request["features"] = "chunks/../features.txt"
    (stage / "request.json").write_text(json.dumps(request), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(HELPER), "--stage-root", str(stage)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "escapes the stage" in result.stderr
    assert not (stage / "response.json").exists()
    linked_stage = tmp_path / "linked-stage"
    linked_stage.symlink_to(stage, target_is_directory=True)
    result = subprocess.run(
        [sys.executable, str(HELPER), "--stage-root", str(linked_stage)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "must not be a symlink" in result.stderr
    assert not (stage / "response.json").exists()
    request["features"] = "features.txt"
    (stage / "latent.tsv").write_text("preexisting\n", encoding="utf-8")
    (stage / ".response.tmp").write_text("preexisting\n", encoding="utf-8")
    (stage / "request.json").write_text(json.dumps(request), encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(HELPER), "--stage-root", str(stage)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert (stage / "latent.tsv").read_text(encoding="utf-8") == "preexisting\n"
    assert (stage / ".response.tmp").read_text(encoding="utf-8") == "preexisting\n"


def test_rejects_oversized_matrix_header_before_reading_values(tmp_path):
    stage = tmp_path / "oversized"
    request = _stage(
        stage, [f"cell{i}" for i in range(4)], [f"peak{i}" for i in range(4)], seed=3
    )
    request.update(
        operation="fit",
        n_layers=1,
        n_hidden=8,
        dropout=0.1,
        max_epochs=1,
        lr=1e-3,
        early_stopping=False,
        seed=3,
    )
    (stage / "request.json").write_text(json.dumps(request), encoding="utf-8")
    (stage / request["chunks"][0]["matrix"]).write_text(
        "%%MatrixMarket matrix coordinate real general\n4 2 100000001\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [sys.executable, str(HELPER), "--stage-root", str(stage)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert result.returncode != 0
    assert "header exceeds typed bounds" in result.stderr
    assert not (stage / "latent.tsv").exists()
    assert not (stage / "response.json").exists()
    assert not (stage / "model_output").exists()
