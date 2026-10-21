"""Dataset size is not a private transport limit; chunk shape still is."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import io as scipy_io
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parents[1] / "inst" / "python"))
from peakvi_private import PrivatePeakVIError, _matrix  # noqa: E402
from scvi_private import _metadata  # noqa: E402


@pytest.mark.parametrize("n_cells,n_features", [(120895, 3), (3, 100001)])
def test_large_axes_load_with_exact_sparse_values_and_metadata(
    tmp_path, n_cells, n_features
):
    cells = [f"cell_{i}" for i in range(n_cells)]
    features = [f"gene_{i}" for i in range(n_features)]
    (tmp_path / "features.txt").write_text("\n".join(features) + "\n")
    chunks = []
    for ordinal, start in enumerate(range(0, n_cells, 4096)):
        selected = cells[start : start + 4096]
        spec = {"cells": f"{ordinal}.cells", "matrix": f"{ordinal}.mtx"}
        (tmp_path / spec["cells"]).write_text("\n".join(selected) + "\n")
        matrix = sparse.csc_matrix(
            (
                np.ones(len(selected)),
                (np.zeros(len(selected), dtype=int), np.arange(len(selected))),
            ),
            shape=(n_features, len(selected)),
        )
        scipy_io.mmwrite(tmp_path / spec["matrix"], matrix)
        chunks.append(spec)
    data, nnz, largest = _matrix(
        tmp_path, {"features": "features.txt", "chunks": chunks}
    )
    assert data.shape == (n_cells, n_features)
    assert list(data.obs_names) == cells
    assert list(data.var_names) == features
    assert nnz == data.X.nnz == n_cells
    assert np.array_equal(data.X[:, 0].toarray().ravel(), np.ones(n_cells))
    assert largest == min(n_cells, 4096)
    for batched in (False, True):
        values = ["sample_b" if i % 2 else "sample_a" for i in range(n_cells)]
        fields = (
            [
                {
                    "parameter": "batch_key",
                    "name": "sample_id",
                    "field_id": "core:sample_id",
                    "value_type": "categorical",
                    "storage_type": "TEXT",
                    "values": values,
                }
            ]
            if batched
            else []
        )
        snapshot = {
            "schema": "scworkbench.ordered_cell_metadata.snapshot.v1",
            "matrix_input": "matrix",
            "matrix_id": "counts",
            "input_descriptor_identity": {},
            "input_content_identity": {},
            "cell_count": n_cells,
            "cell_id": cells,
            "cell_order_identity": {},
            "fields": fields,
        }
        (tmp_path / "metadata.json").write_text(json.dumps(snapshot))
        request = {
            "metadata": "metadata.json",
            "batch_key": "sample_id" if batched else None,
            "categorical_covariates": [],
            "continuous_covariates": [],
        }
        obs = _metadata(tmp_path, request, cells)
        assert list(obs.index) == cells
        if batched:
            assert list(obs["sample_id"]) == values
    (tmp_path / chunks[0]["cells"]).write_text(cells[0] + "\n" + cells[0] + "\n")
    with pytest.raises(PrivatePeakVIError):
        _matrix(tmp_path, {"features": "features.txt", "chunks": chunks})
