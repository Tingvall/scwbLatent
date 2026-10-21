# scwbLatent

This independent R operation pack is the sole current owner of
`scwblatent::peakvi_fit`, `scwblatent::peakvi_project`, and
`scwblatent::scvi_train`, and `scwblatent::model_bound_latent_union`.
The unqualified core implementations have been removed. It retains the
accepted ABI-v3 latent
calculation, but does not rewrite or reinterpret historical PeakVI
runs. Activate the installed pack with
`wb_operation_pack_activate()` after setting `SCWB_LATENT_PYTHON` to a Python
environment with `scvi-tools`, `anndata`, `scipy`, and `torch`. GPU scVI
requires a CUDA-enabled PyTorch build and an allocated, visible device.

`peakvi_fit` consumes one authenticated matrix through ordered, exactly-once
typed slices. It emits a fitted file-backed latent space, contained PeakVI
model, and fit summary in one promotion group. `peakvi_project` consumes a
G1-authenticated model and disjoint query-matrix slices, loads the retained
model without refitting, and emits a new immutable query space plus summary.
Loading scvi's `model.pt` is pickle-capable; callers must explicitly provide
`pickle_execution_trust = "operator_verified_source_v1"` after independently
verifying the model's source. Core byte hashes authenticate the retained tree
but do not make pickle deserialization safe, and this declaration is not a
signature or sandbox.
The caller chooses the query `space_id_out`; collisions reject and overwrite
is prohibited. The core issues output identities, validates source/model
linkage and complete receipts, and owns publication and rollback.

`model_bound_latent_union` reads complete, separately authenticated fitted
and projected space axes. Its producer-lineage read requires a retained
projection summary and query space from the same accepted projection run,
with that run consuming the retained model paired to the fitted space.
The pack combines the two coordinate matrices in the core-sealed G7 output
order and writes a linked space and summary; it performs no fitting or
projection. Staged producer lineage is rejected. The qualified union accepts
PeakVI and scVI model identities, while the current real producer fixture
covers PeakVI. Its old-core comparison policy requires exactly equal ordered
cell IDs, coordinate values, and summary fields, and does not compare legacy
metadata hashes.

The pack stages at most `n_features * chunk_size` R matrix elements per
typed read. It writes Matrix Market chunks to a private temporary directory;
the Python process sees only stage-relative files and no bundle path or
publication handle. The fit summary records measured maximum R chunk size,
cell/feature counts, and the returned latent-file checksum. The private
Python computation materializes its AnnData working matrix, so this pilot is
not a large-data qualification.

The private matrix loader imposes no fixed total cell, feature, chunk-count or
nonzero-count ceiling. Chunk size remains bounded by the declared 4096-cell
maximum; the full sparse AnnData working matrix must fit available memory.
Removing dataset caps does not make scVI training out-of-core. Ordered scVI
metadata currently supports at most 16 selected covariate fields, and contained
models retain their declared 32-file/1-GiB publication limits.
Private Python execution disables bytecode-cache writes, so fitting does not
modify the installed pack or invalidate its installation identity.

PeakVI `batch_key` remains deliberately restricted to null. By contrast,
`scvi_train` receives core-resolved, typed, exactly ordered batch,
categorical and continuous metadata in an authenticated snapshot. Its matrix
grant declares the nonnegative-integer policy and the private scVI helper
cannot access bundle storage. Its fit publishes a latent space, contained
model and summary in one atomic closure. `use_gpu` defaults to `false` and
selects explicit CPU execution; `true` selects one visible CUDA device and
rejects if unavailable, unusable or not observed during training. The summary
records the request, observed training device, completed epochs, training-batch
evidence, and Torch/CUDA/GPU runtime identity. A GPU request alone is not proof
of GPU training. PeakVI remains CPU-only. Scientific
comparison to the old core uses a prospective
numerical policy (ordered cells and latent dimensions equal, finite
coordinates, and pairwise-distance Spearman correlation at least 0.9),
not byte- or hash-identity claims. Its frozen old-core fixture is in
`tests/fixtures/scvi_old_core_reference.json`; the pack comparator exercises
registered dispatch. The relocated PeakVI legacy tests under
`tests/migration_evidence/` document the retired route; they are not a second
production owner or part of the live pack test suite.

For a GPU fit, set `SCWB_LATENT_PYTHON` to a scvi-tools environment whose
PyTorch build reports CUDA availability inside the allocated job. The helper
uses `cuda:0` as the first scheduler-visible logical device, not physical GPU
zero. It checks allocation before loading private matrix chunks and never
changes a GPU request to CPU. A CUDA-capable driver with CPU-only PyTorch is
insufficient. The fit summary keeps `use_gpu` as the requested setting and
records `observed_training_device`, `completed_epochs`,
`training_batches_observed`, `torch_version`, `torch_cuda_version`,
`cuda_visible_devices`, `gpu_name`, and `gpu_compute_capability`. Empty GPU
fields on a CPU fit are intentional; `observed_training_device` is captured
from a real training batch and checked against the trainer's root device.

After installing the pinned core and this pack into the job's R library, run
the small registered-operation GPU check from the pack source root:

```sh
SCWB_LATENT_PYTHON=/path/to/cuda-scvi-env/bin/python \
  Rscript tests/composite/run_scvi_gpu_smoke.R
```

The script requires exactly one scheduler-visible GPU and runs unbatched and
`sample_id`-batched synthetic count fits through the public operation. It
promotes model and summary outputs, finalizes, reopens, and validates the
recipient. It does not run the Honeybee experiment.
