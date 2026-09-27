# scwbLatent (Phase 6 CPU PeakVI and scVI owner)

This independent R operation pack is the sole current owner of
`scwblatent::peakvi_fit`, `scwblatent::peakvi_project`, and
`scwblatent::scvi_train`. The unqualified core implementations have been
removed. It retains the accepted ABI-v3 CPU
pilot calculation, but does not rewrite or reinterpret historical PeakVI
runs. Activate the installed pack with
`wb_operation_pack_activate()` after setting `SCWB_LATENT_PYTHON` to a Python
environment with `scvi-tools`, `anndata`, `scipy`, and CPU `torch`.

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

The pack stages at most `n_features * chunk_size` R matrix elements per
typed read. It writes Matrix Market chunks to a private temporary directory;
the Python process sees only stage-relative files and no bundle path or
publication handle. The fit summary records measured maximum R chunk size,
cell/feature counts, and the returned latent-file checksum. The private
Python computation materializes its AnnData working matrix, so this pilot is
not a large-data or GPU qualification.

PeakVI `batch_key` remains deliberately restricted to null. By contrast,
`scvi_train` receives core-resolved, typed, exactly ordered batch,
categorical and continuous metadata in an authenticated snapshot. Its matrix
grant declares the nonnegative-integer policy and the private CPU scVI helper
cannot access bundle storage. Its fit publishes a latent space, contained
model and summary in one atomic closure. GPU is out of scope. Scientific
comparison to the old core uses a prospective
numerical policy (ordered cells and latent dimensions equal, finite
coordinates, and pairwise-distance Spearman correlation at least 0.9),
not byte- or hash-identity claims. Its frozen old-core fixture is in
`tests/fixtures/scvi_old_core_reference.json`; the pack comparator exercises
registered dispatch. The relocated PeakVI legacy tests under
`tests/migration_evidence/` document the retired route; they are not a second
production owner or part of the live pack test suite.
