# Retired core latent tests

These files are preserved method-owned migration evidence from the accepted
old-core PeakVI implementation. They deliberately are not executed by the
current pack suite: their direct core operation and PyBench-v2 entry points
were removed from the no-pack production catalog. The live pack lifecycle,
missing-pack/dependency, numerical-comparator and private Python contracts
reside in `tests/test_peakvi_pilot_contract.R`,
`tests/composite/old_core_peakvi_oracle.R`, and
`tests/test_peakvi_private.py`.

`legacy_model_bound_latent_union_contract.R` is a byte-for-byte copy of the
old core union contract at scworkbench
`b564d10cdb5aee757a0c8401ab15f69594d12ea0`. It records the frozen
row-bind, summary, disjointness, and dimension expectations. The live
registered union comparison is in `tests/test_latent_union_migration_contract.R`;
the old core test is not run after ownership migration.
