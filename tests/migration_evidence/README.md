# Retired core PeakVI tests

These files are preserved method-owned migration evidence from the accepted
old-core PeakVI implementation. They deliberately are not executed by the
current pack suite: their direct core operation and PyBench-v2 entry points
were removed from the no-pack production catalog. The live pack lifecycle,
missing-pack/dependency, numerical-comparator and private Python contracts
reside in `tests/test_peakvi_pilot_contract.R`,
`tests/composite/old_core_peakvi_oracle.R`, and
`tests/test_peakvi_private.py`.
