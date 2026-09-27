#!/usr/bin/env Rscript

library(scworkbench)
source("tests/test_latent_union_migration_contract.R")
stopifnot(test_latent_union_retained_producer_contract())
cat("PHASE6_LATENT_UNION_OWNERSHIP_OK\n")
