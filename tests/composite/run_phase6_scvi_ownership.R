required <- c("SCWB_LATENT_CORE_SOURCE", "SCWB_LATENT_PACK_SOURCE",
  "SCWB_LATENT_PYBENCH_SOURCE", "SCWB_LATENT_CORE_COMMIT",
  "SCWB_LATENT_PACK_COMMIT", "SCWB_LATENT_PYBENCH_COMMIT",
  "SCWB_LATENT_PYTHON")
values <- Sys.getenv(required, unset = "")
if (any(!nzchar(values))) stop("Missing composite pin(s): ",
  paste(required[!nzchar(values)], collapse = ", "))
names(values) <- required
revision <- function(path) trimws(system2("git",
  c("-C", shQuote(path), "rev-parse", "HEAD"), stdout = TRUE))
stopifnot(identical(revision(values[["SCWB_LATENT_CORE_SOURCE"]]),
    values[["SCWB_LATENT_CORE_COMMIT"]]),
  identical(revision(values[["SCWB_LATENT_PACK_SOURCE"]]),
    values[["SCWB_LATENT_PACK_COMMIT"]]),
  identical(revision(values[["SCWB_LATENT_PYBENCH_SOURCE"]]),
    values[["SCWB_LATENT_PYBENCH_COMMIT"]]))
protocol <- c("operation-pack-v3.schema.json",
  "operation-pack-transition-request-v3.schema.json",
  "operation-pack-transition-receipt-v3.schema.json",
  file.path("fixtures", "operation-pack-v3",
    "phase6-ordered-cell-metadata-identities.json"))
for (name in protocol) {
  r <- file.path(values[["SCWB_LATENT_CORE_SOURCE"]], "inst", "protocol", name)
  py <- file.path(values[["SCWB_LATENT_PYBENCH_SOURCE"]], "protocol", name)
  stopifnot(file.exists(r), file.exists(py),
    identical(unname(tools::md5sum(r)), unname(tools::md5sum(py))))
}
setwd(values[["SCWB_LATENT_PACK_SOURCE"]])
library(scworkbench)
library(scwbLatent)
source(file.path("tests", "test_scvi_train_pack_contract.R"))
tests <- c("test_a_scvi_missing_pack_and_dependency_predispatch_contract",
  "test_b_scvi_historical_pack_free_read_only_inspection_contract",
  "test_scvi_cpu_fit_metadata_and_promotion_contract",
  "test_scvi_separate_categorical_and_continuous_covariate_contract",
  "test_scvi_count_and_covariate_authority_rejection_contract",
  "test_scvi_publication_fault_zero_residue_and_retry_contract",
  "test_scvi_incomplete_model_tree_rejects_promotion_contract",
  "test_scvi_durable_process_replacement_and_model_integrity_contract")
for (name in tests) {
  cat("RUN ", name, "\n", sep = "")
  stopifnot(isTRUE(get(name)()))
  cat("PASS ", name, "\n", sep = "")
}
cat("PHASE6 SCVI OWNERSHIP COMPOSITE passed=", length(tests),
  " failed=0 skipped=0\n", sep = "")
