required <- c("SCWB_LATENT_CORE_SOURCE", "SCWB_LATENT_PACK_SOURCE",
  "SCWB_LATENT_LEGACY_CORE_SOURCE", "SCWB_LATENT_PYBENCH_SOURCE",
  "SCWB_LATENT_CORE_COMMIT", "SCWB_LATENT_PACK_COMMIT",
  "SCWB_LATENT_PYBENCH_COMMIT", "SCWB_LATENT_PYTHON")
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
    values[["SCWB_LATENT_PYBENCH_COMMIT"]]),
  identical(revision(values[["SCWB_LATENT_LEGACY_CORE_SOURCE"]]),
    "d624cfd0cc1b84440681967fc30a214004e2c591"))
setwd(values[["SCWB_LATENT_PACK_SOURCE"]])
library(scworkbench)
library(scwbLatent)
source(file.path(values[["SCWB_LATENT_PACK_SOURCE"]], "tests",
  "test_peakvi_pilot_contract.R"))
tests <- c("test_a_peakvi_missing_dependency_and_pack_predispatch_contract",
  "test_peakvi_cpu_fit_project_contract",
  "test_peakvi_model_linked_projection_widths_contract",
  "test_peakvi_registered_source_authentication_contract",
  "test_peakvi_staged_receipt_closure_rejection_contract",
  "test_peakvi_old_core_numerical_policy_contract")
for (name in tests) {
  cat("RUN ", name, "\n", sep = "")
  stopifnot(isTRUE(get(name)()))
  cat("PASS ", name, "\n", sep = "")
}
cat("PHASE6 OWNERSHIP COMPOSITE passed=", length(tests),
  " failed=0 skipped=0\n", sep = "")
