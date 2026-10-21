.latent_runtime <- new.env(parent = emptyenv())
.latent_runtime$python <- NULL

.latent_python <- function() {
  configured <- Sys.getenv("SCWB_LATENT_PYTHON", unset = "")
  if (!nzchar(configured)) stop(
    "scwbLatent requires SCWB_LATENT_PYTHON before activation.",
    call. = FALSE)
  # Preserve the venv entry path: resolving its python symlink escapes the
  # virtual environment and silently changes the dependency set.
  python <- file.path(normalizePath(dirname(configured), winslash = "/",
    mustWork = TRUE), basename(configured))
  if (!file.exists(python) || dir.exists(python)) stop(
    "scwbLatent Python interpreter is absent.", call. = FALSE)
  probe <- suppressWarnings(system2(python,
    c("-c", shQuote("import scvi, anndata, scipy, torch; assert hasattr(scvi.model, 'PEAKVI') and hasattr(scvi.model, 'SCVI')")),
    stdout = TRUE, stderr = TRUE))
  status <- attr(probe, "status") %||% 0L
  if (!identical(as.integer(status), 0L)) stop(
    "scwbLatent Python dependency probe failed before dispatch: ",
    paste(probe, collapse = "\n"), call. = FALSE)
  python
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Register scwbLatent operations
#'
#' Activation verifies the selected private Python runtime before the
#' operation pack is allowed to enter the catalog.
#' @export
register_operation_pack <- function() {
  .latent_runtime$python <- .latent_python()
  scworkbench::wb_operation_pack_register(
    system.file("scworkbench", "operation-pack.json",
      package = "scwbLatent"),
    list(`scwblatent::scvi_train` = .latent_scvi_train,
      `scwblatent::peakvi_fit` = .latent_peakvi_fit,
      `scwblatent::peakvi_project` = .latent_peakvi_project,
      `scwblatent::model_bound_latent_union` =
        .latent_model_bound_latent_union))
}
