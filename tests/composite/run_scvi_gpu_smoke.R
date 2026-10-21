# Run only in an allocated GPU job with scwbLatent and the pinned core installed.
python <- Sys.getenv("SCWB_LATENT_PYTHON", unset = "")
if (!nzchar(python) || !file.exists(python))
  stop("SCWB_LATENT_PYTHON must name the allocated job's CUDA scVI interpreter")
probe <- suppressWarnings(system2(python, c("-c", shQuote(paste(
  "import torch, scvi;",
  "assert torch.cuda.is_available() and torch.cuda.device_count() == 1,",
  "'expected exactly one scheduler-visible CUDA device';",
  "print(torch.__version__, torch.version.cuda, torch.cuda.get_device_name(0), scvi.__version__)"))),
  stdout = TRUE, stderr = TRUE))
if (!is.null(attr(probe, "status")) && attr(probe, "status") != 0L)
  stop(paste(probe, collapse = "\n"))
cat("CUDA_RUNTIME ", paste(probe, collapse = " "), "\n", sep = "")

library(scworkbench)
library(scwbLatent)
scworkbench::wb_operation_pack_activate(system.file("scworkbench",
  "operation-pack.json", package = "scwbLatent"),
  lib.loc = dirname(find.package("scwbLatent")))

run_scvi_gpu_smoke <- function() {
  root <- tempfile("scvi-gpu-smoke-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  cells <- paste0("cell_", c(2L, 10L, 1L,
    setdiff(seq_len(24L), c(2L, 10L, 1L))))
  features <- paste0("gene_", c(3L, 1L,
    setdiff(seq_len(32L), c(3L, 1L))))
  values <- matrix(as.integer((seq_len(length(cells) *
    length(features)) * 7L) %% 11L), nrow = length(features),
    dimnames = list(features, cells))
  bundle <- scworkbench::wb_init_bundle(root)
  on.exit(try(scworkbench::wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  scworkbench::wb_define_cells(bundle, cells, metadata = data.frame(
    cell_id = cells, sample_id = rep(c("sample_A", "sample_B"),
      each = length(cells) / 2L), stringsAsFactors = FALSE))
  scworkbench::wb_attach_matrix(bundle, "RNA_counts",
    Matrix::Matrix(values, sparse = TRUE), layer = "counts",
    feature_type = "gene", cell_ids = cells)
  registered <- as.matrix(scworkbench::wb_matrix(bundle, "RNA_counts"))
  stopifnot(identical(rownames(registered), features),
    identical(colnames(registered), cells),
    isTRUE(all.equal(unname(registered), unname(values),
      check.attributes = FALSE)))
  model_ids <- character()
  for (case in list(list(label = "unbatched", batch_key = NULL),
      list(label = "sample_id", batch_key = "sample_id"))) {
    fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
      bundle, "scwblatent::scvi_train", inputs = list(
        matrix = "RNA_counts", cell_set = "all", min_cells = 1L,
        batch_key = case$batch_key, use_gpu = TRUE, max_epochs = 2L,
        early_stopping = FALSE, n_latent = 2L, n_layers = 1L,
        n_hidden = 16L, chunk_size = 6L))
    latent <- scworkbench::wb_space(fit$session$staging_bundle,
      fit$run$outputs$space$space_id, cells = cells)
    if (is.data.frame(latent)) latent <- as.matrix(latent[,
      setdiff(names(latent), "cell_id"), drop = FALSE])
    stopifnot(identical(as.integer(dim(latent)), c(length(cells), 2L)),
      all(is.finite(latent)),
      isTRUE(fit$run$inputs$user_inputs$use_gpu))
    promoted <- scworkbench::wb_session_promote_current(bundle,
      fit$session, outputs = c("model", "fit_summary"),
      created_by = "scvi-gpu-smoke")
    stopifnot(all(promoted$promoted$status == "ok"))
    model_id <- promoted$promoted$artifact_id[
      promoted$promoted$output == "model"][[1L]]
    summary_id <- promoted$promoted$artifact_id[
      promoted$promoted$output == "fit_summary"][[1L]]
    summary <- scworkbench::wb_art_load(bundle, summary_id,
      expected_kind = "obs_table")$table
    model <- scworkbench::wb_art_load(bundle, model_id,
      expected_kind = "model_ref")
    stopifnot(isTRUE(summary$use_gpu[[1L]]),
      identical(summary$observed_training_device[[1L]], "cuda:0"),
      as.integer(summary$completed_epochs[[1L]]) >= 1L,
      as.integer(summary$completed_epochs[[1L]]) <= 2L,
      as.integer(summary$training_batches_observed[[1L]]) >= 1L,
      nzchar(summary$torch_cuda_version[[1L]]),
      nzchar(summary$gpu_name[[1L]]),
      identical(as.integer(summary$n_cells[[1L]]), length(cells)),
      identical(as.integer(summary$n_features[[1L]]), length(features)),
      identical(as.character(summary$batch_key[[1L]]),
        if (is.null(case$batch_key)) "" else case$batch_key),
      identical(model$model_linkage$source_matrix_id, "RNA_counts"))
    model_ids <- c(model_ids, model_id)
    cat("GPU_FIT ", case$label, " epochs=",
      summary$completed_epochs[[1L]], " device=",
      summary$observed_training_device[[1L]], "\n", sep = "")
  }
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(bundle)$ready))
  scworkbench::wb_close_bundle(bundle)
  bundle <- scworkbench::wb_open_bundle(root, read_only = TRUE)
  on.exit(try(scworkbench::wb_close_bundle(bundle), silent = TRUE), add = TRUE)
  report <- scworkbench::wb_bundle_validation_report(root, level = "recipient")
  stopifnot(isTRUE(report$ok), length(model_ids) == 2L,
    !anyDuplicated(model_ids))
  for (model_id in model_ids)
    stopifnot(identical(scworkbench::wb_art_load(bundle, model_id,
      expected_kind = "model_ref")$model_type, "scvi_train"))
  cat("SCVI_GPU_SMOKE_PASS\n")
  TRUE
}

stopifnot(run_scvi_gpu_smoke())
