# scVI's method code lives in this pack. Core provides only authenticated
# matrix chunks, ordered cell metadata, and atomic output writers.

.latent_scvi_train <- function(readers, writers, params, context) {
  stage <- tempfile("scwblatent-scvi-private-")
  dir.create(stage, recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  staged <- .latent_stage_matrix(readers$chunks,
    context$typed_slice_reads$chunks, params$chunk_size, stage)
  metadata <- readers$cell_metadata()
  cells <- as.character(unlist(metadata$cell_id, use.names = FALSE))
  if (!identical(metadata$schema,
      "scworkbench.ordered_cell_metadata.snapshot.v1") ||
      !identical(cells, staged$cells) ||
      !identical(as.character(metadata$matrix_id),
        as.character(context$input_ids$matrix)))
    .latent_fail("ordered cell metadata differs from the matrix grant")
  names <- vapply(metadata$fields %||% list(), function(field)
    as.character(field$name)[1], character(1))
  if (anyDuplicated(names) || length(metadata$fields %||% list()) > 16L)
    .latent_fail("ordered cell metadata fields are duplicated or oversized")
  for (field in metadata$fields %||% list()) {
    if (length(field$values %||% list()) != length(cells) ||
        !identical(as.character(field$value_type)[1],
          if (identical(as.character(field$parameter)[1],
              "continuous_covariates")) "continuous" else "categorical"))
      .latent_fail("ordered cell metadata type or cardinality changed")
  }
  metadata_path <- file.path(stage, "cell_metadata.json")
  jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE,
    null = "null", digits = NA)
  staged$digests[["cell_metadata.json"]] <- digest::digest(
    file = metadata_path, algo = "sha256")
  request <- list(schema = "scwblatent.scvi_private.v1",
    operation = "fit", chunks = staged$chunks,
    features = "features.txt", metadata = "cell_metadata.json",
    n_latent = as.integer(params$n_latent),
    n_layers = as.integer(params$n_layers),
    n_hidden = as.integer(params$n_hidden),
    dropout = as.numeric(params$dropout),
    max_epochs = as.integer(params$max_epochs),
    lr = as.numeric(params$lr),
    early_stopping = isTRUE(params$early_stopping),
    seed = as.integer(params$seed),
    batch_key = params$batch_key,
    categorical_covariates = unname(as.list(as.character(
      params$categorical_covariates %||% character()))),
    continuous_covariates = unname(as.list(as.character(
      params$continuous_covariates %||% character()))))
  computed <- .latent_private_compute(request, stage, staged,
    script_name = "scvi_private.py")
  response <- computed$response
  if (!isTRUE(response$trained) ||
      !identical(as.character(response$operation)[1], "fit") ||
      !identical(as.character(response$model_output_dir)[1],
        "model_output") ||
      !identical(as.integer(response$matrix_nnz) >= 0L, TRUE) ||
      !identical(as.integer(response$effective_params$seed),
        as.integer(params$seed)) ||
      !identical(as.integer(response$effective_params$n_latent),
        as.integer(params$n_latent)) ||
      !identical(as.integer(response$effective_params$n_hidden),
        as.integer(params$n_hidden)) ||
      !identical(as.integer(response$effective_params$n_layers),
        as.integer(params$n_layers)) ||
      !identical(as.integer(response$effective_params$max_epochs),
        as.integer(params$max_epochs)) ||
      !isTRUE(all.equal(as.numeric(response$effective_params$dropout),
        as.numeric(params$dropout))) ||
      !isTRUE(all.equal(as.numeric(response$effective_params$lr),
        as.numeric(params$lr))) ||
      !identical(isTRUE(response$effective_params$early_stopping),
        isTRUE(params$early_stopping)) ||
      !identical(as.character(response$effective_params$batch_key %||% ""),
        as.character(params$batch_key %||% "")) ||
      !identical(as.character(response$effective_params$
        categorical_covariates %||% character()),
        as.character(params$categorical_covariates %||% character())) ||
      !identical(as.character(response$effective_params$
        continuous_covariates %||% character()),
        as.character(params$continuous_covariates %||% character())))
    .latent_fail("scVI response changed the effective training parameters")
  space <- writers$write_space("space", .latent_space_primitive(
    computed$values))
  writers$write_terminal_primitive("fit_summary",
    scworkbench::wb_prim_obs_table(data.frame(
      source_matrix_id = as.character(context$input_ids$matrix)[1],
      n_cells = as.integer(length(staged$cells)),
      n_features = as.integer(length(staged$features)),
      n_latent = as.integer(ncol(computed$values)),
      max_chunk_cells = as.integer(params$chunk_size),
      max_chunk_elements = as.numeric(staged$max_chunk_elements),
      latent_sha256 = computed$latent_sha256,
      metadata_sha256 = staged$digests[["cell_metadata.json"]],
      batch_key = as.character(params$batch_key %||% "")[1],
      categorical_covariates_json = jsonlite::toJSON(
        as.character(unlist(request$categorical_covariates,
          use.names = FALSE)), auto_unbox = FALSE),
      continuous_covariates_json = jsonlite::toJSON(
        as.character(unlist(request$continuous_covariates,
          use.names = FALSE)), auto_unbox = FALSE),
      max_epochs = as.integer(params$max_epochs),
      n_hidden = as.integer(params$n_hidden),
      n_layers = as.integer(params$n_layers),
      dropout = as.numeric(params$dropout),
      learning_rate = as.numeric(params$lr),
      early_stopping = isTRUE(params$early_stopping),
      seed = as.integer(params$seed),
      use_gpu = FALSE,
      scvi_version = as.character(response$scvi_version),
      stringsAsFactors = FALSE), table_type = "scvi_fit_summary"))
  root <- file.path(stage, "model_output")
  if (!dir.exists(root) || nzchar(Sys.readlink(root)))
    .latent_fail("scVI model tree is absent or linked")
  relative <- sort(list.files(root, recursive = TRUE, all.files = TRUE,
    no.. = TRUE), method = "radix")
  if (!length(relative) || !"model.pt" %in% relative)
    .latent_fail("scVI model tree lacks model.pt")
  tree <- writers$open_model_tree("model")
  for (entry in relative) {
    source <- file.path(root, entry)
    if (!file.exists(source) || dir.exists(source) ||
        nzchar(Sys.readlink(source)) ||
        any(strsplit(entry, "/", fixed = TRUE)[[1]] %in%
          c("", ".", "..")))
      .latent_fail("scVI model tree contains an unsafe entry")
    destination <- tree$path(entry)
    if (!isTRUE(file.copy(source, destination, overwrite = FALSE)) ||
        !identical(digest::digest(file = source, algo = "sha256"),
          digest::digest(file = destination, algo = "sha256")))
      .latent_fail("scVI model tree changed during contained publication")
  }
  model <- scworkbench::wb_prim_model_ref(space$allocated_id,
    space$allocated_id, engine = "scvi")
  model$model_type <- "scvi_train"
  tree$finish(model)
  invisible(list(status = "ok"))
}
