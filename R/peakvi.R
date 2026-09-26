.latent_space_primitive <- function(values, space_id = NULL) {
  value <- list(kind = "space", values = values,
    dimensions = as.integer(dim(values)), cell_id = rownames(values),
    storage_profile = "space_file_backed_hdf5_v1")
  if (!is.null(space_id)) value$space_id <- as.character(space_id)[1]
  structure(value, class = "wb_primitive")
}

.latent_fit_request <- function(params, staged) list(
  schema = "scwblatent.peakvi_private.v1", operation = "fit",
  chunks = staged$chunks, features = "features.txt",
  n_latent = as.integer(params$n_latent),
  n_layers = as.integer(params$n_layers),
  n_hidden = as.integer(params$n_hidden),
  dropout = as.numeric(params$dropout),
  max_epochs = as.integer(params$max_epochs),
  lr = as.numeric(params$lr),
  early_stopping = isTRUE(params$early_stopping),
  seed = as.integer(params$seed))

.latent_peakvi_fit <- function(readers, writers, params, context) {
  stage <- tempfile("scwblatent-private-")
  dir.create(stage, recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  staged <- .latent_stage_matrix(readers$chunks,
    context$typed_slice_reads$chunks, params$chunk_size, stage)
  computed <- .latent_private_compute(.latent_fit_request(params, staged),
    stage, staged)
  if (!isTRUE(computed$response$trained) ||
      !identical(as.character(computed$response$operation)[1], "fit") ||
      !identical(as.character(computed$response$model_output_dir)[1],
        "model_output"))
    .latent_fail("fit response did not attest training and contained model output")
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
      max_chunk_bytes = as.numeric(staged$max_chunk_bytes),
      total_matrix_elements = as.numeric(staged$total_elements),
      latent_sha256 = computed$latent_sha256,
      scvi_version = as.character(computed$response$scvi_version),
      stringsAsFactors = FALSE), table_type = "peakvi_fit_summary"))
  root <- file.path(stage, "model_output")
  if (!dir.exists(root) || nzchar(Sys.readlink(root)))
    .latent_fail("fit model tree is absent or linked")
  relative <- sort(list.files(root, recursive = TRUE, all.files = TRUE,
    no.. = TRUE), method = "radix")
  if (!length(relative) || !"model.pt" %in% relative)
    .latent_fail("fit model tree lacks its primary model.pt")
  tree <- writers$open_model_tree("model")
  for (entry in relative) {
    source <- file.path(root, entry)
    if (!file.exists(source) || dir.exists(source) ||
        nzchar(Sys.readlink(source)) ||
        any(strsplit(entry, "/", fixed = TRUE)[[1]] %in% c("", ".", "..")))
      .latent_fail("private model tree contains an unsafe entry")
    destination <- tree$path(entry)
    if (!isTRUE(file.copy(source, destination, overwrite = FALSE)))
      .latent_fail("could not copy private model entry to issued writer")
    if (!identical(digest::digest(file = source, algo = "sha256"),
        digest::digest(file = destination, algo = "sha256")))
      .latent_fail("model tree changed while closing core writer")
  }
  model <- scworkbench::wb_prim_model_ref(space$allocated_id,
    space$allocated_id, engine = "peakvi")
  model$model_type <- "peakvi"
  tree$finish(model)
  invisible(list(status = "ok"))
}

.latent_peakvi_project <- function(readers, writers, params, context) {
  model <- readers$fitted_model()
  if (is.null(model$fitted_cells) ||
      is.null(model$fitted_features) ||
      !identical(as.character(model$model_source_linkage$schema)[1],
        "scworkbench.operation_pack.matrix_model_source_linkage.v1"))
    .latent_fail("project requires a matrix-trained authenticated G1 model")
  capability <- context$typed_slice_reads$query_chunks
  query_cells <- as.character(capability$cells)
  if (length(intersect(query_cells, as.character(model$fitted_cells))))
    .latent_fail("project query overlaps the fitted reference axis")
  if (!identical(as.character(capability$features),
      as.character(model$fitted_features)))
    .latent_fail("project query feature axis differs from the fitted model")
  stage <- tempfile("scwblatent-private-")
  dir.create(stage, recursive = TRUE)
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  staged <- .latent_stage_matrix(readers$query_chunks, capability,
    params$chunk_size, stage)
  .latent_stage_model_read(model, stage)
  request <- list(schema = "scwblatent.peakvi_private.v1",
    operation = "project", chunks = staged$chunks, features = "features.txt",
    n_latent = as.integer(model$latent_dimension_count),
    model_input_dir = "model_input")
  computed <- .latent_private_compute(request, stage, staged)
  if (isTRUE(computed$response$trained) ||
      !identical(as.character(computed$response$operation)[1], "project") ||
      !identical(as.integer(ncol(computed$values)),
        as.integer(model$latent_dimension_count)))
    .latent_fail("projection response attempted a refit or changed latent dimensions")
  writers$write_space("space", .latent_space_primitive(computed$values,
    as.character(params$space_id_out)[1]))
  writers$write_terminal_primitive("projection_summary",
    scworkbench::wb_prim_obs_table(data.frame(
      source_matrix_id = as.character(context$input_ids$query_matrix)[1],
      fitted_matrix_id = as.character(
        model$model_source_linkage$source_matrix_id)[1],
      n_query_cells = as.integer(length(query_cells)),
      n_latent = as.integer(ncol(computed$values)),
      transformation = "transfer_without_refit",
      model_load_method = as.character(
        computed$response$model_load_method)[1],
      latent_sha256 = computed$latent_sha256,
      stringsAsFactors = FALSE), table_type = "peakvi_project_summary"))
  invisible(list(status = "ok"))
}
