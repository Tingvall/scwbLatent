.latent_union_read_space <- function(read, capability) {
  cells <- as.character(capability$cells)
  width <- as.integer(capability$dimension_count)
  chunk_size <- as.integer(capability$requested_chunk_size)
  if (!length(cells) || width < 1L || chunk_size < 1L ||
      anyNA(cells) || anyDuplicated(cells))
    .latent_fail("union typed-space grant has an invalid full axis")
  starts <- seq.int(1L, length(cells), by = chunk_size)
  blocks <- lapply(starts, function(start) {
    selected <- cells[seq.int(start,
      min(length(cells), start + chunk_size - 1L))]
    chunk <- read(cells = selected)
    values <- as.matrix(chunk$values)
    if (!identical(as.character(chunk$cell_id), selected) ||
        !identical(as.integer(dim(values)),
          c(as.integer(length(selected)), width)) ||
        !all(is.finite(values)))
      .latent_fail("union typed-space chunk differs from its issued axis")
    values
  })
  values <- do.call(rbind, blocks)
  rownames(values) <- cells
  values
}

.latent_model_bound_latent_union <- function(readers, writers, params,
    context) {
  lineage <- readers$projection_lineage()
  model <- readers$fitted_model()
  projection <- readers$projection_summary()
  reference_grant <- context$typed_space_reads$reference_chunks
  query_grant <- context$typed_space_reads$query_chunks
  reference <- .latent_union_read_space(readers$reference_chunks,
    reference_grant)
  query <- .latent_union_read_space(readers$query_chunks, query_grant)
  linkage <- model$model_source_linkage %||% list()
  engine <- tolower(as.character(linkage$engine %||% "")[1])
  model_id <- as.character(lineage$model$artifact_id)[1]
  reference_id <- as.character(lineage$fitted_space$space_id)[1]
  query_id <- as.character(lineage$query_space$space_id)[1]
  if (!engine %in% c("peakvi", "scvi") ||
      !identical(as.character(linkage$paired_fitted_space_id %||%
        "")[1], reference_id) ||
      !nzchar(model_id) ||
      !identical(as.character(context$input_ids$reference_space)[1],
        reference_id) || !identical(as.character(context$input_ids$query_space)[1],
          query_id) || ncol(reference) != ncol(query) ||
      ncol(reference) < 2L || length(intersect(rownames(reference),
        rownames(query))))
    .latent_fail("union model, space axes, or dimensions disagree")
  table <- projection$table %||% data.frame()
  expected_type <- paste0(engine, "_project_summary")
  no_refit <- identical(as.character(table$projection_mode %||% "")[1],
      "transform_without_refit") ||
    identical(as.character(table$transformation %||% "")[1],
      "transfer_without_refit")
  declared_dimensions <- suppressWarnings(as.integer(table$n_latent %||%
    NA_integer_))[1]
  if (!is.data.frame(table) || nrow(table) != 1L ||
      !identical(as.character(projection$table_type %||% "")[1],
        expected_type) || !no_refit || is.na(declared_dimensions) ||
      declared_dimensions != ncol(query) ||
      (!is.null(table$engine) && !identical(tolower(as.character(
        table$engine)[1]), engine)) ||
      (!is.null(table$output_space_id) && !identical(as.character(
        table$output_space_id)[1], query_id)) ||
      (!is.null(table$model_artifact_id) && !identical(as.character(
        table$model_artifact_id)[1], model_id)) ||
      (!is.null(table$model_space_id) && !identical(as.character(
        table$model_space_id)[1], reference_id)))
    .latent_fail("projection summary does not attest the linked no-refit transfer")
  source_cells <- c(rownames(reference), rownames(query))
  final_cells <- as.character(context$effective_cells)
  if (length(final_cells) != length(source_cells) ||
      anyDuplicated(final_cells) || !setequal(final_cells, source_cells))
    .latent_fail("effective cells are not the complete source-space union")
  colnames(query) <- colnames(reference)
  combined <- rbind(reference, query)
  combined <- combined[match(final_cells, rownames(combined)), , drop = FALSE]
  output_id <- as.character(params$space_id_out)[1]
  writers$write_space("space", .latent_space_primitive(combined, output_id))
  summary <- data.frame(engine = engine, model_artifact_id = model_id,
    model_space_id = reference_id, reference_space_id = reference_id,
    query_space_id = query_id, output_space_id = output_id,
    n_reference_cells = as.integer(nrow(reference)),
    n_query_cells = as.integer(nrow(query)),
    n_cells = as.integer(length(final_cells)),
    n_dimensions = as.integer(ncol(combined)), cells_disjoint = TRUE,
    coordinate_operation = "row_bind_without_transformation",
    projection_receipt_type = expected_type,
    projection_mode = "transform_without_refit",
    stringsAsFactors = FALSE)
  writers$write_terminal_primitive("union_summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "model_bound_latent_union_summary",
      entity_kind = "space", meta = list(engine = engine,
        model_artifact_id = model_id, model_space_id = reference_id,
        coordinate_operation = "row_bind_without_transformation",
        projection_receipt_type = expected_type,
        projection_mode = "transform_without_refit")))
  declaration <- scworkbench::wb_derived_cell_output_axis_declaration(
    context$derived_cell_output_axis_grants[[1L]], final_cells,
    member_selections = list(rownames(reference), rownames(query)))
  list(status = "ok", derived_cell_output_axis_declarations =
    list(declaration))
}
