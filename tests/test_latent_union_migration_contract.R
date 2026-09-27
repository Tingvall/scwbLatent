peakvi_test_path <- if (file.exists(file.path("tests",
    "test_peakvi_pilot_contract.R"))) file.path("tests",
  "test_peakvi_pilot_contract.R") else "test_peakvi_pilot_contract.R"
source(peakvi_test_path, local = TRUE)

test_latent_union_retained_producer_contract <- function() {
  state <- .latent_test_fixture()
  on.exit(.latent_test_cleanup(state), add = TRUE)
  b <- state$bundle
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(), b,
    "scwblatent::peakvi_fit", inputs = list(matrix = "ATAC_reference",
      n_latent = 2L, n_layers = 1L, n_hidden = 16L, max_epochs = 1L,
      chunk_size = 5L, cell_set = "all", min_cells = 1L))
  fit_promoted <- scworkbench::wb_session_promote_current(b, fit$session,
    outputs = c("model", "fit_summary"), created_by = "union-contract")
  stopifnot(all(fit_promoted$promoted$status == "ok"))
  model_id <- fit_promoted$promoted$artifact_id[
    fit_promoted$promoted$output == "model"][[1L]]
  fit_summary_id <- fit_promoted$promoted$artifact_id[
    fit_promoted$promoted$output == "fit_summary"][[1L]]
  fit <- scworkbench::wb_session_discard(fit$session)
  model <- scworkbench::wb_art_load(b, model_id,
    expected_kind = "model_ref")
  project <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(), b,
    "scwblatent::peakvi_project", inputs = list(
      model = scworkbench::wb_ref_artifact(model_id, "model"),
      query_matrix = "ATAC_query", space_id_out = "union_query",
      chunk_size = 4L, cell_set = "all", min_cells = 1L,
      pickle_execution_trust = "operator_verified_source_v1"))
  stopifnot(!"union_query" %in% scworkbench::wb_spaces(b)$space_id,
    nrow(scworkbench::wb_space(scworkbench::wb_session_view(
      project$session), "union_query")) == 8L)
  scworkbench:::.wb_operation_pack_v3_full_space_assert_current(
    scworkbench::wb_session_view(project$session), "union_query")
  staged_lineage <- tryCatch(scworkbench::wb_session_run_op(
    project$session, b, "scwblatent::model_bound_latent_union",
    inputs = list(reference_space = model$model_space_id,
      query_space = "union_query",
      model = scworkbench::wb_ref_artifact(model_id, "model"),
      projection_summary = scworkbench::wb_ref_run(
        project$run$run_id, "projection_summary"),
      cell_set = "all", min_cells = 1L,
      space_id_out = "staged_joint", chunk_size = 4L)),
    error = identity)
  stopifnot(inherits(staged_lineage, "error"),
    !"staged_joint" %in% scworkbench::wb_spaces(b)$space_id)
  project_promoted <- scworkbench::wb_session_promote_current(b,
    project$session, outputs = "projection_summary",
    created_by = "union-contract")
  stopifnot(all(project_promoted$promoted$status == "ok"))
  summary_id <- project_promoted$promoted$artifact_id[
    project_promoted$promoted$output == "projection_summary"][[1L]]
  staged_summary <- project$run$outputs$projection_summary
  project <- scworkbench::wb_session_discard(project$session)
  inputs <- function(output = "union_joint", summary =
      scworkbench::wb_ref_artifact(summary_id, "projection_summary"),
      query = "union_query") list(reference_space = model$model_space_id,
        query_space = query,
        model = scworkbench::wb_ref_artifact(model_id, "model"),
        projection_summary = summary, cell_set = "all", min_cells = 1L,
        space_id_out = output, chunk_size = 4L)
  run <- function(...) scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), b,
    "scwblatent::model_bound_latent_union", inputs = inputs(...))
  inventory <- function() list(spaces = scworkbench::wb_spaces(b)$space_id,
    artifacts = nrow(scworkbench::wb_art_list(b)),
    heads = nrow(scworkbench::wb_prov_list_heads(b)))
  baseline <- inventory()
  stage_probe <- run(output = "stage_probe_joint")
  stage_bundle <- scworkbench::wb_session_view(stage_probe$session)
  fixture_export <- Sys.getenv("SCWB_LATENT_UNION_FIXTURE_EXPORT", unset = "")
  if (nzchar(fixture_export)) {
    dir.create(fixture_export, recursive = TRUE, showWarnings = FALSE)
    transition_id <- stage_probe$run$operation_pack_transition$transition_id
    transition_root <- file.path(stage_bundle$dir, "operation_pack_v2",
      transition_id)
    stopifnot(all(file.copy(file.path(transition_root,
      c("request.json", "receipt.json")), fixture_export,
      overwrite = TRUE)))
  }
  stage_payload <- DBI::dbGetQuery(stage_bundle$con,
    "SELECT path FROM space_payloads WHERE space_id='union_query';")
  stopifnot(nrow(stage_payload) == 1L)
  stage_path <- file.path(stage_bundle$dir, stage_payload$path[[1L]])
  committed_path <- file.path(b$dir, stage_payload$path[[1L]])
  committed_hash <- unname(tools::md5sum(committed_path))
  stopifnot(nzchar(Sys.readlink(stage_path)))
  private_copy <- tempfile("union-stage-space-", fileext = ".h5")
  stopifnot(file.copy(stage_path, private_copy), unlink(stage_path) == 0L,
    file.copy(private_copy, stage_path))
  unlink(private_copy)
  stage_values <- rhdf5::h5read(stage_path, "coordinates")
  stage_values[1L, 1L] <- stage_values[1L, 1L] + 1
  rhdf5::h5write(stage_values, stage_path, "coordinates")
  stopifnot(identical(unname(tools::md5sum(committed_path)),
    committed_hash))
  changed_stage <- tryCatch(scworkbench::wb_session_run_op(
    stage_probe$session, b, "scwblatent::model_bound_latent_union",
    inputs = inputs(output = "changed_stage_joint")), error = identity)
  stopifnot(inherits(changed_stage, "error"),
    grepl("typed-space|space source|content identity",
      conditionMessage(changed_stage), ignore.case = TRUE),
    identical(inventory(), baseline))
  DBI::dbExecute(stage_bundle$con,
    "UPDATE space_payloads SET content_sha256=? WHERE space_id='union_query';",
    params = list(scworkbench:::.wb_prov_sha256_file(stage_path)))
  rehashed_stage <- tryCatch(scworkbench::wb_session_run_op(
    stage_probe$session, b, "scwblatent::model_bound_latent_union",
    inputs = inputs(output = "rehashed_stage_joint")), error = identity)
  stopifnot(inherits(rehashed_stage, "error"),
    identical(inventory(), baseline))
  stage_probe <- scworkbench::wb_session_discard(stage_probe$session)
  forged <- scworkbench::wb_prim_obs_table(data.frame(
    engine = "peakvi", projection_mode = "transform_without_refit",
    output_space_id = "union_query", model_space_id = model$model_space_id,
    model_artifact_id = model_id, n_latent = 2L),
    table_type = "peakvi_project_summary", entity_kind = "space")
  for (candidate in list(forged, staged_summary)) {
    rejected <- tryCatch(run(output = "forged_joint", summary = candidate),
      error = identity)
    stopifnot(inherits(rejected, "error"),
      identical(inventory(), baseline))
  }
  options(scworkbench.operation_pack_v2_test_hook = function(event, ...)
    if (identical(event, "before_publication"))
      stop("union publication fault"))
  on.exit(options(scworkbench.operation_pack_v2_test_hook = NULL),
    add = TRUE)
  fault <- tryCatch(run(output = "fault_joint"), error = identity)
  options(scworkbench.operation_pack_v2_test_hook = NULL)
  stopifnot(inherits(fault, "error"),
    identical(inventory(), baseline))
  union <- run()
  stopifnot(identical(sort(names(union$run$outputs)),
    c("space", "union_summary")))
  reference <- scworkbench::wb_space(b, model$model_space_id)
  query <- scworkbench::wb_space(b, "union_query")
  expected_cells <- c(rownames(reference), rownames(query))
  expected_cells <- scworkbench::wb_cells(b, cell_set = "all")[
    scworkbench::wb_cells(b, cell_set = "all") %in% expected_cells]
  expected <- rbind(as.matrix(reference[, -1L, drop = FALSE]),
    as.matrix(query[, -1L, drop = FALSE]))
  expected <- expected[match(expected_cells, rownames(expected)), ,
    drop = FALSE]
  actual <- scworkbench::wb_space(scworkbench::wb_session_view(
    union$session), "union_joint")
  stopifnot(identical(rownames(actual), expected_cells),
    isTRUE(all.equal(unname(as.matrix(actual[, -1L, drop = FALSE])),
      unname(expected), tolerance = 0)))
  summary <- union$run$outputs$union_summary$table
  stopifnot(identical(summary$engine, "peakvi"),
    identical(summary$model_artifact_id, model_id),
    identical(summary$model_space_id, model$model_space_id),
    identical(summary$query_space_id, "union_query"),
    identical(summary$output_space_id, "union_joint"),
    identical(summary$n_reference_cells, 24L),
    identical(summary$n_query_cells, 8L),
    identical(summary$n_cells, 32L),
    identical(summary$n_dimensions, 2L),
    identical(summary$cells_disjoint, TRUE),
    identical(summary$coordinate_operation,
      "row_bind_without_transformation"),
    identical(summary$projection_receipt_type,
      "peakvi_project_summary"),
    identical(summary$projection_mode, "transform_without_refit"))
  promoted <- scworkbench::wb_session_promote_current(b,
    union$session, outputs = "union_summary", created_by = "union-contract")
  stopifnot(all(promoted$promoted$status == "ok"),
    identical(rownames(scworkbench::wb_space(b, "union_joint")),
      expected_cells),
    isTRUE(scworkbench::wb_finalize_recipient_bundle(b)$ready))
  scworkbench::wb_close_bundle(b)
  state$bundle <- NULL
  relocated <- tempfile("union-relocated-")
  fs::dir_copy(state$root, relocated)
  on.exit(unlink(relocated, recursive = TRUE, force = TRUE), add = TRUE)
  reopened <- scworkbench::wb_open_bundle(relocated, read_only = TRUE)
  on.exit(if (!is.null(reopened))
    try(scworkbench::wb_close_bundle(reopened), silent = TRUE), add = TRUE)
  stopifnot(isTRUE(scworkbench::wb_bundle_validation_report(relocated,
      level = "recipient")$ok),
    identical(rownames(scworkbench::wb_space(reopened, "union_joint")),
      expected_cells))
  scworkbench::wb_close_bundle(reopened)
  reopened <- NULL
  history_probe <- tempfile("union-pack-free-", fileext = ".R")
  on.exit(unlink(history_probe), add = TRUE)
  writeLines(c("library(scworkbench)",
    "root <- commandArgs(trailingOnly = TRUE)[[1L]]",
    "stopifnot(!'scwblatent::model_bound_latent_union' %in% wb_list_ops())",
    "bundle <- wb_open_bundle(root, read_only = TRUE)",
    "stopifnot(isTRUE(wb_bundle_validation_report(root, level = 'recipient')$ok),",
    "  nrow(wb_space(bundle, 'union_joint')) == 32L,",
    "  nrow(wb_prov_list_heads(bundle)) > 0L)",
    "wb_close_bundle(bundle)"), history_probe)
  probe <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(history_probe), shQuote(relocated)),
    stdout = TRUE, stderr = TRUE,
    env = paste0("R_LIBS=", paste(.libPaths(),
      collapse = .Platform$path.sep))))
  probe_status <- attr(probe, "status")
  if (is.null(probe_status)) probe_status <- 0L
  stopifnot(identical(as.integer(probe_status), 0L))
  state$bundle <- scworkbench::wb_open_bundle(state$root)
  b <- state$bundle
  tampered_root <- tempfile("union-tampered-")
  fs::dir_copy(state$root, tampered_root)
  on.exit(unlink(tampered_root, recursive = TRUE, force = TRUE),
    add = TRUE)
  tampered <- scworkbench::wb_open_bundle(tampered_root)
  payload <- DBI::dbGetQuery(tampered$con,
    "SELECT path FROM space_payloads WHERE space_id='union_query';")
  stopifnot(nrow(payload) == 1L)
  scworkbench::wb_close_bundle(tampered)
  h5 <- file.path(tampered_root, payload$path[[1L]])
  values <- rhdf5::h5read(h5, "coordinates")
  values[1L, 1L] <- values[1L, 1L] + 1
  rhdf5::h5write(values, h5, "coordinates")
  tampered <- scworkbench::wb_open_bundle(tampered_root)
  b <- tampered
  tampered_inventory <- inventory()
  changed_bytes <- tryCatch(run(output = "changed_joint"),
    error = identity)
  stopifnot(inherits(changed_bytes, "error"),
    identical(inventory(), tampered_inventory))
  scworkbench::wb_close_bundle(tampered)
  b <- state$bundle
  scworkbench::wb_entity_supersede(b, "artifact", fit_summary_id,
    summary_id, reason = "source supersession contract")
  after_supersede <- inventory()
  stale <- tryCatch(run(output = "stale_joint"), error = identity)
  stopifnot(inherits(stale, "error"),
    identical(inventory(), after_supersede),
    !isTRUE(scworkbench::wb_bundle_validation_report(state$root,
      level = "recipient")$ok))
  TRUE
}
