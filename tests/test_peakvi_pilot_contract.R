.latent_test_fixture <- function() {
  stopifnot(!any(c("peakvi_fit", "peakvi_project") %in%
    scworkbench::wb_list_ops()))
  library <- getOption("scwblatent.phase6.fixture_library", NULL)
  if (is.null(library)) {
    installed <- tryCatch(find.package("scwbLatent"),
      error = function(e) "")
    if (nzchar(installed)) library <- dirname(installed)
  }
  if (is.null(library) || !dir.exists(library)) {
    library <- tempfile("latent-test-lib-"); dir.create(library)
    package_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
    log <- tempfile("latent-install-", fileext = ".log")
    status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL",
      paste0("--library=", shQuote(library)), shQuote(package_root)),
      stdout = log, stderr = log)
    if (!identical(as.integer(status), 0L))
      stop(paste(readLines(log, warn = FALSE), collapse = "\n"))
    options(scwblatent.phase6.fixture_library = library)
  }
  old <- .libPaths(); .libPaths(c(library, old))
  if (inherits(tryCatch(scworkbench::wb_operation_pack_inspect(
      "scwblatent"), error = identity), "error")) {
    manifest <- system.file("scworkbench", "operation-pack.json",
      package = "scwbLatent", lib.loc = library)
    scworkbench::wb_operation_pack_activate(manifest, lib.loc = library)
  }
  registered <- scworkbench::wb_list_ops()
  stopifnot(sum(registered == "scwblatent::peakvi_fit") == 1L,
    sum(registered == "scwblatent::peakvi_project") == 1L)
  root <- tempfile("latent-bundle-")
  bundle <- scworkbench::wb_init_bundle(root)
  reference <- paste0("ref_", seq_len(24L))
  query <- paste0("query_", seq_len(8L))
  features <- paste0("peak_", seq_len(32L))
  scworkbench::wb_define_cells(bundle, c(reference, query))
  values <- matrix((seq_len(length(features) * length(reference)) %% 5L) + 1L,
    nrow = length(features), dimnames = list(features, reference))
  scworkbench::wb_attach_matrix(bundle, "ATAC_reference",
    Matrix::Matrix(values, sparse = TRUE), feature_type = "peak",
    cell_ids = reference)
  query_values <- matrix((seq_len(length(features) * length(query)) %% 7L) + 1L,
    nrow = length(features), dimnames = list(features, query))
  scworkbench::wb_attach_matrix(bundle, "ATAC_query",
    Matrix::Matrix(query_values, sparse = TRUE), feature_type = "peak",
    cell_ids = query)
  list(library = library, old = old, root = root, bundle = bundle,
    reference = reference, query = query, features = features)
}

.latent_test_cleanup <- function(state) {
  try(scworkbench::wb_close_bundle(state$bundle), silent = TRUE)
  # Once executed, the package remains active for other tests in this process.
  .libPaths(state$old)
  unlink(state$root, recursive = TRUE, force = TRUE)
}

test_peakvi_cpu_fit_project_contract <- function() {
  state <- .latent_test_fixture()
  on.exit(.latent_test_cleanup(state), add = TRUE)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    state$bundle, "scwblatent::peakvi_fit", inputs = list(
      matrix = "ATAC_reference", n_latent = 2L, n_layers = 1L,
      n_hidden = 16L, max_epochs = 1L, chunk_size = 5L,
      cell_set = "all", min_cells = 1L))
  stopifnot(identical(sort(names(fit$run$outputs)),
    c("fit_summary", "model", "space")))
  promoted <- scworkbench::wb_session_promote_current(state$bundle,
    fit$session, outputs = c("model", "fit_summary"),
    created_by = "phase6-pilot")
  stopifnot(all(promoted$promoted$status == "ok"))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  fit$session <- scworkbench::wb_session_discard(fit$session)
  model <- scworkbench::wb_art_load(state$bundle, model_id,
    expected_kind = "model_ref")
  stopifnot(identical(model$model_linkage$source_matrix_id,
    "ATAC_reference"), identical(as.integer(model$model_linkage$
    latent_dimension_count), 2L))
  fit_summary_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "fit_summary"][[1L]]
  fit_summary <- scworkbench::wb_art_load(state$bundle, fit_summary_id,
    expected_kind = "obs_table")$table
  stopifnot(identical(as.integer(fit_summary$max_chunk_elements), 160L),
    as.integer(fit_summary$max_chunk_cells) == 5L,
    as.integer(fit_summary$n_cells) == 24L)
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    state$bundle)$ready))
  relocated <- tempfile("latent-relocated-")
  on.exit(unlink(relocated, recursive = TRUE, force = TRUE), add = TRUE)
  scworkbench::wb_close_bundle(state$bundle)
  fs::dir_copy(state$root, relocated)
  state$bundle <- scworkbench::wb_open_bundle(relocated, read_only = TRUE)
  recipient <- scworkbench::wb_bundle_validation_report(relocated,
    level = "recipient")
  if (!isTRUE(recipient$ok))
    stop(paste(recipient$diagnostics$code, collapse = ", "))
  stopifnot(isTRUE(recipient$ok))
  reopened <- scworkbench::wb_art_load(state$bundle, model_id,
    expected_kind = "model_ref")
  stopifnot(identical(reopened$model_linkage,
    model$model_linkage))
  scworkbench::wb_close_bundle(state$bundle)
  state$bundle <- scworkbench::wb_open_bundle(relocated)
  project_once <- function(query_matrix = "ATAC_query",
      trust = "operator_verified_source_v1")
    scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    state$bundle, "scwblatent::peakvi_project", inputs = list(
      model = scworkbench::wb_ref_artifact(model_id, "model"),
      query_matrix = query_matrix, space_id_out = "peakvi_query",
      chunk_size = 4L, cell_set = "all", min_cells = 1L,
      pickle_execution_trust = trust))
  baseline <- list(artifacts = nrow(scworkbench::wb_art_list(state$bundle)),
    spaces = scworkbench::wb_spaces(state$bundle)$space_id)
  for (trust in list(NULL, "unverified_source")) {
    rejected <- tryCatch(project_once(trust = trust), error = identity)
    stopifnot(inherits(rejected, "error"),
      identical(nrow(scworkbench::wb_art_list(state$bundle)),
        baseline$artifacts),
      identical(scworkbench::wb_spaces(state$bundle)$space_id,
        baseline$spaces))
  }
  overlapping <- tryCatch(project_once("ATAC_reference"), error = identity)
  stopifnot(inherits(overlapping, "error"),
    grepl("overlaps the fitted reference axis",
      conditionMessage(overlapping), fixed = TRUE),
    identical(nrow(scworkbench::wb_art_list(state$bundle)),
      baseline$artifacts),
    identical(scworkbench::wb_spaces(state$bundle)$space_id,
      baseline$spaces))
  options(scworkbench.operation_pack_v2_test_hook = function(event, ...)
    if (identical(event, "before_publication"))
      stop("phase6 pilot publication fault"))
  on.exit(options(scworkbench.operation_pack_v2_test_hook = NULL),
    add = TRUE)
  failed <- tryCatch(project_once(), error = identity)
  options(scworkbench.operation_pack_v2_test_hook = NULL)
  stopifnot(inherits(failed, "error"),
    identical(nrow(scworkbench::wb_art_list(state$bundle)),
      baseline$artifacts),
    identical(scworkbench::wb_spaces(state$bundle)$space_id,
      baseline$spaces))
  project <- project_once()
  stopifnot(identical(sort(names(project$run$outputs)),
    c("projection_summary", "space")))
  projected <- scworkbench::wb_session_promote_current(state$bundle,
    project$session, outputs = "projection_summary",
    created_by = "phase6-pilot")
  stopifnot(all(projected$promoted$status == "ok"),
    "peakvi_query" %in% scworkbench::wb_spaces(state$bundle)$space_id,
    identical(as.integer(dim(scworkbench::wb_space(state$bundle,
      "peakvi_query"))), c(8L, 3L)))
  project$session <- scworkbench::wb_session_discard(project$session)
  collision <- tryCatch(project_once(), error = identity)
  stopifnot(inherits(collision, "error"),
    sum(scworkbench::wb_spaces(state$bundle)$space_id ==
      "peakvi_query") == 1L)
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    state$bundle)$ready))
  scworkbench::wb_close_bundle(state$bundle)
  state$bundle <- scworkbench::wb_open_bundle(relocated, read_only = TRUE)
  stopifnot(isTRUE(scworkbench::wb_bundle_validation_report(relocated,
    level = "recipient")$ok))
  TRUE
}

test_peakvi_registered_source_authentication_contract <- function() {
  state <- .latent_test_fixture()
  on.exit(.latent_test_cleanup(state), add = TRUE)
  manifest <- jsonlite::read_json(system.file("scworkbench",
    "operation-pack.json", package = "scwbLatent", lib.loc = state$library),
    simplifyVector = FALSE)
  fit_definition <- manifest$operations[[which(vapply(manifest$operations,
    function(x) identical(x$name, "peakvi_fit"), logical(1)))]]
  project_definition <- manifest$operations[[which(vapply(manifest$operations,
    function(x) identical(x$name, "peakvi_project"), logical(1)))]]
  model_output <- fit_definition$outputs[[which(vapply(
    fit_definition$outputs, function(x) identical(x$name, "model"),
    logical(1)))]]
  stopifnot(identical(model_output$
      matrix_trained_contained_model_closure_policy$policy,
      "matrix_trained_contained_model_closure.v1"),
    identical(unlist(model_output$
      matrix_trained_contained_model_closure_policy$
      required_sibling_terminals, use.names = FALSE),
      c("space", "fit_summary")),
    identical(fit_definition$promotion_closure$model,
      list("fit_summary", "model", "space")),
    identical(project_definition$self_linked_contained_model_reads[[1L]]$
      linkage_policy, "authenticated_existing_source_and_cell_linkage_v1"))
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    state$bundle, "scwblatent::peakvi_fit", inputs = list(
      matrix = "ATAC_reference", n_latent = 2L, n_layers = 1L,
      n_hidden = 16L, max_epochs = 1L, chunk_size = 5L,
      cell_set = "all", min_cells = 1L))
  promoted <- scworkbench::wb_session_promote_current(state$bundle,
    fit$session, outputs = c("model", "fit_summary"),
    created_by = "phase6-source-authentication")
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  scworkbench::wb_session_discard(fit$session)
  scworkbench::wb_close_bundle(state$bundle)
  state$bundle <- NULL
  for (attack in c("matrix_bytes", "model_tree", "linkage")) {
    copy <- tempfile(paste0("latent-auth-", attack, "-"))
    fs::dir_copy(state$root, copy)
    on.exit(unlink(copy, recursive = TRUE, force = TRUE), add = TRUE)
    damaged <- scworkbench::wb_open_bundle(copy)
    if (identical(attack, "matrix_bytes")) {
      fragment <- DBI::dbGetQuery(damaged$con,
        "SELECT path FROM matrix_fragments WHERE matrix_id='ATAC_reference';")
      stopifnot(nrow(fragment) == 1L)
      target <- file.path(copy, fragment$path[[1L]], "val")
      stopifnot(file.exists(target))
      connection <- file(target, open = "ab")
      writeBin(as.raw(1L), connection); close(connection)
    } else if (identical(attack, "model_tree")) {
      model <- scworkbench::wb_art_load(damaged, model_id,
        expected_kind = "model_ref")
      target <- file.path(model$model_tree_root, model$model_primary_path)
      stopifnot(file.exists(target))
      Sys.chmod(target, mode = "0644")
      connection <- file(target, open = "ab")
      writeBin(as.raw(1L), connection); close(connection)
    } else {
      stopifnot(DBI::dbExecute(damaged$con,
        "UPDATE artifact_model_ref SET linkage_json='{}' WHERE artifact_id=?;",
        params = list(model_id)) == 1L)
    }
    inventory <- function() {
      files <- list.files(copy, recursive = TRUE, all.files = TRUE,
        full.names = TRUE, no.. = TRUE)
      files <- sort(files[file.exists(files) & !dir.exists(files)],
        method = "radix")
      # A failed preflight may refresh the source manifest or checkpoint
      # SQLite's WAL. Compare logical table counts separately below.
      files <- files[!substring(files, nchar(copy) + 2L) %in% c(
        "manifest.json", "db/app.sqlite", "db/app.sqlite-wal",
        "db/app.sqlite-shm")]
      tables <- sort(DBI::dbListTables(damaged$con), method = "radix")
      row_counts <- stats::setNames(vapply(tables, function(table)
        as.integer(DBI::dbGetQuery(damaged$con,
          paste0('SELECT COUNT(*) AS n FROM "', table, '";'))$n[[1L]]),
        integer(1)), tables)
      list(artifacts = nrow(scworkbench::wb_art_list(damaged)),
        spaces = scworkbench::wb_spaces(damaged)$space_id,
        heads = nrow(scworkbench::wb_prov_list_heads(damaged)),
        row_counts = row_counts,
        files = stats::setNames(unname(tools::md5sum(files)),
          substring(files, nchar(copy) + 2L)))
    }
    before <- inventory()
    rejected <- tryCatch(scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), damaged,
      "scwblatent::peakvi_project", inputs = list(
        model = scworkbench::wb_ref_artifact(model_id, "model"),
        query_matrix = "ATAC_query", space_id_out = "peakvi_query",
        chunk_size = 4L, cell_set = "all", min_cells = 1L,
        pickle_execution_trust = "operator_verified_source_v1")),
      error = identity)
    after <- inventory()
    if (!identical(after, before)) {
      changed <- union(names(before$files)[before$files != after$files[names(before$files)]],
        setdiff(names(after$files), names(before$files)))
      stop("PeakVI ", attack, " rejection left changed files: ",
        paste(changed, collapse = ", "), "; output counts: ",
        paste(unlist(before[1:3]), collapse = "/"), " -> ",
        paste(unlist(after[1:3]), collapse = "/"))
    }
    stopifnot(inherits(rejected, "error"))
    scworkbench::wb_close_bundle(damaged)
  }
  TRUE
}

test_peakvi_staged_receipt_closure_rejection_contract <- function() {
  for (attack in c("receipt_identity", "missing_sibling")) {
    state <- .latent_test_fixture()
    inventory <- function() {
      files <- list.files(state$root, recursive = TRUE, all.files = TRUE,
        full.names = TRUE, no.. = TRUE)
      files <- sort(files[file.exists(files) & !dir.exists(files)],
        method = "radix")
      relative <- substring(files, nchar(state$root) + 2L)
      files <- files[!relative %in% c("manifest.json", "db/app.sqlite",
        "db/app.sqlite-wal", "db/app.sqlite-shm")]
      tables <- sort(DBI::dbListTables(state$bundle$con), method = "radix")
      list(artifacts = scworkbench::wb_art_list(state$bundle)$artifact_id,
        spaces = scworkbench::wb_spaces(state$bundle)$space_id,
        heads = scworkbench::wb_prov_list_heads(state$bundle),
        row_counts = stats::setNames(vapply(tables, function(table)
          as.integer(DBI::dbGetQuery(state$bundle$con,
            paste0('SELECT COUNT(*) AS n FROM "', table, '";'))$n[[1L]]),
          integer(1)), tables),
        files = stats::setNames(unname(tools::md5sum(files)),
          substring(files, nchar(state$root) + 2L)))
    }
    before <- inventory()
    fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
      state$bundle, "scwblatent::peakvi_fit", inputs = list(
        matrix = "ATAC_reference", n_latent = 2L, n_layers = 1L,
        n_hidden = 16L, max_epochs = 1L, chunk_size = 5L,
        cell_set = "all", min_cells = 1L))
    record <- fit$run$operation_pack_transition
    path <- file.path(scworkbench:::.wb_operation_pack_v2_store(
      fit$session$staging_bundle), record$transition_id, "receipt.json")
    receipt <- jsonlite::read_json(path, simplifyVector = FALSE)
    if (identical(attack, "receipt_identity")) {
      receipt$receipt_identity <- paste0("sha256:", strrep("0", 64L))
    } else {
      receipt$outputs <- Filter(function(output)
        !identical(output$name, "fit_summary"), receipt$outputs)
      stopifnot(length(receipt$outputs) == 2L)
    }
    jsonlite::write_json(receipt, path, auto_unbox = TRUE, null = "null")
    rejected <- tryCatch(scworkbench::wb_session_promote_current(
      state$bundle, fit$session, outputs = c("model", "fit_summary"),
      created_by = "phase6-receipt-fault"), error = identity)
    stopifnot(inherits(rejected, "error"), identical(inventory(), before))
    scworkbench::wb_session_discard(fit$session)
    .latent_test_cleanup(state)
  }
  TRUE
}

# Prospective numerical policy: exact ordered cells/dimensions, finite
# nonconstant latent coordinates, finite nondegenerate pairwise geometry, and
# pairwise-distance Spearman >= 0.95 between old-core and pack embeddings.
# Rotations, signs, coordinates, model bytes and hashes are not equated.
test_peakvi_old_core_numerical_policy_contract <- function() {
  source <- Sys.getenv("SCWB_LATENT_LEGACY_CORE_SOURCE", unset = "")
  if (!nzchar(source) || !dir.exists(source))
    stop("Set SCWB_LATENT_LEGACY_CORE_SOURCE to the accepted pilot core checkout.")
  revision <- trimws(system2("git", c("-C", shQuote(source), "rev-parse",
    "HEAD"), stdout = TRUE))
  stopifnot(identical(revision,
    "d624cfd0cc1b84440681967fc30a214004e2c591"))
  python <- Sys.getenv("SCWB_LATENT_PYTHON", unset = "")
  stopifnot(nzchar(python), file.exists(python))
  oracle <- tempfile("old-core-peakvi-", fileext = ".rds")
  log <- tempfile("old-core-peakvi-", fileext = ".log")
  on.exit(unlink(c(oracle, log), force = TRUE), add = TRUE)
  script <- normalizePath("tests/composite/old_core_peakvi_oracle.R",
    winslash = "/", mustWork = TRUE)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), shQuote(source), shQuote(oracle)),
    stdout = log, stderr = log,
    env = paste0("SCWORKBENCH_PYBENCH_PYTHON=", shQuote(python)))
  if (!identical(as.integer(status), 0L))
    stop(paste(readLines(log, warn = FALSE), collapse = "\n"))
  legacy <- readRDS(oracle)
  state <- .latent_test_fixture()
  on.exit(.latent_test_cleanup(state), add = TRUE)
  modern <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    state$bundle, "scwblatent::peakvi_fit", inputs = list(
      matrix = "ATAC_reference", n_latent = 2L, n_layers = 1L,
      n_hidden = 16L, max_epochs = 1L, chunk_size = 5L,
      cell_set = "all", min_cells = 1L))
  new_values <- as.matrix(scworkbench::wb_space(
    scworkbench::wb_session_view(modern$session),
    modern$run$outputs$space$space_id)[-1L])
  old_values <- legacy$coordinates
  stopifnot(identical(dim(old_values), c(24L, 2L)),
    identical(dim(new_values), dim(old_values)),
    identical(legacy$cells, state$reference),
    identical(rownames(old_values), state$reference),
    identical(rownames(new_values), state$reference),
    all(is.finite(old_values)), all(is.finite(new_values)),
    all(apply(old_values, 2L, stats::sd) > 0),
    all(apply(new_values, 2L, stats::sd) > 0),
    all(is.finite(stats::dist(old_values))),
    all(is.finite(stats::dist(new_values))),
    as.numeric(stats::median(stats::dist(old_values))) > 0,
    as.numeric(stats::median(stats::dist(new_values))) > 0)
  distance_rank_correlation <- suppressWarnings(stats::cor(
    as.numeric(stats::dist(old_values)),
    as.numeric(stats::dist(new_values)), method = "spearman"))
  message("PeakVI old/new pairwise-distance Spearman = ",
    format(distance_rank_correlation, digits = 5L))
  stopifnot(is.finite(distance_rank_correlation),
    distance_rank_correlation >= 0.95)
  scworkbench::wb_session_discard(modern$session)
  TRUE
}

test_a_peakvi_missing_dependency_and_pack_predispatch_contract <- function() {
  state <- .latent_test_fixture()
  on.exit(.latent_test_cleanup(state), add = TRUE)
  payload_inventory <- function() {
    paths <- list.files(state$root, recursive = TRUE, all.files = TRUE,
      full.names = TRUE, no.. = TRUE)
    paths <- sort(paths[file.exists(paths) & !dir.exists(paths)],
      method = "radix")
    stats::setNames(unname(tools::md5sum(paths)),
      substring(paths, nchar(state$root) + 2L))
  }
  before <- list(artifacts = nrow(scworkbench::wb_art_list(state$bundle)),
    spaces = scworkbench::wb_spaces(state$bundle)$space_id,
    heads = nrow(scworkbench::wb_prov_list_heads(state$bundle)),
    files = payload_inventory())
  scworkbench::wb_operation_pack_unload("scwblatent")
  stopifnot(!any(c("peakvi_fit", "peakvi_project",
    "scwblatent::peakvi_fit", "scwblatent::peakvi_project") %in%
    scworkbench::wb_list_ops()))
  absent <- tryCatch(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), state$bundle,
    "scwblatent::peakvi_fit", inputs = list(matrix = "ATAC_reference",
      n_latent = 2L, cell_set = "all", min_cells = 1L)), error = identity)
  stopifnot(inherits(absent, "error"))
  absent_project <- tryCatch(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), state$bundle,
    "scwblatent::peakvi_project", inputs = list()), error = identity)
  stopifnot(inherits(absent_project, "error"))
  original <- Sys.getenv("SCWB_LATENT_PYTHON", unset = NA_character_)
  on.exit(if (is.na(original)) Sys.unsetenv("SCWB_LATENT_PYTHON")
    else Sys.setenv(SCWB_LATENT_PYTHON = original), add = TRUE)
  Sys.setenv(SCWB_LATENT_PYTHON = file.path(state$root,
    "missing-python"))
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbLatent", lib.loc = state$library)
  missing <- tryCatch(scworkbench::wb_operation_pack_activate(manifest,
    lib.loc = state$library), error = identity)
  stopifnot(inherits(missing, "error"),
    !"scwblatent::peakvi_fit" %in% as.character(
      unlist(scworkbench::wb_list_ops(), use.names = FALSE)),
    !"scwblatent::peakvi_project" %in% scworkbench::wb_list_ops(),
    identical(nrow(scworkbench::wb_art_list(state$bundle)),
      before$artifacts),
    identical(scworkbench::wb_spaces(state$bundle)$space_id,
      before$spaces),
    identical(nrow(scworkbench::wb_prov_list_heads(state$bundle)),
      before$heads),
    identical(payload_inventory(), before$files))
  TRUE
}
