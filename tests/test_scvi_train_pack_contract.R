.latent_scvi_fixture <- function() {
  library <- getOption("scwblatent.phase6.fixture_library", NULL)
  if (is.null(library)) library <- dirname(find.package("scwbLatent"))
  if (inherits(tryCatch(scworkbench::wb_operation_pack_inspect(
      "scwblatent"), error = identity), "error")) {
    manifest <- system.file("scworkbench", "operation-pack.json",
      package = "scwbLatent", lib.loc = library)
    scworkbench::wb_operation_pack_activate(manifest, lib.loc = library)
  }
  root <- tempfile("scvi-pack-bundle-")
  bundle <- scworkbench::wb_init_bundle(root)
  cells <- paste0("rna_", seq_len(24L))
  features <- paste0("gene_", seq_len(32L))
  metadata <- data.frame(cell_id = cells,
    batch = rep(c("batch_b", "batch_a"), each = 12L),
    condition = rep(c("stim", "ctrl"), 12L),
    quality = rep(seq(0.1, 0.8, length.out = 8L), 3L),
    stringsAsFactors = FALSE)
  scworkbench::wb_define_cells(bundle, cells, metadata = metadata)
  values <- matrix((seq_len(length(features) * length(cells)) %% 11L),
    nrow = length(features), dimnames = list(features, cells))
  scworkbench::wb_attach_matrix(bundle, "RNA_counts",
    Matrix::Matrix(values, sparse = TRUE), layer = "counts",
    feature_type = "gene", cell_ids = cells)
  list(bundle = bundle, root = root, cells = cells, features = features,
    metadata = metadata)
}

.latent_scvi_fixture_close <- function(fixture) {
  try(scworkbench::wb_close_bundle(fixture$bundle), silent = TRUE)
  unlink(fixture$root, recursive = TRUE, force = TRUE)
}

.latent_scvi_source_only <- function(bundle) {
  count <- function(table) as.integer(DBI::dbGetQuery(bundle$con,
    paste0("SELECT COUNT(*) AS n FROM ", table, ";"))$n[[1L]])
  list(artifacts = nrow(scworkbench::wb_art_list(bundle)),
    spaces = as.character(scworkbench::wb_spaces(bundle)$space_id),
    runs = count("prov_runs"),
    output_bindings = count("prov_run_outputs"),
    materializations = count("prov_materializations"),
    models = count("artifact_model_ref"),
    space_payloads = count("space_payloads"),
    transitions = if (dir.exists(file.path(bundle$dir,
      "operation_pack_v2"))) sort(list.files(file.path(bundle$dir,
        "operation_pack_v2"), recursive = TRUE, all.files = TRUE,
        no.. = TRUE)) else character())
}

test_a_scvi_missing_pack_and_dependency_predispatch_contract <- function() {
  # This test runs first in a fresh package process, before activation.
  stopifnot(!"scwblatent::scvi_train" %in% scworkbench::wb_list_ops())
  root <- tempfile("scvi-missing-pack-")
  bundle <- scworkbench::wb_init_bundle(root)
  on.exit({try(scworkbench::wb_close_bundle(bundle), silent = TRUE);
    unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
  baseline <- .latent_scvi_source_only(bundle)
  missing <- tryCatch(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), bundle, "scwblatent::scvi_train",
    inputs = list(matrix = "absent", cell_set = "all", min_cells = 1L)),
    error = identity)
  stopifnot(inherits(missing, "error"),
    identical(.latent_scvi_source_only(bundle), baseline))
  original <- Sys.getenv("SCWB_LATENT_PYTHON", unset = NA_character_)
  on.exit(if (is.na(original)) Sys.unsetenv("SCWB_LATENT_PYTHON") else
    Sys.setenv(SCWB_LATENT_PYTHON = original), add = TRUE)
  Sys.setenv(SCWB_LATENT_PYTHON = file.path(root, "absent-python"))
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbLatent")
  absent_python <- tryCatch(scworkbench::wb_operation_pack_activate(manifest,
    lib.loc = dirname(find.package("scwbLatent"))), error = identity)
  stopifnot(inherits(absent_python, "error"),
    !"scwblatent::scvi_train" %in% scworkbench::wb_list_ops(),
    identical(.latent_scvi_source_only(bundle), baseline))
  TRUE
}

test_b_scvi_historical_pack_free_read_only_inspection_contract <- function() {
  stopifnot(!"scwblatent::scvi_train" %in% scworkbench::wb_list_ops())
  archive <- file.path("tests", "fixtures",
    "phase6_scvi_old_core_historical.bundle.tar.gz")
  stopifnot(identical(unname(digest::digest(file = archive,
    algo = "sha256")),
    "be1114ecf1be9415f2744c963a9fcd1582f967c1f3c2bf8789ccffc278920a28"))
  root <- tempfile("scvi-historical-read-only-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(archive, exdir = root)
  path <- file.path(root, "phase6_scvi_old_core_historical.bundle")
  fingerprint <- function() {
    entries <- list.files(path, recursive = TRUE, full.names = TRUE,
      all.files = TRUE, no.. = TRUE)
    files <- sort(entries[!dir.exists(entries)], method = "radix")
    stats::setNames(vapply(files, function(file) digest::digest(
      file = file, algo = "sha256"), character(1)),
      substring(files, nchar(path) + 2L))
  }
  before <- fingerprint()
  b <- scworkbench::wb_open_bundle(path, read_only = TRUE)
  on.exit(try(scworkbench::wb_close_bundle(b), silent = TRUE), add = TRUE)
  artifacts <- scworkbench::wb_art_list(b)
  stopifnot(nrow(artifacts) == 1L)
  artifact_id <- as.character(artifacts$artifact_id[[1L]])
  model <- scworkbench::wb_art_load(b, artifact_id,
    expected_kind = "model_ref")
  cert <- scworkbench::wb_entity_certificate(b, "artifact", artifact_id)
  dag <- scworkbench::wb_prov_bundle_dag(b)
  stopifnot(identical(as.character(cert$producer$operation$name),
      "scvi_train"),
    identical(model$model_type, "scvi_train"),
    file.exists(file.path(model$model_path, "model.txt")),
    nrow(dag$nodes) > 0L)
  scworkbench::wb_close_bundle(b)
  b <- NULL
  report <- scworkbench::wb_bundle_validation_report(path, "recipient")
  stopifnot(isTRUE(report$ok),
    identical(report$pre_validation_fingerprint,
      report$post_validation_fingerprint),
    identical(before, fingerprint()),
    !"scwblatent::scvi_train" %in% scworkbench::wb_list_ops())
  TRUE
}

test_scvi_cpu_fit_metadata_and_promotion_contract <- function() {
  fixture <- .latent_scvi_fixture()
  on.exit(.latent_scvi_fixture_close(fixture), add = TRUE)
  stopifnot(sum(scworkbench::wb_list_ops() ==
    "scwblatent::scvi_train") == 1L)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    fixture$bundle, "scwblatent::scvi_train", inputs = list(
      matrix = "RNA_counts", n_latent = 2L, n_layers = 1L,
      n_hidden = 16L, max_epochs = 1L, chunk_size = 6L,
      batch_key = "batch", categorical_covariates = "condition",
      continuous_covariates = "quality",
      cell_set = "all", min_cells = 1L))
  stopifnot(identical(sort(names(fit$run$outputs)),
    c("fit_summary", "model", "space")))
  reference <- jsonlite::read_json(file.path("tests", "fixtures",
    "scvi_old_core_reference.json"), simplifyVector = TRUE)
  observed <- scworkbench::wb_space(fit$session$staging_bundle,
    fit$run$outputs$space$space_id, cells = fixture$cells)
  if (is.data.frame(observed)) {
    stopifnot("cell_id" %in% names(observed),
      identical(as.character(observed$cell_id), fixture$cells))
    observed <- as.matrix(observed[, setdiff(names(observed), "cell_id"),
      drop = FALSE])
  }
  reference_latent <- as.matrix(reference$latent)
  stopifnot(identical(as.character(reference$cell_id), fixture$cells),
    identical(as.integer(dim(observed)), as.integer(dim(reference_latent))),
    all(is.finite(observed)))
  distance_rank <- suppressWarnings(stats::cor(
    as.numeric(stats::dist(observed)),
    as.numeric(stats::dist(reference_latent)), method = "spearman"))
  stopifnot(is.finite(distance_rank), distance_rank >=
    as.numeric(reference$policy$distance_rank_spearman_minimum))
  cat("SCVI OLD-CORE DISTANCE-RANK SPEARMAN ",
    format(distance_rank, digits = 6), "\n", sep = "")
  transition_id <- as.character(fit$run$operation_pack_transition$
    transition_id)[1]
  promoted <- scworkbench::wb_session_promote_current(fixture$bundle,
    fit$session, outputs = c("model", "fit_summary"),
    created_by = "phase6-scvi-pack")
  stopifnot(all(promoted$promoted$status == "ok"))
  model_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "model"][[1L]]
  summary_id <- promoted$promoted$artifact_id[
    promoted$promoted$output == "fit_summary"][[1L]]
  summary <- scworkbench::wb_art_load(fixture$bundle, summary_id,
    expected_kind = "obs_table")$table
  stopifnot(identical(summary$source_matrix_id[[1]], "RNA_counts"),
    identical(summary$batch_key[[1]], "batch"),
    identical(as.integer(summary$seed[[1]]), 1L),
    identical(as.integer(summary$n_cells[[1]]), 24L),
    identical(as.integer(summary$n_features[[1]]), 32L),
    identical(as.integer(summary$n_latent[[1]]), 2L),
    identical(as.integer(summary$n_hidden[[1]]), 16L),
    identical(as.integer(summary$n_layers[[1]]), 1L),
    identical(as.integer(summary$max_epochs[[1]]), 1L),
    identical(as.integer(summary$max_chunk_cells[[1]]), 6L),
    as.numeric(summary$max_chunk_elements[[1]]) <= 32 * 6,
    isTRUE(all.equal(as.numeric(summary$dropout[[1]]), 0.1)),
    isTRUE(all.equal(as.numeric(summary$learning_rate[[1]]), 0.001)),
    isTRUE(summary$early_stopping[[1]]),
    isFALSE(summary$use_gpu[[1]]),
    identical(as.character(summary$categorical_covariates_json[[1]]),
      '["condition"]'),
    identical(as.character(summary$continuous_covariates_json[[1]]),
      '["quality"]'),
    nzchar(as.character(summary$scvi_version[[1]])),
    grepl("^[0-9a-f]{64}$", summary$latent_sha256[[1]]),
    grepl("^[0-9a-f]{64}$", summary$metadata_sha256[[1]]))
  model <- scworkbench::wb_art_load(fixture$bundle, model_id,
    expected_kind = "model_ref")
  stopifnot(identical(model$engine, "scvi"),
    identical(model$model_type, "scvi_train"),
    identical(as.character(model$model_linkage$source_matrix_id),
      "RNA_counts"))
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    fixture$bundle)$ready))
  relocated <- tempfile("scvi-pack-relocated-")
  on.exit(unlink(relocated, recursive = TRUE, force = TRUE), add = TRUE)
  scworkbench::wb_close_bundle(fixture$bundle)
  fs::dir_copy(fixture$root, relocated)
  fixture$bundle <- scworkbench::wb_open_bundle(relocated,
    read_only = TRUE)
  recipient <- scworkbench::wb_bundle_validation_report(relocated,
    level = "recipient")
  if (!isTRUE(recipient$ok)) stop(paste(recipient$diagnostics$code,
    collapse = ", "))
  stopifnot(isTRUE(recipient$ok),
    identical(scworkbench::wb_art_load(fixture$bundle, model_id,
      expected_kind = "model_ref")$model_linkage,
      model$model_linkage))
  retained <- scworkbench::wb_operation_pack_transition_registry_reads(
    fixture$bundle, transition_id)
  snapshot <- retained$cell_metadata$snapshot
  stopifnot(inherits(snapshot, "wb_operation_pack_ordered_cell_metadata"),
    identical(as.character(unlist(snapshot$cell_id, use.names = FALSE)),
      fixture$cells),
    identical(vapply(snapshot$fields, `[[`, character(1), "name"),
      c("batch", "condition", "quality")),
    identical(as.character(unlist(snapshot$fields[[1]]$values)),
      fixture$metadata$batch),
    identical(as.character(unlist(snapshot$fields[[2]]$values)),
      fixture$metadata$condition),
    isTRUE(all.equal(as.numeric(unlist(snapshot$fields[[3]]$values)),
      fixture$metadata$quality)),
    identical(as.integer(retained$cell_metadata$observation$use_count), 1L))
  TRUE
}

test_scvi_separate_categorical_and_continuous_covariate_contract <- function() {
  fixture <- .latent_scvi_fixture()
  on.exit(.latent_scvi_fixture_close(fixture), add = TRUE)
  allocated <- character()
  for (case in list(
      list(categorical_covariates = "condition",
        continuous_covariates = character()),
      list(categorical_covariates = character(),
        continuous_covariates = "quality"))) {
    fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
      fixture$bundle, "scwblatent::scvi_train", inputs = c(list(
        matrix = "RNA_counts", n_latent = 2L, n_layers = 1L,
        n_hidden = 16L, max_epochs = 1L, chunk_size = 6L,
        batch_key = "batch", cell_set = "all", min_cells = 1L), case))
    stopifnot(identical(sort(names(fit$run$outputs)),
      c("fit_summary", "model", "space")),
      identical(as.character(fit$run$inputs$user_inputs$batch_key),
        "batch"))
    allocated <- c(allocated,
      as.character(fit$run$outputs$space$space_id))
    fit$session <- scworkbench::wb_session_discard(fit$session)
  }
  stopifnot(length(allocated) == 2L, !anyDuplicated(allocated))
  TRUE
}

test_scvi_count_and_covariate_authority_rejection_contract <- function() {
  fixture <- .latent_scvi_fixture()
  on.exit(.latent_scvi_fixture_close(fixture), add = TRUE)
  values <- matrix(seq_len(length(fixture$features) * length(
    fixture$cells)) / 7, nrow = length(fixture$features),
    dimnames = list(fixture$features, fixture$cells))
  scworkbench::wb_attach_matrix(fixture$bundle, "RNA_fractional",
    Matrix::Matrix(values, sparse = TRUE), layer = "counts",
    feature_type = "gene")
  scworkbench::wb_attach_matrix(fixture$bundle, "RNA_log_norm",
    Matrix::Matrix(round(values), sparse = TRUE), layer = "lognorm",
    feature_type = "gene")
  negative <- matrix(1, nrow = length(fixture$features),
    ncol = length(fixture$cells), dimnames = list(fixture$features,
      fixture$cells))
  negative[1L, 1L] <- -1
  scworkbench::wb_attach_matrix(fixture$bundle, "RNA_negative",
    Matrix::Matrix(negative, sparse = TRUE), layer = "counts",
    feature_type = "gene")
  baseline <- .latent_scvi_source_only(fixture$bundle)
  run <- function(matrix = "RNA_counts", batch_key = "batch",
      categorical_covariates = "condition",
      continuous_covariates = "quality")
    scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
      fixture$bundle, "scwblatent::scvi_train", inputs = list(
        matrix = matrix, n_latent = 2L, n_layers = 1L,
        n_hidden = 16L, max_epochs = 1L, chunk_size = 6L,
        batch_key = batch_key,
        categorical_covariates = categorical_covariates,
        continuous_covariates = continuous_covariates,
        cell_set = "all", min_cells = 1L))
  for (case in list(
      list(matrix = "RNA_log_norm"), list(matrix = "RNA_fractional"),
      list(matrix = "RNA_negative"),
      list(batch_key = "missing_field"),
      list(continuous_covariates = "condition"),
      list(categorical_covariates = "batch"))) {
    failure <- tryCatch(do.call(run, case), error = identity)
    stopifnot(inherits(failure, "error"),
      identical(.latent_scvi_source_only(fixture$bundle), baseline))
  }
  DBI::dbExecute(fixture$bundle$con,
    "UPDATE cells SET quality=NULL WHERE cell_id=?;",
    params = list(fixture$cells[[1L]]))
  missing_value <- tryCatch(run(), error = identity)
  stopifnot(inherits(missing_value, "error"),
    identical(.latent_scvi_source_only(fixture$bundle), baseline))
  TRUE
}

test_scvi_publication_fault_zero_residue_and_retry_contract <- function() {
  fixture <- .latent_scvi_fixture()
  on.exit({options(scworkbench.operation_pack_v2_test_hook = NULL);
    .latent_scvi_fixture_close(fixture)}, add = TRUE)
  baseline <- .latent_scvi_source_only(fixture$bundle)
  run <- function() scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), fixture$bundle,
    "scwblatent::scvi_train", inputs = list(
      matrix = "RNA_counts", n_latent = 2L, n_layers = 1L,
      n_hidden = 16L, max_epochs = 1L, chunk_size = 6L,
      batch_key = "batch", categorical_covariates = "condition",
      continuous_covariates = "quality", cell_set = "all",
      min_cells = 1L))
  options(scworkbench.operation_pack_v2_test_hook = function(event, ...)
    if (identical(event, "before_publication"))
      stop("scVI injected publication failure"))
  failure <- tryCatch(run(), error = identity)
  options(scworkbench.operation_pack_v2_test_hook = NULL)
  stopifnot(inherits(failure, "error"),
    identical(.latent_scvi_source_only(fixture$bundle), baseline))
  retry <- run()
  stopifnot(identical(sort(names(retry$run$outputs)),
    c("fit_summary", "model", "space")))
  retry$session <- scworkbench::wb_session_discard(retry$session)
  stopifnot(identical(.latent_scvi_source_only(fixture$bundle), baseline))
  TRUE
}

test_scvi_incomplete_model_tree_rejects_promotion_contract <- function() {
  fixture <- .latent_scvi_fixture()
  on.exit(.latent_scvi_fixture_close(fixture), add = TRUE)
  baseline <- .latent_scvi_source_only(fixture$bundle)
  fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    fixture$bundle, "scwblatent::scvi_train", inputs = list(
      matrix = "RNA_counts", n_latent = 2L, n_layers = 1L,
      n_hidden = 16L, max_epochs = 1L, chunk_size = 6L,
      batch_key = "batch", cell_set = "all", min_cells = 1L))
  model_file <- file.path(fit$session$staging_bundle$dir,
    fit$run$outputs$model$model_path, "model.pt")
  stopifnot(file.exists(model_file))
  unlink(model_file)
  failure <- tryCatch(scworkbench::wb_session_promote_current(
    fixture$bundle, fit$session, outputs = "model",
    created_by = "scvi-incomplete-tree-contract"), error = identity)
  stopifnot(inherits(failure, "error"),
    identical(.latent_scvi_source_only(fixture$bundle), baseline))
  fit$session <- scworkbench::wb_session_discard(fit$session)
  TRUE
}

test_scvi_durable_process_replacement_and_model_integrity_contract <- function() {
  fixture <- .latent_scvi_fixture()
  on.exit(.latent_scvi_fixture_close(fixture), add = TRUE)
  workspace_root <- tempfile("scvi-pack-durable-")
  dir.create(workspace_root)
  on.exit(unlink(workspace_root, recursive = TRUE, force = TRUE), add = TRUE)
  session <- scworkbench::wb_session_create(fixture$bundle,
    workspace_root, "scvi_pack_recovery")
  fit <- scworkbench::wb_session_run_op(session, fixture$bundle,
    "scwblatent::scvi_train", inputs = list(
      matrix = "RNA_counts", n_latent = 2L, n_layers = 1L,
      n_hidden = 16L, max_epochs = 1L, chunk_size = 6L,
      batch_key = "batch", categorical_covariates = "condition",
      continuous_covariates = "quality", cell_set = "all",
      min_cells = 1L))
  run_id <- as.character(fit$run$run_id)
  expected_space_id <- as.character(fit$run$outputs$space$space_id)
  workspace <- fit$session$durable$workspace
  expected <- fit$session$durable$payload_identities[[run_id]]
  stopifnot(all(c("space", "model", "fit_summary") %in% names(expected)))
  scworkbench::wb_close_bundle(fixture$bundle)
  fixture$bundle <- NULL
  rm(session, fit)
  gc()
  result <- tempfile("scvi-pack-recovered-", fileext = ".rds")
  log <- tempfile("scvi-pack-recovered-", fileext = ".log")
  on.exit(unlink(c(result, log), force = TRUE), add = TRUE)
  code <- paste0(
    "library(scworkbench);",
    "manifest<-system.file('scworkbench','operation-pack.json',package='scwbLatent');",
    "wb_operation_pack_activate(manifest,lib.loc=dirname(find.package('scwbLatent')));",
    "b<-wb_open_bundle(", deparse(fixture$root), ");",
    "s<-wb_session_resume(b,", deparse(workspace), ");",
    "stopifnot(identical(names(s$runs),", deparse(run_id), "));",
    "p<-wb_session_promote_current(b,s,outputs=c('model','fit_summary'),",
    "created_by='scvi-pack-recovery');",
    "stopifnot(all(p$promoted$status=='ok'),",
    "identical(p$session$durable$state,'promoted'));",
    "f<-wb_finalize_recipient_bundle(b);stopifnot(isTRUE(f$ready));",
    "wb_close_bundle(b);",
    "v<-wb_bundle_validation_report(", deparse(fixture$root), ",'recipient');",
    "stopifnot(isTRUE(v$ok),identical(v$pre_validation_fingerprint,",
    "v$post_validation_fingerprint));",
    "saveRDS(p$promoted,", deparse(result), ")")
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "-e", shQuote(code)), stdout = log, stderr = log)
  if (!identical(as.integer(status), 0L) || !file.exists(result))
    stop(paste(readLines(log, warn = FALSE), collapse = "\n"),
      call. = FALSE)
  promoted <- readRDS(result)
  fixture$bundle <- scworkbench::wb_open_bundle(fixture$root,
    read_only = TRUE)
  model_id <- promoted$artifact_id[promoted$output == "model"][[1L]]
  summary_id <- promoted$artifact_id[
    promoted$output == "fit_summary"][[1L]]
  model <- scworkbench::wb_art_load(fixture$bundle, model_id,
    expected_kind = "model_ref")
  summary <- scworkbench::wb_art_load(fixture$bundle, summary_id,
    expected_kind = "obs_table")
  latent <- scworkbench::wb_space(fixture$bundle,
    expected_space_id, cells = fixture$cells)
  stopifnot(nrow(latent) == length(fixture$cells),
    identical(expected$model$content_sha256,
      scworkbench:::.wb_prov_payload_identity(fixture$bundle,
        model)$content_sha256),
    identical(expected$fit_summary$content_sha256,
      scworkbench:::.wb_prov_payload_identity(fixture$bundle,
        summary)$content_sha256))
  entity_id <- DBI::dbGetQuery(fixture$bundle$con, "
    SELECT binding.entity_id FROM prov_run_outputs output
    JOIN prov_run_entity_bindings binding
      ON binding.snapshot_id=output.snapshot_id
     AND binding.run_id=output.run_id
     AND binding.binding_name=output.output_name
    WHERE output.artifact_id=? AND binding.entity_kind='model_ref';",
    params = list(model_id))$entity_id[[1L]]
  presence <- scworkbench:::.wb_model_ref_physical_presence(
    fixture$bundle, entity_id)
  stopifnot(isTRUE(presence$present),
    file.exists(file.path(presence$absolute_model_path, "model.pt")))
  damaged <- tempfile("scvi-pack-missing-model-")
  on.exit(unlink(damaged, recursive = TRUE, force = TRUE), add = TRUE)
  fs::dir_copy(fixture$root, damaged)
  unlink(file.path(damaged, presence$model_path, "model.pt"))
  invalid <- scworkbench::wb_bundle_validation_report(damaged,
    level = "content")
  if (!isFALSE(invalid$ok) || !any(grepl("MODEL|PAYLOAD",
      invalid$diagnostics$code)))
    stop(paste("damaged model diagnostic:", paste(
      invalid$diagnostics$code, collapse = ", ")))
  TRUE
}
