.latent_fail <- function(message) stop(paste0("scwbLatent: ", message),
  call. = FALSE)

.latent_stage_matrix <- function(reader, capability, chunk_size, root) {
  cells <- as.character(capability$cells)
  features <- as.character(capability$features)
  chunk_size <- as.integer(chunk_size)[1]
  if (!length(cells) || !length(features) || is.na(chunk_size) ||
      chunk_size < 1L || anyDuplicated(cells) || anyDuplicated(features))
    .latent_fail("invalid authenticated matrix axes or chunk size")
  dir.create(file.path(root, "chunks"), recursive = TRUE,
    showWarnings = FALSE)
  writeLines(features, file.path(root, "features.txt"), useBytes = TRUE)
  chunks <- unname(split(cells, ceiling(seq_along(cells) / chunk_size)))
  entries <- vector("list", length(chunks))
  digests <- list()
  maximum_bytes <- 0
  maximum_elements <- 0
  for (i in seq_along(chunks)) {
    selected <- as.character(chunks[[i]])
    piece <- reader(features = features, cells = selected)
    if (!identical(as.character(piece$feature_id), features) ||
        !identical(as.character(piece$cell_id), selected) ||
        !identical(as.integer(dim(piece$values)),
          as.integer(c(length(features), length(selected)))))
      .latent_fail("typed matrix chunk differs from authenticated axes")
    maximum_bytes <- max(maximum_bytes,
      as.numeric(utils::object.size(piece$values)))
    maximum_elements <- max(maximum_elements,
      as.double(length(features)) * as.double(length(selected)))
    prefix <- sprintf("chunks/%06d", i)
    matrix_relative <- paste0(prefix, ".mtx")
    cells_relative <- paste0(prefix, ".cells.txt")
    Matrix::writeMM(Matrix::Matrix(piece$values, sparse = TRUE),
      file.path(root, matrix_relative))
    writeLines(selected, file.path(root, cells_relative), useBytes = TRUE)
    entries[[i]] <- list(matrix = matrix_relative, cells = cells_relative)
    digests[[matrix_relative]] <- digest::digest(
      file = file.path(root, matrix_relative), algo = "sha256")
    digests[[cells_relative]] <- digest::digest(
      file = file.path(root, cells_relative), algo = "sha256")
  }
  digests[["features.txt"]] <- digest::digest(
    file = file.path(root, "features.txt"), algo = "sha256")
  if (maximum_elements > as.double(length(features)) * chunk_size)
    .latent_fail("matrix staging exceeded the declared chunk bound")
  list(cells = cells, features = features, chunks = entries,
    digests = digests, max_chunk_bytes = maximum_bytes,
    max_chunk_elements = maximum_elements,
    total_elements = as.double(length(cells)) * as.double(length(features)))
}

.latent_stage_model_read <- function(payload, root) {
  if (!is.list(payload) || !is.function(payload$runtime_private_path) ||
      !length(payload$entries))
    .latent_fail("G1 model read did not provide an authenticated private tree")
  target <- file.path(root, "model_input")
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  for (entry in payload$entries) {
    relative <- as.character(entry$path)[1]
    if (!nzchar(relative) || grepl("\\\\", relative) ||
        startsWith(relative, "/") ||
        any(strsplit(relative, "/", fixed = TRUE)[[1]] %in%
          c("", ".", "..")))
      .latent_fail("G1 model tree contains an unsafe relative path")
    source <- payload$runtime_private_path(relative)
    destination <- file.path(target, relative)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(file.copy(source, destination, overwrite = FALSE)))
      .latent_fail("could not copy authenticated G1 model entry privately")
    if (!identical(digest::digest(file = destination, algo = "sha256"),
        as.character(entry$sha256)[1]))
      .latent_fail("private model copy changed authenticated bytes")
  }
  invisible("model_input")
}

.latent_private_compute <- function(request, stage, staged) {
  python <- .latent_runtime$python
  if (!nzchar(as.character(python %||% "")[1]) ||
      !file.exists(python)) .latent_fail(
        "CPU PeakVI interpreter disappeared before execution")
  script <- system.file("python", "peakvi_private.py", package = "scwbLatent")
  if (!nzchar(script) || !file.exists(script))
    .latent_fail("private PeakVI compute helper is absent")
  jsonlite::write_json(request, file.path(stage, "request.json"),
    auto_unbox = TRUE, null = "null", digits = NA)
  before <- staged$digests
  status <- suppressWarnings(system2(python,
    c(shQuote(script), "--stage-root", shQuote(stage)),
    stdout = file.path(stage, "stdout.log"),
    stderr = file.path(stage, "stderr.log")))
  if (!identical(as.integer(status), 0L))
    .latent_fail(paste0("private PeakVI compute failed: ", paste(
      readLines(file.path(stage, "stderr.log"), warn = FALSE),
      collapse = "\n")))
  for (relative in names(before))
    if (!identical(digest::digest(file = file.path(stage, relative),
        algo = "sha256"), before[[relative]]))
      .latent_fail("private input bytes changed during compute")
  response_path <- file.path(stage, "response.json")
  latent_path <- file.path(stage, "latent.tsv")
  if (!file.exists(response_path) || !file.exists(latent_path) ||
      nzchar(Sys.readlink(response_path)) || nzchar(Sys.readlink(latent_path)))
    .latent_fail("private compute response is missing or linked")
  response <- jsonlite::read_json(response_path, simplifyVector = TRUE)
  if (!identical(as.character(response$status)[1], "ok"))
    .latent_fail("private compute did not return success")
  latent_sha256 <- digest::digest(file = latent_path, algo = "sha256")
  if (!identical(as.character(response$latent_sha256)[1], latent_sha256))
    .latent_fail("private latent output checksum differs from response")
  latent <- utils::read.delim(latent_path, check.names = FALSE,
    stringsAsFactors = FALSE)
  expected <- c("cell_id", paste0("latent_", seq_len(as.integer(
    request$n_latent))))
  if (!identical(names(latent), expected) ||
      !identical(as.character(latent$cell_id), staged$cells) ||
      !identical(as.integer(response$latent_rows),
        as.integer(length(staged$cells))) ||
      !identical(as.integer(response$latent_columns),
        as.integer(request$n_latent)))
    .latent_fail("private latent output differs from authenticated axes")
  values <- as.matrix(latent[, -1L, drop = FALSE])
  storage.mode(values) <- "double"
  if (any(!is.finite(values)))
    .latent_fail("private latent output contains non-finite values")
  dimnames(values) <- list(staged$cells, expected[-1L])
  list(values = values, response = response,
    latent_sha256 = latent_sha256)
}
