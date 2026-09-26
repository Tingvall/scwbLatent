args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Expected old-core source and output RDS paths.")
run_oracle <- function(args) {
source <- normalizePath(args[[1L]], mustWork = TRUE)
output <- args[[2L]]
pkgload::load_all(source, quiet = TRUE)
root <- tempfile("peakvi-old-core-oracle-")
b <- scworkbench::wb_init_bundle(root)
on.exit({try(scworkbench::wb_close_bundle(b), silent = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
reference <- paste0("ref_", seq_len(24L))
query <- paste0("query_", seq_len(8L))
features <- paste0("peak_", seq_len(32L))
scworkbench::wb_define_cells(b, c(reference, query))
values <- matrix((seq_len(length(features) * length(reference)) %% 5L) + 1L,
  nrow = length(features), dimnames = list(features, reference))
scworkbench::wb_attach_matrix(b, "ATAC_reference",
  Matrix::Matrix(values, sparse = TRUE), feature_type = "peak",
  cell_ids = reference)
query_values <- matrix((seq_len(length(features) * length(query)) %% 7L) + 1L,
  nrow = length(features), dimnames = list(features, query))
scworkbench::wb_attach_matrix(b, "ATAC_query",
  Matrix::Matrix(query_values, sparse = TRUE), feature_type = "peak",
  cell_ids = query)
fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(), b,
  "peakvi_fit", inputs = list(matrix = "ATAC_reference", n_latent = 2L,
    n_layers = 1L, n_hidden = 16L, max_epochs = 1L,
    cell_set = "all", min_cells = 1L,
    space_id_out = "legacy_peakvi"))
space <- scworkbench::wb_space(scworkbench::wb_session_view(fit$session),
  "legacy_peakvi")
coords <- as.matrix(space[-1L])
saveRDS(list(cells = rownames(coords), coordinates = coords,
  core_revision = "d624cfd0cc1b84440681967fc30a214004e2c591"), output)
scworkbench::wb_session_discard(fit$session)
}
run_oracle(args)
