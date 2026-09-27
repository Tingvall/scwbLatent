args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("Expected accepted old-core source and destination bundle paths.")
source <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
destination <- args[[2L]]
revision <- trimws(system2("git", c("-C", shQuote(source),
  "rev-parse", "HEAD"), stdout = TRUE))
if (!identical(revision,
    "b564d10cdb5aee757a0c8401ab15f69594d12ea0"))
  stop("Historical scVI fixture requires the exact accepted old core.")
pkgload::load_all(source, quiet = TRUE)
b <- scworkbench::wb_init_bundle(destination)
cells <- c("cell_2", "cell_10", "cell_1")
scworkbench::wb_define_cells(b, cells)
scworkbench::wb_attach_space(b, "historical_scvi_space",
  data.frame(cell_id = cells, z1 = c(1, 2, 3), z2 = c(3, 2, 1)),
  source = "historical inspection fixture")
model_dir <- file.path(destination, "artifact_models", "historical_scvi")
dir.create(model_dir, recursive = TRUE)
writeLines("retained historical scVI model fixture",
  file.path(model_dir, "model.txt"))
model <- scworkbench::wb_prim_model_ref("historical_scvi_space",
  "historical_scvi_space", engine = "scvi")
model$model_path <- model_dir
model$model_type <- "scvi_train"
run <- scworkbench::wb_new_run("scvi_train", inputs = list(
  cell_set = "all", parents = character()), outputs = list(model = model))
session <- scworkbench::wb_session_add_run(scworkbench::wb_session_new(),
  run)
scworkbench::wb_session_promote_current(b, session, outputs = "model",
  created_by = "accepted-old-core-scvi-fixture")
scworkbench::wb_finalize_recipient_bundle(b)
scworkbench::wb_close_bundle(b)
