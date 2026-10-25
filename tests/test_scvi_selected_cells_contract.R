# Real one-epoch CPU acceptance. Source the existing pack helper contract first.
test_scvi_selected_cells_contract <- function() {
  if (!"scwblatent::scvi_train" %in% scworkbench::wb_list_ops()) scworkbench::wb_operation_pack_load("scwbLatent")
  root <- tempfile("scvi-selected-"); b <- scworkbench::wb_init_bundle(root)
  on.exit({try(scworkbench::wb_close_bundle(b),silent=TRUE);unlink(root,recursive=TRUE)},add=TRUE)
  cells <- paste0("cell_",c(10L,2L,1L,setdiff(1:24,c(10L,2L,1L))))
  features <- paste0("gene_",c(3L,1L,setdiff(1:32,c(3L,1L))))
  md <- data.frame(cell_id=cells,batch=rep(c("B","A"),12L),
    condition=rep(c("ctrl","stim","stim"),8L),quality=seq_along(cells)/25,
    include=rep(c("yes","yes","no"),8L),stringsAsFactors=FALSE)
  scworkbench::wb_define_cells(b,cells,metadata=md)
  x <- matrix(as.integer((seq_len(32L*24L)*7L)%%11L),32L,dimnames=list(features,cells))
  scworkbench::wb_attach_matrix(b,"RNA_counts",Matrix::Matrix(x,sparse=TRUE),layer="counts",feature_type="gene")
  selected <- as.character(scworkbench::wb_sel_meta(scworkbench::wb_sel_all(b),"include",values="yes"))
  stopifnot(length(selected)==16L,length(unique(md$batch[match(selected,cells)]))==2L)
  scworkbench::wb_cell_set_replace(b,"set:selected",selected)
  settings <- list(matrix="RNA_counts",batch_key="batch",categorical_covariates="condition",
    continuous_covariates="quality",n_latent=2L,n_layers=1L,n_hidden=16L,
    max_epochs=1L,early_stopping=FALSE,use_gpu=FALSE,chunk_size=4L,min_cells=1L)
  baseline <- .latent_scvi_source_only(b)
  for (bad in list(c(selected,selected[[1]]),c(selected[-1L],"foreign"))) {
    failed <- tryCatch(scworkbench::wb_session_run_op(scworkbench::wb_session_new(),b,
      "scwblatent::scvi_train",inputs=c(settings,list(cell_ids=bad))),error=identity)
    stopifnot(inherits(failed,"error"),identical(.latent_scvi_source_only(b),baseline))
  }
  records <- list()
  for (route in c("named_set","metadata_selector","all")) {
    ids <- if(route=="all") cells else selected
    choice <- if(route=="named_set") list(cell_set="set:selected") else if(route=="metadata_selector")
      list(cell_ids=as.character(scworkbench::wb_sel_meta(scworkbench::wb_sel_all(b),"include",values="yes"))) else list(cell_set="all")
    start <- proc.time()[["elapsed"]]
    fit <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),b,"scwblatent::scvi_train",inputs=c(settings,choice))
    space_id <- fit$run$outputs$space$space_id
    latent <- scworkbench::wb_space(fit$session$staging_bundle,space_id)
    if(is.data.frame(latent)) {
      stopifnot(identical(as.character(latent$cell_id),ids))
      latent <- as.matrix(latent[,setdiff(names(latent),"cell_id"),drop=FALSE])
    } else stopifnot(identical(rownames(latent),ids))
    stopifnot(identical(as.integer(dim(latent)),c(length(ids),2L)),all(is.finite(latent)))
    promoted <- scworkbench::wb_session_promote_current(b,fit$session,outputs=c("model","fit_summary"))
    stopifnot(all(promoted$promoted$status=="ok"))
    transition <- as.character(fit$run$operation_pack_transition$transition_id)
    reads <- scworkbench::wb_operation_pack_transition_registry_reads(b,transition)
    snapshot <- reads$cell_metadata$snapshot
    stopifnot(identical(as.character(unlist(snapshot$cell_id,use.names=FALSE)),ids),
      snapshot$cell_count==length(ids),
      identical(vapply(snapshot$fields,`[[`,"","name"),c("batch","condition","quality")))
    index <- match(ids,cells)
    for(j in seq_along(snapshot$fields)) stopifnot(isTRUE(all.equal(
      unlist(snapshot$fields[[j]]$values,use.names=FALSE),md[[c("batch","condition","quality")[[j]]]][index])))
    # Inspect the authenticated retained request/receipt; no training is mocked.
    record <- scworkbench:::.wb_operation_pack_v2_record_validate(file.path(root,"operation_pack_v2",transition),load_payloads=TRUE)
    stopifnot(record$ok)
    descriptor <- record$request$inputs[[1]]$descriptor
    grant <- record$request$typed_slice_grants[[1]]
    axis <- jsonlite::read_json(file.path(root,"operation_pack_v2",transition,grant$cells_axis$relative_path),simplifyVector=TRUE)
    stopifnot(identical(as.character(axis$ids),ids), identical(as.integer(descriptor$dimensions),c(32L,length(ids))))
    stopifnot(length(record$receipt$typed_slice_usage)==1L,
      record$receipt$typed_slice_usage[[1]]$coverage_status=="exact",
      sum(vapply(record$receipt$typed_slice_usage[[1]]$calls,function(call)as.integer(call$cell_count),integer(1)))==length(ids),
      record$receipt$typed_slice_usage[[1]]$call_count==ceiling(length(ids)/4))
    model_id <- promoted$promoted$artifact_id[promoted$promoted$output=="model"][[1]]
    summary_id <- promoted$promoted$artifact_id[promoted$promoted$output=="fit_summary"][[1]]
    summary <- scworkbench::wb_art_load(b,summary_id,expected_kind="obs_table")$table
    stopifnot(summary$n_cells[[1]]==length(ids),summary$observed_training_device[[1]]=="cpu",
      summary$completed_epochs[[1]]==1L,summary$training_batches_observed[[1]]>=1L,
      summary$max_chunk_elements[[1]]<=32L*4L)
    records[[route]] <- list(space_id=space_id,model_id=model_id,summary_id=summary_id,cells=ids)
    cat("REAL_CPU_FIT",route,"cells",length(ids),"batches",length(unique(md$batch[index])),
      "chunk_calls",record$receipt$typed_slice_usage[[1]]$call_count,"seconds",proc.time()[["elapsed"]]-start,"\n")
  }
  scworkbench::wb_finalize_recipient_bundle(b);scworkbench::wb_close_bundle(b)
  expected <- file.path(dirname(root),paste0(basename(root),"-expected.rds"));saveRDS(records,expected)
  on.exit(unlink(expected),add=TRUE)
  code <- paste0(".libPaths(c(",deparse(dirname(find.package("scworkbench"))),",.libPaths()));library(scworkbench);",
    "b<-wb_open_bundle(",deparse(root),",read_only=TRUE);x<-readRDS(",deparse(expected),");",
    "for(r in x){v<-wb_space(b,r$space_id);if(is.data.frame(v)){stopifnot(identical(v$cell_id,r$cells));v<-as.matrix(v[,-1,drop=FALSE])}else stopifnot(identical(rownames(v),r$cells));stopifnot(all(is.finite(v)),nrow(v)==length(r$cells),wb_art_load(b,r$model_id,expected_kind='model_ref')$engine=='scvi',wb_art_load(b,r$summary_id,expected_kind='obs_table')$table$n_cells[[1]]==length(r$cells))};wb_close_bundle(b);cat('FRESH_PROCESS_REOPEN_PASS\\n')")
  child <- processx::run(file.path(R.home("bin"),"Rscript"),c("-e",code),error_on_status=FALSE,timeout=120000)
  if(child$status!=0L) stop(child$stderr,child$stdout)
  cat(child$stdout)
  stopifnot(scworkbench::wb_bundle_validation_report(root,"recipient")$ok)
  TRUE
}
