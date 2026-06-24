suppressPackageStartupMessages({
  user_lib <- "/home/dt2024/dt2024020307/R/x86_64-pc-linux-gnu-library/4.2"
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
  library(Seurat)
  library(CellChat)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 50 * 1024^3)

input_dir <- Sys.getenv(
  "SCPCP_INPUT_DIR",
  "/home/dt2024/share/project/scpca_cellchat/scpca/SCPCP000004/rds"
)
result_dir <- Sys.getenv(
  "SCPCP_RESULT_DIR",
  "/home/dt2024/dt2024020307/cellchat/results/SCPCP000004_official_cellchat_R"
)

nboot <- as.integer(Sys.getenv("SCPCP_NBOOT", "100"))
seed_use <- as.integer(Sys.getenv("SCPCP_SEED", "1"))
min_cells <- as.integer(Sys.getenv("SCPCP_MIN_CELLS", "10"))
limit <- as.integer(Sys.getenv("SCPCP_LIMIT", "0"))
overwrite <- tolower(Sys.getenv("SCPCP_OVERWRITE", "FALSE")) %in% c("true", "1", "yes", "y")
workers <- as.integer(Sys.getenv("SCPCP_WORKERS", "12"))

if (requireNamespace("future", quietly = TRUE)) {
  if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::multisession, workers = workers)
  }
}

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
out_rds_dir <- file.path(result_dir, "cellchat_obj_rds")
out_comm_dir <- file.path(result_dir, "communication_csv")
out_group_dir <- file.path(result_dir, "group_counts")
dir.create(out_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_comm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_group_dir, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(result_dir, "run_summary.csv")
timing_csv <- file.path(result_dir, "all_step_timings.csv")
step_total_csv <- file.path(result_dir, "step_total_by_step.csv")

run_started <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
message("SCPCP000004 official CellChat R run started: ", run_started)
message("Input dir: ", input_dir)
message("Result dir: ", result_dir)
message("nboot=", nboot, " seed=", seed_use, " min_cells=", min_cells, " workers=", workers)

select_group_col <- function(meta) {
  preferred <- c(
    "openscpca_celltype_annotation",
    "consensus_celltype_annotation",
    "cellassign_celltype_annotation",
    "singler_celltype_annotation",
    "scimilarity_celltype_annotation",
    "cell_type",
    "celltype",
    "annotation",
    "seurat_clusters",
    "cluster"
  )
  for (nm in preferred) {
    if (nm %in% colnames(meta)) return(nm)
  }
  stop("No supported annotation/group column found in Seurat metadata.", call. = FALSE)
}

invalid_label <- function(x) {
  y <- tolower(trimws(as.character(x)))
  is.na(x) | y == "" | y %in% c("openscpca-excluded", "openscpca excluded", "excluded")
}

write_table <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE)
}

append_step <- function(store, sample_id, step, start_time, end_time) {
  store[[length(store) + 1L]] <- data.frame(
    sample_id = sample_id,
    step = step,
    start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
    end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
    elapsed_sec = as.numeric(difftime(end_time, start_time, units = "secs")),
    stringsAsFactors = FALSE
  )
  store
}

time_step <- function(sample_id, step, expr, timing_store) {
  message("[", sample_id, "] ", step, " ...")
  start_time <- Sys.time()
  value <- force(expr)
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  message("[", sample_id, "] ", step, " done in ", sprintf("%.2f", elapsed), " sec")
  list(value = value, timing_store = append_step(timing_store, sample_id, step, start_time, end_time))
}

timed <- function(sample_id, step, timing_store, code) {
  time_step(sample_id, step, eval.parent(substitute(code)), timing_store)
}

summary_rows <- list()
timing_rows <- list()

flush_outputs <- function() {
  summary_df <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
  timing_df <- if (length(timing_rows)) do.call(rbind, timing_rows) else data.frame()
  write_table(summary_df, summary_csv)
  write_table(timing_df, timing_csv)
  if (nrow(timing_df)) {
    step_total <- aggregate(elapsed_sec ~ step, timing_df, sum)
    step_total <- step_total[order(step_total$elapsed_sec, decreasing = TRUE), , drop = FALSE]
    write_table(step_total, step_total_csv)
  }
}

files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("No .rds files found in input_dir: ", input_dir, call. = FALSE)
if (is.finite(limit) && limit > 0L) files <- head(files, limit)

for (file in files) {
  sample_id <- sub("_processed_seurat\\.rds$", "", basename(file))
  sample_start <- Sys.time()
  out_rds <- file.path(out_rds_dir, paste0(sample_id, "_cellchat_official_R.rds"))
  status <- "success"
  message_text <- ""
  group_col <- NA_character_
  assay_use <- NA_character_
  n_cells_input <- NA_integer_
  n_cells_keep <- NA_integer_
  n_groups <- NA_integer_
  n_lr <- NA_integer_
  n_pathways <- NA_integer_

  if (file.exists(out_rds) && !overwrite) {
    message("[", sample_id, "] existing result found, skipping: ", out_rds)
    sample_end <- Sys.time()
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      sample_id = sample_id,
      input_file = file,
      output_rds = out_rds,
      status = "skipped_existing",
      message = "Existing result found; set SCPCP_OVERWRITE=TRUE to rerun.",
      group_col = NA_character_,
      assay = NA_character_,
      n_cells_input = NA_integer_,
      n_cells_used = NA_integer_,
      n_groups = NA_integer_,
      n_lr = NA_integer_,
      n_pathways = NA_integer_,
      start_time = format(sample_start, "%Y-%m-%d %H:%M:%S"),
      end_time = format(sample_end, "%Y-%m-%d %H:%M:%S"),
      elapsed_sec = as.numeric(difftime(sample_end, sample_start, units = "secs")),
      stringsAsFactors = FALSE
    )
    flush_outputs()
    next
  }

  message("============================================================")
  message("[", sample_id, "] input: ", file)

  tryCatch({
    res <- timed(sample_id, "readRDS", timing_rows, {
      readRDS(file)
    })
    obj <- res$value
    timing_rows <- res$timing_store

    if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object.", call. = FALSE)
    n_cells_input <- ncol(obj)
    group_col <- select_group_col(obj@meta.data)
    assay_use <- if ("RNA" %in% names(obj@assays)) "RNA" else Seurat::DefaultAssay(obj)
    Seurat::DefaultAssay(obj) <- assay_use

    res <- timed(sample_id, "filter_cells_and_groups", timing_rows, {
      meta <- obj@meta.data
      keep <- rep(TRUE, nrow(meta))
      if ("scpca_filter" %in% colnames(meta)) {
        keep <- keep & !is.na(meta$scpca_filter) & as.character(meta$scpca_filter) == "Keep"
      }
      labels <- as.character(meta[[group_col]])
      keep <- keep & !invalid_label(labels)
      cells_keep <- rownames(meta)[keep]
      obj2 <- subset(obj, cells = cells_keep)
      labels2 <- trimws(as.character(obj2@meta.data[[group_col]]))
      obj2@meta.data[[group_col]] <- labels2
      group_counts <- as.data.frame(table(obj2@meta.data[[group_col]]), stringsAsFactors = FALSE)
      colnames(group_counts) <- c("group", "n_cells")
      valid_groups <- group_counts$group[group_counts$n_cells >= min_cells]
      obj2 <- subset(obj2, cells = rownames(obj2@meta.data)[obj2@meta.data[[group_col]] %in% valid_groups])
      obj2@meta.data[[group_col]] <- factor(as.character(obj2@meta.data[[group_col]]), levels = valid_groups)
      group_counts <- as.data.frame(table(obj2@meta.data[[group_col]]), stringsAsFactors = FALSE)
      colnames(group_counts) <- c("group", "n_cells")
      group_counts <- group_counts[order(group_counts$n_cells, decreasing = TRUE), , drop = FALSE]
      group_counts$sample_id <- sample_id
      group_counts <- group_counts[, c("sample_id", "group", "n_cells")]
      write_table(group_counts, file.path(out_group_dir, paste0(sample_id, "_group_counts.csv")))
      list(obj = obj2, group_counts = group_counts)
    })
    obj <- res$value$obj
    group_counts <- res$value$group_counts
    timing_rows <- res$timing_store
    n_cells_keep <- ncol(obj)
    n_groups <- nrow(group_counts)

    if (n_cells_keep < 2L) stop("Fewer than 2 cells remain after filtering.", call. = FALSE)
    if (n_groups < 2L) stop("Fewer than 2 annotation groups remain after filtering.", call. = FALSE)

    res <- timed(sample_id, "check_or_normalize_data_slot", timing_rows, {
      data_slot <- Seurat::GetAssayData(obj, assay = assay_use, slot = "data")
      has_positive <- if (inherits(data_slot, "sparseMatrix")) {
        length(data_slot@x) > 0L && any(data_slot@x > 0)
      } else {
        length(data_slot) > 0L && any(data_slot > 0)
      }
      if (!has_positive) {
        obj <- Seurat::NormalizeData(obj, assay = assay_use, verbose = FALSE)
      }
      obj
    })
    obj <- res$value
    timing_rows <- res$timing_store

    res <- timed(sample_id, "createCellChat", timing_rows, {
      CellChat::createCellChat(object = obj, group.by = group_col, assay = assay_use)
    })
    cellchat <- res$value
    timing_rows <- res$timing_store
    cellchat@DB <- CellChat::CellChatDB.human

    res <- timed(sample_id, "subsetData", timing_rows, {
      CellChat::subsetData(cellchat)
    })
    cellchat <- res$value
    timing_rows <- res$timing_store

    res <- timed(sample_id, "identifyOverExpressedGenes", timing_rows, {
      CellChat::identifyOverExpressedGenes(cellchat)
    })
    cellchat <- res$value
    timing_rows <- res$timing_store

    res <- timed(sample_id, "identifyOverExpressedInteractions", timing_rows, {
      CellChat::identifyOverExpressedInteractions(cellchat)
    })
    cellchat <- res$value
    timing_rows <- res$timing_store

    if (is.null(cellchat@LR$LRsig) || nrow(cellchat@LR$LRsig) == 0L) {
      stop("No significant LR pairs after identifyOverExpressedInteractions().", call. = FALSE)
    }

    res <- timed(sample_id, "computeCommunProb_official_R", timing_rows, {
      CellChat::computeCommunProb(
        cellchat,
        type = "triMean",
        nboot = nboot,
        seed.use = seed_use,
        population.size = FALSE
      )
    })
    cellchat <- res$value
    timing_rows <- res$timing_store

    res <- timed(sample_id, "filterCommunication", timing_rows, {
      CellChat::filterCommunication(cellchat, min.cells = min_cells)
    })
    cellchat <- res$value
    timing_rows <- res$timing_store

    res <- timed(sample_id, "computeCommunProbPathway_official_R", timing_rows, {
      CellChat::computeCommunProbPathway(cellchat)
    })
    cellchat <- res$value
    timing_rows <- res$timing_store

    res <- timed(sample_id, "aggregateNet_official_R", timing_rows, {
      CellChat::aggregateNet(cellchat)
    })
    cellchat <- res$value
    timing_rows <- res$timing_store

    n_lr <- if (!is.null(cellchat@net$prob)) dim(cellchat@net$prob)[3] else NA_integer_
    n_pathways <- if (!is.null(cellchat@netP$pathways)) length(cellchat@netP$pathways) else NA_integer_

    res <- timed(sample_id, "write_outputs", timing_rows, {
      saveRDS(cellchat, out_rds, compress = FALSE)
      comm <- tryCatch(CellChat::subsetCommunication(cellchat), error = function(e) NULL)
      if (!is.null(comm)) {
        write_table(comm, file.path(out_comm_dir, paste0(sample_id, "_communication.csv")))
      }
      TRUE
    })
    timing_rows <- res$timing_store

    rm(obj, cellchat)
    gc()
  }, error = function(e) {
    status <<- "failed"
    message_text <<- conditionMessage(e)
    message("[", sample_id, "] FAILED: ", message_text)
  })

  sample_end <- Sys.time()
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    sample_id = sample_id,
    input_file = file,
    output_rds = out_rds,
    status = status,
    message = message_text,
    group_col = group_col,
    assay = assay_use,
    n_cells_input = n_cells_input,
    n_cells_used = n_cells_keep,
    n_groups = n_groups,
    n_lr = n_lr,
    n_pathways = n_pathways,
    start_time = format(sample_start, "%Y-%m-%d %H:%M:%S"),
    end_time = format(sample_end, "%Y-%m-%d %H:%M:%S"),
    elapsed_sec = as.numeric(difftime(sample_end, sample_start, units = "secs")),
    stringsAsFactors = FALSE
  )
  flush_outputs()
}

run_ended <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
message("SCPCP000004 official CellChat R run finished: ", run_ended)
message("Summary: ", summary_csv)
message("Timings: ", timing_csv)
