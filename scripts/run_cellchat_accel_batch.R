#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(CellChatAccelRcpp)
})

parse_args <- function(args) {
  out <- list(
    input_dir = NULL,
    output_dir = "cellchat_accel_results",
    pattern = "\\.rds$",
    recursive = FALSE,
    group_col = "openscpca_celltype_annotation",
    assay = "RNA",
    species = "human",
    min_cells = 10L,
    nboot = 100L,
    seed = 1L,
    exclude_label = "openscpca-excluded",
    skip_existing = TRUE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unknown argument: ", key, call. = FALSE)
    key <- sub("^--", "", key)
    if (!key %in% names(out)) stop("Unknown option: --", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for --", key, call. = FALSE)
    val <- args[[i + 1L]]
    out[[key]] <- val
    i <- i + 2L
  }
  out$recursive <- tolower(as.character(out$recursive)) %in% c("1", "true", "yes", "y")
  out$skip_existing <- tolower(as.character(out$skip_existing)) %in% c("1", "true", "yes", "y")
  out$min_cells <- as.integer(out$min_cells)
  out$nboot <- as.integer(out$nboot)
  out$seed <- as.integer(out$seed)
  out
}

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/run_cellchat_accel_batch.R \\\n",
    "    --input_dir /path/to/rds_dir \\\n",
    "    --output_dir /path/to/results \\\n",
    "    --group_col openscpca_celltype_annotation \\\n",
    "    --pattern '\\\\.rds$' \\\n",
    "    --nboot 100 \\\n",
    "    --min_cells 10\n\n",
    "Required:\n",
    "  --input_dir   Directory containing Seurat .rds files.\n\n",
    "Common options:\n",
    "  --output_dir      Output directory. Default: cellchat_accel_results\n",
    "  --group_col       Seurat metadata column used as cell group labels.\n",
    "  --pattern         Regex for input files. Default: \\\\.rds$\n",
    "  --recursive       true/false. Default: false\n",
    "  --assay           Seurat assay. Default: RNA\n",
    "  --species         human or mouse. Default: human\n",
    "  --min_cells       Minimum cells per group. Default: 10\n",
    "  --nboot           Number of bootstrap permutations. Default: 100\n",
    "  --seed            Random seed. Default: 1\n",
    "  --exclude_label   Metadata label to exclude. Default: openscpca-excluded\n",
    "  --skip_existing   true/false. Default: true\n",
    sep = ""
  )
}

record <- function(rows, sample_id, step, expr) {
  gc()
  t0 <- proc.time()
  w0 <- Sys.time()
  status <- "ok"
  error <- NA_character_
  value <- tryCatch(force(expr), error = function(e) {
    status <<- "error"
    error <<- conditionMessage(e)
    NULL
  })
  elapsed <- proc.time() - t0
  row <- data.frame(
    sample_id = sample_id,
    step = step,
    status = status,
    elapsed_sec = unname(elapsed[["elapsed"]]),
    user_sec = unname(elapsed[["user.self"]]),
    system_sec = unname(elapsed[["sys.self"]]),
    wall_sec = as.numeric(difftime(Sys.time(), w0, units = "secs")),
    error = error,
    stringsAsFactors = FALSE
  )
  cat(format(Sys.time(), "%F %T"), sample_id, step, status,
      round(row$elapsed_sec, 3), "sec\n")
  list(value = value, rows = c(rows, list(row)), ok = identical(status, "ok"))
}

filter_seurat_object <- function(obj, group_col, min_cells, exclude_label) {
  if (!group_col %in% colnames(obj@meta.data)) {
    stop("Missing group column in Seurat metadata: ", group_col, call. = FALSE)
  }
  group <- as.character(obj@meta.data[[group_col]])
  keep <- !is.na(group) & nzchar(group)
  if (nzchar(exclude_label)) keep <- keep & group != exclude_label
  obj <- subset(obj, cells = colnames(obj)[keep])

  group <- as.character(obj@meta.data[[group_col]])
  tab <- table(group)
  valid <- names(tab)[tab >= min_cells]
  if (!length(valid)) {
    stop("No cell groups have at least ", min_cells, " cells after filtering.", call. = FALSE)
  }
  subset(obj, cells = colnames(obj)[group %in% valid])
}

select_cellchat_db <- function(species) {
  species <- tolower(species)
  if (identical(species, "human")) return(CellChatDB.human)
  if (identical(species, "mouse")) return(CellChatDB.mouse)
  stop("--species must be 'human' or 'mouse'.", call. = FALSE)
}

process_one <- function(path, args, db) {
  sample_id <- sub("\\.rds$", "", basename(path), ignore.case = TRUE)
  out_rds <- file.path(args$output_dir, paste0(sample_id, "_cellchat_accel.rds"))
  out_lr <- file.path(args$output_dir, paste0(sample_id, "_LR_detail.csv"))
  rows <- list()

  if (args$skip_existing && file.exists(out_rds)) {
    cat("Skipping existing output:", out_rds, "\n")
    return(data.frame(
      sample_id = sample_id,
      input_file = path,
      output_file = out_rds,
      status = "skipped_existing",
      stringsAsFactors = FALSE
    ))
  }

  rr <- record(rows, sample_id, "01_readRDS", readRDS(path))
  rows <- rr$rows
  if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
  obj <- rr$value

  if (inherits(obj, "CellChat")) {
    cc <- obj
    cc@DB <- db
  } else {
    rr <- record(rows, sample_id, "02_filter_seurat", {
      Seurat::DefaultAssay(obj) <- args$assay
      filter_seurat_object(obj, args$group_col, args$min_cells, args$exclude_label)
    })
    rows <- rr$rows
    if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
    obj <- rr$value

    rr <- record(rows, sample_id, "03_createCellChat",
                 CellChat::createCellChat(object = obj, group.by = args$group_col, assay = args$assay))
    rows <- rr$rows
    if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
    cc <- rr$value
    rm(obj)
    gc()
    cc@DB <- db
  }

  rr <- record(rows, sample_id, "04_subsetData", CellChat::subsetData(cc))
  rows <- rr$rows
  if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
  cc <- rr$value

  rr <- record(rows, sample_id, "05_identifyOverExpressedGenes",
               CellChat::identifyOverExpressedGenes(cc))
  rows <- rr$rows
  if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
  cc <- rr$value

  rr <- record(rows, sample_id, "06_identifyOverExpressedInteractions",
               CellChat::identifyOverExpressedInteractions(cc))
  rows <- rr$rows
  if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
  cc <- rr$value

  rr <- record(rows, sample_id, "07_computeCommunProbAccelRcpp",
               CellChatAccelRcpp::computeCommunProbAccelRcpp(cc, nboot = args$nboot, seed.use = args$seed))
  rows <- rr$rows
  if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
  cc <- rr$value

  rr <- record(rows, sample_id, "08_filterCommunication",
               CellChat::filterCommunication(cc, min.cells = args$min_cells))
  rows <- rr$rows
  if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
  cc <- rr$value

  rr <- record(rows, sample_id, "09_computeCommunProbPathwayAccelRcpp",
               CellChatAccelRcpp::computeCommunProbPathwayAccelRcpp(cc))
  rows <- rr$rows
  if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
  cc <- rr$value

  rr <- record(rows, sample_id, "10_aggregateNetAccelRcpp",
               CellChatAccelRcpp::aggregateNetAccelRcpp(cc))
  rows <- rr$rows
  if (!rr$ok) return(write_sample_failure(args$output_dir, sample_id, path, out_rds, rows))
  cc <- rr$value

  lr <- tryCatch(CellChat::subsetCommunication(cc), error = function(e) data.frame(error = conditionMessage(e)))
  utils::write.csv(lr, out_lr, row.names = FALSE)
  saveRDS(cc, out_rds, compress = FALSE)
  write_sample_timings(args$output_dir, sample_id, rows)

  data.frame(
    sample_id = sample_id,
    input_file = path,
    output_file = out_rds,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

write_sample_timings <- function(output_dir, sample_id, rows) {
  timing <- do.call(rbind, rows)
  utils::write.csv(timing, file.path(output_dir, paste0(sample_id, "_step_timings.csv")), row.names = FALSE)
}

write_sample_failure <- function(output_dir, sample_id, input_file, output_file, rows) {
  write_sample_timings(output_dir, sample_id, rows)
  last <- tail(do.call(rbind, rows), 1)
  data.frame(
    sample_id = sample_id,
    input_file = input_file,
    output_file = output_file,
    status = "error",
    error_step = last$step,
    error = last$error,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$input_dir) || identical(args$input_dir, "")) {
    usage()
    stop("--input_dir is required.", call. = FALSE)
  }
  if (!dir.exists(args$input_dir)) stop("Input directory does not exist: ", args$input_dir, call. = FALSE)
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(args$input_dir, pattern = args$pattern, recursive = args$recursive,
                      full.names = TRUE)
  files <- sort(files)
  if (!length(files)) stop("No input files matched pattern: ", args$pattern, call. = FALSE)

  db <- select_cellchat_db(args$species)
  cat("Input files:", length(files), "\n")
  cat("Output dir:", normalizePath(args$output_dir), "\n")
  cat("Group column:", args$group_col, "\n")
  cat("nboot:", args$nboot, "\n")

  summaries <- vector("list", length(files))
  for (i in seq_along(files)) {
    cat("\n=== [", i, "/", length(files), "] ", basename(files[[i]]), " ===\n", sep = "")
    summaries[[i]] <- process_one(files[[i]], args, db)
    gc()
  }

  summary_df <- do.call(rbind, summaries)
  utils::write.csv(summary_df, file.path(args$output_dir, "batch_summary.csv"), row.names = FALSE)

  timing_files <- list.files(args$output_dir, pattern = "_step_timings\\.csv$", full.names = TRUE)
  if (length(timing_files)) {
    timings <- do.call(rbind, lapply(timing_files, utils::read.csv))
    utils::write.csv(timings, file.path(args$output_dir, "all_step_timings.csv"), row.names = FALSE)
    ok <- timings[timings$status == "ok", , drop = FALSE]
    if (nrow(ok)) {
      by_step <- aggregate(elapsed_sec ~ step, ok, sum)
      by_step <- by_step[order(-by_step$elapsed_sec), ]
      utils::write.csv(by_step, file.path(args$output_dir, "step_total_by_step.csv"), row.names = FALSE)
    }
  }

  print(summary_df)
}

main()
