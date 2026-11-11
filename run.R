#!/usr/bin/env Rscript
# scDblFinder doublet detection module.
#
# runs on the per-dataset h5ad and uses scDblFinder's own `samples` argument so doublets are
# searched for independently within each pool/sample 

suppressPackageStartupMessages({
  library(scDblFinder)
  library(SingleCellExperiment)
  library(anndataR)
  library(data.table)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("scDblFinder doublet detection module")
p <- add_base_args(p)              # --output_dir, --name
p <- add_stage_args(p, "doublets") # --dataset_h5ad
p <- add_argument(p, "--dbr_sd", type = "numeric", default = 0.015,
                   help = "scDblFinder dbr.sd: uncertainty band around the expected doublet rate")
args <- parse_args(p, argv = strip_output_args())

cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  cat("  reading dataset h5ad ..\n")
  sce <- read_h5ad(args$dataset_h5ad, as = "SingleCellExperiment", x_mapping = "counts")

  if (!"pool_id" %in% colnames(colData(sce))) {
    stop("colData(sce) is missing 'pool_id")
  }

  cat("  running scDblFinder..\n")
  dbl <- scDblFinder(sce, samples = "pool_id", returnType = "table",
    dbr.sd = args$dbr_sd, verbose = FALSE)

  dbl_dt <- as.data.table(as.data.frame(dbl), keep.rownames = "cell_id")
  dbl_dt <- dbl_dt[type == "real"]

  out_dt <- data.table(
    pool_id = colData(sce)$pool_id[match(dbl_dt$cell_id, colnames(sce))],
    cell_id = dbl_dt$cell_id,
    scdbl_score = dbl_dt$score,
    scdbl_class = dbl_dt$class,
    multiplet = dbl_dt$class == "doublet"
  )
  cat(sprintf("  %d / %d cells called doublet\n",
    sum(out_dt$multiplet), nrow(out_dt)))

  out_path <- file.path(args$output_dir, paste0(args$name, "_doublets.tsv"))
  fwrite(out_dt, out_path, sep = "\t", quote = FALSE)
  cat(sprintf("  wrote: %s\n", out_path))
}

if (sys.nframe() == 0L) {
  main()
}
