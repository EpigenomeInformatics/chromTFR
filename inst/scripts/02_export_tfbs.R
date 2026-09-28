#!/usr/bin/env Rscript

#####################################################################
# 02_export_tfbs.R
# Binding sites of the selected motifs, resized to the footprint window
# and optionally restricted to a region set.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(logger)
  library(GenomicRanges)
  library(methylTFRAnnotationHg38)
  library(chromTFR)
})

args <- commandArgs(trailingOnly = TRUE)
motifSet <- if (length(args) > 0) args[1] else "jaspar2020"
# The distal set shares its binding sites with the genome wide one
tfSet <- sub("_distal$", "", motifSet)
motifs <- c("CTCF", "FOS::JUN", "IRF1")

out.dir <- "chromTFR/tfbs"
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
# Only needed for a distal motif set
distal.file <- "path/to/distal_regions.RDS"

log_info("Motif set: ", motifSet)
out.file <- file.path(out.dir, paste0("tfbs_", motifSet, ".RDS"))
if (file.exists(out.file)) {
  stop("Binding sites already exist, delete them to re-export: ", out.file)
}

tf_bindsites <- getTFbindsites(motifSet = tfSet)
absent <- setdiff(motifs, names(tf_bindsites))
if (length(absent) > 0) {
  log_warn(paste(absent, collapse = ", "), " not in ", tfSet, ", skipped")
  motifs <- setdiff(motifs, absent)
}
if (length(motifs) == 0) stop("None of the motifs are in the annotation")

enhancer <- NULL
if (grepl("_distal$", motifSet)) {
  if (!file.exists(distal.file)) {
    stop("Distal regions file does not exist: ", distal.file)
  }
  enhancer <- readRDS(distal.file)
  log_info("Loaded ", length(enhancer), " distal regions")
}

tfbs <- lapply(motifs, function(m) {
  gr <- prepareTFBS(tf_bindsites[[m]], enhancer)
  log_info(m, ": ", length(gr), " sites of width ", width(gr)[1])
  gr
})
names(tfbs) <- motifs
saveRDS(tfbs, out.file)
log_success("Wrote ", out.file)
