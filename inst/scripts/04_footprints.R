#!/usr/bin/env Rscript

#####################################################################
# 04_footprints.R
# Footprints of the selected motifs, two groups compared as aggregates.
# The expected profile is either the GC bin or the k-mer model.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(logger)
  library(GenomicRanges)
  library(methylTFRAnnotationHg38)
  library(chromTFR)
})

motifSet <- "jaspar2020_distal"
tfSet <- "jaspar2020"
motifs <- c("IRF1", "IRF2", "IRF5", "FOS::JUN", "BATF::JUN")
expected.model <- "kmer"

GRP1 <- "Tororo"
GRP2 <- "Jinja"
group_colors <- c("Tororo" = "#B23A48", "Jinja" = "#3B6EA5")

plot.window <- 200L
flank.norm <- 30L
smooth.window <- 5L

ins.dir <- file.path("/scratch/icbb/igunduz/chromTFR", "insertions")
tfbs.dir <- file.path("/scratch/icbb/igunduz/chromTFR", "tfbs")
fig.dir <- file.path("/scratch/icbb/igunduz/chromTFR", "atac")
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)
distal.file <- file.path(
  "/scratch/icbb/igunduz/methylTFR_manuscript/github",
  "methylTFRAnnotationHg38_old/inst/extdata/distal_regions.RDS"
)

manifest <- fread(file.path(ins.dir, "group_manifest.tsv"))
manifest <- manifest[group %in% c(GRP1, GRP2)]
groups <- manifest$group
ins <- lapply(manifest$file, readRDS)
names(ins) <- groups
log_info("Groups: ", paste(groups, manifest$n, sep = ": n = ",
  collapse = ", "
))

if (expected.model == "kmer") {
  tfbs <- readRDS(file.path(tfbs.dir, paste0("tfbs_", motifSet, ".RDS")))
  expected <- readRDS(file.path(ins.dir,
    paste0("kmer_expected_", motifSet, ".RDS")
  ))
  motifs <- intersect(motifs, names(tfbs))
} else {
  tf_bindsites <- getTFbindsites(motifSet = tfSet)
  gcfreqs <- getGCfreq(motifSet = motifSet)
  gc_dist <- getGenomeGC()
  motifs <- intersect(motifs, intersect(names(tf_bindsites), names(gcfreqs)))
  enhancer <- if (motifSet == "jaspar2020_distal") {
    readRDS(distal.file)
  } else {
    NULL
  }
  if (!is.null(enhancer)) {
    gc_dist <- IRanges::subsetByOverlaps(gc_dist, enhancer)
  }
  tfbs <- lapply(motifs, function(m) {
    prepareTFBS(tf_bindsites[[m]], enhancer)
  })
  names(tfbs) <- motifs
  binIns <- lapply(ins, addGCBintoAccessome, gcdist = gc_dist)
  expected <- lapply(groups, function(g) {
    prof <- lapply(motifs, function(m) {
      computeAccExpectations(binIns[[g]], gcfreqs[[m]])
    })
    names(prof) <- motifs
    prof
  })
  names(expected) <- groups
}
if (length(motifs) == 0) stop("None of the motifs are in the annotation")

profiles <- lapply(motifs, function(m) {
  exp.m <- lapply(groups, function(g) expected[[g]][[m]])
  names(exp.m) <- groups
  accFootprintData(tfbs[[m]], ins, exp.m,
    smooth = smooth.window, flankNorm = flank.norm
  )
})
names(profiles) <- motifs
profiles <- Filter(Negate(is.null), profiles)
if (length(profiles) == 0) stop("No footprint could be drawn")

p <- plotAccFootprintGrid(profiles, tfbs, group_colors,
  method = "division", flankNorm = flank.norm, plotWindow = plot.window,
  corrected = paste(expected.model, "corrected")
)

file <- file.path(fig.dir, paste0(motifSet, "_", expected.model,
  "_footprints.pdf"
))
ggsave(file, p,
  width = 6 * length(profiles), height = 9, bg = "transparent"
)
log_success("Wrote ", file)
