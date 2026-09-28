#!/usr/bin/env Rscript

#####################################################################
# 05_plot_heatmap.R
# Bias corrected deviations per sample and motif, row wise z-scored.
# The Tn5 model is fitted per sample, so the correction survives the
# z-scoring.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(logger)
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ComplexHeatmap)
  library(methylTFRAnnotationHg38)
  library(chromTFR)
})

ins.dir <- "chromTFR/insertions"
fig.dir <- "chromTFR/figures"
if (!dir.exists(fig.dir)) dir.create(fig.dir, recursive = TRUE)

motifSet <- "altius"
motifs <- c("ap1_1", "ctcf", "irf_3")
GRP <- "group"
ann.cols <- c("group")

kmer <- 6L
max.ins <- 1e6
max.sites <- 25000
genome <- BSgenome.Hsapiens.UCSC.hg38

manifest <- fread(file.path(ins.dir, "manifest.tsv"))
if (!all(file.exists(manifest$file))) {
  stop("Insertions missing, run 01_prep_insertions.R")
}
samples <- manifest$sample
ann <- as.data.frame(readRDS(file.path(ins.dir, "sample_ann.RDS")))
regions <- readRDS(file.path(ins.dir, "regions.RDS"))

tf_bindsites <- getTFbindsites(motifSet = motifSet)
gcfreqs <- getGCfreq(motifSet = motifSet)
motifs <- intersect(motifs, intersect(names(tf_bindsites), names(gcfreqs)))
if (length(motifs) == 0) stop("None of the motifs are in the annotation")
tfbs <- lapply(motifs, function(m) prepareTFBS(tf_bindsites[[m]]))
names(tfbs) <- motifs

bg.file <- file.path(ins.dir, paste0("kmer_bg_", kmer, "mer.RDS"))
if (file.exists(bg.file)) {
  bg <- readRDS(bg.file)
} else {
  bg <- kmerBackground(regions, genome, k = kmer)
  saveRDS(bg, bg.file)
}

bias.file <- file.path(ins.dir,
  paste0("kmer_bias_persample_", kmer, "mer.RDS")
)
bias <- if (file.exists(bias.file)) readRDS(bias.file) else list()
for (s in setdiff(samples, names(bias))) {
  log_info("Fitting the Tn5 model of ", s)
  bias[[s]] <- computeKmerBias(
    readRDS(manifest$file[match(s, samples)]), regions, genome,
    k = kmer, max.ins = max.ins, bg = bg
  )
  gc()
}
saveRDS(bias, bias.file)
bias <- bias[samples]

exp.file <- file.path(ins.dir, paste0("exp_dev_", motifSet, ".RDS"))
expDev <- if (file.exists(exp.file)) readRDS(exp.file) else NULL
if (!is.null(expDev) && !identical(colnames(expDev), samples)) expDev <- NULL
todo <- setdiff(motifs, rownames(expDev))
if (length(todo) > 0) {
  expDev <- rbind(expDev, accExpectedMatrix(tfbs[todo], bias, genome,
    k = kmer, max.sites = max.sites
  ))
  saveRDS(expDev, exp.file)
}
expDev <- expDev[motifs, , drop = FALSE]

obs.file <- file.path(ins.dir, paste0("obs_dev_", motifSet, ".RDS"))
obsDev <- if (file.exists(obs.file)) readRDS(obs.file) else NULL
if (!is.null(obsDev) && !identical(colnames(obsDev), samples)) obsDev <- NULL
todo <- setdiff(motifs, rownames(obsDev))
if (length(todo) > 0) {
  log_info("Computing ", length(todo), " motifs, ", nrow(obsDev),
    " from the cache"
  )
  obsDev <- rbind(obsDev, accObservedMatrix(manifest$file, tfbs[todo],
    samples
  ))
  saveRDS(obsDev, obs.file)
}
obsDev <- obsDev[motifs, , drop = FALSE]

devs <- obsDev - expDev
devs <- devs[rowSums(is.na(devs)) == 0, , drop = FALSE]
z <- computeRowZScore(devs)
log_info("Heatmap of ", nrow(z), " motifs over ", ncol(z), " samples")

ht <- plotAccDeviationHeatmap(z, ann, ann.cols, split = GRP)
file <- file.path(fig.dir, paste0(motifSet, "_deviation_heatmap.pdf"))
pdf(file, width = 8, height = 0.3 * nrow(z) + 3)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
invisible(dev.off())
log_success("Wrote ", file)

fwrite(data.table(motif = rownames(z), z),
  file.path(ins.dir, paste0("deviation_zscores_", motifSet, ".csv"))
)
