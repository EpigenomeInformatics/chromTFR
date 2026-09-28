#!/usr/bin/env Rscript

#####################################################################
# 03_prep_kmer_bias.R
# Tn5 k-mer preference per group and the expected profile per motif.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(logger)
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(chromTFR)
})

args <- commandArgs(trailingOnly = TRUE)
motifSet <- if (length(args) > 0) args[1] else "jaspar2020_distal"

ins.dir <- file.path("/scratch/icbb/igunduz/chromTFR", "insertions")
tfbs.dir <- file.path("/scratch/icbb/igunduz/chromTFR", "tfbs")
genome <- BSgenome.Hsapiens.UCSC.hg38
kmer <- 6L
max.ins <- 2e6
max.sites <- 25000

log_info("Motif set: ", motifSet)
exp.file <- file.path(ins.dir, paste0("kmer_expected_", motifSet, ".RDS"))
if (file.exists(exp.file)) {
  stop("Expected profiles already exist, delete them to recompute: ", exp.file)
}

regions <- readRDS(file.path(ins.dir, "regions.RDS"))
manifest <- fread(file.path(ins.dir, "group_manifest.tsv"))
tfbs <- readRDS(file.path(tfbs.dir, paste0("tfbs_", motifSet, ".RDS")))

bias.file <- file.path(ins.dir, paste0("kmer_bias_", kmer, "mer.RDS"))
if (file.exists(bias.file)) {
  bias <- readRDS(bias.file)
} else {
  bg <- kmerBackground(regions, genome, k = kmer)
  bias <- lapply(manifest$file, function(f) {
    computeKmerBias(readRDS(f), regions, genome,
      k = kmer, max.ins = max.ins, bg = bg
    )
  })
  names(bias) <- manifest$group
  saveRDS(bias, bias.file)
  log_success("Wrote ", bias.file)
}

expected <- lapply(manifest$group, function(g) {
  prof <- lapply(names(tfbs), function(m) {
    log_info("Bias profile of ", m, " in ", g)
    kmerBiasProfile(tfbs[[m]], genome, bias[[g]],
      k = kmer, max.sites = max.sites
    )
  })
  names(prof) <- names(tfbs)
  prof
})
names(expected) <- manifest$group
saveRDS(expected, exp.file)
log_success("Wrote ", exp.file)
