#!/usr/bin/env Rscript

#####################################################################
# run_footprints_vd2.R
# One script for the vd2 year 0.5 dataset: insertions, group
# aggregates, Tn5 k-mer bias and the footprints of the Altius
# archetypes, low against high aEIR. Every step caches.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(logger)
  library(GenomicRanges)
  library(ChrAccR)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(methylTFRAnnotationHg38)
  library(chromTFR)
})

dsa.dir <- paste0(
  "/icbb/projects/mmaran/malaria/dsATAC/longSaturationAnalysis/",
  "vd2.year0.5/data/dsATAC_filtered"
)
out.dir <- file.path("/scratch/icbb/igunduz/chromTFR", "vd2_year0.5")
fig.dir <- file.path("/scratch/icbb/igunduz/chromTFR", "atac")
for (d in c(out.dir, fig.dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

motifSet <- "altius"
motifs <- c("ap1_1", "ets_2", "ewsr1_fli1", "tead")

# Responders, marked in red in the chromVAR heatmap
keep.samples <- c(
  "VD2_2016_d3340_88", "VD2_2016_d3381_164", "VD2_2016_d3401_216",
  "VD2_2016_d3278_60", "VD2_2016_d3280_168", "VD2_2016_d3369_112",
  "VD2_2013_d3278_52", "VD2_2013_d3340_80", "VD2_2014_d3381_160",
  "VD2_2013_d3401_208", "VD2_2013_d3280_44", "VD2_2013_d3369_104",
  "VD2_2014_d3369_108", "VD2_2014_d3398_200"
)

GRP <- "EirGrp"
GRP1 <- "low"
GRP2 <- "high"
group_colors <- c("low" = "#3B6EA5", "high" = "#B23A48")

region.type <- ".peaks.cons"
region.extend <- 500L
tn5.shift <- c(4L, -5L)

kmer <- 6L
max.ins <- 2e6
max.sites <- 25000
genome <- BSgenome.Hsapiens.UCSC.hg38

plot.window <- 200L
flank.norm <- 30L
smooth.window <- 5L

#####################################################################
# Insertion sites and group aggregates
#####################################################################

dsa <- loadAccDataset(dsa.dir)
samples <- getAccSamples(dsa)
ann <- getAccSampleAnnotation(dsa)

absent.samples <- setdiff(keep.samples, samples)
if (length(absent.samples) > 0) {
  stop("Samples not in the dataset: ", paste(absent.samples, collapse = ", "))
}
idx <- match(keep.samples, samples)
samples <- samples[idx]
ann <- ann[idx, , drop = FALSE]

if (!GRP %in% colnames(ann)) {
  stop("Column '", GRP, "' not found, available: ",
    paste(colnames(ann), collapse = ", ")
  )
}
grp <- as.character(ann[[GRP]])
if (!all(c(GRP1, GRP2) %in% grp)) {
  stop("Groups not found in ", GRP, ", available: ",
    paste(unique(grp), collapse = ", ")
  )
}
log_info("Samples: ", length(samples), ", ", GRP1, " n = ", sum(grp == GRP1),
  ", ", GRP2, " n = ", sum(grp == GRP2)
)

regions.file <- file.path(out.dir, "regions.RDS")
if (file.exists(regions.file)) {
  regions <- readRDS(regions.file)
} else {
  regions <- getAccRegions(dsa, regionType = region.type,
    extend = region.extend
  )
  saveRDS(regions, regions.file)
}

tags <- gsub("[^A-Za-z0-9_.-]", "_", samples)
sample.files <- file.path(out.dir, paste0("ins_", tags, ".RDS"))
for (i in seq_along(samples)) {
  if (file.exists(sample.files[i])) next
  saveRDS(getInsertionSites(dsa, samples[i],
    regions = regions, shift = tn5.shift, normalize = TRUE
  ), sample.files[i])
  gc()
  log_success("Wrote ", sample.files[i])
}

groups <- c(GRP1, GRP2)
group.files <- file.path(out.dir, paste0("ins_group_", groups, ".RDS"))
names(group.files) <- groups
for (g in groups) {
  if (file.exists(group.files[g])) next
  log_info("Merging ", sum(grp == g), " samples of ", g)
  saveRDS(mergeInsertionSites(sample.files[grp == g]), group.files[g])
  gc()
}
ins <- lapply(groups, function(g) readRDS(group.files[g]))
names(ins) <- groups

#####################################################################
# Binding sites and the Tn5 k-mer background
#####################################################################

tf_bindsites <- getTFbindsites(motifSet = motifSet)
gcfreqs <- getGCfreq(motifSet = motifSet)
absent <- setdiff(motifs, intersect(names(tf_bindsites), names(gcfreqs)))
if (length(absent) > 0) {
  log_warn(paste(absent, collapse = ", "), " not in ", motifSet, ", skipped")
  motifs <- setdiff(motifs, absent)
}
if (length(motifs) == 0) stop("None of the motifs are in the annotation")

tfbs <- lapply(motifs, function(m) prepareTFBS(tf_bindsites[[m]]))
names(tfbs) <- motifs

bias.file <- file.path(out.dir, paste0("kmer_bias_", kmer, "mer.RDS"))
if (file.exists(bias.file)) {
  bias <- readRDS(bias.file)
} else {
  bg <- kmerBackground(regions, genome, k = kmer)
  bias <- lapply(groups, function(g) {
    computeKmerBias(ins[[g]], regions, genome,
      k = kmer, max.ins = max.ins, bg = bg
    )
  })
  names(bias) <- groups
  saveRDS(bias, bias.file)
  log_success("Wrote ", bias.file)
}

exp.file <- file.path(out.dir, paste0("kmer_expected_", motifSet, ".RDS"))
if (file.exists(exp.file)) {
  expected <- readRDS(exp.file)
} else {
  expected <- lapply(groups, function(g) {
    prof <- lapply(motifs, function(m) {
      log_info("Bias profile of ", m, " in ", g)
      kmerBiasProfile(tfbs[[m]], genome, bias[[g]],
        k = kmer, max.sites = max.sites
      )
    })
    names(prof) <- motifs
    prof
  })
  names(expected) <- groups
  saveRDS(expected, exp.file)
  log_success("Wrote ", exp.file)
}

#####################################################################
# Footprints
#####################################################################

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
  corrected = "k-mer corrected"
)
file <- file.path(fig.dir, paste0("vd2_", motifSet, "_aEIR_footprints.pdf"))
ggsave(file, p, width = 6 * length(profiles), height = 9,
  bg = "transparent"
)
log_success("Wrote ", file)
