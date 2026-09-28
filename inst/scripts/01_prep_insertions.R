#!/usr/bin/env Rscript

#####################################################################
# 01_prep_insertions.R
# Insertion sites per sample and per group of a ChrAccR dataset.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(logger)
  library(GenomicRanges)
  library(ChrAccR)
  library(chromTFR)
})

dsa.dir <- "/icbb/projects/mmaran/malaria/pilot/longPeaksYear0.5/"
out.dir <- file.path("/scratch/icbb/igunduz/chromTFR", "insertions")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

region.type <- ".peaks.cons"
region.extend <- 500L
# Set to c(0L, 0L) if the fragments are already Tn5 shifted
tn5.shift <- c(4L, -5L)

GRP <- "region"
GRP1 <- "Tororo"
GRP2 <- "Jinja"

dsa <- loadAccDataset(dsa.dir)
samples <- getAccSamples(dsa)
ann <- getAccSampleAnnotation(dsa)
if (!GRP %in% colnames(ann)) {
  stop("Column '", GRP, "' not found, available: ",
    paste(colnames(ann), collapse = ", ")
  )
}
grp <- as.character(ann[[GRP]])
log_info("Samples: ", length(samples))

regions <- getAccRegions(dsa, regionType = region.type,
  extend = region.extend
)
saveRDS(regions, file.path(out.dir, "regions.RDS"))
saveRDS(ann, file.path(out.dir, "sample_ann.RDS"))

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

fwrite(data.table(sample = samples, group = grp, file = sample.files),
  file.path(out.dir, "manifest.tsv"), sep = "\t"
)

groups <- c(GRP1, GRP2)
group.files <- file.path(out.dir, paste0("ins_group_", groups, ".RDS"))
for (i in seq_along(groups)) {
  if (file.exists(group.files[i])) next
  files <- sample.files[grp == groups[i]]
  if (length(files) == 0) stop("No sample of the group ", groups[i])
  log_info("Merging ", length(files), " samples of ", groups[i])
  saveRDS(mergeInsertionSites(files), group.files[i])
  gc()
}

fwrite(data.table(
  group = groups,
  n = c(sum(grp == GRP1), sum(grp == GRP2)),
  file = group.files
), file.path(out.dir, "group_manifest.tsv"), sep = "\t")
log_success("Wrote the manifests")
