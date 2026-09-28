#' Deviation scores of a ChrAccR dataset
#'
#' The ATAC counterpart of \code{run_methyltfr}. The methylTFR annotation
#' is used unchanged, methylation levels are replaced by insertion
#' densities per base and per site, and the expected profile comes from
#' the GC bin composition of the binding sites.
#'
#' @param dsa A \code{DsATAC} object or the path of a ChrAccR dataset.
#' @param tf_bindsites A \code{GRangesList} of binding site positions.
#' @param gcfreqs A \code{list} of GC bin frequency tables.
#' @param gc_dist A \code{GRanges} of the genome wide GC distribution.
#' @param sampleIds Samples to process, all of them by default.
#' @param motifs Motifs to process, all of them by default.
#' @param regions Optional \code{GRanges} restricting the insertions.
#' @param enhancer Optional \code{GRanges} restricting the binding sites.
#' @param threads Thread count for parallel processing.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @param shift Offsets applied to the forward and reverse fragment ends.
#' @param normalize If TRUE, counts are scaled to insertions per million.
#' @param sample_ann Optional \code{data.frame} of sample annotation.
#' @return A \code{SummarizedExperiment} with the deviations, their
#' row-wise z-scores and the expected deviations.
#' @examples
#' \dontrun{
#' devs <- runChromTFR(dsa, tf_bindsites, gcfreqs, gc_dist)
#' }
#' @importFrom BiocParallel bplapply
#' @importFrom IRanges subsetByOverlaps
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom logger log_info log_success
#' @export
runChromTFR <- function(
    dsa, tf_bindsites = NULL, gcfreqs = NULL, gc_dist = NULL,
    sampleIds = NULL, motifs = NULL, regions = NULL, enhancer = NULL,
    threads = 1, ignoreStrand = TRUE, shift = c(4L, -5L),
    normalize = TRUE, sample_ann = NULL
) {
    if (is.character(dsa)) {
        dsa <- loadAccDataset(dsa)
    }
    if (any(vapply(
        list(tf_bindsites, gcfreqs, gc_dist), is.null, logical(1)
    ))) {
        stop("Please load the annotation objects for the given genome")
    }
    if (is.null(sampleIds)) {
        sampleIds <- getAccSamples(dsa)
    }
    if (is.null(sample_ann)) {
        sample_ann <- getAccSampleAnnotation(dsa)
        sample_ann <- sample_ann[match(
            sampleIds, getAccSamples(dsa)
        ), , drop = FALSE]
    }
    if (is.null(motifs)) {
        motifs <- validAccMotifs(tf_bindsites, gcfreqs)
    }
    BPPARAM <- bpparamFromThreads(threads)
    if (!is.null(enhancer)) {
        gc_dist <- subsetByOverlaps(gc_dist, enhancer,
            ignore.strand = ignoreStrand
        )
    }
    readIns <- makeInsertionReader(dsa, sampleIds,
        regions = regions, shift = shift, normalize = normalize
    )

    dev <- matrix(NA_real_,
        nrow = length(motifs), ncol = length(sampleIds),
        dimnames = list(motifs, sampleIds)
    )
    expDev <- dev
    for (i in seq_along(sampleIds)) {
        log_info("Processing ", sampleIds[i])
        ins <- readIns(i)
        binIns <- addGCBintoAccessome(ins, gc_dist, ignoreStrand)
        res <- bplapply(motifs, computeAccDeviation,
            ins = ins, tf_bindsites = tf_bindsites, gcfreqs = gcfreqs,
            binIns = binIns, enhancer = enhancer,
            ignoreStrand = ignoreStrand, BPPARAM = BPPARAM
        )
        dev[, i] <- vapply(res, function(x) x$dev[1], numeric(1))
        expDev[, i] <- vapply(res, function(x) x$exp_dev[1], numeric(1))
        rm(ins, binIns, res)
        gc()
    }
    log_success("Computed all deviations successfully")

    SummarizedExperiment(
        assays = list(
            deviations = dev,
            z = computeRowZScore(dev),
            expected = expDev
        ),
        colData = DataFrame(sample_ann),
        rowData = DataFrame(motifs = motifs)
    )
}
