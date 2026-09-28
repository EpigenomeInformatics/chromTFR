#' Check that ChrAccR is installed
#'
#' @return Invisible TRUE, called for the error it raises.
#' @keywords internal
.requireChrAccR <- function() {
    if (!requireNamespace("ChrAccR", quietly = TRUE)) {
        stop("The ChrAccR package is required for this function")
    }
    invisible(TRUE)
}


#' Load a ChrAccR dataset
#'
#' @param path Path of the ChrAccR dataset directory.
#' @return A \code{DsATAC} object.
#' @examples
#' \dontrun{
#' dsa <- loadAccDataset("dsATAC_filtered")
#' }
#' @importFrom logger log_success
#' @export
loadAccDataset <- function(path) {
    .requireChrAccR()
    if (is.null(path) || !is.character(path) || !dir.exists(path)) {
        stop("Dataset directory does not exist, please check the path")
    }
    dsa <- ChrAccR::loadDsAcc(path)
    log_success("Loaded the dataset from ", path)
    return(dsa)
}


#' Sample identifiers of a ChrAccR dataset
#'
#' @param dsa A \code{DsATAC} object.
#' @return A character vector of sample identifiers.
#' @examples
#' \dontrun{
#' getAccSamples(dsa)
#' }
#' @export
getAccSamples <- function(dsa) {
    .requireChrAccR()
    ids <- ChrAccR::getSamples(dsa)
    if (length(ids) == 0) {
        stop("No samples found in the dataset")
    }
    return(as.character(ids))
}


#' Sample annotation of a ChrAccR dataset
#'
#' @param dsa A \code{DsATAC} object.
#' @return A \code{data.frame} with one row per sample.
#' @examples
#' \dontrun{
#' getAccSampleAnnotation(dsa)
#' }
#' @export
getAccSampleAnnotation <- function(dsa) {
    .requireChrAccR()
    ann <- as.data.frame(ChrAccR::getSampleAnnot(dsa),
        stringsAsFactors = FALSE
    )
    ids <- getAccSamples(dsa)
    if (nrow(ann) != length(ids)) {
        stop("Sample annotation must have one row per sample")
    }
    if (!"sampleName" %in% colnames(ann)) {
        ann$sampleName <- ids
    }
    return(ann)
}


#' Region set of a ChrAccR dataset
#'
#' Regions are extended so that they still cover the motif flanks.
#'
#' @param dsa A \code{DsATAC} object.
#' @param regionType Region type stored in the dataset.
#' @param extend Number of bases added on each side.
#' @return A \code{GRanges} object.
#' @examples
#' \dontrun{
#' getAccRegions(dsa, regionType = ".peaks.cons", extend = 500)
#' }
#' @importFrom GenomicRanges resize width
#' @export
getAccRegions <- function(dsa, regionType = "peaks", extend = 500) {
    .requireChrAccR()
    rts <- ChrAccR::getRegionTypes(dsa)
    if (!regionType %in% rts) {
        stop(sprintf(
            "Region type '%s' was not found, available: %s",
            regionType, paste(rts, collapse = ", ")
        ))
    }
    gr <- ChrAccR::getCoord(dsa, regionType)
    if (is.numeric(extend) && extend > 0) {
        gr <- resize(gr, width(gr) + 2 * extend, fix = "center")
    }
    return(gr)
}


#' Tn5 insertion sites of one sample
#'
#' Both fragment ends are kept and shifted to the Tn5 cut position.
#' Duplicated positions are collapsed into a per base count, held in the
#' \code{score} column so that the object can be used like a methylTFR
#' methylome.
#'
#' @param dsa A \code{DsATAC} object.
#' @param sampleId Sample identifier.
#' @param regions Optional \code{GRanges} restricting the insertions.
#' @param shift Offsets applied to the forward and reverse fragment ends.
#' @param normalize If TRUE, counts are scaled to insertions per million.
#' @return A \code{GRanges} of single base insertion sites.
#' @examples
#' \dontrun{
#' ins <- getInsertionSites(dsa, "sample_1", regions = peaks)
#' }
#' @importFrom GenomicRanges GRanges seqnames start end
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom logger log_info
#' @import data.table
#' @export
getInsertionSites <- function(
    dsa, sampleId, regions = NULL, shift = c(4L, -5L), normalize = TRUE
) {
    .requireChrAccR()
    frags <- ChrAccR::getFragmentGr(dsa, sampleId)
    if (is.null(frags) || length(frags) == 0) {
        stop("No fragment data found for sample ", sampleId)
    }
    ins <- c(
        GRanges(seqnames(frags), IRanges(start(frags) + shift[1], width = 1)),
        GRanges(seqnames(frags), IRanges(end(frags) + shift[2], width = 1))
    )
    total <- length(ins)
    if (!is.null(regions)) {
        ins <- subsetByOverlaps(ins, regions, ignore.strand = TRUE)
    }
    dt <- data.table(
        chr = as.character(seqnames(ins)),
        pos = start(ins)
    )[, .(count = .N), by = .(chr, pos)]

    sf <- if (isTRUE(normalize)) 1e6 / total else 1
    gr <- GRanges(dt$chr, IRanges(dt$pos, width = 1),
        score = dt$count * sf,
        coverage = dt$count
    )
    log_info(sampleId, ": ", total, " insertions, ", length(gr), " positions")
    return(sort(gr))
}


#' Pool the insertion sites of several samples
#'
#' The files are read one at a time and collapsed after every sample, so
#' a group costs the memory of its largest sample. Scores are averaged
#' over the samples, counts are summed.
#'
#' @param files Paths of the per sample insertion site files.
#' @return A \code{GRanges} of single base insertion sites.
#' @examples
#' \dontrun{
#' pooled <- mergeInsertionSites(c("ins_a.RDS", "ins_b.RDS"))
#' }
#' @importFrom GenomicRanges GRanges seqnames start
#' @importFrom IRanges IRanges
#' @import data.table
#' @export
mergeInsertionSites <- function(files) {
    if (length(files) == 0) {
        stop("No insertion site files to merge")
    }
    acc <- NULL
    for (f in files) {
        gr <- readRDS(f)
        dt <- data.table(
            chr = as.character(seqnames(gr)),
            pos = start(gr),
            score = gr$score,
            coverage = gr$coverage
        )
        acc <- rbind(acc, dt)[, .(
            score = sum(score), coverage = sum(coverage)
        ), by = .(chr, pos)]
        rm(gr, dt)
        gc()
    }
    acc[, score := score / length(files)]
    gr <- GRanges(acc$chr, IRanges(acc$pos, width = 1),
        score = acc$score,
        coverage = acc$coverage
    )
    return(sort(gr))
}


#' Per sample reader for runChromTFR
#'
#' @param dsa A \code{DsATAC} object.
#' @param sampleIds A character vector of sample identifiers.
#' @param regions Optional \code{GRanges} restricting the insertions.
#' @param shift Offsets applied to the forward and reverse fragment ends.
#' @param normalize If TRUE, counts are scaled to insertions per million.
#' @return A function of a single sample index.
#' @examples
#' \dontrun{
#' reader <- makeInsertionReader(dsa, getAccSamples(dsa))
#' }
#' @export
makeInsertionReader <- function(
    dsa, sampleIds, regions = NULL, shift = c(4L, -5L), normalize = TRUE
) {
    function(i) {
        getInsertionSites(dsa, sampleIds[i],
            regions = regions, shift = shift, normalize = normalize
        )
    }
}
