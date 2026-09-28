#' Prepare the binding sites of one motif
#'
#' @param tfbs A \code{GRanges} of binding sites of one motif.
#' @param enhancer Optional \code{GRanges} restricting the sites.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @return A \code{GRanges} object.
#' @importFrom GenomicRanges resize width
#' @importFrom IRanges subsetByOverlaps
#' @examples
#' \dontrun{
#' prepareTFBS(tf_bindsites[["CTCF"]])
#' }
#' @export
prepareTFBS <- function(tfbs, enhancer = NULL, ignoreStrand = TRUE) {
    tfbs <- resize(tfbs, width(tfbs)[1] + 130, fix = "center")
    if (!is.null(enhancer)) {
        tfbs <- subsetByOverlaps(tfbs, enhancer, ignore.strand = ignoreStrand)
    }
    return(tfbs)
}


#' Insertion profile around a set of binding sites
#'
#' Positions without an insertion are kept as zeros, so the profile is an
#' insertion density per base and per site rather than a mean over the
#' covered positions only.
#'
#' @param ins A \code{GRanges} of insertion sites with a \code{score}
#' column.
#' @param tfbs A \code{GRanges} of binding sites, prepared with
#' \code{prepareTFBS}.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @return A \code{data.table} with the columns \code{x} and
#' \code{avg_ins}, or NULL if no insertion falls in the sites.
#' @examples
#' \dontrun{
#' accProfile(ins, prepareTFBS(tf_bindsites[["CTCF"]]))
#' }
#' @importFrom GenomicRanges findOverlaps start end width
#' @importFrom S4Vectors mcols
#' @import data.table
#' @export
accProfile <- function(ins, tfbs, ignoreStrand = TRUE) {
    hits <- findOverlaps(ins, tfbs,
        type = "within", ignore.strand = ignoreStrand
    )
    if (length(hits@from) == 0) {
        return(NULL)
    }
    S4Vectors::mcols(tfbs)$mid_point <- round(end(tfbs) +
        ((start(tfbs) - end(tfbs)) / 2))
    prof <- data.table(
        x = start(ins[hits@from]) - tfbs[hits@to]$mid_point,
        ins = ins[hits@from]$score
    )[, .(ins = sum(ins)), by = x]

    flank <- floor(width(tfbs)[1] / 2)
    grid <- data.table(x = seq(-flank, flank))
    prof <- merge(grid, prof, by = "x", all.x = TRUE)
    prof[is.na(ins), ins := 0]
    prof[, avg_ins := ins / length(tfbs)]
    return(prof[, .(x, avg_ins)])
}


#' Central over flanking insertion density
#'
#' Each interval is divided by its width, so windows of different size
#' stay comparable.
#'
#' @param data A \code{data.table} with the columns \code{x} and
#' \code{avg_ins}.
#' @param breaks Interval borders of the footprint windows.
#' @return A numeric value.
#' @examples
#' prof <- data.table::data.table(x = -250:250, avg_ins = 1)
#' accDeviationScore(prof)
#' @importFrom stats na.omit
#' @import data.table
#' @export
accDeviationScore <- function(
    data, breaks = c(-250, -200, -25, 25, 200, 250)
) {
    data <- copy(data)
    data[, cuts := cut(x, breaks)]
    ivl <- data[, .(dens = sum(avg_ins)), by = cuts]
    ivl <- na.omit(ivl[order(cuts)])
    ivl[, dens := dens / diff(breaks)[as.integer(cuts)]]
    n <- nrow(ivl)
    if (n == 0) {
        return(NA_real_)
    }
    avg_first_last <- (ivl$dens[1] + ivl$dens[n]) / 2
    return(ivl$dens[(n + 1) %/% 2] / avg_first_last)
}


#' Insertion density of each GC bin
#'
#' The total width of the bin is the denominator, so bases without an
#' insertion count as zero.
#'
#' @param ins A \code{GRanges} of insertion sites.
#' @param gcdist A \code{GRanges} of the genome wide GC distribution.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @return A \code{matrix} of GC bins and their insertion density.
#' @examples
#' \dontrun{
#' addGCBintoAccessome(ins, gc_dist)
#' }
#' @importFrom GenomicRanges findOverlaps width
#' @importFrom methods is
#' @import data.table
#' @export
addGCBintoAccessome <- function(ins, gcdist, ignoreStrand = TRUE) {
    if (is.null(ins) || !is(ins, "GRanges")) {
        stop("Please provide valid insertion sites as GRanges")
    }
    if (is.null(gcdist) || !is(gcdist, "GRanges")) {
        stop("Please provide a valid GC distribution")
    }
    bp <- data.table(
        gcbin = gcdist$GC_bin, w = width(gcdist)
    )[, .(bp = sum(w)), by = gcbin]

    hits <- findOverlaps(ins, gcdist, ignore.strand = ignoreStrand)
    if (length(hits@from) == 0) {
        stop("No insertion sites found in the GC distribution")
    }
    obs <- data.table(
        ins = ins[hits@from]$score,
        gcbin = gcdist[hits@to]$GC_bin
    )[, .(total = sum(ins)), by = gcbin]

    acc <- merge(bp, obs, by = "gcbin", all.x = TRUE)
    acc[is.na(total), total := 0]
    acc[, avg_ins := total / bp]
    acc <- acc[order(gcbin)]
    return(as.matrix(acc[, .(gcbin, avg_ins)]))
}


#' Expected insertion profile from the GC composition
#'
#' @param binIns GC bin densities from \code{\link{addGCBintoAccessome}}.
#' @param gcfreq A \code{matrix} of GC bin frequencies of one motif.
#' @return A \code{data.table} with the columns \code{x} and
#' \code{avg_ins}.
#' @examples
#' \dontrun{
#' computeAccExpectations(binIns, gcfreqs[["CTCF"]])
#' }
#' @import data.table
#' @export
computeAccExpectations <- function(binIns, gcfreq) {
    if (!is.matrix(binIns) || !is.matrix(gcfreq)) {
        stop("Please provide the GC bin tables as matrices")
    }
    if (nrow(gcfreq) != nrow(binIns)) {
        stop("The GC bins of the annotation and the dataset do not match")
    }
    expData <- t(gcfreq) %*% binIns[, 2]
    mpos <- round(seq(-floor(length(expData) / 2),
        floor(length(expData) / 2),
        length.out = length(expData)
    ))
    return(data.table(x = mpos, avg_ins = as.numeric(expData)))
}


#' Bias corrected deviation of one motif
#'
#' @param motif Motif name.
#' @param ins A \code{GRanges} of insertion sites.
#' @param tf_bindsites A \code{GRangesList} of binding site positions.
#' @param gcfreqs A \code{list} of GC bin frequency tables.
#' @param binIns GC bin densities of the sample.
#' @param enhancer Optional \code{GRanges} restricting the sites.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @return A \code{data.table} with the deviation and the expected
#' deviation.
#' @examples
#' \dontrun{
#' computeAccDeviation("CTCF", ins, tf_bindsites, gcfreqs, binIns)
#' }
#' @import data.table
#' @export
computeAccDeviation <- function(
    motif, ins, tf_bindsites, gcfreqs, binIns, enhancer = NULL,
    ignoreStrand = TRUE
) {
    if (is.null(motif) || !is.character(motif)) {
        stop("Please provide a valid motif name")
    }
    tfbs <- prepareTFBS(tf_bindsites[[motif]], enhancer, ignoreStrand)
    obs <- accProfile(ins, tfbs, ignoreStrand)
    if (is.null(obs) || length(tfbs) == 0) {
        return(data.table(dev = NA_real_, exp_dev = NA_real_))
    }
    expProf <- computeAccExpectations(binIns, gcfreqs[[motif]])
    obsDev <- accDeviationScore(obs)
    expDev <- accDeviationScore(expProf)
    return(data.table(dev = obsDev - expDev, exp_dev = expDev))
}


#' Drop motifs without binding sites or frequency table
#'
#' @param tf_bindsites A \code{GRangesList} of binding site positions.
#' @param gcfreqs A \code{list} of GC bin frequency tables.
#' @return A character vector of motif names.
#' @importFrom logger log_info
#' @keywords internal
validAccMotifs <- function(tf_bindsites, gcfreqs) {
    motifs <- names(gcfreqs)
    keep <- vapply(motifs, function(m) {
        !is.null(tf_bindsites[[m]]) && length(tf_bindsites[[m]]) > 0 &&
            !is.null(gcfreqs[[m]])
    }, logical(1))
    if (any(!keep)) {
        log_info("Discarding ", sum(!keep), " motifs without annotation")
        motifs <- motifs[keep]
    }
    if (length(motifs) == 0) {
        stop("No valid motifs remaining after validation")
    }
    return(motifs)
}


#' Row-wise z-score of a matrix
#'
#' @param mat A \code{matrix}.
#' @return A \code{matrix} of row-wise z-scores.
#' @examples
#' computeRowZScore(matrix(c(1, 2, 3, 4), nrow = 2))
#' @importFrom matrixStats rowMeans2 rowSds
#' @export
computeRowZScore <- function(mat) {
    mat <- (mat - rowMeans2(mat, na.rm = TRUE)) /
        rowSds(mat, na.rm = TRUE)
    mat[is.nan(mat) | is.na(mat)] <- 0
    return(mat)
}


#' Parallel back-end from a thread count
#'
#' @param threads Thread count for parallel processing.
#' @return A \code{BiocParallelParam} object.
#' @importFrom BiocParallel SerialParam MulticoreParam SnowParam
#' @keywords internal
bpparamFromThreads <- function(threads) {
    if (is.null(threads) || !is.numeric(threads) || threads <= 1) {
        return(SerialParam())
    }
    if (.Platform$OS.type == "windows") {
        return(SnowParam(workers = threads))
    }
    return(MulticoreParam(workers = threads))
}
