#' Running mean of an insertion profile
#'
#' @param data A \code{data.table} with the columns \code{x} and
#' \code{avg_ins}.
#' @param window Width of the running mean, 1 to switch it off.
#' @return A \code{data.table}.
#' @importFrom stats filter na.omit
#' @import data.table
#' @keywords internal
smoothProfile <- function(data, window = 5L) {
    if (is.null(window) || window <= 1) {
        return(data)
    }
    data <- data[order(x)]
    data[, avg_ins := as.numeric(stats::filter(
        avg_ins, rep(1 / window, window), sides = 2
    ))]
    return(na.omit(data))
}


#' Deviation of an observed and expected profile pair
#'
#' @param plot_data A \code{data.table} with the columns \code{x},
#' \code{avg_ins} and \code{type}.
#' @return A numeric value.
#' @examples
#' \dontrun{
#' accDeviation(plotData)
#' }
#' @import data.table
#' @export
accDeviation <- function(plot_data) {
    obs <- accDeviationScore(plot_data[type == "Observed", .(x, avg_ins)])
    exp <- accDeviationScore(plot_data[type == "Expected", .(x, avg_ins)])
    return(obs - exp)
}


#' Observed and expected profile of one motif in one sample
#'
#' @param motif Motif name.
#' @param tf_bindsites A \code{GRangesList} of binding site positions.
#' @param ins A \code{GRanges} of insertion sites.
#' @param sample_name Optional sample label.
#' @param gc_dist A \code{GRanges} of the genome wide GC distribution.
#' @param gcfreqs A \code{list} of GC bin frequency tables.
#' @param enhancer Optional \code{GRanges} restricting the sites.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @param smooth Width of the running mean applied to the profiles.
#' @param returnPlotData If TRUE, the plot data is returned as well.
#' @return A \code{ggplot} object, or a list with the plot and its data.
#' @examples
#' \dontrun{
#' plotExpectedAccFootprint("CTCF", tf_bindsites, ins,
#'     gc_dist = gc_dist, gcfreqs = gcfreqs)
#' }
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab ggtitle theme
#' @importFrom ggplot2 theme_classic scale_color_manual
#' @importFrom IRanges subsetByOverlaps
#' @import data.table
#' @export
plotExpectedAccFootprint <- function(
    motif, tf_bindsites, ins, sample_name = NULL, gc_dist, gcfreqs,
    enhancer = NULL, ignoreStrand = TRUE, smooth = 5L,
    returnPlotData = FALSE
) {
    if (is.null(sample_name)) {
        sample_name <- "sample"
    }
    tfbs <- prepareTFBS(tf_bindsites[[motif]], enhancer, ignoreStrand)
    obs <- accProfile(ins, tfbs, ignoreStrand)
    if (is.null(obs)) {
        stop("No insertion sites found in the ", motif, " binding sites")
    }
    gcd <- if (is.null(enhancer)) {
        gc_dist
    } else {
        subsetByOverlaps(gc_dist, enhancer, ignore.strand = ignoreStrand)
    }
    binIns <- addGCBintoAccessome(ins, gcd, ignoreStrand)
    expProf <- computeAccExpectations(binIns, gcfreqs[[motif]])

    obs[, type := "Observed"]
    expProf[, type := "Expected"]
    plot_data <- rbindlist(list(
        smoothProfile(obs, smooth), smoothProfile(expProf, smooth)
    ))

    p1 <- ggplot(plot_data, aes(x = x, y = avg_ins, color = type)) +
        geom_line() +
        xlab("Distance from motif centre") +
        ylab("Insertions per bp and site") +
        theme_classic() +
        ggtitle(paste("TF footprint for", motif, "in", sample_name)) +
        scale_color_manual(values = c(
            "Expected" = "blue", "Observed" = "red"
        )) +
        theme(legend.position = "bottom")

    if (returnPlotData) {
        return(list(plot = p1, plotDF = plot_data))
    }
    return(p1)
}


#' Combine the observed and expected profiles
#'
#' @param plot_data A \code{data.table} with the columns \code{x},
#' \code{avg_ins} and \code{type}.
#' @param method Either "substraction" or "division".
#' @param sample_name Sample label used in the curve label.
#' @return A \code{list} with the corrected data and the axis label.
#' @import data.table
#' @keywords internal
accFootprintDifference <- function(plot_data, method, sample_name) {
    if (method == "substraction") {
        d <- plot_data[, .(avg_ins = avg_ins[type == "Observed"] -
            avg_ins[type == "Expected"]), by = x]
        d[, type := paste("Obs. sub. Exp.", sample_name)]
        lab <- "(Observed - Expected)"
    } else {
        d <- plot_data[, .(avg_ins = avg_ins[type == "Observed"] /
            avg_ins[type == "Expected"]), by = x]
        d[, type := paste("Obs. div. Exp.", sample_name)]
        lab <- "(Observed / Expected)"
    }
    return(list(data = d, lab = lab))
}


#' Set the outer flanks of a corrected profile to the baseline
#'
#' @param data The corrected profile as a \code{data.table}.
#' @param method Either "substraction" or "division".
#' @param flankNorm Width of the flanking window used for normalisation.
#' @return A \code{data.table}.
#' @import data.table
#' @keywords internal
normaliseAccFlank <- function(data, method, flankNorm) {
    if (is.null(flankNorm) || flankNorm <= 0) {
        return(data)
    }
    flank <- max(abs(data$x), na.rm = TRUE)
    idx <- abs(data$x) >= flank - flankNorm
    normFactor <- mean(data$avg_ins[idx], na.rm = TRUE)
    if (method == "substraction") {
        data[, avg_ins := avg_ins - normFactor]
    } else {
        data[, avg_ins := avg_ins / normFactor]
    }
    return(data)
}


#' Bias corrected footprint of one motif
#'
#' @param motif Motif name.
#' @param tf_bindsites A \code{GRangesList} of binding site positions.
#' @param ins A \code{GRanges} of insertion sites.
#' @param sample_name Optional sample label.
#' @param gc_dist A \code{GRanges} of the genome wide GC distribution.
#' @param gcfreqs A \code{list} of GC bin frequency tables.
#' @param enhancer Optional \code{GRanges} restricting the sites.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @param method Either "substraction" or "division".
#' @param flankNorm Width of the flanking window used for normalisation.
#' @param smooth Width of the running mean applied to the profiles.
#' @param plotWindow Half width of the plotted window.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' plotAccMotifFootprint("CTCF", tf_bindsites, ins,
#'     gc_dist = gc_dist, gcfreqs = gcfreqs)
#' }
#' @importFrom ggplot2 ggplot aes geom_line geom_hline xlab ylab ggtitle
#' @importFrom ggplot2 theme_classic theme coord_cartesian
#' @import data.table
#' @export
plotAccMotifFootprint <- function(
    motif, tf_bindsites, ins, sample_name = NULL, gc_dist, gcfreqs,
    enhancer = NULL, ignoreStrand = TRUE, method = "division",
    flankNorm = 50, smooth = 5L, plotWindow = 200L
) {
    if (is.null(method) || !method %in% c("substraction", "division")) {
        method <- "division"
        warning("method is not provided, using the default division")
    }
    plot_data <- plotExpectedAccFootprint(
        motif, tf_bindsites, ins,
        sample_name = sample_name, gc_dist = gc_dist, gcfreqs = gcfreqs,
        enhancer = enhancer, ignoreStrand = ignoreStrand, smooth = smooth,
        returnPlotData = TRUE
    )$plotDF

    dev <- accDeviation(plot_data)
    diff <- accFootprintDifference(plot_data, method, sample_name)
    data <- normaliseAccFlank(diff$data, method, flankNorm)

    ggplot(data, aes(x = x, y = avg_ins, color = type)) +
        geom_hline(
            yintercept = if (method == "substraction") 0 else 1,
            linetype = "dotted", colour = "grey55"
        ) +
        geom_line() +
        xlab("Distance from motif centre") +
        ylab(paste0("Insertion difference ", diff$lab)) +
        theme_classic() +
        ggtitle(paste0(
            "TF footprint difference for ", motif, " in ", sample_name,
            " (dev = ", format(round(dev, 2), nsmall = 2), ")"
        )) +
        theme(legend.position = "bottom") +
        coord_cartesian(xlim = c(-plotWindow, plotWindow))
}
