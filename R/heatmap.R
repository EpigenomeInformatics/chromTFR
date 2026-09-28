#' Observed deviation of every motif in every sample
#'
#' The insertions of one sample are read at a time, so the memory cost is
#' that of a single sample.
#'
#' @param files Paths of the per sample insertion site files.
#' @param tfbs A named \code{list} of prepared binding sites.
#' @param samples Sample identifiers, used as column names.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @return A \code{matrix} of motifs by samples.
#' @examples
#' \dontrun{
#' accObservedMatrix(files, tfbs, samples)
#' }
#' @importFrom logger log_info
#' @export
accObservedMatrix <- function(files, tfbs, samples, ignoreStrand = TRUE) {
    motifs <- names(tfbs)
    obs <- matrix(NA_real_,
        nrow = length(motifs), ncol = length(samples),
        dimnames = list(motifs, samples)
    )
    for (i in seq_along(samples)) {
        log_info("Processing ", samples[i])
        ins <- readRDS(files[i])
        for (m in motifs) {
            p <- accProfile(ins, tfbs[[m]], ignoreStrand)
            if (is.null(p)) next
            obs[m, i] <- accDeviationScore(p)
        }
        rm(ins)
        gc()
    }
    return(obs)
}


#' Expected deviation of every motif in every sample
#'
#' The sequence of each motif is read once and scored against the k-mer
#' weights of every sample, so the correction differs by sample and motif
#' and survives a row-wise z-scoring.
#'
#' @param tfbs A named \code{list} of prepared binding sites.
#' @param bias A named \code{list} of per sample k-mer weights.
#' @param genome A \code{BSgenome} object.
#' @param k Length of the k-mer.
#' @param max.sites Number of binding sites drawn per motif.
#' @return A \code{matrix} of motifs by samples.
#' @examples
#' \dontrun{
#' accExpectedMatrix(tfbs, bias, Hsapiens)
#' }
#' @importFrom logger log_info
#' @import data.table
#' @export
accExpectedMatrix <- function(tfbs, bias, genome, k = 6L, max.sites = 25000) {
    motifs <- names(tfbs)
    samples <- names(bias)
    kmers <- names(bias[[1]])
    expDev <- matrix(NA_real_,
        nrow = length(motifs), ncol = length(samples),
        dimnames = list(motifs, samples)
    )
    for (m in motifs) {
        log_info("Expected profile of ", m)
        prof <- kmerIndexProfile(tfbs[[m]], genome, kmers,
            k = k, max.sites = max.sites
        )
        expDev[m, ] <- vapply(samples, function(s) {
            accDeviationScore(profileFromBias(prof, bias[[s]]))
        }, numeric(1))
        rm(prof)
        gc()
    }
    return(expDev)
}


#' Heatmap of bias corrected deviation z-scores
#'
#' @param z A \code{matrix} of motifs by samples, row-wise z-scored.
#' @param ann Optional \code{data.frame} of sample annotation, one row
#' per column of \code{z}.
#' @param annCols Columns of \code{ann} drawn as annotation bars.
#' @param split Optional column of \code{ann} splitting the columns.
#' @param limits Range covered by the colour scale.
#' @param name Legend title.
#' @return A \code{Heatmap} object.
#' @examples
#' \dontrun{
#' plotAccDeviationHeatmap(z, ann, c("region"), split = "region")
#' }
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @export
plotAccDeviationHeatmap <- function(
    z, ann = NULL, annCols = character(), split = NULL, limits = c(-2, 2),
    name = "Bias corrected\ndeviation z-score"
) {
    colFun <- colorRamp2(
        seq(limits[1], limits[2], length.out = 5),
        c("#4575B4", "#91BFDB", "#FFFFBF", "#FC8D59", "#D73027")
    )
    annCols <- intersect(annCols, colnames(ann))
    topAnn <- NULL
    if (!is.null(ann) && length(annCols) > 0) {
        topAnn <- HeatmapAnnotation(
            df = ann[, annCols, drop = FALSE],
            annotation_name_side = "right"
        )
    }
    Heatmap(z,
        name = name,
        col = colFun,
        top_annotation = topAnn,
        column_split = if (!is.null(split) && split %in% colnames(ann)) {
            ann[[split]]
        } else {
            NULL
        },
        cluster_rows = nrow(z) > 2,
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 9),
        heatmap_legend_param = list(direction = "horizontal")
    )
}
