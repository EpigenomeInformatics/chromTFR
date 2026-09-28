#' Count the k-mers of the accessible genome
#'
#' The background depends on the regions only, so it is counted once and
#' handed to \code{\link{computeKmerBias}} for every sample.
#'
#' @param regions A \code{GRanges} of the accessible regions.
#' @param genome A \code{BSgenome} object.
#' @param k Length of the k-mer.
#' @param chunk Number of regions read from the genome at a time.
#' @return A named numeric vector of k-mer counts.
#' @examples
#' \dontrun{
#' bg <- kmerBackground(peaks, BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
#' }
#' @importFrom logger log_info
#' @export
kmerBackground <- function(regions, genome, k = 6L, chunk = 5000) {
    log_info("Counting background k-mers in ", length(regions), " regions")
    bg <- NULL
    rid <- split(seq_along(regions), ceiling(seq_along(regions) / chunk))
    for (ii in rid) {
        f <- Biostrings::oligonucleotideFrequency(
            Biostrings::getSeq(genome, regions[ii]),
            width = k, step = 1, simplify.as = "collapsed"
        )
        bg <- if (is.null(bg)) f else bg + f
    }
    return(bg)
}


#' Tn5 sequence preference of one sample
#'
#' The enrichment of every k-mer at the observed cut sites over its
#' frequency in the accessible genome. The k-mer is centred on the cut
#' and the weights are symmetrised over the two strands, because both
#' fragment ends are counted as insertions.
#'
#' @param ins A \code{GRanges} of insertion sites with a \code{score}
#' column.
#' @param regions A \code{GRanges} of the accessible regions.
#' @param genome A \code{BSgenome} object.
#' @param k Length of the k-mer.
#' @param max.ins Number of insertion sites drawn for the estimate.
#' @param chunk Number of sites read from the genome at a time.
#' @param bg Background counts from \code{\link{kmerBackground}}, or NULL
#' to count them here.
#' @return A named numeric vector of k-mer weights.
#' @examples
#' \dontrun{
#' bias <- computeKmerBias(ins, peaks, Hsapiens, bg = bg)
#' }
#' @importFrom GenomicRanges GRanges seqnames start
#' @importFrom IRanges IRanges
#' @importFrom logger log_info
#' @import data.table
#' @export
computeKmerBias <- function(
    ins, regions, genome, k = 6L, max.ins = 2e6, chunk = 1e5, bg = NULL
) {
    flank <- floor(k / 2)
    win <- GRanges(seqnames(ins), IRanges(start(ins) - flank, width = k))
    win$score <- ins$score
    win <- win[start(win) > 0]
    if (length(win) > max.ins) {
        win <- win[sample(length(win), max.ins)]
    }
    log_info("Counting k-mers at ", length(win), " insertion sites")

    obs <- NULL
    idx <- split(seq_along(win), ceiling(seq_along(win) / chunk))
    for (ii in idx) {
        km <- as.character(Biostrings::getSeq(genome, win[ii]))
        obs <- rbind(obs, data.table(kmer = km, n = win$score[ii]))[
            , .(n = sum(n)), by = kmer
        ]
    }
    if (is.null(bg)) {
        bg <- kmerBackground(regions, genome, k = k)
    }

    o <- obs$n[match(names(bg), obs$kmer)]
    o[is.na(o)] <- 0
    rc <- as.character(Biostrings::reverseComplement(
        Biostrings::DNAStringSet(names(bg))
    ))
    j <- match(rc, names(bg))
    o <- (o + o[j]) / 2
    b <- (bg + bg[j]) / 2

    w <- ((o + 1) / sum(o + 1)) / ((b + 1) / sum(b + 1))
    names(w) <- names(bg)
    log_info("k-mer weights from ", round(min(w), 2), " to ", round(max(w), 2))
    return(w)
}


#' Expected cut propensity along a set of binding sites
#'
#' @param tfbs A \code{GRanges} of binding sites, prepared with
#' \code{prepareTFBS}.
#' @param genome A \code{BSgenome} object.
#' @param bias Named k-mer weights from \code{\link{computeKmerBias}}.
#' @param k Length of the k-mer.
#' @param max.sites Number of binding sites drawn for the estimate.
#' @return A \code{data.table} with the columns \code{x} and \code{w}.
#' @examples
#' \dontrun{
#' kmerBiasProfile(tfbs, Hsapiens, bias)
#' }
#' @importFrom GenomicRanges width strand strand<-
#' @importFrom logger log_info
#' @import data.table
#' @export
kmerBiasProfile <- function(
    tfbs, genome, bias, k = 6L, max.sites = 25000
) {
    if (length(tfbs) > max.sites) {
        tfbs <- tfbs[sample(length(tfbs), max.sites)]
    }
    strand(tfbs) <- "*"
    seqs <- Biostrings::getSeq(genome, tfbs)
    w <- width(tfbs)[1]
    flank <- floor(k / 2)
    npos <- w - k + 1
    log_info("Bias profile over ", length(tfbs), " sites")

    val <- vapply(seq_len(npos), function(p) {
        km <- as.character(Biostrings::subseq(seqs, p, p + k - 1))
        mean(bias[km], na.rm = TRUE)
    }, numeric(1))

    x <- seq_len(npos) + flank - 1 - round((w - 1) / 2)
    return(data.table(x = x, w = val))
}


#' k-mer identities along a set of binding sites
#'
#' The sequence is read once and reused for every sample, so fitting the
#' Tn5 model per sample costs a table lookup rather than a second pass
#' over the genome.
#'
#' @param tfbs A \code{GRanges} of binding sites, prepared with
#' \code{prepareTFBS}.
#' @param genome A \code{BSgenome} object.
#' @param kmers The k-mer vocabulary, the names of a bias vector.
#' @param k Length of the k-mer.
#' @param max.sites Number of binding sites drawn for the estimate.
#' @return A \code{list} with the offsets \code{x} and the index matrix
#' \code{idx} of sites by positions.
#' @examples
#' \dontrun{
#' prof <- kmerIndexProfile(tfbs, Hsapiens, names(bias))
#' }
#' @importFrom GenomicRanges width strand strand<-
#' @importFrom logger log_info
#' @export
kmerIndexProfile <- function(
    tfbs, genome, kmers, k = 6L, max.sites = 25000
) {
    if (length(tfbs) > max.sites) {
        tfbs <- tfbs[sample(length(tfbs), max.sites)]
    }
    strand(tfbs) <- "*"
    seqs <- Biostrings::getSeq(genome, tfbs)
    w <- width(tfbs)[1]
    flank <- floor(k / 2)
    npos <- w - k + 1
    log_info("Reading ", npos, " positions of ", length(tfbs), " sites")

    idx <- vapply(seq_len(npos), function(p) {
        match(as.character(Biostrings::subseq(seqs, p, p + k - 1)), kmers)
    }, integer(length(seqs)))

    x <- seq_len(npos) + flank - 1 - round((w - 1) / 2)
    return(list(x = x, idx = idx))
}


#' Expected profile of one sample from stored k-mer identities
#'
#' @param prof Output of \code{\link{kmerIndexProfile}}.
#' @param bias Named k-mer weights of the sample.
#' @return A \code{data.table} with the columns \code{x} and
#' \code{avg_ins}.
#' @examples
#' \dontrun{
#' profileFromBias(prof, bias[["sample_1"]])
#' }
#' @import data.table
#' @export
profileFromBias <- function(prof, bias) {
    w <- matrix(bias[prof$idx], nrow = nrow(prof$idx))
    return(data.table(x = prof$x, avg_ins = colMeans(w, na.rm = TRUE)))
}


#' Match an expected profile to the flanks of an observed one
#'
#' Only the plotted curves depend on this. The deviation is a ratio of
#' the centre over the flanks, so it is unchanged by the scale.
#'
#' @param exp_profile A \code{data.table} with the columns \code{x} and
#' \code{w}.
#' @param obs_profile A \code{data.table} with the columns \code{x} and
#' \code{avg_ins}.
#' @param flankNorm Width of the flanking window used for matching.
#' @return A \code{data.table} with the columns \code{x} and
#' \code{avg_ins}.
#' @examples
#' \dontrun{
#' scaleToFlank(expProfile, obsProfile)
#' }
#' @import data.table
#' @export
scaleToFlank <- function(exp_profile, obs_profile, flankNorm = 50) {
    exp_profile <- copy(exp_profile)
    if (!"w" %in% names(exp_profile)) {
        exp_profile[, w := avg_ins]
    }
    fo <- max(abs(obs_profile$x), na.rm = TRUE)
    fe <- max(abs(exp_profile$x), na.rm = TRUE)
    mo <- mean(obs_profile$avg_ins[abs(obs_profile$x) >= fo - flankNorm],
        na.rm = TRUE
    )
    me <- mean(exp_profile$w[abs(exp_profile$x) >= fe - flankNorm],
        na.rm = TRUE
    )
    exp_profile[, avg_ins := w * (mo / me)]
    return(exp_profile[, .(x, avg_ins)])
}
