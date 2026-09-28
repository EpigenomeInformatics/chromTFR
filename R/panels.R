#' Footprint data of one motif across groups
#'
#' @param tfbs A \code{GRanges} of binding sites, prepared with
#' \code{prepareTFBS}.
#' @param ins A named \code{list} of insertion site \code{GRanges}.
#' @param expected A named \code{list} of expected profiles with the
#' columns \code{x} and \code{w}, or a single profile used for all.
#' @param smooth Width of the running mean applied to the profiles.
#' @param flankNorm Width of the flanking window used for matching.
#' @param ignoreStrand If TRUE, strand information is ignored.
#' @return A \code{data.table} with the columns \code{x},
#' \code{avg_ins}, \code{type} and \code{group}, or NULL.
#' @examples
#' \dontrun{
#' accFootprintData(tfbs, insList, expectedList)
#' }
#' @importFrom logger log_warn
#' @import data.table
#' @export
accFootprintData <- function(
    tfbs, ins, expected, smooth = 5L, flankNorm = 50, ignoreStrand = TRUE
) {
    groups <- names(ins)
    df <- rbindlist(Filter(Negate(is.null), lapply(groups, function(g) {
        obs <- accProfile(ins[[g]], tfbs, ignoreStrand)
        if (is.null(obs)) {
            log_warn(g, ": no insertions in the binding sites")
            return(NULL)
        }
        e <- if (is.data.frame(expected)) expected else expected[[g]]
        e <- scaleToFlank(e, obs, flankNorm)
        obs[, type := "Observed"]
        e[, type := "Expected"]
        d <- rbindlist(list(smoothProfile(obs, smooth),
            smoothProfile(e, smooth)))
        d[, group := g]
        d
    })))
    if (nrow(df) == 0) {
        return(NULL)
    }
    df[, group := factor(group, levels = intersect(groups, unique(group)))]
    return(df[])
}


#' Theme shared by the footprint panels
#'
#' @param base_size Base font size.
#' @return A \code{ggplot2} theme.
#' @importFrom ggplot2 theme_classic theme element_text
#' @importFrom grid unit
#' @keywords internal
accPanelTheme <- function(base_size = 10) {
    theme_classic(base_size = base_size) +
        theme(
            plot.title = element_text(hjust = 0, face = "plain",
                size = base_size),
            legend.position = "bottom",
            legend.key.size = unit(3.5, "mm")
        )
}


#' Observed against expected footprint panel
#'
#' @param df Footprint data from \code{\link{accFootprintData}}.
#' @param title Panel title.
#' @param colors Named vector of group colours.
#' @param plotWindow Half width of the plotted window.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' accRawPanel(df, "CTCF", c(a = "red", b = "blue"))
#' }
#' @importFrom ggplot2 ggplot aes geom_line scale_colour_manual labs
#' @importFrom ggplot2 scale_linetype_manual coord_cartesian
#' @export
accRawPanel <- function(df, title, colors, plotWindow = 200L) {
    ggplot(df, aes(x = x, y = avg_ins, colour = group, linetype = type)) +
        geom_line(linewidth = 0.4) +
        scale_colour_manual(values = colors[levels(df$group)], name = NULL) +
        scale_linetype_manual(
            values = c("Observed" = "solid", "Expected" = "dashed"),
            name = NULL
        ) +
        coord_cartesian(xlim = c(-plotWindow, plotWindow)) +
        labs(
            title = title,
            x = "Distance from motif centre",
            y = "Insertions per bp and site"
        ) +
        accPanelTheme()
}


#' Bias corrected footprint panel
#'
#' Each group is labelled with its deviation score.
#'
#' @param df Footprint data from \code{\link{accFootprintData}}.
#' @param title Panel title.
#' @param colors Named vector of group colours.
#' @param method Either "substraction" or "division".
#' @param flankNorm Width of the flanking window used for normalisation.
#' @param plotWindow Half width of the plotted window.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' accDiffPanel(df, "CTCF", c(a = "red", b = "blue"))
#' }
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs
#' @importFrom ggplot2 scale_colour_manual coord_cartesian
#' @importFrom stats setNames
#' @importFrom logger log_info
#' @import data.table
#' @export
accDiffPanel <- function(
    df, title, colors, method = "division", flankNorm = 30L,
    plotWindow = 200L
) {
    d <- df[, .(avg_ins = if (method == "substraction") {
        avg_ins[type == "Observed"] - avg_ins[type == "Expected"]
    } else {
        avg_ins[type == "Observed"] / avg_ins[type == "Expected"]
    }), by = .(x, group)]

    flank <- max(abs(d$x), na.rm = TRUE)
    base <- if (method == "substraction") 0 else 1
    d[, avg_ins := if (method == "substraction") {
        avg_ins - mean(avg_ins[abs(x) >= flank - flankNorm], na.rm = TRUE)
    } else {
        avg_ins / mean(avg_ins[abs(x) >= flank - flankNorm], na.rm = TRUE)
    }, by = group]

    devs <- df[, .(dev = accDeviation(.SD)),
        by = group, .SDcols = c("x", "avg_ins", "type")
    ]
    log_info(title, " deviations: ", paste(devs$group, round(devs$dev, 3),
        sep = " = ", collapse = ", "
    ))

    present <- levels(droplevels(d$group))
    labels <- paste0(present, " (", format(
        round(devs$dev[match(present, devs$group)], 2), nsmall = 2
    ), ")")
    d[, group := factor(labels[match(as.character(group), present)],
        levels = labels
    )]
    cols <- setNames(unname(colors[present]), labels)

    ylab <- if (method == "substraction") {
        "Insertion difference (Observed - Expected)"
    } else {
        "Insertion ratio (Observed / Expected)"
    }
    ggplot(d, aes(x = x, y = avg_ins, colour = group)) +
        geom_hline(yintercept = base, linetype = "dotted",
            colour = "grey55") +
        geom_line(linewidth = 0.4) +
        scale_colour_manual(values = cols, name = NULL) +
        coord_cartesian(xlim = c(-plotWindow, plotWindow)) +
        labs(title = title, x = "Distance from motif centre", y = ylab) +
        accPanelTheme()
}


#' Grid of footprint panels
#'
#' One column per motif, the observed profiles on top and the corrected
#' ones below.
#'
#' @param profiles A named \code{list} of footprint data tables.
#' @param tfbs A named \code{list} of prepared binding sites, used for
#' the site counts in the titles.
#' @param colors Named vector of group colours.
#' @param method Either "substraction" or "division".
#' @param flankNorm Width of the flanking window used for normalisation.
#' @param plotWindow Half width of the plotted window.
#' @param corrected Label of the corrected row.
#' @return A \code{patchwork} object.
#' @examples
#' \dontrun{
#' plotAccFootprintGrid(profiles, tfbs, colors)
#' }
#' @importFrom patchwork wrap_plots
#' @export
plotAccFootprintGrid <- function(
    profiles, tfbs, colors, method = "division", flankNorm = 30L,
    plotWindow = 200L, corrected = "bias corrected"
) {
    motifs <- names(profiles)
    panels <- c(
        lapply(motifs, function(m) {
            accRawPanel(profiles[[m]],
                paste0(m, " (", length(tfbs[[m]]), " sites)"),
                colors, plotWindow
            )
        }),
        lapply(motifs, function(m) {
            accDiffPanel(profiles[[m]], paste0(m, ", ", corrected),
                colors, method, flankNorm, plotWindow
            )
        })
    )
    wrap_plots(panels, ncol = length(motifs))
}
