#' Add Ideogram Annotation to Coverage Plot.
#'
#' @param genome Single specie names, which must be one of the result from \code{ucscGenomes()$db}.
#' If missing, will invoke a menu for users to choose from. Default: hg19.
#' @param mark.color The color to mark plot region on ideogram. Default: "red".
#' @param mark.alpha The alpha to mark plot region on ideogram. Default: 0.7.
#' @param mark.line.size The line size to mark plot region on ideogram. Default: 1.
#' @param add.shadow Logical value, whether to add shadow polygon. Default: TRUE.
#' @param shadow.color The color to fill shadow polygon. Default: "grey".
#' @param shadow.alpha The alpha of shadow polygon. Default: 0.7.
#' @param shadow.line.size The line size of shadow polygon. Default: 1.
#' @param highlight.centromere Logical value, whether to highlight centromere region. Default: FALSE.
#' @param highlight.color The color to mark centromere region on ideogram. Default: "green".
#' @param highlight.alpha The alpha to mark centromere region on ideogram. Default: 0.7.
#' @param highlight.line.size The line size to mark centromere region on ideogram. Default: 1.
#' @param highlight.shadow.color The color to fill centromere shadow polygon. Default: "black".
#' @param highlight.shadow.alpha The alpha of centromere shadow polygon. Default: 0.7.
#' @param highlight.shadow.line.size The line size of centromere shadow polygon. Default: 1.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of ideogram annotation to coverage plot. Default: 0.2.
#'
#' @return Plot.
#' @importFrom magrittr %>%
#' @importFrom ggbio layout_karyogram
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom ggplot2 ggplot_add ggplot geom_rect aes_string geom_polygon theme_classic theme element_blank
#' scale_x_continuous scale_y_continuous margin
#' @importFrom patchwork wrap_plots
#' @importFrom methods extends
#' @importFrom utils menu
#' @importFrom GenomeInfoDb seqlengths seqlengths<- seqnames
#' @importFrom GenomicRanges trim GRanges
#' @importFrom S4Vectors values<-
#' @export
#'
geom_ideogram <- function(genome = "hg19", mark.color = "red", mark.alpha = 0.7, mark.line.size = 1,
                          add.shadow = TRUE, shadow.color = "grey", shadow.alpha = 0.7, shadow.line.size = 1,
                          highlight.centromere = FALSE, highlight.color = "green", highlight.alpha = 0.7, highlight.line.size = 1,
                          highlight.shadow.color = "black", highlight.shadow.alpha = 0.7, highlight.shadow.line.size = 1,
                          plot.space = 0.1, plot.height = 0.1) {
  structure(list(
    genome = genome, mark.color = mark.color, mark.alpha = mark.alpha, mark.line.size = mark.line.size,
    add.shadow = add.shadow, shadow.color = shadow.color, shadow.alpha = shadow.alpha, shadow.line.size = shadow.line.size,
    highlight.centromere = highlight.centromere, highlight.color = highlight.color, highlight.alpha = highlight.alpha,
    highlight.line.size = highlight.line.size, highlight.shadow.color = highlight.shadow.color, highlight.shadow.alpha = highlight.shadow.alpha,
    highlight.shadow.line.size = highlight.shadow.line.size, plot.space = plot.space, plot.height = plot.height
  ),
  class = "ideogram"
  )
}

#' @export
ggplot_add.ideogram <- function(object, plot, object_name) {
  # if (length(plot$layers) == 0) {
  #   # geom_base
  #   # get plot data
  #   plot.data <- plot[[1]]$layers[[1]]$data
  #   # prepare plot range
  #   plot.chr <- as.character(plot.data[1, "seqnames"])
  #   plot.region.start <- plot.data[1, "start"]
  #   plot.region.end <- plot.data[nrow(plot.data), "end"]
  # } else {
  #   # get plot data
  #   plot.data <- plot$layers[[1]]$data
  #   # prepare plot range
  #   plot.chr <- as.character(plot.data[1, "seqnames"])
  #   plot.region.start <- plot$coordinates$limits$x[1]
  #   plot.region.end <- plot$coordinates$limits$x[2]
  # }
  # get plot data, plot data should contain bins
  if ("patchwork" %in% class(plot)) {
    plot.data <- plot[[1]]$layers[[1]]$data
  } else {
    plot.data <- plot$layers[[1]]$data
  }
  plot.chr <- as.character(plot.data[1, "seqnames"])
  # plot.region.start <- plot.data[1, "start"]
  plot.region.start <- min(plot.data[, "start"])
  # plot.region.end <- plot.data[nrow(plot.data), "end"]
  plot.region.end <- max(plot.data[, "end"])
  plot.region.gr <- GenomicRanges::GRanges(
    seqnames = plot.chr,
    ranges = IRanges::IRanges(plot.region.start, plot.region.end)
  )
  # get parameters
  genome <- object$genome
  mark.color <- object$mark.color
  mark.alpha <- object$mark.alpha
  mark.line.size <- object$mark.line.size
  add.shadow <- object$add.shadow
  shadow.color <- object$shadow.color
  shadow.alpha <- object$shadow.alpha
  shadow.line.size <- object$shadow.line.size
  highlight.centromere <- object$highlight.centromere
  highlight.color <- object$highlight.color
  highlight.alpha <- object$highlight.alpha
  highlight.line.size <- object$highlight.line.size
  highlight.shadow.color <- object$highlight.shadow.color
  highlight.shadow.alpha <- object$highlight.shadow.alpha
  highlight.shadow.line.size <- object$highlight.shadow.line.size
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # get genome and chr ideogram
  genome.info <- suppressWarnings(getIdeogram(genome = genome, subchr = plot.chr, cytobands = TRUE))
  genome.info.df <- genome.info %>% as.data.frame()
  # get genome length
  genome.length <- genome.info.df[nrow(genome.info.df), "end"]

  # create basic plot
  ideogram.plot <- ggplot() +
    ggbio::layout_karyogram(genome.info, cytobands = TRUE, geom = NULL, rect.height = 10)

  # mark plot region
  zoom.df <- data.frame(
    xmin = plot.region.start,
    xmax = plot.region.end,
    ymin = 0 - 0.2,
    ymax = 10 + 0.2,
    seqnames = plot.chr
  )
  ideogram.plot <- ideogram.plot +
    geom_rect(
      data = zoom.df,
      aes_string(xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax"),
      color = mark.color, fill = mark.color, size = mark.line.size,
      alpha = mark.alpha
    )

  # highlight centromere
  centromere.df <- genome.info.df[genome.info.df$gieStain == "acen", ]
  if (highlight.centromere) {
    if (nrow(centromere.df) < 1) {
      warning("There is no centromere region!")
    } else {
      # check overlap
      centromere.gr <- GenomicRanges::makeGRangesFromDataFrame(df = centromere.df, keep.extra.columns = TRUE)
      centromere.valid.df <- as.data.frame(IRanges::subsetByOverlaps(x = centromere.gr, ranges = plot.region.gr))
      if (nrow(centromere.valid.df) >= 1) {
        # get centromere start site
        centromere.start <- centromere.df[1, "start"]
        if (centromere.start < plot.region.start) {
          centromere.start <- plot.region.start
        }
        # # get centromere end site
        centromere.end <- centromere.df[nrow(centromere.df), "end"]
        if (centromere.end > plot.region.end) {
          centromere.end <- plot.region.end
        }
        # create hightlight region
        highlight.region <- data.frame(
          xmin = centromere.start,
          xmax = centromere.end,
          ymin = 0 - 0.2,
          ymax = 10 + 0.2,
          seqnames = plot.chr
        )
        ideogram.plot <- ideogram.plot +
          geom_rect(
            data = highlight.region,
            aes_string(xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax"),
            color = highlight.color, fill = highlight.color, size = highlight.line.size,
            alpha = highlight.alpha
          )
        # add shadow for plot region
        centromere.relative.start <- ((centromere.start - plot.region.start) / (plot.region.end - plot.region.start)) * genome.length
        centromere.relative.end <- ((centromere.end - plot.region.start) / (plot.region.end - plot.region.start)) * genome.length
        centromere.shadow.df <- data.frame(
          x = c(centromere.start, centromere.relative.start, centromere.relative.end, centromere.end),
          y = c(10.2, 16, 16, 10.2)
        )
        ideogram.plot <- ideogram.plot +
          geom_polygon(
            data = centromere.shadow.df,
            aes_string(x = "x", y = "y"),
            color = highlight.shadow.color, fill = highlight.shadow.color, size = highlight.shadow.line.size,
            alpha = highlight.shadow.alpha
          )
      } else {
        warning("Centromere region is not in plot region!")
      }
    }
  }

  if (add.shadow) {
    # add shadow for plot region
    shadow.df <- data.frame(
      x = c(plot.region.start, 0, genome.length, plot.region.end),
      y = c(10.2, 16, 16, 10.2)
    )
    ideogram.plot <- ideogram.plot +
      geom_polygon(
        data = shadow.df,
        aes_string(x = "x", y = "y"),
        color = shadow.color, fill = shadow.color, size = shadow.line.size,
        alpha = shadow.alpha
      )
  }

  # add theme
  ideogram.plot <- ideogram.plot + theme_ideogram()
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    ideogram.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
