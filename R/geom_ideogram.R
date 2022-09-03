#' Add Ideogram Annotation to Coverage Plot.
#'
#' @param genome Single specie names, which must be one of the result from \code{ucscGenomes()$db}.
#' If missing, will invoke a menu for users to choose from. Default: hg19.
#' @param mark.color The color to mark plot region on ideogram. Default: "red".
#' @param mark.alpha The alpha to mark plot region on ideogram. Default: 0.7.
#' @param mark.line.size The line size to mark plot region on ideogram. Default: 1.
#' @param add.shadow Logical value, whether to add shadow polygon. Default: TRUE.
#' @param shadow.color The color to fill shadow polygon. Default: grey.
#' @param shadow.alpha The alpha of shadow polygon. Default: 0.7.
#' @param shadow.line.size The line size of shadow polygon. Default: 1.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of ideogram annotation to coverage plot. Default: 0.2.
#'
#' @return Plot.
#' @importFrom magrittr %>%
#' @importFrom ggbio layout_karyogram
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
                          plot.space = 0.1, plot.height = 0.1) {
  structure(list(
    genome = genome, mark.color = mark.color, mark.alpha = mark.alpha, mark.line.size = mark.line.size,
    add.shadow = add.shadow, shadow.color = shadow.color, shadow.alpha = shadow.alpha,
    shadow.line.size = shadow.line.size, plot.space = plot.space, plot.height = plot.height
  ),
  class = "ideogram"
  )
}

#' @export
ggplot_add.ideogram <- function(object, plot, object_name) {
  if (length(plot$layers) == 0) {
    # geom_base
    # get plot data
    plot.data <- plot[[1]]$layers[[1]]$data
    # prepare plot range
    plot.chr <- as.character(plot.data[1, "seqnames"])
    plot.region.start <- plot.data[1, "start"]
    plot.region.end <- plot.data[nrow(plot.data), "end"]
  } else {
    # get plot data
    plot.data <- plot$layers[[1]]$data
    # prepare plot range
    plot.chr <- as.character(plot.data[1, "seqnames"])
    plot.region.start <- plot$coordinates$limits$x[1]
    plot.region.end <- plot$coordinates$limits$x[2]
  }

  # get parameters
  genome <- object$genome
  mark.color <- object$mark.color
  mark.alpha <- object$mark.alpha
  mark.line.size <- object$mark.line.size
  add.shadow <- object$add.shadow
  shadow.color <- object$shadow.color
  shadow.alpha <- object$shadow.alpha
  shadow.line.size <- object$shadow.line.size
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # get genome and chr ideogram
  genome.info <- suppressWarnings(getIdeogram(genome = genome, subchr = plot.chr, cytobands = TRUE))
  genome.info.df <- genome.info %>% as.data.frame()
  # get genome length
  genome.length <- genome.info.df[nrow(genome.info.df), "end"]

  # create basci plot
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
  if (add.shadow) {
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
