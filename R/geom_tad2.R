#' Add Contact Map to Coverage Plot.
#'
#' @param matrix Matrix (n x n) contains contact map information.
#' @param granges The rownames and colnames information of matrix.
#' @param zlim The maximum value of color to plot. Larger values will be truncated to this value.
#' @param color.bias Bias parameter for color palette, it is the same as \code{bias} parameter of \code{\link{colorRampPalette}}. Default: 1.
#' @param color.ramp A name of a color palette, choose from rownames of \code{\link{brewer.pal.info}}. Default: YlOrRd.
#' @param color.palette A vector of colors. Overwrites the palette for \code{color.ramp}. Default: NULL.
#' @param smooth.func An optional smoothing function for the matrix. Default: NULL.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of contact map to coverage plot. Default: 1.
#' @param top Logical value, whether to place this plot on the coverage plot. Default: TRUE.
#' @param show.rect Logical value, whether to add rect border to the plot. Default: FALSE.
#'
#' @return Plot.
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeMatrix matrixPlotter
#' @importFrom ggplot2 ggplot_add ggplot labs theme_classic theme element_blank element_rect
#' element_text margin scale_y_continuous scale_x_continuous coord_cartesian
#' @importFrom patchwork wrap_plots
#'
#'
#' @examples
#' library(ggcoverage)
#' library(HiCDataHumanIMR90)
#' data(Dixon2012_IMR90, package = "HiCDataHumanIMR90")
#' mat <- as.matrix(hic_imr90_40@.Data[[1]]@intdata)
#' granges <- hic_imr90_40@.Data[[1]]@xgi
#' # prepare coverage dataframe
#' df <- data.frame(
#'   seqnames = "chr1", start = seq(from = 50000000, to = 59999000, by = 1000),
#'   end = seq(from = 50001000, to = 60000000, by = 1000), score = sample(1:100, 10000, replace = TRUE),
#'   Type = "Example", Group = "Example"
#' )
#' # create plot
#' ggcoverage(
#'   data = df, color = "grey", region = "chr1:50000000-56000000",
#'   mark.region = NULL, range.position = "out"
#' ) +
#'   geom_tad2(matrix = log2(mat + 1), granges = granges, zlim = 5, color.palette = c("blue", "red"))
geom_tad2 <- function(matrix, granges, zlim = NULL, color.bias = 1, color.ramp = "YlOrRd",
                      color.palette = NULL, smooth.func = NULL, plot.space = 0.1, plot.height = 1, top = TRUE, show.rect = FALSE) {
  structure(list(
    matrix = matrix, granges = granges, zlim = zlim, color.bias = color.bias,
    color.ramp = color.ramp, color.palette = color.palette, smooth.func = smooth.func,
    plot.space = plot.space, plot.height = plot.height, top = top, show.rect = show.rect
  ),
  class = "tad2"
  )
}

#' @export
ggplot_add.tad2 <- function(object, plot, object_name) {
  # get plot data
  # track.data <- plot$layers[[1]]$data
  # get plot data, plot data should contain bins
  if ("patchwork" %in% class(plot)) {
    track.data <- plot[[1]]$layers[[1]]$data
  } else {
    track.data <- plot$layers[[1]]$data
  }
  # prepare plot range
  # the plot region are not normal, so start is minimum value
  plot.range.chr <- track.data[1, "seqnames"]
  # plot.range.start <- track.data[1, "start"]
  plot.range.start <- min(track.data[, "start"])
  # plot.range.end <- track.data[nrow(track.data), "end"]
  plot.range.end <- max(track.data[, "end"])
  plot.range.gr <- GenomicRanges::GRanges(
    seqnames = plot.range.chr,
    ranges = IRanges::IRanges(plot.range.start, plot.range.end)
  )
  # get parameters
  matrix <- object$matrix
  granges <- object$granges
  zlim <- object$zlim
  color.bias <- object$color.bias
  color.ramp <- object$color.ramp
  color.palette <- object$color.palette
  smooth.func <- object$smooth.func
  plot.space <- object$plot.space
  plot.height <- object$plot.height
  top <- object$top
  show.rect <- object$show.rect

  # filter range and mat
  used.range.index <- as.data.frame(IRanges::findOverlaps(granges, plot.range.gr))[, "queryHits"]
  matrix.used <- matrix[used.range.index, used.range.index]
  used.range.gr <- IRanges::subsetByOverlaps(x = granges, ranges = plot.range.gr)

  # plot
  tad.plot <- GenomeMatrix::matrixPlotter(
    mat = matrix.used, granges = used.range.gr, plotGR = plot.range.gr,
    extend = 0, heightProp = 1 / 2, zlim = zlim, colorBias = color.bias,
    highlight = NULL, colRamp = color.ramp, colPalette = color.palette,
    smoothFilt = smooth.func
  )
  # add theme
  tad.final.plot <- tad.plot +
    labs(y = "ContactMap") +
    theme_tad2(x.range = c(plot.range.start, plot.range.end), margin.len = plot.space, show.rect = show.rect)
  # assemble plot
  if (top) {
    patchwork::wrap_plots(tad.final.plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
      plot,
      ncol = 1, heights = c(plot.height, 1)
    )
  } else {
    patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
      tad.final.plot,
      ncol = 1, heights = c(1, plot.height)
    )
  }
}
