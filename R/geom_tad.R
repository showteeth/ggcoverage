#' Add Contact Map to Coverage Plot.
#'
#' @param matrix Matrix (n x n) contains contact map information.
#' @param granges The rownames and colnames information of matrix.
#' @param color.palette One of the RColorbrewer or viridis colour palettes.
#' Parameter of \code{\link{Brick_vizart_plot_heatmap}}. Default: NULL.
#' @param value.cut If present, values beyond a certain quantile will be capped to that quantile.
#' Parameter of \code{\link{Brick_vizart_plot_heatmap}}. Default: NULL.
#' @param transform.fun If any sort of transformations should be applied to the data before plotting.
#' Parameter of \code{\link{Brick_vizart_plot_heatmap}}. Default: NULL.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of contact map to coverage plot. Default: 1.
#' @param top Logical value, whether to place this plot on the coverage plot. Default: TRUE.
#' @param show.rect Logical value, whether to add rect border to the plot. Default: FALSE.
#'
#' @return Plot.
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges findOverlaps subsetByOverlaps
#' @import HiCBricks
#' @importFrom ggplot2 ggplot_add ggplot labs theme_classic theme element_blank element_rect
#' element_text margin scale_y_continuous scale_x_continuous
#' @importFrom patchwork wrap_plots
#' @importFrom utils write.table
#'
#' @examples
#' library(ggcoverage)
#' library(GenomicRanges)
#'
#' # prepare track dataframe
#' track.file = system.file("extdata", "HiC", "H3K36me3.bw", package = "ggcoverage")
#' track.df = LoadTrackFile(track.file = track.file, format = "bw",
#'                          region = "chr2L:8050000-8300000", extend = 0)
#' track.df$score = ifelse(track.df$score <0, 0, track.df$score)
#' # check the data
#' head(track.df)
#'
#' # Load Hi-C data
#' hic.mat.file = system.file("extdata", "HiC", "HiC_mat.txt", package = "ggcoverage")
#' hic.mat = read.table(file = hic.mat.file, sep = "\t")
#' hic.mat = as.matrix(hic.mat)
#'
#' # bin data
#' hic.bin.file = system.file("extdata", "HiC", "HiC_bin.txt", package = "ggcoverage")
#' hic.bin = read.table(file = hic.bin.file, sep = "\t")
#' colnames(hic.bin) = c("chr", "start", "end")
#' hic.bin.gr = GenomicRanges::makeGRangesFromDataFrame(df = hic.bin)
#'
#' # transfrom function
#' FailSafe_log10 <- function(x){
#'   x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
#'   return(log10(x+1))
#' }
#'
#' # load link data: prepare arcs
#' link.file = system.file("extdata", "HiC", "HiC_link.bedpe", package = "ggcoverage")
#'
#' # basic coverage
#' basic.coverage = ggcoverage(
#'   data = track.df, color = "grey",
#'   mark.region = NULL, range.position = "out"
#' )
#'
#' # add annotations
#' basic.coverage +
#'   geom_tad(matrix = hic.mat, granges = hic.bin.gr, value.cut = 0.99,
#'            color.palette = "viridis", transform.fun = FailSafe_log10,
#'            top = FALSE, show.rect = TRUE) +
#'   geom_link(link.file = link.file, file.type = "bedpe", show.rect = TRUE)
#'
#' @export
geom_tad <- function(matrix, granges, color.palette = NULL, value.cut = NULL,
                     transform.fun = NULL, plot.space = 0.1, plot.height = 1,
                     top = TRUE, show.rect = FALSE) {
  structure(list(
    matrix = matrix, granges = granges, color.palette = color.palette,
    value.cut = value.cut, transform.fun = transform.fun,
    plot.space = plot.space, plot.height = plot.height, top = top, show.rect = show.rect
  ),
  class = "tad"
  )
}
#' @export
ggplot_add.tad <- function(object, plot, object_name) {
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
  color.palette <- object$color.palette
  value.cut <- object$value.cut
  transform.fun <- object$transform.fun
  plot.space <- object$plot.space
  plot.height <- object$plot.height
  top <- object$top
  show.rect <- object$show.rect

  # create tmp directory
  out.dir <- file.path(tempdir(), "ggcoverage")
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }

  # filter used matrix and bin
  used.range.index <- as.data.frame(IRanges::findOverlaps(granges, plot.range.gr))[, "queryHits"]
  matrix.used <- matrix[used.range.index, used.range.index]
  used.range.df <- as.data.frame(IRanges::subsetByOverlaps(x = granges, ranges = plot.range.gr))

  # write used mat and bin to file
  write.table(
    x = used.range.df[c(1, 2, 3)], file = file.path(out.dir, "used.range.txt"),
    sep = "\t", quote = F, row.names = F, col.names = F
  )
  write.table(
    x = as.matrix(matrix.used), file = file.path(out.dir, "used.matrix.txt"),
    sep = "\t", quote = F, row.names = F, col.names = F
  )

  # create brick container
  BrickContainer <- HiCBricks::Create_many_Bricks(
    BinTable = file.path(out.dir, "used.range.txt"),
    bin_delim = "\t", output_directory = out.dir, file_prefix = "ggcoverage",
    experiment_name = "ggcoverage", resolution = 100000,
    remove_existing = TRUE, impose_discontinuity = FALSE
  )

  # load matrix
  Load.matrix <- HiCBricks::Brick_load_matrix(
    Brick = BrickContainer, chr1 = as.character(plot.range.chr),
    chr2 = as.character(plot.range.chr),
    matrix_file = file.path(out.dir, "used.matrix.txt"), delim = "\t",
    remove_prior = TRUE, resolution = 100000
  )

  # create plot
  tad.plot <-
    HiCBricks::Brick_vizart_plot_heatmap(
      File = file.path(out.dir, "ggcoverage_hic_plot.pdf"),
      Bricks = list(BrickContainer),
      x_coords = paste(plot.range.chr, plot.range.start, plot.range.end, sep = ":"),
      y_coords = paste(plot.range.chr, plot.range.start, plot.range.end, sep = ":"),
      resolution = 100000, FUN = transform.fun, value_cap = value.cut,
      distance = NULL, legend_title = "", palette = color.palette,
      width = 10, height = 10, rotate = TRUE, return_object = TRUE
    )
  # add theme
  tad.final.plot <- tad.plot +
    labs(y = "ContactMap") +
    theme_tad(margin.len = plot.space, show.rect = show.rect)
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
