#' Add Peak Annotation to Coverage Plot.
#'
#' @param bed.file The path to consensus peaks file. Default: NULL.
#' @param peak.df The dataframe contains consensus peaks. Default: NULL.
#' @param peak.color Peak colors. Default: NULL.
#' @param peak.size The line size of peak. Default: 5.
#' @param color.by Name of optional column in bed file/data frame which is used for coloring. Default: NULL.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of peak annotation to coverage plot. Default: 0.1.
#'
#' @return Plot.
#' @importFrom utils read.table
#' @importFrom dplyr arrange
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot_add ggplot geom_segment aes_string theme_classic
#'   theme element_blank element_text element_rect margin scale_x_continuous
#'   scale_y_continuous scale_color_continuous coord_cartesian
#' @export
#'
#' @examples
#' # load metadata
#' sample_meta <- data.frame(
#'   SampleName = c(
#'     "Chr18_MCF7_ER_1",
#'     "Chr18_MCF7_ER_2",
#'     "Chr18_MCF7_ER_3",
#'     "Chr18_MCF7_input"
#'   ),
#'   Type = c("MCF7_ER_1", "MCF7_ER_2", "MCF7_ER_3", "MCF7_input"),
#'   Group = c("IP", "IP", "IP", "Input")
#' )
#'
#' # import coverage track
#' track_folder <- system.file("extdata", "ChIP-seq", package = "ggcoverage")
#' track_df <- LoadTrackFile(
#'   track.folder = track_folder,
#'   format = "bw",
#'   region = "chr18:76822285-76900000",
#'   meta.info = sample_meta
#' )
#'
#' # create mock peak file
#' df_peaks <- data.frame(
#'   seqnames = c("chr18", "chr18", "chr18"),
#'   start = c(76822533, 76846900, 76880000),
#'   end = c(76836900, 76860000, 76887000),
#'   score = c(4, 6, 13)
#' )
#'
#' # plot with default color
#' ggcoverage(data = track_df) +
#'   geom_peak(peak.df = df_peaks, peak.size = 3)
#'
#' # plot with color by 'score' variable
#' ggcoverage(data = track_df) +
#'   geom_peak(peak.df = df_peaks, peak.size = 3, color.by = "score")
#'
#' # plot with color by 'score' variable and custom color scale
#' ggcoverage(data = track_df) +
#'   geom_peak(peak.df = df_peaks, peak.size = 3, color.by = "score", peak.color = rainbow(5))
#'
geom_peak <- function(bed.file = NULL, peak.df = NULL, peak.color = NULL, peak.size = 5,
                      color.by = NULL, plot.space = 0.1, plot.height = 0.1) {
  structure(
    list(
      bed.file = bed.file, peak.df = peak.df, peak.color = peak.color, peak.size = peak.size,
      color.by = color.by, plot.space = plot.space, plot.height = plot.height
    ),
    class = "peak"
  )
}

#' @export
ggplot_add.peak <- function(object, plot, object_name) {
  if ("patchwork" %in% class(plot)) {
    plot.data <- plot[[1]]$layers[[1]]$data
  } else {
    plot.data <- plot$layers[[1]]$data
  }
  plot.chr <- as.character(plot.data[1, "seqnames"])
  plot.region.start <- min(plot.data[, "start"])
  plot.region.end <- max(plot.data[, "end"])

  # get parameters
  bed.file <- object$bed.file
  peak.df <- object$peak.df
  peak.color <- object$peak.color
  peak.size <- object$peak.size
  color.by <- object$color.by
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # prepare peak dataframe
  if (!is.null(bed.file)) {
    bed.info <- utils::read.table(file = bed.file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  } else if (!is.null(peak.df)) {
    bed.info <- peak.df
  }
  colnames(bed.info)[c(1, 2, 3)] <- c("seqnames", "start", "end")
  bed.info$start <- as.numeric(bed.info$start) + 1
  valid.bed <- GetRegion(chr = plot.chr, df = bed.info, start = plot.region.start, end = plot.region.end)

  # color management
  if (!is.null(color.by)) {
    if (is.numeric(valid.bed[[color.by]])) {
      if (!is.null(peak.color)) {
        scale_colors <- scale_color_gradientn(colours = peak.color)
      } else {
        scale_colors <- scale_color_continuous()
      }
    } else {
      if (length(peak.color) < length(unique(valid.bed[[color.by]]))) {
        warning("Fewer colors provided than there are groups in ", color.by, " variable, falling back to default colors")
        auto_colors <- AutoColor(data = valid.bed[[color.by]], pal = "Set3")
        scale_colors <- scale_color_manual(values = auto_colors)
      } else {
        scale_colors <- scale_color_manual(values = peak.color)
      }
    }
  } else {
    color.by <- "color"
    if (is.null(peak.color)) {
      peak.color <- "black"
    }
    valid.bed$color <- peak.color[1]
    scale_colors <- scale_color_manual(values = peak.color)
  }

  peak.plot <- ggplot() +
    geom_segment(
      data = valid.bed,
      mapping = aes_string(
        x = "start",
        y = "1",
        xend = "end",
        yend = "1",
        color = color.by
      ),
      size = peak.size,
    ) +
    scale_colors +
    labs(y = "Peak")

  # add theme
  peak.plot <- peak.plot +
    theme_peak(margin.len = plot.space, x.range = c(plot.region.start, plot.region.end)) +
    theme(legend.position = "none")
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    peak.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
