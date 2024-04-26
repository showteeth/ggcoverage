#' Create Coverage Plot.
#'
#' @param data Coverage dataframe loaded by \code{\link{LoadTrackFile}}.
#' @param single.nuc Logical value, whether to visualize at single nucleotide level. Default: FALSE.
#' @param mapping Set of aesthetic mappings created by \code{aes} or \code{aes_}. Default: NULL.
#' @param color Track color. Default: NULL (select automatically).
#' @param rect.color The color of every bin. Default: NA.
#' @param plot.type The type of the plot, choose from facet (separate plot for every sample) and
#' joint (combine all sample in a single plot). Default: facet.
#' @param facet.key Sample type key to create coverage plot. Default: Type.
#' @param joint.avg Logical value, whether to show average coverage across \code{group.key}. Default: FALSE.
#' @param facet.order The order of Coverage plot. Default: NULL.
#' @param facet.color The color of sample text. Default: NULL (select automatically).
#' @param facet.y.scale The shared type of y-axis scales across facets, choose from free (facets have different y-axis scales),
#' fixed (facets have same y-axis scales). Default: free.
#' @param group.key Group of samples. Default: NULL.
#' @param range.size The label size of range text, used when \code{range.position} is in. Default: 3.
#' @param range.position The position of y axis range, chosen from in (move y axis in the plot) and
#' out (normal y axis). Default: in.
#' @param plot.space The space between every facet. Default: 0.2.
#' @param mark.region Mark region on the plot. Default: NULL.
#' @param mark.color The color of marked region. Default: "grey".
#' @param mark.alpha The alpha of marked region. Default: 0.5.
#' @param show.mark.label Logical value, whether to show mark label (use label column in \code{mark.region}). Default: TRUE.
#' @param mark.label.size The label size of mark label. Default: 4.
#'
#' @return A ggplot2 object.
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual geom_rect geom_text aes theme_classic theme unit
#' element_blank annotate rel scale_y_continuous expansion scale_x_continuous coord_cartesian geom_step
#' @importFrom scales comma
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang as_label .data
#' @importFrom stats as.formula
#' @importFrom ggh4x facet_wrap2 strip_themed
#' @importFrom dplyr group_by summarise
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' library(ggcoverage)
#' library(rtracklayer)
#' library(ggplot2)
#'
#' # import track data
#' meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
#' sample.meta <- read.csv(meta.file)
#' track.folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
#'
#' track.df <- LoadTrackFile(
#'   track.folder = track.folder, format = "bw",
#'   region = "chr14:21,677,306-21,737,601",
#'   extend = 2000, meta.info = sample.meta
#' )
#'
#' gtf.file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
#' gtf.gr <- rtracklayer::import.gff(con = gtf.file, format = "gtf")
#'
#' # plot tracks with coloring by 'Group' variable
#' ggcoverage(data = track.df, facet.key = "Type", group.key = "Group")
#'
#' # plot tracks without coloring by any group
#' ggcoverage(data = track.df, facet.key = "Type", group.key = NULL)
#'
#' # plot tracks with coloring each facet differently (facet.key == group.key)
#' ggcoverage(data = track.df, facet.key = "Type", group.key = "Type")
#'
#' # supply your own colors
#' ggcoverage(
#'   data = track.df, facet.key = "Type",
#'   group.key = "Type", color = 1:4,
#'   facet.color = 1:4
#' )
#'
#' # plot tracks together in one panel instead of separately;
#' # 'facet.key' is not needed
#' ggcoverage(
#'   data = track.df, group.key = "Type",
#'   plot.type = "joint"
#' )
#'
#' # use a custom theme
#' ggcoverage(data = track.df, facet.key = "Type") +
#'   theme_bw()
#'
#' # mark a region
#' ggcoverage(
#'   data = track.df, facet.key = "Type",
#'   mark.region = data.frame(
#'     start = c(21678900,21732001,21737590),
#'     end = c(21679900,21732400,21737650),
#'     label=c("M1", "M2", "M3")),
#'   mark.color = grey(0.4)
#' )
#'
#' # position range labels outside of tracks
#' ggcoverage(data = track.df, facet.key = "Type", range.position = "out")
#'
ggcoverage <- function(data, single.nuc = FALSE, mapping = NULL, color = NULL,
                       rect.color = NA, plot.type = c("facet", "joint"), facet.key = "Type", joint.avg = FALSE,
                       facet.order = NULL, facet.color = NULL, facet.y.scale = c("free", "fixed"),
                       group.key = "Group", range.size = 3, range.position = c("in", "out"), plot.space = 0.2,
                       mark.region = NULL, mark.color = "grey", mark.alpha = 0.5, show.mark.label = TRUE, mark.label.size = 4) {
  # check parameters
  range.position <- match.arg(arg = range.position)

  # create plot
  coverage.plot <- ggplot() +
    geom_coverage(
      data = data, mapping = mapping, color = color, rect.color = rect.color,
      single.nuc = single.nuc, plot.type = plot.type,
      facet.key = facet.key, joint.avg = joint.avg, facet.order = facet.order,
      facet.color = facet.color, facet.y.scale = facet.y.scale,
      group.key = group.key, range.size = range.size, range.position = range.position,
      mark.region = mark.region, mark.color = mark.color, mark.alpha = mark.alpha, show.mark.label = show.mark.label,
      mark.label.size = mark.label.size
    )
  # add theme
  if (range.position == "in") {
    coverage.plot +
      theme_coverage(space = plot.space)
  } else if (range.position == "out") {
    coverage.plot +
      theme_coverage2(space = plot.space)
  }
}
