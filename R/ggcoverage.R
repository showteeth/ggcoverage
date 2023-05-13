#' Create Coverage Plot.
#'
#' @param data Coverage dataframe loaded by \code{\link{LoadTrackFile}}.
#' @param single.nuc Logical value, whether to visualize at single nucleotide level. Default: FALSE.
#' @param mapping Set of aesthetic mappings created by \code{aes} or \code{aes_}. Default: NULL.
#' @param color Track color. Default: NULL (select automatically).
#' @param rect.color The color of every bin. Default: NA.
#' @param facet.key Sample type key to create coverage plot. Default: Type.
#' @param facet.order The order of Coverage plot. Default: NULL.
#' @param facet.color The color of sample text. Default: NULL (select automatically).
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
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual geom_rect geom_text aes theme_classic theme unit
#' element_blank annotate rel scale_y_continuous expansion scale_x_continuous coord_cartesian
#' @importFrom scales comma
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang as_label
#' @importFrom stats as.formula
#' @importFrom ggh4x facet_wrap2 strip_themed
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_text_repel
#' @export
#'
#' @examples
#' # library(ggcoverage)
#' # library(utils)
#' # library(rtracklayer)
#' # meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
#' # sample.meta <- utils::read.csv(meta.file)
#' # track folder
#' # track.folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
#' # load bigwig file
#' # track.df <- LoadTrackFile(track.folder = track.folder, format = "bw",region = "chr14:21,677,306-21,737,601",
#' #                           extend = 2000, meta.info = sample.meta)
#' # gtf.file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
#' # gtf.gr <- rtracklayer::import.gff(con = gtf.file, format = "gtf")
#' # ggcoverage(data = track.df, color = "auto", range.position = "out")
ggcoverage <- function(data, single.nuc = FALSE, mapping = NULL, color = NULL,
                       rect.color = NA, facet.key = "Type", facet.order = NULL, facet.color = NULL,
                       group.key = "Group", range.size = 3, range.position = c("in", "out"), plot.space = 0.2,
                       mark.region = NULL, mark.color = "grey", mark.alpha = 0.5, show.mark.label = TRUE, mark.label.size = 4) {
  # check parameters
  range.position <- match.arg(arg = range.position)

  # create plot
  coverage.plot <- ggplot() +
    geom_coverage(
      data = data, mapping = mapping, color = color, rect.color = rect.color,
      single.nuc = single.nuc, facet.key = facet.key, facet.order = facet.order, facet.color = facet.color,
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
