#' Create Coverage Plot.
#'
#' @param data Coverage dataframe loaded by \code{\link{LoadTrackFile}}.
#' @param region Region used to create coverage plot, eg: chr14:21,677,306-21,737,601 or chr14:21,677,306.
#' Default: NULL.
#' @param gtf.gr Granges object of GTF, created with \code{\link{import.gff}}. Default: NULL.
#' @param extend Extend length of \code{region}. Default: 2000.
#' @param gene.name The name of gene. Default: HNRNPC.
#' @param gene.name.type Gene name type (filed of \code{gtf.gr}), chosen from gene_name and gene_id.
#' Default: gene_name.
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
#' @importFrom dplyr filter arrange
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
#' # track.df <- LoadTrackFile(
#' #   track.folder = track.folder, format = "bw",
#' #   meta.info = sample.meta
#' # )
#' # gtf.file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
#' # gtf.gr <- rtracklayer::import.gff(con = gtf.file, format = "gtf")
#' # ggcoverage(data = track.df, color = "auto", range.position = "out")
ggcoverage <- function(data, region = "chr14:21,677,306-21,737,601", gtf.gr = NULL, extend = 2000,
                       gene.name = "HNRNPC", gene.name.type = c("gene_name", "gene_id"), single.nuc = FALSE,
                       mapping = NULL, color = NULL, rect.color = NA, facet.key = "Type", facet.order = NULL, facet.color = NULL,
                       group.key = "Group", range.size = 3, range.position = c("in", "out"), plot.space = 0.2,
                       mark.region = NULL, mark.color = "grey", mark.alpha = 0.5, show.mark.label = TRUE, mark.label.size = 4) {
  # check parameters
  gene.name.type <- match.arg(arg = gene.name.type)
  range.position <- match.arg(arg = range.position)

  if (single.nuc) {
    data <- data
  } else {
    # formating data
    data <- FormatTrack(
      data = data, region = region, gtf.gr = gtf.gr, extend = extend,
      gene.name = gene.name, gene.name.type = gene.name.type
    )
  }
  plot.range.start <- data[1, "start"]
  plot.range.end <- data[nrow(data), "end"]

  # create plot
  coverage.plot <- ggplot() +
    geom_coverage(
      data = data, mapping = mapping,
      color = color, rect.color = rect.color, facet.key = facet.key, facet.order = facet.color, facet.color = facet.color,
      group.key = group.key, range.size = range.size, range.position = range.position,
      mark.region = mark.region, mark.color = mark.color, mark.alpha = mark.alpha, show.mark.label = show.mark.label,
      mark.label.size = mark.label.size
    )
  # add theme
  if (range.position == "in") {
    coverage.plot +
      theme_coverage(space = plot.space, x.range = c(plot.range.start, plot.range.end))
  } else if (range.position == "out") {
    coverage.plot +
      theme_coverage2(space = plot.space, x.range = c(plot.range.start, plot.range.end))
  }
}
