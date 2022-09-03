#' Add Peak Annotation to Coverage Plot.
#'
#' @param bed.file The path to bed file.
#' @param peak.color Peak color. Default: black.
#' @param peak.size The line size of peak. Default: 5.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of peak annotation to coverage plot. Default: 0.2.
#'
#' @return Plot.
#' @importFrom utils read.table
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot_add ggplot geom_segment aes_string theme_classic theme element_blank element_text
#' element_rect margin scale_x_continuous scale_y_continuous coord_cartesian
#' @export
#'
#' @examples
#' # library(ggcoverage)
#' # library(rtracklayer)
#' # sample.meta <- data.frame(
#' #   SampleName = c("Chr18_MCF7_ER_1", "Chr18_MCF7_ER_2", "Chr18_MCF7_ER_3", "Chr18_MCF7_input"),
#' #   Type = c("MCF7_ER_1", "MCF7_ER_2", "MCF7_ER_3", "MCF7_input"),
#' #   Group = c("IP", "IP", "IP", "Input")
#' # )
#' # track folder
#' # track.folder <- system.file("extdata", "ChIP-seq", package = "ggcoverage")
#' # load bigwig file
#' # track.df <- LoadTrackFile(
#' #   track.folder = track.folder, format = "bw",
#' #   meta.info = sample.meta
#' # )
#' # gtf.file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
#' # gtf.gr <- rtracklayer::import.gff(con = gtf.file, format = "gtf")
#' # create mark region
#' # mark.region <- data.frame(start = c(76822533), end = c(76823743), label = c("Promoter"))
#' # basic.coverage <- ggcoverage(
#' #   data = track.df, color = "auto", region = "chr18:76822285-76900000",
#' #   mark.region = mark.region, show.mark.label = FALSE
#' # )
#' # get consensus peak file
#' # peak.file <- system.file("extdata", "ChIP-seq", "consensus.peak", package = "ggcoverage")
#' # basic.coverage + geom_gene(gtf.gr = gtf.gr) + geom_peak(bed.file = peak.file)
geom_peak <- function(bed.file, peak.color = "black", peak.size = 5,
                      plot.space = 0.1, plot.height = 0.1) {
  structure(list(
    bed.file = bed.file, peak.color = peak.color, peak.size = peak.size,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "peak"
  )
}

#' @export
ggplot_add.peak <- function(object, plot, object_name) {
  # get plot data
  plot.data <- plot$layers[[1]]$data
  # prepare plot range
  plot.chr <- as.character(plot.data[1, "seqnames"])
  plot.region.start <- plot$coordinates$limits$x[1]
  plot.region.end <- plot$coordinates$limits$x[2]

  # get parameters
  bed.file <- object$bed.file
  peak.color <- object$peak.color
  peak.size <- object$peak.size
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # prepare bed file
  bed.info <- utils::read.table(file = bed.file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  bed.info <- bed.info[c(1, 2, 3)]
  colnames(bed.info) <- c("seqnames", "start", "end")
  # convert to 1-based
  bed.info$start <- as.numeric(bed.info$start) + 1

  # get valid bed
  valid.bed <- GetRegion(chr = plot.chr, df = bed.info, start = plot.region.start, end = plot.region.end)

  peak.plot <- ggplot() +
    geom_segment(
      data = valid.bed,
      mapping = aes_string(
        x = "start",
        y = "1",
        xend = "end",
        yend = "1"
      ),
      size = peak.size,
      color = peak.color
    ) +
    labs(y = "Peak")

  # add theme
  peak.plot <- peak.plot + theme_peak(margin.len = plot.space, x.range = c(plot.region.start, plot.region.end))
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    peak.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
