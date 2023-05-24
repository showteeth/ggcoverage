#' Add Feature Annotation to Coverage Plot.
#'
#' @param feature.file The path to feature file, should contain four columns. Default: NULL.
#' @param feature.df The dataframe contains feature informatin, should contain four columns. Default: NULL.
#' @param feature.color Feature color. Default: black.
#' @param feature.size The line size of peak. Default: 5.
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
#' # coverage.file <- system.file("extdata", "Proteomics", "MS_BSA_coverage.xlsx", package = "ggcoverage")
#' # fasta.file <- system.file("extdata", "Proteomics", "MS_BSA_coverage.fasta", package = "ggcoverage")
#' # protein.id = "sp|P02769|ALBU_BOVIN"
#' # protein.coverage = ggprotein(coverage.file = coverage.file, fasta.file = fasta.file, protein.id = protein.id)
#' # feature.df = data.frame(ProteinID = protein.id, start = c(1, 19, 25), end = c(18, 24, 607),
#' #                         Type = c("Signal", "Propeptide", "Chain"))
#' # protein.coverage +
#' #    geom_feature(feature.df = feature.df, feature.color = c("#4d81be","#173b5e","#6a521d"))
geom_feature <- function(feature.file = NULL, feature.df = NULL, feature.color = "black", feature.size = 5,
                         plot.space = 0.1, plot.height = 0.1) {
  structure(list(
    feature.file = feature.file, feature.df = feature.df, feature.color = feature.color, feature.size = feature.size,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "feature"
  )
}

#' @export
ggplot_add.feature <- function(object, plot, object_name) {
  # get plot data
  if ("patchwork" %in% class(plot)) {
    plot.data <- plot[[1]]$layers[[1]]$data
  } else {
    plot.data <- plot$layers[[1]]$data
  }
  # prepare plot range
  if ("ProteinID" %in% colnames(plot.data)) {
    plot.chr <- as.character(plot.data[1, "ProteinID"])
    plot.region.start <- min(plot.data[, "start"])
    plot.region.end <- max(plot.data[, "end"])
  } else {
    plot.chr <- as.character(plot.data[1, "seqnames"])
    plot.region.start <- min(plot.data[, "start"])
    plot.region.end <- max(plot.data[, "end"])
  }

  # get parameters
  feature.file <- object$feature.file
  feature.df <- object$feature.df
  feature.color <- object$feature.color
  feature.size <- object$feature.size
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # prepare peak dataframe
  if (!is.null(feature.file)) {
    feature.info <- utils::read.table(file = feature.file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  } else if (!is.null(feature.df)) {
    feature.info <- feature.df
  }
  feature.info <- feature.info[c(1, 2, 3, 4)]
  colnames(feature.info) <- c("seqnames", "start", "end", "Type")

  # get valid feature
  valid.feature <- GetRegion(chr = plot.chr, df = feature.info, start = plot.region.start, end = plot.region.end)

  # prepare fill color
  feature.type <- valid.feature$Type %>% unique()
  if (length(feature.type) < length(feature.color)) {
    used.feature.color <- feature.color[1:length(feature.type)]
    if (is.null(names(used.feature.color))) {
      names(used.feature.color) <- feature.type
    }
  } else {
    warning("The color you provided is smaller than Type column in data, select automatically!")
    used.feature.color <- AutoColor(data = valid.feature, n = 9, name = "Set1", key = "Type")
  }

  # create plot
  feature.plot <- ggplot() +
    geom_segment(
      data = valid.feature,
      mapping = aes_string(
        x = "start",
        y = "1",
        xend = "end",
        yend = "1",
        color = "Type"
      ),
      size = feature.size
    ) +
    labs(y = "Feature")

  # add theme
  feature.plot <- feature.plot + theme_feature(
    margin.len = plot.space, x.range = c(plot.region.start, plot.region.end),
    fill.color = used.feature.color
  )
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    feature.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
