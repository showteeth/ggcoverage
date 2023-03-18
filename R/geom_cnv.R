#' Add CNV Annotation to Coverage Plot.
#'
#' @param cnv.df Dataframe contains copy number information, should contains at least
#' three columns (chr, position, bin value, copy number).
#' @param bin.col Column index for bin value (point). Default: 3.
#' @param cn.col Column index for copy number (line). Default: 4.
#' @param ref.cn Reference copy number (ploidy). Default: 2.
#' @param bin.point.color Point color of bin value. Default: "grey".
#' @param bin.point.alpha Point alpha of bin value. Default: 0.6.
#' @param cn.line.color Line color of copy number. Default: "red".
#' @param ref.line.color Line color of reference copy number (ploidy). Default: "black".
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of contact map to coverage plot. Default: 0.2.
#'
#' @return Plot.
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame start end
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom ggplot2 ggplot_add ggplot geom_point geom_line geom_hline aes_string labs theme_classic theme element_blank element_rect
#' element_text margin scale_x_continuous coord_cartesian
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
#' # library(ggcoverage)
#' # library(utils)
#' # library("BSgenome.Hsapiens.UCSC.hg19")
#' # # prepare files
#' # cnv.file <- system.file("extdata", "DNA-seq", "SRR054616_copynumber.txt", package = "ggcoverage")
#' # track.file <- system.file("extdata", "DNA-seq", "SRR054616.bw", package = "ggcoverage")
#' # # read CNV
#' # cnv.df = read.table(file = cnv.file, sep = "\t", header = TRUE)
#' # # load track
#' # track.df = LoadTrackFile(track.file = track.file, format = "bw")
#' # track.df$seqnames = paste0("chr", track.df$seqnames)
#' # # plot
#' # ggcoverage(data = track.df, color = "grey", region = "chr4:1-160000000",
#' #            mark.region = NULL, range.position = "out") +
#' #   geom_gc(bs.fa.seq=BSgenome.Hsapiens.UCSC.hg19) +
#' #   geom_cnv(cnv.df = cnv.df, bin.col = 3, cn.col = 4) +
#' #   geom_ideogram(genome = "hg19",plot.space = 0, highlight.centromere = TRUE)
geom_cnv <- function(cnv.df, bin.col = 3, cn.col = 4, ref.cn = 2,
                     bin.point.color = "grey", bin.point.alpha = 0.6, cn.line.color = "red",
                     ref.line.color = "black", plot.space = 0.1, plot.height = 0.2) {
  structure(list(
    cnv.df = cnv.df, bin.col = bin.col, cn.col = cn.col, ref.cn = ref.cn,
    bin.point.color = bin.point.color, bin.point.alpha = bin.point.alpha, cn.line.color = cn.line.color,
    ref.line.color = ref.line.color, plot.space = plot.space, plot.height = plot.height
  ),
  class = "cnv"
  )
}


#' @export
ggplot_add.cnv <- function(object, plot, object_name) {
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
  cnv.df <- object$cnv.df
  bin.col <- object$bin.col
  cn.col <- object$cn.col
  ref.cn <- object$ref.cn
  bin.point.color <- object$bin.point.color
  bin.point.alpha <- object$bin.point.alpha
  cn.line.color <- object$cn.line.color
  ref.line.color <- object$ref.line.color
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # subset dataframe
  used.cols <- c(1, 2)
  df.names <- c("chr", "start")
  if (!is.null(bin.col)) {
    used.cols <- c(used.cols, bin.col)
    df.names <- c(df.names, "binValue")
  }
  if (!is.null(cn.col)) {
    used.cols <- c(used.cols, cn.col)
    df.names <- c(df.names, "copyNumber")
  }
  if (is.null(bin.col) & is.null(cn.col)) {
    stop("Please provide at least one of bin.col and cn.col!")
  }
  cnv.df <- cnv.df[used.cols]
  colnames(cnv.df) <- df.names

  # create gr object
  cnv.df$end <- cnv.df$start + 1
  cnv.gr <- GenomicRanges::makeGRangesFromDataFrame(df = cnv.df, keep.extra.columns = TRUE)
  # filter CNV gr
  cnv.valid.df <- as.data.frame(IRanges::subsetByOverlaps(x = cnv.gr, ranges = plot.range.gr))
  # remove cnv outside region
  cnv.valid.df <- cnv.valid.df[cnv.valid.df$start >= GenomicRanges::start(x = plot.range.gr) &
    cnv.valid.df$end <= GenomicRanges::end(x = plot.range.gr), ]
  rownames(cnv.valid.df) <- 1:nrow(cnv.valid.df)

  # create empty plot
  cnv.basic.plot <- ggplot()
  # create basic plot
  if ("binValue" %in% colnames(cnv.valid.df)) {
    cnv.basic.plot <- cnv.basic.plot +
      geom_point(
        data = cnv.valid.df, mapping = aes_string(x = "start", y = "binValue"),
        color = bin.point.color, alpha = bin.point.alpha
      )
  }
  if ("copyNumber" %in% colnames(cnv.valid.df)) {
    cnv.basic.plot <- cnv.basic.plot +
      geom_line(
        data = cnv.valid.df, mapping = aes_string(x = "start", y = "copyNumber"),
        color = cn.line.color
      )
  }
  # add reference ploidy
  if (!is.null(ref.cn)) {
    cnv.basic.plot <- cnv.basic.plot +
      geom_hline(yintercept = ref.cn, color = ref.line.color)
  }

  # create final plot
  cnv.plot <-
    cnv.basic.plot +
    labs(y = "CopyNumber") +
    theme_cnv(x.range = c(plot.range.start, plot.range.end), margin.len = plot.space)
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    cnv.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
