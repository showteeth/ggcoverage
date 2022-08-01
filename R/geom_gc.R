#' Add GC Content Annotation to Coverage Plot.
#'
#' @param fa.file Genome fasta file. Default: NULL.
#' @param bs.fa.seq BSgenome for species. Default: NULL.
#' @param chr.split Split between chromosome name and description in \code{fa.file}. Default: "[[:space:]]".
#' @param guide.line GC content guide line. Default: NULL (use mean GC content).
#' @param line.color GC line color. Default: "black".
#' @param guide.line.color The color of guide line. Default: "red".
#' @param guide.line.type The line type of guide line. Default: "dashed".
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of GC content annotation to coverage plot. Default: 0.2.
#'
#' @return Plot.
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom Biostrings readDNAStringSet letterFrequency
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom BSgenome getSeq
#' @importFrom ggplot2 ggplot_add ggplot geom_line aes_string geom_hline labs theme_classic theme element_blank
#' element_text element_rect margin scale_x_continuous scale_y_continuous coord_cartesian
#' @export
#'
#' @examples
#' library(ggcoverage)
#' library(utils)
#' library(rtracklayer)
#' library("BSgenome.Hsapiens.UCSC.hg19")
#' # track folder
#' track.file <- system.file("extdata", "DNA-seq", "CNV_example.txt", package = "ggcoverage")
#' track.df <- utils::read.table(track.file, header = TRUE)
#' gtf.file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
#' gtf.gr <- rtracklayer::import.gff(con = gtf.file, format = "gtf")
#' basic.coverage <- ggcoverage(
#'   data = track.df, color = NULL, mark.region = NULL,
#'   region = "chr4:61750000-62,700,000", range.position = "out"
#' )
#' basic.coverage + geom_gc(bs.fa.seq = BSgenome.Hsapiens.UCSC.hg19)
geom_gc <- function(fa.file = NULL, bs.fa.seq = NULL, chr.split = "[[:space:]]", guide.line = NULL,
                    line.color = "black", guide.line.color = "red", guide.line.type = "dashed",
                    plot.space = 0.1, plot.height = 0.2) {
  structure(list(
    fa.file = fa.file, bs.fa.seq = bs.fa.seq, chr.split = chr.split, guide.line = guide.line,
    line.color = line.color, guide.line.color = guide.line.color, guide.line.type = guide.line.type,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "gc"
  )
}

#' @export
ggplot_add.gc <- function(object, plot, object_name) {
  # get plot data, plot data should contain bins
  plot.data <- plot$layers[[1]]$data
  # prepare plot range
  plot.chr <- as.character(plot.data[1, "seqnames"])
  plot.region.start <- plot$coordinates$limits$x[1]
  plot.region.end <- plot$coordinates$limits$x[2]
  # filter with first sample
  plot.data <- plot.data %>% dplyr::filter(Type == plot.data$Type[1])

  # get parameters
  fa.file <- object$fa.file
  bs.fa.seq <- object$bs.fa.seq
  chr.split <- object$chr.split
  guide.line <- object$guide.line
  line.color <- object$line.color
  guide.line.color <- object$guide.line.color
  guide.line.type <- object$guide.line.type
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # get sequence
  if (is.null(bs.fa.seq)) {
    if (is.null(fa.file)) {
      stop("Please provide either fa.seq or fa.file!")
    } else {
      fa.seq <- Biostrings::readDNAStringSet(fa.file)
      # change fasta name
      names(fa.seq) <- sapply(strsplit(x = names(fa.seq), split = chr.split), "[", 1)
      # filter sequence by chromosome
      fa.seq.selected <- fa.seq[names(fa.seq) == plot.chr]
    }
  } else {
    fa.seq.selected <- bs.fa.seq
  }

  # convert plot data to GRanges
  plot.data.gr <- GenomicRanges::makeGRangesFromDataFrame(plot.data,
    ignore.strand = TRUE,
    seqnames.field = c("seqnames"),
    start.field = "start",
    end.field = "end"
  )
  # get GRanges' sequence
  range.seqs <- BSgenome::getSeq(fa.seq.selected, plot.data.gr)
  # calculate GC content
  plot.data$GC <- as.numeric(Biostrings::letterFrequency(x = range.seqs, letters = "GC", as.prob = TRUE))

  # prepare x position
  plot.data$Middle <- (plot.data$start + plot.data$end) / 2
  # create plot
  gc.plot <- ggplot() +
    geom_line(data = plot.data, aes_string(x = "Middle", y = "GC", group = "1"), color = line.color)
  # add guide line
  if (is.null(guide.line)) {
    # use mean as guide line
    guide.line <- mean(plot.data$GC)
  }
  gc.plot <- gc.plot +
    geom_hline(yintercept = guide.line, color = guide.line.color, linetype = guide.line.type) +
    labs(y = "GC")

  # add theme
  gc.plot <- gc.plot + theme_gc(margin.len = plot.space, x.range = c(plot.region.start, plot.region.end))
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    gc.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
