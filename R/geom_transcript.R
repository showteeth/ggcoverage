#' Add Transcript Annotation to Coverage Plot.
#'
#' @param gtf.gr Granges object of GTF, created with \code{\link{import.gff}}. Default: NULL.
#' @param gene.name Gene name of all transcripts. Default: HNRNPC.
#' @param overlap.tx.gap The gap between transcript groups. Default: 0.1.
#' @param overlap.style The style of transcript groups, choose from loose (each transcript occupies single line)
#' and tight (place non-overlap transcripts in one line). Default: loose.
#' @param tx.size The line size of transcript. Default: 1.
#' @param utr.size The line size of UTR. Default: 2.
#' @param exon.size The line size of exon. Default: 3.
#' @param arrow.size The line size of arrow. Default: 1.5.
#' @param arrow.gap The gap distance between intermittent arrows. Default: NULL.
#'   Set arrow.num and arrow.gap to NULL to suppress intermittent arrows.
#' @param arrow.num Total number of intermittent arrows over whole region. Default: 50.
#'   Set arrow.num and arrow.gap to NULL to suppress intermittent arrows.
#' @param color.by Color the line by. Default: strand.
#' @param fill.color Color used for \code{color.by}.
#' Default: blue for - (minus strand), green for + (plus strand).
#' @param label.size The size of transcript label. Default: 3.
#' @param label.vjust The vjust of transcript label. Default: 2.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of transcript annotation to coverage plot. Default: 0.2.
#'
#' @return Plot.
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame setdiff
#' @importFrom IRanges IRanges subsetByOverlaps findOverlaps
#' @importFrom dplyr filter select arrange
#' @importFrom ggplot2 ggplot_add ggplot geom_segment aes_string arrow unit geom_text labs theme_classic theme element_blank
#' element_text element_rect margin scale_y_continuous scale_color_manual scale_x_continuous coord_cartesian
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
#' library(ggcoverage)
#' library(utils)
#' library(rtracklayer)
#'
#' # load metadata
#' meta_file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
#' sample_meta <- read.csv(meta_file)
#'
#' # track folder
#' track_folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
#'
#' # load bigwig file
#' track_df <- LoadTrackFile(
#'   track.folder = track_folder,
#'   format = "bw",
#'   region = "chr14:21,677,306-21,737,601",
#'   extend = 2000,
#'   meta.info = sample_meta
#' )
#'
#' # load GTF file
#' gtf_file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
#' gtf_gr <- rtracklayer::import.gff(con = gtf_file, format = "gtf")
#'
#' # plot coverage and gene annotation
#' basic.coverage <- ggcoverage(data = track_df, range.position = "out")
#' basic.coverage +
#'   geom_transcript(gtf.gr = gtf_gr, label.vjust = 1.5)
geom_transcript <-
  function(gtf.gr,
           gene.name = "HNRNPC",
           overlap.tx.gap = 0.1,
           overlap.style = "loose",
           tx.size = 1,
           utr.size = 2,
           exon.size = 3,
           arrow.size = 3,
           arrow.gap = NULL,
           arrow.num = 50,
           color.by = "strand",
           fill.color = c(
             "-" = "cornflowerblue",
             "+" = "darkolivegreen3"
           ),
           label.size = 3,
           label.vjust = 2,
           plot.space = 0.1,
           plot.height = 1) {
    structure(
      list(
        gtf.gr = gtf.gr,
        gene.name = gene.name,
        overlap.tx.gap = overlap.tx.gap,
        overlap.style = overlap.style,
        tx.size = tx.size,
        utr.size = utr.size,
        exon.size = exon.size,
        arrow.size = arrow.size,
        arrow.gap = arrow.gap,
        arrow.num = arrow.num,
        color.by = color.by,
        fill.color = fill.color,
        label.size = label.size,
        label.vjust = label.vjust,
        plot.space = plot.space,
        plot.height = plot.height
      ),
      class = "transcript"
    )
  }

#' @export
ggplot_add.transcript <- function(object, plot, object_name) {
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
  gtf.gr <- object$gtf.gr
  gene.name <- object$gene.name
  overlap.tx.gap <- object$overlap.tx.gap
  overlap.style <- object$overlap.style
  tx.size <- object$tx.size
  utr.size <- object$utr.size
  exon.size <- object$exon.size
  arrow.size <- object$arrow.size
  color.by <- object$color.by
  fill.color <- object$fill.color
  arrow.gap <- object$arrow.gap
  arrow.num <- object$arrow.num
  label.size <- object$label.size
  label.vjust <- object$label.vjust
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # filter gene in region
  gtf.df.used <- IRanges::subsetByOverlaps(x = gtf.gr, ranges = plot.range.gr) %>% as.data.frame()
  # check information
  used.gtf.columns <- c("seqnames", "start", "end", "strand", "type", "transcript_name")
  if (!all(used.gtf.columns %in% colnames(gtf.df.used))) {
    used.unique <- setdiff(used.gtf.columns, colnames(gtf.df.used))
    stop(paste0(paste0(used.unique, collapse = ", "), " is not in provided GTF file!"))
  }
  # select used features
  gene.tx.df <- gtf.df.used %>%
    dplyr::filter(gene_name == gene.name) %>%
    dplyr::filter(type %in% c("transcript", "exon", "UTR")) %>%
    dplyr::select(c("seqnames", "start", "end", "strand", "type", "transcript_name"))
  # modify region
  gene.tx.df[gene.tx.df$start <= plot.range.start, "start"] <- plot.range.start
  gene.tx.df[gene.tx.df$end >= plot.range.end, "end"] <- plot.range.end

  # prepare plot dataframe
  # seperate: transcript, UTR, exon
  gene.tx.df.tx <- gene.tx.df %>% dplyr::filter(type == "transcript")
  # convert dataframe to GR
  tx.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.tx.df.tx,
    ignore.strand = TRUE,
    seqnames.field = c("seqnames"),
    start.field = "start",
    end.field = "end",
    strand.field = "strand"
  )
  # divide genes to non-overlap groups
  if (overlap.style == "loose") {
    tx.group.idx <- GetGeneGroup(tx.gr, overlap.gene.gap = overlap.tx.gap)
  } else if (overlap.style == "tight") {
    tx.group.idx <- GetGeneGroupTight(tx.gr, overlap.gene.gap = overlap.tx.gap)
  }
  group.num <- length(unique(tx.group.idx))
  gene.tx.df.tx$group <- tx.group.idx
  gene.tx.df.tx$start <- as.numeric(gene.tx.df.tx$start)
  gene.tx.df.tx$end <- as.numeric(gene.tx.df.tx$end)

  # get exon region
  gene.tx.df <- merge(gene.tx.df, gene.tx.df.tx[c("transcript_name", "group")], by = "transcript_name")
  gene.tx.df.exon <- gene.tx.df %>% dplyr::filter(type == "exon")
  gene.tx.df.exon$start <- as.numeric(gene.tx.df.exon$start)
  gene.tx.df.exon$end <- as.numeric(gene.tx.df.exon$end)

  # get utr region
  gene.tx.df.utr <- gene.tx.df %>% dplyr::filter(type == "UTR")
  gene.tx.df.utr$start <- as.numeric(gene.tx.df.utr$start)
  gene.tx.df.utr$end <- as.numeric(gene.tx.df.utr$end)

  # substract UTR from exon
  tx.exon.utr <- SplitTxExonUTR(exon.df = gene.tx.df.exon, utr.df = gene.tx.df.utr)
  gene.tx.df.exon <- tx.exon.utr$exon
  gene.tx.df.utr <- tx.exon.utr$utr

  # create basic plot
  tx.plot <- ggplot() +
    geom_arrows(gene.tx.df.tx, color.by, tx.size, arrow.size)

  # deal with missing UTR
  if (is.null(gene.tx.df.utr)) {
    warning("No UTR detected in provided GTF!")
  } else {
    tx.plot <- tx.plot +
      geom_arrows(gene.tx.df.utr, color.by, utr.size, arrow.size)
  }
  tx.plot <- tx.plot +
    geom_arrows(gene.tx.df.exon, color.by, exon.size, arrow.size) +
    theme_classic()

  if (is.null(arrow.gap)) {
    if (is.null(arrow.num)) {
      stop("Please provide either arrow.num or arrow.gap!")
    } else {
      arrow.gap <- (plot.range.end - plot.range.start) / arrow.num
    }
  }
  arrow.list <- list()
  # create arrow based on gene
  for (i in 1:nrow(gene.tx.df.tx)) {
    tx.seq <- as.character(gene.tx.df.tx[i, "seqnames"])
    tx.start <- as.numeric(gene.tx.df.tx[i, "start"])
    tx.end <- as.numeric(gene.tx.df.tx[i, "end"])
    tx.strand <- as.character(gene.tx.df.tx[i, "strand"])
    tx.type <- as.character(gene.tx.df.tx[i, "type"])
    tx.name <- as.character(gene.tx.df.tx[i, "transcript_name"])
    tx.group <- as.numeric(gene.tx.df.tx[i, "group"])
    tx.gap <- tx.end - tx.start
    if (tx.gap <= arrow.gap) {
      # create only one arrow
      arrow.pos <- floor((tx.end + tx.start) / 2)
      arrow.list[[tx.name]] <- c(
        tx.seq, arrow.pos, arrow.pos + 1, tx.strand,
        tx.type, tx.name, tx.group
      )
    } else {
      tx.arrow.num <- floor(tx.gap / arrow.gap)
      tx.arrow.start <- (arrow.gap * 0:tx.arrow.num) + tx.start
      tx.arrow.end <- tx.arrow.start + 1
      for (trn in 1:length(tx.arrow.start)) {
        arrow.list[[paste(tx.name, trn, sep = "_")]] <-
          c(
            tx.seq, tx.arrow.start[trn], tx.arrow.end[trn], tx.strand,
            tx.type, tx.name, tx.group
          )
      }
    }
  }
  arrow.df <- do.call(rbind, arrow.list) %>% as.data.frame()
  colnames(arrow.df) <- c("seqnames", "start", "end", "strand", "type", "transcript_name", "group")
  arrow.df$start <- as.numeric(arrow.df$start)
  arrow.df$end <- as.numeric(arrow.df$end)
  arrow.df$group <- as.numeric(arrow.df$group)
  # add arrow
  tx.arrow.plot <- tx.plot +
    geom_arrows(arrow.df, color.by, tx.size / 2, arrow.size, 35, TRUE)

  # prepare label dataframe
  label.df <- data.frame(
    pos = (gene.tx.df.tx$start + gene.tx.df.tx$end) / 2,
    group = gene.tx.df.tx$group,
    gene = gene.tx.df.tx$transcript_name
  )
  # add label to plot
  tx.final.plot <- tx.arrow.plot +
    geom_text(
      data = label.df,
      mapping = aes_string(x = "pos", y = "group", label = "gene"),
      vjust = label.vjust, hjust = 0.5, size = label.size
    ) +
    labs(y = "Transcript") +
    theme_transcript(
      overlap.tx.gap = overlap.tx.gap, group.num = group.num,
      fill.color = fill.color, x.range = c(plot.range.start, plot.range.end),
      margin.len = plot.space
    )
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    tx.final.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
