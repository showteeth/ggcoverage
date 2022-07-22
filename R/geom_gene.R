#' Add Gene Annotation to Coverage Plot.
#'
#' @param gtf.gr Granges object of GTF, created with \code{\link{import.gff}}. Default: NULL.
#' @param overlap.gene.gap The gap between gene groups. Default: 0.1.
#' @param gene.size The line size of gene. Default: 1.
#' @param utr.size The line size of UTR. Default: 2.
#' @param exon.size The line size of exon. Default: 4.
#' @param arrow.size The line size of arrow. Default: 1.
#' @param color.by Color the line by. Default: strand.
#' @param fill.color Color used for \code{color.by}.
#' Default: darkblue for - (minus strand), darkgreen for + (plus strand).
#' @param show.utr Logical value, whether to show UTR. Default: TRUE.
#' @param arrow.gap The gap distance between arrow. Default: NULL.
#' @param arrow.num Total arrow num of whole region. Default: 50.
#' @param arrow.length The length of arrow. Default: 0.06.
#' @param label.size The size of gene label. Default: 3.
#' @param label.vjust The vjust of gene label. Default: 2.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of gene annotation to coverage plot. Default: 0.2.
#'
#' @return Plot.
#' @importFrom magrittr %>%
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
#' meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
#' sample.meta <- utils::read.csv(meta.file)
#' # track folder
#' track.folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
#' # load bigwig file
#' track.df <- LoadTrackFile(
#'   track.folder = track.folder, format = "bw",
#'   meta.info = sample.meta
#' )
#' gtf.file <- system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
#' gtf.gr <- rtracklayer::import.gff(con = gtf.file, format = "gtf")
#' basic.coverage <- ggcoverage(data = track.df, color = "auto", range.position = "out")
#' basic.coverage + geom_gene(gtf.gr = gtf.gr)
geom_gene <- function(gtf.gr, overlap.gene.gap = 0.1, gene.size = 1,
                      utr.size = 2, exon.size = 4, arrow.size = 1, color.by = "strand",
                      fill.color = c("-" = "darkblue", "+" = "darkgreen"), show.utr = TRUE,
                      arrow.gap = NULL, arrow.num = 50, arrow.length = 0.06,
                      label.size = 3, label.vjust = 2, plot.space = 0.1, plot.height = 0.2) {
  structure(list(
    gtf.gr = gtf.gr,
    overlap.gene.gap = overlap.gene.gap, gene.size = gene.size,
    utr.size = utr.size, exon.size = exon.size, arrow.size = arrow.size, color.by = color.by,
    fill.color = fill.color, show.utr = show.utr, arrow.gap = arrow.gap, arrow.num = arrow.num,
    arrow.length = arrow.length, label.size = label.size, label.vjust = label.vjust,
    plot.space = plot.space, plot.height = plot.height
  ),
  class = "gene"
  )
}

#' @export
ggplot_add.gene <- function(object, plot, object_name) {
  # get plot data
  track.data <- plot$layers[[1]]$data
  # prepare plot range
  plot.range.chr <- track.data[1, "seqnames"]
  plot.range.start <- track.data[1, "start"]
  plot.range.end <- track.data[nrow(track.data), "end"]
  plot.range.gr <- GenomicRanges::GRanges(
    seqnames = plot.range.chr,
    ranges = IRanges::IRanges(plot.range.start, plot.range.end)
  )
  # get parameters
  gtf.gr <- object$gtf.gr
  overlap.gene.gap <- object$overlap.gene.gap
  gene.size <- object$gene.size
  utr.size <- object$utr.size
  exon.size <- object$exon.size
  arrow.size <- object$arrow.size
  color.by <- object$color.by
  fill.color <- object$fill.color
  show.utr <- object$show.utr
  arrow.gap <- object$arrow.gap
  arrow.num <- object$arrow.num
  arrow.length <- object$arrow.length
  label.size <- object$label.size
  label.vjust <- object$label.vjust
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # process
  # get gene in region
  # gtf.gr <- rtracklayer::import.gff(gtf.file,format = 'gtf')
  gtf.df.used <- IRanges::subsetByOverlaps(x = gtf.gr, ranges = plot.range.gr) %>% as.data.frame()
  # select used features
  gene.info.used <- gtf.df.used %>%
    dplyr::filter(type %in% c("gene", "exon", "UTR")) %>%
    dplyr::select(c("seqnames", "start", "end", "strand", "type", "gene_type", "gene_name"))
  # modify region
  gene.info.used[gene.info.used$start <= plot.range.start, "start"] <- plot.range.start
  gene.info.used[gene.info.used$end >= plot.range.end, "end"] <- plot.range.end
  gene.info.used.gene <- gene.info.used %>% dplyr::filter(type == "gene")
  # convert dataframe to GR
  used.gene.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.info.used.gene,
    ignore.strand = TRUE,
    seqnames.field = c("seqnames"),
    start.field = "start",
    end.field = "end",
    strand.field = "strand"
  )
  # divide genes to non-overlap groups
  gene.group.idx <- GetGeneGroup(used.gene.gr, overlap.gene.gap = overlap.gene.gap)
  group.num <- length(unique(gene.group.idx))
  gene.info.used.gene$group <- gene.group.idx
  gene.info.used.gene$start <- as.numeric(gene.info.used.gene$start)
  gene.info.used.gene$end <- as.numeric(gene.info.used.gene$end)

  # get exon region
  gene.info.used <- merge(gene.info.used, gene.info.used.gene[c("gene_name", "group")], by = "gene_name")
  gene.info.used.exon <- gene.info.used %>% dplyr::filter(type == "exon")
  gene.info.used.exon$start <- as.numeric(gene.info.used.exon$start)
  gene.info.used.exon$end <- as.numeric(gene.info.used.exon$end)

  # create plot without arrow
  if (show.utr) {
    # get utr region
    gene.info.used.utr <- gene.info.used %>% dplyr::filter(type == "UTR")
    gene.info.used.utr$start <- as.numeric(gene.info.used.utr$start)
    gene.info.used.utr$end <- as.numeric(gene.info.used.utr$end)
    # substract UTR from exon
    gene.exon.utr <- SplitExonUTR(exon.df = gene.info.used.exon, utr.df = gene.info.used.utr)
    gene.info.used.exon <- gene.exon.utr$exon
    gene.info.used.utr <- gene.exon.utr$utr
    gene.plot <- ggplot() +
      geom_segment(
        data = gene.info.used.gene,
        mapping = aes_string(
          x = "start",
          y = "group",
          xend = "end",
          yend = "group",
          color = color.by
        ),
        show.legend = FALSE,
        size = gene.size
      ) +
      geom_segment(
        data = gene.info.used.utr,
        mapping = aes_string(
          x = "start",
          y = "group",
          xend = "end",
          yend = "group",
          color = color.by
        ),
        show.legend = FALSE,
        size = utr.size
      ) +
      geom_segment(
        data = gene.info.used.exon,
        mapping = aes_string(
          x = "start",
          y = "group",
          xend = "end",
          yend = "group",
          color = color.by
        ),
        show.legend = FALSE,
        size = exon.size
      )
  } else {
    gene.plot <- ggplot() +
      geom_segment(
        data = gene.info.used.gene,
        mapping = aes_string(
          x = "start",
          y = "group",
          xend = "end",
          yend = "group",
          color = color.by
        ),
        show.legend = FALSE,
        size = gene.size
      ) +
      geom_segment(
        data = gene.info.used.exon,
        mapping = aes_string(
          x = "start",
          y = "group",
          xend = "end",
          yend = "group",
          color = color.by
        ),
        show.legend = FALSE,
        size = exon.size
      )
  }

  if (is.null(arrow.gap)) {
    if (is.null(arrow.num)) {
      stop("Please provide either arrow.num or arrow.gap!")
    } else {
      arrow.gap <- (plot.range.end - plot.range.start) / arrow.num
    }
  }
  arrow.list <- list()
  # create arrow based on gene
  for (i in 1:nrow(gene.info.used.gene)) {
    gene.seq <- as.character(gene.info.used.gene[i, "seqnames"])
    gene.start <- as.numeric(gene.info.used.gene[i, "start"])
    gene.end <- as.numeric(gene.info.used.gene[i, "end"])
    gene.strand <- as.character(gene.info.used.gene[i, "strand"])
    gene.type <- as.character(gene.info.used.gene[i, "type"])
    gene.gene_type <- as.character(gene.info.used.gene[i, "gene_type"])
    gene.name <- as.character(gene.info.used.gene[i, "gene_name"])
    gene.group <- as.numeric(gene.info.used.gene[i, "group"])
    gene.gap <- gene.end - gene.start
    if (gene.gap <= arrow.gap) {
      # create only one arrow
      arrow.pos <- floor((gene.end + gene.start) / 2)
      arrow.list[[gene.name]] <- c(
        gene.seq, arrow.pos, arrow.pos + 1, gene.strand,
        gene.type, gene.gene_type, gene.name, gene.group
      )
    } else {
      gene.arrow.num <- floor(gene.gap / arrow.gap)
      gene.arrow.start <- (arrow.gap * 0:gene.arrow.num) + gene.start
      gene.arrow.end <- gene.arrow.start + 1
      for (grn in 1:length(gene.arrow.start)) {
        arrow.list[[paste(gene.name, grn, sep = "_")]] <-
          c(
            gene.seq, gene.arrow.start[grn], gene.arrow.end[grn], gene.strand,
            gene.type, gene.gene_type, gene.name, gene.group
          )
      }
    }
  }
  arrow.df <- do.call(rbind, arrow.list) %>% as.data.frame()
  colnames(arrow.df) <- c("seqnames", "start", "end", "strand", "type", "gene_type", "gene_name", "group")
  arrow.df$start <- as.numeric(arrow.df$start)
  arrow.df$end <- as.numeric(arrow.df$end)
  arrow.df$group <- as.numeric(arrow.df$group)

  gene.arrow.plot <- gene.plot + geom_segment(
    data = arrow.df,
    mapping = aes_string(
      x = "start",
      y = "group",
      xend = "end",
      yend = "group",
      color = color.by
    ),
    arrow = arrow(
      ends = ifelse(arrow.df$strand == "-", "first", "last"),
      type = "open",
      angle = 45,
      length = unit(x = arrow.length, units = "inches")
    ),
    show.legend = FALSE,
    size = arrow.size
  )

  label.df <- data.frame(
    pos = (gene.info.used.gene$start + gene.info.used.gene$end) / 2,
    group = gene.info.used.gene$group,
    gene = gene.info.used.gene$gene_name
  )

  gene.final.plot <- gene.arrow.plot +
    geom_text(
      data = label.df,
      mapping = aes_string(x = "pos", y = "group", label = "gene"),
      vjust = label.vjust, size = label.size
    ) +
    labs(y = "Gene") +
    theme_gene(
      overlap.gene.gap = overlap.gene.gap, group.num = group.num,
      fill.color = fill.color, x.range = c(plot.range.start, plot.range.end),
      margin.len = plot.space
    )
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    gene.final.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
