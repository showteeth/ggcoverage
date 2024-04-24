#' Add Gene Annotation to Coverage Plot.
#'
#' @param gtf.gr Granges object of GTF, created with \code{\link{import.gff}}. Default: NULL.
#' @param overlap.gene.gap The gap between gene groups. Default: 0.1.
#' @param overlap.style The style of gene groups, choose from loose (each gene occupies single line)
#' and tight (place non-overlap genes in one line). Default: loose.
#' @param gene.size The line size of gene. Default: 1.
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
#' @param show.utr Logical value, whether to show UTR. Default: TRUE.
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
#'   geom_gene(gtf.gr = gtf_gr)
geom_gene <- function(gtf.gr,
                      overlap.gene.gap = 0.1,
                      overlap.style = "loose",
                      gene.size = 1,
                      utr.size = 2,
                      exon.size = 3,
                      arrow.size = 1.5,
                      arrow.gap = NULL,
                      arrow.num = 50,
                      color.by = "strand",
                      fill.color = c("-" = "cornflowerblue",
                                     "+" = "darkolivegreen3"),
                      show.utr = FALSE,
                      label.size = 3,
                      label.vjust = 2,
                      plot.space = 0.1,
                      plot.height = 0.2) {
  structure(
    list(
      gtf.gr = gtf.gr,
      overlap.gene.gap = overlap.gene.gap,
      overlap.style = overlap.style,
      gene.size = gene.size,
      utr.size = utr.size,
      exon.size = exon.size,
      arrow.size = arrow.size,
      arrow.gap = arrow.gap,
      arrow.num = arrow.num,
      color.by = color.by,
      fill.color = fill.color,
      show.utr = show.utr,
      label.size = label.size,
      label.vjust = label.vjust,
      plot.space = plot.space,
      plot.height = plot.height
    ),
    class = "gene"
  )
}

#' @export
ggplot_add.gene <- function(object, plot, object_name) {
  # get plot data
  if ("patchwork" %in% class(plot)) {
    track.data <- plot[[1]]$layers[[1]]$data
  } else {
    track.data <- plot$layers[[1]]$data
  }
  # prepare plot range
  # the plot region are not normal, so start is minimum value
  plot.range.chr <- track.data[1, "seqnames"]
  plot.range.start <- min(track.data[, "start"])
  plot.range.end <- max(track.data[, "end"])
  plot.range.gr <- GenomicRanges::GRanges(
    seqnames = plot.range.chr,
    ranges = IRanges::IRanges(plot.range.start, plot.range.end)
  )
  # get parameters
  gtf.gr <- object$gtf.gr
  overlap.gene.gap <- object$overlap.gene.gap
  overlap.style <- object$overlap.style
  gene.size <- object$gene.size
  utr.size <- object$utr.size
  exon.size <- object$exon.size
  arrow.size <- object$arrow.size
  color.by <- object$color.by
  fill.color <- object$fill.color
  show.utr <- object$show.utr
  arrow.gap <- object$arrow.gap
  arrow.num <- object$arrow.num
  label.size <- object$label.size
  label.vjust <- object$label.vjust
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # get gene in region
  gtf.df.used <- IRanges::subsetByOverlaps(x = gtf.gr, ranges = plot.range.gr) %>% as.data.frame()
  # check information
  used.gtf.columns <- c("seqnames", "start", "end", "strand", "type", "gene_name")
  gene.type.names <- c("gene_type", "gene_biotype")
  if (all(used.gtf.columns %in% colnames(gtf.df.used))) {
    used.gene.type.name <- intersect(gene.type.names, colnames(gtf.df.used))
    if (length(used.gene.type.name) < 1) {
      stop("gene_type/gene_biotype is not in provided GTF file!")
    } else {
      # select used features
      gene.info.used <- gtf.df.used %>%
        dplyr::filter(type %in% c("gene", "exon", "UTR")) %>%
        dplyr::select(c("seqnames", "start", "end", "strand", "type", used.gene.type.name, "gene_name"))
    }
    colnames(gene.info.used) <- c("seqnames", "start", "end", "strand", "type", "gene_type", "gene_name")
  } else {
    used.unique <- setdiff(c(used.gtf.columns, gene.type.names), colnames(gtf.df.used))
    stop(paste0(paste0(used.unique, collapse = ", "), " is not in provided GTF file!"))
  }
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
  if (overlap.style == "loose") {
    gene.group.idx <- GetGeneGroup(used.gene.gr, overlap.gene.gap = overlap.gene.gap)
  } else if (overlap.style == "tight") {
    gene.group.idx <- GetGeneGroupTight(used.gene.gr, overlap.gene.gap = overlap.gene.gap)
  }
  group.num <- length(unique(gene.group.idx))
  gene.info.used.gene$group <- gene.group.idx
  gene.info.used.gene$start <- as.numeric(gene.info.used.gene$start)
  gene.info.used.gene$end <- as.numeric(gene.info.used.gene$end)

  # get exon region
  gene.info.used <- merge(gene.info.used, gene.info.used.gene[c("gene_name", "group")], by = "gene_name")
  gene.info.used.exon <- gene.info.used %>% dplyr::filter(type == "exon")
  gene.info.used.exon$start <- as.numeric(gene.info.used.exon$start)
  gene.info.used.exon$end <- as.numeric(gene.info.used.exon$end)
  # get utr region
  gene.info.used.utr <- gene.info.used %>% dplyr::filter(type == "UTR")
  gene.info.used.utr$start <- as.numeric(gene.info.used.utr$start)
  gene.info.used.utr$end <- as.numeric(gene.info.used.utr$end)
  # change UTR
  if (show.utr & nrow(gene.info.used.utr) == 0) {
    warning("No UTR detected in provided GTF, omitting plotting UTRs.")
    show.utr <- FALSE
  }
  # plot genomic features with arrow at the end
  if (show.utr) {
    # substract UTR from exon
    gene.exon.utr <- SplitExonUTR(exon.df = gene.info.used.exon, utr.df = gene.info.used.utr)
    gene.info.used.exon <- gene.exon.utr$exon
    gene.info.used.utr <- gene.exon.utr$utr
  }
  gene.plot <- ggplot() +
    geom_arrows(gene.info.used.gene, color.by, gene.size, arrow.size) +
    geom_arrows(gene.info.used.exon, color.by, exon.size, arrow.size)
  if (show.utr) {
    gene.plot <- gene.plot +
      geom_arrows(gene.info.used.utr, color.by, utr.size, arrow.size)
  }

  if (!is.null(arrow.gap) || !is.null(arrow.num)) {
    if (!is.null(arrow.num)) {
      arrow.gap <- (plot.range.end - plot.range.start) / arrow.num
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
    gene.plot <- gene.plot +
      geom_arrows(arrow.df, color.by, gene.size / 2, arrow.size, 35, TRUE)
  }

  label.df <- data.frame(
    pos = (gene.info.used.gene$start + gene.info.used.gene$end) / 2,
    group = gene.info.used.gene$group,
    gene = gene.info.used.gene$gene_name
  )

  gene.final.plot <- gene.plot +
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
