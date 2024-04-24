# prepare GR
PrepareRegion <- function(region = NULL,
                          gtf.gr = NULL, gene.name = "HNRNPC", gene.name.type = c("gene_name", "gene_id"),
                          extend = 2000) {
  # check parameters
  gene.name.type <- match.arg(arg = gene.name.type)

  if (!is.null(region)) {
    region.split <- unlist(strsplit(x = region, split = ":"))
    region.chr <- region.split[1]
    region.se <- unlist(strsplit(x = region.split[2], split = "-"))
    if (length(region.se) == 1) {
      region.start <- as.numeric(gsub(pattern = ",", replacement = "", x = region.se[1]))
      region.end <- region.start
    } else if (length(region.se) == 2) {
      region.start <- as.numeric(gsub(pattern = ",", replacement = "", x = region.se[1]))
      region.end <- as.numeric(gsub(pattern = ",", replacement = "", x = region.se[2]))
    }
  } else if (!is.null(gtf.gr)) {
    # get gene related gtf
    gene.gtf.info <- gtf.gr %>%
      as.data.frame() %>%
      dplyr::filter(type == "gene")
    # get specific gene
    gene.gtf.info.used <- gene.gtf.info[gene.gtf.info[, gene.name.type] == gene.name, ]
    # get position
    region.chr <- as.character(gene.gtf.info.used$seqnames)
    region.start <- gene.gtf.info.used$start
    region.end <- gene.gtf.info.used$end
  }

  # extend region
  if (!(is.null(extend))) {
    # extend start
    region.start <- region.start - extend
    # avoid invalid start region
    if (region.start < 1) {
      region.start <- 1
    }
    region.end <- region.end + extend
  }

  if (region.start == region.end) {
    # avoid invalid start region
    message("The start position is same as end position, and the extension is 0. Automatically increase by 1!")
    region.end <- region.start + 1
  }

  # prepare gr
  gr <- GenomicRanges::GRanges(
    seqnames = region.chr,
    ranges = IRanges::IRanges(region.start, region.end)
  )
  return(gr)
}

# select color automatically
AutoColor <- function(data, n, name, key) {
  getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n, name))
  # sample group with same color
  group.info <- unique(data[, key])
  fill.color <- getPalette(length(group.info))
  names(fill.color) <- group.info
  return(fill.color)
}

# create aa plot dataframe with padding offset
AAPadding <- function(len, offset = 0, aa.seq) {
  aa.seq.full <- rep(aa.seq, each = 3)
  aa.seq.full.len <- length(aa.seq.full)
  aa.padding.len <- len - aa.seq.full.len - offset
  final.aa.seq <- c(rep("B", offset), aa.seq.full, rep("B", aa.padding.len))
  aa.anno.seq <- unlist(strsplit(x = gsub(pattern = "(.*)", replacement = "X\\1X", x = aa.seq), split = ""))
  aa.anno.seq <- gsub(pattern = "X", replacement = "", x = aa.anno.seq)
  final.anno.seq <- c(rep("", offset), aa.anno.seq, rep("", aa.padding.len))
  final.aa.df <- data.frame(aa = final.aa.seq, anno = final.anno.seq)
  return(final.aa.df)
}

# divide genes to non-overlap groups (V2)
GetGeneGroup <- function(gene.gr, fc = "queryHits", sc = "subjectHits", overlap.gene.gap = 1) {
  overlap.df <- IRanges::findOverlaps(gene.gr, gene.gr, ignore.strand = TRUE) %>%
    as.data.frame() %>%
    dplyr::arrange(.data[[fc]], .data[[sc]])
  overlap.list <- lapply(split(overlap.df, overlap.df[, fc]), function(x) {
    x[, sc]
  })
  overlap.list <- overlap.list[order(sapply(overlap.list, length), decreasing = TRUE)]
  group.idx <- rep(1, length(unique(overlap.df[, fc])))
  for (i in names(overlap.list)) {
    overlap.vec <- overlap.list[[i]]
    i <- as.numeric(i)
    curr.group <- group.idx[i]
    remain.vec <- setdiff(overlap.vec, i)
    for (j in remain.vec) {
      if (group.idx[j] == curr.group) {
        group.idx[j] <- curr.group + overlap.gene.gap
      }
    }
  }
  return(group.idx)
}

# place non-overlapping transcripts together
GetGeneGroupTight <- function(gene.gr, overlap.gene.gap = 1) {
  # convert to dataframe
  gene.gr.df <- as.data.frame(gene.gr)
  gene.gr.df$ID <- seq_len(nrow(gene.gr.df))
  # split to group
  group.flag <- 1
  group.list <- list()
  while (nrow(gene.gr.df) > 0) {
    A <- gene.gr.df[1, ]
    vec <- c(A$ID)
    if (nrow(gene.gr.df) >= 2) {
      for (i in 2:nrow(gene.gr.df)) {
        if (gene.gr.df[i, "start"] > A[1, "end"]) {
          vec <- c(vec, gene.gr.df[i, "ID"])
          A$end <- gene.gr.df[i, "end"]
        }
      }
    }
    group.list[[paste0("G", group.flag)]] <- vec
    gene.gr.df <- gene.gr.df %>% dplyr::filter(!.data$ID %in% vec)
    group.flag <- group.flag + 1
  }
  # get group index
  group.index <- c()
  for (g in seq_along(group.list)) {
    g.index <- 1 + overlap.gene.gap * (g - 1)
    g.index.vec <- rep(g.index, length(group.list[[g]]))
    names(g.index.vec) <- group.list[[g]]
    group.index <- c(group.index, g.index.vec)
  }
  group.index <- group.index[order(as.numeric(names(group.index)))]
  return(group.index)
}

# substract UTR from exon
SplitExonUTR <- function(exon.df, utr.df) {
  # get metadata
  exon.meta <- exon.df[c("gene_name", "type", "gene_type", "group")] %>% unique()
  # get list
  exon.list <- split(exon.df, exon.df$gene_name)
  utr.list <- split(utr.df, utr.df$gene_name)
  # all genes
  gene.process <- union(names(exon.list), names(utr.list))
  # results list
  final.utr.list <- list()
  final.exon.list <- list()
  # process
  for (gene in gene.process) {
    gene.exon <- exon.list[[gene]]
    gene.utr <- utr.list[[gene]]
    if (is.null(gene.exon)) {
      final.utr.list[[gene]] <- gene.utr
    } else if (is.null(gene.utr)) {
      final.exon.list[[gene]] <- gene.exon
    } else {
      exon.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.exon,
        ignore.strand = FALSE, keep.extra.columns = TRUE,
        seqnames.field = c("seqnames"),
        start.field = "start",
        end.field = "end",
        strand.field = "strand"
      )
      utr.gr <- GenomicRanges::makeGRangesFromDataFrame(gene.utr,
        ignore.strand = FALSE, keep.extra.columns = TRUE,
        seqnames.field = c("seqnames"),
        start.field = "start",
        end.field = "end",
        strand.field = "strand"
      )
      pure.exon <- GenomicRanges::setdiff(exon.gr, utr.gr) %>% as.data.frame()
      pure.exon$gene_name <- gene
      pure.exon <- merge(pure.exon, exon.meta, by = "gene_name") %>% dplyr::select(-"width")
      final.exon.list[[gene]] <- pure.exon
      final.utr.list[[gene]] <- gene.utr
    }
  }
  final.utr <- do.call(rbind, final.utr.list)
  rownames(final.utr) <- NULL
  final.exon <- do.call(rbind, final.exon.list)
  rownames(final.exon) <- NULL
  return(list(exon = final.exon, utr = final.utr))
}

# substract UTR from exon
SplitTxExonUTR <- function(exon.df, utr.df) {
  # get metadata
  exon.meta <- exon.df[c("transcript_name", "type", "group")] %>% unique()
  # get list
  exon.list <- split(exon.df, exon.df$transcript_name)
  utr.list <- split(utr.df, utr.df$transcript_name)
  # all transcripts
  tx.process <- union(names(exon.list), names(utr.list))
  # results list
  final.utr.list <- list()
  final.exon.list <- list()
  # process
  for (tx in tx.process) {
    tx.exon <- exon.list[[tx]]
    tx.utr <- utr.list[[tx]]
    if (is.null(tx.exon)) {
      final.utr.list[[tx]] <- tx.utr
    } else if (is.null(tx.utr)) {
      final.exon.list[[tx]] <- tx.exon
    } else {
      exon.gr <- GenomicRanges::makeGRangesFromDataFrame(tx.exon,
        ignore.strand = FALSE, keep.extra.columns = TRUE,
        seqnames.field = c("seqnames"),
        start.field = "start",
        end.field = "end",
        strand.field = "strand"
      )
      utr.gr <- GenomicRanges::makeGRangesFromDataFrame(tx.utr,
        ignore.strand = FALSE, keep.extra.columns = TRUE,
        seqnames.field = c("seqnames"),
        start.field = "start",
        end.field = "end",
        strand.field = "strand"
      )
      pure.exon <- GenomicRanges::setdiff(exon.gr, utr.gr) %>% as.data.frame()
      pure.exon$transcript_name <- tx
      pure.exon <- merge(pure.exon, exon.meta, by = "transcript_name") %>% dplyr::select(-"width")
      final.exon.list[[tx]] <- pure.exon
      final.utr.list[[tx]] <- tx.utr
    }
  }
  final.utr <- do.call(rbind, final.utr.list)
  rownames(final.utr) <- NULL
  final.exon <- do.call(rbind, final.exon.list)
  rownames(final.exon) <- NULL
  return(list(exon = final.exon, utr = final.utr))
}

# From: https://github.com/jorainer/biovizBase/blob/master/R/ideogram.R
# Fix bug: the names on the supplied 'seqlengths' vector must be
# identical to the seqnames
getIdeogram <- function(genome, subchr = NULL, cytobands = TRUE) {
  .gnm <- genome
  lst <- lapply(.gnm, function(genome) {
    if (!(exists("session") && extends(class(session), "BrowserSession"))) {
      session <- rtracklayer::browserSession()
    }
    if (missing(genome)) {
      choices <- rtracklayer::ucscGenomes()[, 1]
      res <- utils::menu(choices, title = "Please specify genome")
      genome <- as.character(choices[res])
    }
    if (cytobands) {
      message("Loading ideogram...")
      tryres <- try(query <-
        rtracklayer::ucscTableQuery(session, "cytoBand", genome))
      if (!inherits(tryres, "try-error")) {
        rtracklayer::tableName(query) <- "cytoBand"
        df <- rtracklayer::getTable(query)
        gr <- GenomicRanges::GRanges(
          seqnames = df$chrom,
          IRanges(start = df$chromStart, end = df$chromEnd)
        )
        S4Vectors::values(gr) <- df[, c("name", "gieStain")]
        message("Loading ranges...")

        gr.r <- rtracklayer::GRangesForUCSCGenome(genome)
        # fix bug
        new.seqlength <- GenomeInfoDb::seqlengths(gr.r)[names(GenomeInfoDb::seqlengths(gr))]
        names(new.seqlength) <- names(GenomeInfoDb::seqlengths(gr))
        suppressWarnings(GenomeInfoDb::seqlengths(gr) <- new.seqlength)
        gr <- GenomicRanges::trim(gr)
      } else {
        message("cytoBand informatin is not available, only get ranges.")
        message("Loading ranges...")
        gr <- rtracklayer::GRangesForUCSCGenome(genome)
        message("Done")
      }
    } else {
      message("Loading...")
      gr <- rtracklayer::GRangesForUCSCGenome(genome)
      message("Done")
    }
    if (length(subchr)) {
      gr <- gr[GenomeInfoDb::seqnames(gr) == subchr]
    }

    gr <- sort(gr)
    gr
  })
  names(lst) <- .gnm
  if (length(lst) == 1) {
    res <- lst[[1]]
  } else {
    res <- lst
  }
  res
}

# used in geom_base
PrepareRect <- function(df, y.center = -0.2) {
  valid.df <- df[df$aa == "B" | df$anno != "", ]
  rect.li <- lapply(seq_len(nrow(valid.df)), function(x) {
    row.info <- valid.df[x, ]
    if (row.info$aa == "B") {
      c(row.info$Pos - 0.5, row.info$Pos + 0.5, row.info$aa)
    } else if (row.info$anno != "") {
      c(row.info$Pos - 1.5, row.info$Pos + 1.5, row.info$aa)
    }
  })
  rect.df <- as.data.frame(t(as.data.frame(rect.li)))
  colnames(rect.df) <- c("xmin", "xmax", "aa")
  rownames(rect.df) <- NULL
  rect.df$xmin <- as.numeric(rect.df$xmin)
  rect.df$xmax <- as.numeric(rect.df$xmax)
  rect.df$ymin <- y.center - 0.1
  rect.df$ymax <- y.center + 0.1
  return(rect.df)
}

#' Plot Object Generated by ggcoverage.
#'
#' @param plot Plot.
#' @param layer.num The number of layers. Default: 1.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' # cov.plot = ggcoverage(data = track.df, color = "auto", region = "chr18:76822285-76900000",
#' #                       range.position = "out", mark.region=mark.region, show.mark.label = TRUE)
#' # plot.data = GetPlotData(plot = cov.plot, layer.num=1)
GetPlotData <- function(plot, layer.num = 1) {
  plot.str <- deparse(substitute(plot))
  if (layer.num == 1) {
    plot.layer.str <- plot.str
  } else if (layer.num >= 2) {
    plot.layer.str <- paste0(plot.str, paste0(rep("[[1]]", layer.num - 1), collapse = ""))
  }
  plot.data <- eval(parse(text = paste0(plot.layer.str, "$layers[[1]]$data")))
  return(plot.data)
}

#' Plot genomic features as arrows.
#' @description
#' This function is a variation of geom_segment to plot (gene) features
#' as arrows. Mainly meant for internal use, not to be called directly.
#'
#' @param data data frame describing arrow position, with columns
#'   start, end, group, and a custom 'color' column
#' @param color name of the color column in the data frame
#' @param line_width line_width of the (arrow) segment
#' @param arrow_size size of the arrow
#' @param arrow_angle angle of the arrow. Default: 35°
#' @param intermittent If TRUE, arrows are only drawn intermittently in
#'   half-transparent white color. Default: FALSE.
#' @importFrom grDevices grey
#' @return A geom layer for ggplot2 objects.
#' @export
geom_arrows <-
  function(data,
           color,
           line_width,
           arrow_size,
           arrow_angle = 35,
           intermittent = FALSE) {
    if (nrow(data)) {
      if (!"strand" %in% colnames(data)) {
        data$strand <- "+"
      }
      if (!intermittent) {
        geom_segment(
          data = data,
          mapping = aes_string(
            x = "start",
            y = "group",
            xend = "end",
            yend = "group",
            color = color
          ),
          arrow = arrow(
            ends = ifelse(data$strand == "+", "last", "first"),
            angle = arrow_angle,
            length = unit(arrow_size, "points"),
            type = "open"
          ),
          lineend = "butt",
          linejoin = "mitre",
          show.legend = FALSE,
          linewidth = line_width
        )
      } else {
        geom_segment(
          data = data,
          mapping = aes_string(
            x = "start",
            y = "group",
            xend = "end",
            yend = "group"
          ),
          arrow = arrow(
            ends = ifelse(data$strand == "+", "last", "first"),
            angle = arrow_angle,
            length = unit(arrow_size, "points"),
            type = "closed"
          ),
          lineend = "butt",
          linejoin = "mitre",
          show.legend = FALSE,
          linewidth = line_width,
          color = grDevices::grey(1, alpha = 0.5)
        )
      }
    }
  }
