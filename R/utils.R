# select color automatically
AutoColor <- function(data, n, name, key) {
  getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n, name))
  # sample group with same color
  group.info <- unique(data[, key])
  fill.color <- getPalette(length(group.info))
  names(fill.color) <- group.info
  return(fill.color)
}

# ceiling for number bigger than zero, floor for number smaller than zero
CeilingNumber <- function(x, digits = 2) {
  # mark number
  if (x < 0) {
    flag <- -1
    x <- abs(x)
  } else {
    flag <- 1
  }
  # transfrom
  if (x > 1) {
    x.ceiling <- round(x + 5 * 10^(-digits - 1), digits)
  } else if (x > 0) {
    x.split <- unlist(strsplit(formatC(x, format = "e"), "e"))
    num.part <- as.numeric(x.split[1])
    sci.part <- as.numeric(x.split[2])
    valid.digits <- digits - 1
    x.ceiling <- round(num.part + 5 * 10^(-valid.digits - 1), valid.digits) * 10^(sci.part)
  }
  # final number
  x.final <- x.ceiling * flag
  return(x.final)
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

# get gene and transcript group
# divide genes to non-overlap groups
# GetGeneGroup = function(gene.gr, fc = "queryHits",sc= "subjectHits", overlap.gene.gap=1){
#   overlap.df = IRanges::findOverlaps(gene.gr,gene.gr, ignore.strand=TRUE) %>% as.data.frame()
#   overlap.list = lapply(split(overlap.df,overlap.df[,fc]), function(x){
#     x[, sc]
#   })
#   group.idx = rep(1, length(unique(overlap.df[,fc])))
#   for (i in seq_along(overlap.list)) {
#     ovrlap.vec <- overlap.list[[i]]
#     curr.group = group.idx[i]
#     remain.vec = setdiff(ovrlap.vec, i)
#     for (j in remain.vec) {
#       if(group.idx[j] == curr.group){
#         group.idx[j] = curr.group + overlap.gene.gap
#       }
#     }
#   }
#   return(group.idx)
# }
GetGeneGroup <- function(gene.gr, fc = "queryHits", sc = "subjectHits", overlap.gene.gap = 1) {
  overlap.df <- IRanges::findOverlaps(gene.gr, gene.gr, ignore.strand = TRUE) %>%
    as.data.frame() %>%
    dplyr::arrange(.data[[fc]], .data[[sc]])
  overlap.list <- lapply(split(overlap.df, overlap.df[, fc]), function(x) {
    x[, sc]
  })
  overlap.list <- overlap.list[order(sapply(overlap.list, length))]
  group.idx <- rep(1, length(unique(overlap.df[, fc])))
  for (i in names(overlap.list)) {
    i <- as.numeric(i)
    ovrlap.vec <- overlap.list[[i]]
    curr.group <- group.idx[i]
    remain.vec <- setdiff(ovrlap.vec, i)
    for (j in remain.vec) {
      if (group.idx[j] == curr.group) {
        group.idx[j] <- curr.group + overlap.gene.gap
      }
    }
  }
  return(group.idx)
}

# SplitExonUTR(exon.df = gene.info.used.exon, utr.df = gene.info.used.utr)
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
# Fix bug: the names on the supplied 'seqlengths' vector must be identical to the seqnames
getIdeogram <- function(genome, subchr = NULL, cytobands = TRUE) {
  .gnm <- genome
  lst <- lapply(.gnm, function(genome) {
    # print(genome)
    ## to remove the "heavy dependency" we put require here.
    # require(rtracklayer)
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
