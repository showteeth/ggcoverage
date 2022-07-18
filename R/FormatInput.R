# select df with specify start and end position
# non-overlap region and
GetRegion <- function(df, chr, start, end = NULL) {
  # subset used chromosome
  df <- df[df$seqnames == chr, ] %>% dplyr::arrange(start)
  rownames(df) <- NULL

  # # get start loci
  # strat.index <- max(which(df$start <= start))
  # if (is.null(end)) {
  #   end.index <- nrow(df)
  # } else {
  #   end.index <- min(which(df$end >= end))
  # }
  # # select possible region
  # df.select <- df[strat.index:end.index, ]
  # # filter
  # df.select <- df.select[df.select$end >= start & df.select$start <= end, ]
  df.select <- df[df$end >= start & df$start <= end, ]
  init.start <- df.select[1, "start"]
  if (init.start < start) {
    df.select[1, "start"] <- start
  }
  if (!is.null(end)) {
    final.end <- df.select[nrow(df.select), "end"]
    if (final.end > end) {
      df.select[nrow(df.select), "end"] <- end
    }
  }
  return(df.select)
}

#' Prepare Input for Creating Coverage Plot.
#'
#' @param data Track dataframe loaded by \code{\link{LoadTrackFile}}.
#' @param region Region used to create coverage plot, eg: chr14:21,677,306-21,737,601 or chr14:21,677,306.
#' Default: NULL.
#' @param gtf.gr Granges object of GTF, created with \code{\link{import.gff}}. Default: NULL.
#' @param gene.name The name of gene. Default: HNRNPC.
#' @param gene.name.type Gene name type (filed of \code{gtf.gr}), chosen from gene_name and gene_id.
#' Default: gene_name.
#' @param extend Extend length of \code{region}. Default: 2000.
#'
#' @return A dataframe.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange
#'
#' @export
#'
FormatTrack <- function(data, region = "chr14:21,677,306-21,737,601",
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
      region.end <- NULL
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
    if (!is.null(region.start)) {
      region.start <- region.start - extend
    }
    # extend end
    if (!is.null(region.start)) {
      region.end <- region.end + extend
    }
  }
  # filter track dataframe
  track.df.chr <- data %>% dplyr::filter(seqnames == region.chr)
  track.chr.split <- split(x = track.df.chr, f = track.df.chr$Type)
  track.used.list <- lapply(track.chr.split, function(x) {
    single.used.track <- GetRegion(x, chr = region.chr, start = region.start, end = region.end)
  })

  # merge dataframe
  track.used.df <- do.call(rbind, track.used.list)
  rownames(track.used.df) <- NULL

  return(track.used.df)
}
