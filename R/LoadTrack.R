#' Load Track File to Dataframe.
#'
#' @param track.file Track file, when \code{track.folder} is not NULL, determined by \code{track.folder}.
#' @param track.folder Track file folder. Default: NULL.
#' @param format Track file format, chosen from bam, wig, bw(bigwig), bedgraph(bedGraph) and txt.
#' @param region Region to extract coverage for, eg: chr14:21,677,306-21,737,601 or chr14:21,677,306.
#'   Default: NULL, coverage is extracted from the first annotated chromosome/sequence.
#' @param extend Extend length of \code{region}. Default: 2000.
#' @param gtf.gr Granges object of GTF, created with \code{\link{import.gff}}. Default: NULL.
#' @param gene.name The name of gene. Default: HNRNPC.
#' @param gene.name.type Gene name type (filed of \code{gtf.gr}), chosen from gene_name and gene_id.
#'   Default: gene_name.
#' @param meta.info Track file metadata. The columns should be: SampleName (\code{track.file} without suffix),
#'   Type (sample with replicates information), Group (sample group). when \code{meta.file} is not NULL,
#'   determined by \code{meta.file}.Default: NULL.
#' @param meta.file File contains track file metadata. Default: "".
#' @param bamcoverage.path The path to \code{bamCoverage}, used when \code{format} is bam. Default: NULL (auto-detect).
#' @param norm.method Methods to normalize the number of reads per bin, chosen from "RPKM", "CPM", "BPM", "RPGC", "None".
#'   Default: RPKM.
#' @param single.nuc Logical value, whether to visualize at single nucleotide level. Default: FALSE.
#' @param single.nuc.region Region for \code{single.nuc}. Default: NULL
#' @param bin.size Size of the bins, in bases. Default: 10. Only used for BAM files, ignored for Wig, Bigwig, etc.
#'   Set to NULL to turn binning off.
#' @param bc.extra.para Extra parameters for \code{bamCoverage}, eg: "--effectiveGenomeSize 2700000000 --ignoreForNormalization chrX"
#' @param n.cores The number of cores to be used for this job. Default: 1.
#'
#' @return A dataframe.
#' @importFrom rtracklayer import
#' @importFrom Rsamtools indexBam ScanBamParam
#' @importFrom utils read.csv
#' @importFrom GenomicAlignments alphabetFrequencyFromBam readGAlignments coverage
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom dplyr %>%
#' @importFrom dplyr select filter mutate all_of group_by summarize
#' @importFrom BiocParallel register MulticoreParam bplapply
#' @importFrom ggplot2 cut_width
#' @export
#'
#' @examples
#' library(ggcoverage)
#' library(utils)
#'
#' meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
#' sample.meta <- utils::read.csv(meta.file)
#'
#' # track folder
#' track.folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
#'
#' # load bigwig file
#' track.df <- LoadTrackFile(
#'   track.folder = track.folder,
#'   format = "bw",
#'   region = "chr14:21,677,306-21,737,601",
#'   extend = 2000,
#'   meta.info = sample.meta
#' )
LoadTrackFile <- function(
    track.file, track.folder = NULL,
    format = c("bam", "wig", "bw", "bedgraph", "txt"),
    region = NULL, extend = 2000,
    gtf.gr = NULL, gene.name = "HNRNPC",
    gene.name.type = c("gene_name", "gene_id"),
    meta.info = NULL, meta.file = "",
    bamcoverage.path = NULL,
    norm.method = c("RPKM", "CPM", "BPM", "RPGC", "None"),
    single.nuc = FALSE, single.nuc.region = NULL,
    bin.size = 10, bc.extra.para = NULL, n.cores = 1) {
  # check parameters
  format <- match.arg(arg = format)
  gene.name.type <- match.arg(arg = gene.name.type)
  norm.method <- match.arg(arg = norm.method)

  # prepare track files
  if (!is.null(track.folder)) {
    track.file <- list.files(path = track.folder, full.names = TRUE, pattern = paste0(format, "$"))
  }

  # get genomic region if supplied, else it is guessed from input
  if (is.null(region)) {
    message("No 'region' specified; extracting coverage for an example range\n(<=100,000 bases, first annotated sequence)")
    if (format == "bam") {
      seqnames <- Rsamtools::scanBamHeader(track.file[1]) %>%
        lapply(function(x) x$targets) %>%
        unname() %>%
        unlist()
      gr <- GenomicRanges::GRanges(
        seqnames = names(seqnames[1]),
        IRanges(start = 1, end = min(100000, seqnames[1]))
      )
    } else if (format %in% c("wig", "bw", "bedgraph")) {
      gr <- range(rtracklayer::import(track.file[1]))
      seqnames <- as.character(seqnames(gr))
      if (GenomicRanges::width(gr) <= 100000) {
        gr <- GenomicRanges::resize(gr, width = 100000)
      }
    }
    message(paste0("Coverage extracted from sequence/chromosome: ", names(seqnames[1])))
  } else {
    gr <- PrepareRegion(
      region = region,
      gtf.gr = gtf.gr,
      gene.name = gene.name,
      gene.name.type = gene.name.type,
      extend = extend
    )
  }

  # get track dataframe
  if (format %in% c("wig", "bw", "bedgraph")) {
    if (single.nuc) {
      stop("To visualize single nucleotide resolution, please use bam file!")
    } else {
      if (is.null(n.cores) || n.cores == 1) {
        track.list <- lapply(track.file, import_bw, gr)
      } else {
        BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
        track.list <- BiocParallel::bplapply(track.file, BPPARAM = BiocParallel::MulticoreParam(), FUN = import_bw, gr)
      }
    }
  } else if (format == "bam") {
    # create index
    if (is.null(n.cores) || n.cores == 1) {
      lapply(track.file, index_bam)
    } else {
      BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
      BiocParallel::bplapply(track.file, BPPARAM = BiocParallel::MulticoreParam(), FUN = index_bam)
    }
    if (single.nuc) {
      if (is.null(n.cores) || n.cores == 1) {
        track.list <- lapply(
          track.file,
          single_nuc_cov,
          single.nuc.region
        )
      } else {
        track.list <- BiocParallel::bplapply(
          track.file,
          BPPARAM = BiocParallel::MulticoreParam(),
          FUN = single_nuc_cov,
          single.nuc.region
        )
      }
    } else {
      if (norm.method == "None") {
        message("Calculating coverage with GenomicAlignments when 'norm.method = None'")
        if (is.null(n.cores) || n.cores == 1) {
          track.list <- lapply(
            track.file, import_bam_ga, gr, bin.size
          )
        } else {
          track.list <- BiocParallel::bplapply(
            track.file,
            BPPARAM = BiocParallel::MulticoreParam(),
            FUN = import_bam_ga, gr, bin.size
          )
        }
      } else {
        message("Calculate coverage with bamCoverage when 'norm.method != None'")
        # require deeptools
        if (is.null(bamcoverage.path)) {
          bamcoverage.path <- Sys.which("bamCoverage")
          if (bamcoverage.path == "") {
            stop("Can not find bamCoverage automatically, please specify 'bamcoverage.path'")
          }
        } else {
          bamcoverage.path <- bamcoverage.path
        }
        if (is.null(n.cores) || n.cores == 1) {
          track.list <- lapply(
            track.file,
            bam_coverage,
            bamcoverage.path,
            bin.size,
            norm.method,
            bc.extra.para, gr
          )
        } else {
          track.list <- BiocParallel::bplapply(
            track.file,
            BPPARAM = BiocParallel::MulticoreParam(),
            FUN = bam_coverage,
            bamcoverage.path,
            bin.size,
            norm.method,
            bc.extra.para, gr
          )
        }
      }
    }
  } else if (format == "txt") {
    if (single.nuc) {
      stop("To visualize single nucleotide, please use bam file!")
    } else {
      # read track file
      if (is.null(n.cores) || n.cores == 1) {
        track.list <- lapply(track.file, import_txt)
      } else {
        BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
        track.list <- BiocParallel::bplapply(track.file, BPPARAM = BiocParallel::MulticoreParam(), FUN = import_txt)
      }
    }
  }
  # get track dataframe
  track.df <- do.call(rbind, track.list)

  # get metadata
  if (file.exists(meta.file)) {
    meta.info.used <- utils::read.csv(meta.file)
  } else if (!is.null(meta.info)) {
    meta.info.used <- meta.info
  } else {
    message("No metadata provided, returning coverage as is.")
    meta.info.used <- NULL
  }

  # merge metadata
  if (is.null(meta.info.used)) {
    # create pseudo-group
    track.df$Type <- track.df$TrackFile
    track.df$Group <- track.df$TrackFile
    track.df$TrackFile <- NULL
  } else {
    meta.info.used$SampleName <- paste(meta.info.used$SampleName, format, sep = ".")
    track.df <- merge(track.df, meta.info.used, by.x = "TrackFile", by.y = "SampleName")
    track.df$TrackFile <- NULL
  }

  # subset txt
  if (format == "txt") {
    track.df <- FormatTrack(
      data = track.df, region = region, gtf.gr = gtf.gr, extend = extend,
      gene.name = gene.name, gene.name.type = gene.name.type
    )
  }
  # return final df
  return(track.df)
}

import_bw <- function(x, gr) {
  single.track.df <- as.data.frame(rtracklayer::import(x, which = gr))
  single.track.df$TrackFile <- basename(x)
  return(single.track.df)
}

import_txt <- function(x) {
  single.track.df <- utils::read.table(x, header = TRUE)
  single.track.df$TrackFile <- basename(x)
  return(single.track.df)
}

import_bam_ga <- function(x, gr, bin.size) {
  # get basename
  track.file.base <- basename(x)
  # load track
  param <- Rsamtools::ScanBamParam(which = gr)
  ga <- GenomicAlignments::readGAlignments(x, param = param)
  ga.cov <- GenomicAlignments::coverage(ga)
  ga.cov.gr <- GenomicRanges::GRanges(ga.cov)
  ga.cov.df <- IRanges::subsetByOverlaps(ga.cov.gr, gr) %>%
    as.data.frame()
  # valid the region
  gr.df <- as.data.frame(gr)
  ga.cov.df[1, "start"] <- gr.df[1, "start"]
  ga.cov.df[nrow(ga.cov.df), "end"] <- gr.df[1, "end"]
  # optional binning
  ga.cov.df <- bin_coverage(ga.cov.df, bin.size)
  # add track file
  ga.cov.df$TrackFile <- track.file.base
  return(ga.cov.df)
}

index_bam <- function(x) {
  bam.index.file <- paste(x, "bai", sep = ".")
  if (!file.exists(bam.index.file)) {
    message("Create index file for: ", basename(x))
    Rsamtools::indexBam(x)
  }
}

single_nuc_cov <- function(x, single.nuc.region) {
  if (!is.null(single.nuc.region)) {
    single.nuc.region <- gsub(pattern = ",", replacement = "", x = single.nuc.region)
    single.nuc.region.chr <- unlist(strsplit(x = single.nuc.region, split = ":"))[1]
    single.nuc.region.se <- unlist(strsplit(x = single.nuc.region, split = ":"))[2]
    single.nuc.region.start <- unlist(strsplit(x = single.nuc.region.se, split = "-"))[1]
    single.nuc.region.end <- unlist(strsplit(x = single.nuc.region.se, split = "-"))[2]
  } else {
    stop("Please provide region for visualizing single nucleotide!")
  }
  single.track.df <- GenomicAlignments::alphabetFrequencyFromBam(x, param = single.nuc.region, baseOnly = TRUE) %>%
    as.data.frame()
  single.track.df <- single.track.df[, c("A", "G", "C", "T")]
  single.track.df$score <- rowSums(single.track.df)
  single.track.df$seqnames <- single.nuc.region.chr
  single.track.df$start <- single.nuc.region.start:single.nuc.region.end
  single.track.df$end <- single.track.df$start + 1
  single.track.df$width <- 1
  single.track.df$strand <- "*"
  single.track.df <- single.track.df %>% dplyr::select(-c("A", "G", "C", "T"))
  # get basename
  track.file.base <- basename(x)
  single.track.df$TrackFile <- track.file.base
  single.track.df <- single.track.df[c(
    "seqnames", "start", "end", "width",
    "strand", "score", "TrackFile"
  )]
  return(single.track.df)
}

bam_coverage <- function(
    x, bamcoverage.path, bin.size, norm.method, bc.extra.para, gr) {
  # bigwig file
  out.bw.file <- tempfile(fileext = c(".bw"))
  # prepare bamCoverage cmd
  bamcoverage.cmd <- paste(
    bamcoverage.path, "-b", x, "-o", out.bw.file,
    "--binSize", bin.size, "--normalizeUsing", norm.method, bc.extra.para
  )
  # run command
  message(paste("Calling bamCoverage: ", bamcoverage.cmd))
  bamcoverage.status <- system(bamcoverage.cmd, intern = TRUE)
  bamcoverage.status.code <- attr(bamcoverage.status, "status")
  if (!is.null(bamcoverage.status.code)) {
    stop("Run bamCoverage error.")
  }
  # import wig, bigwig and bedgraph file
  single.track.df <- as.data.frame(rtracklayer::import(out.bw.file, which = gr))
  single.track.df$TrackFile <- basename(x)
  return(single.track.df)
}

bin_coverage <- function(df, bin.size = 10) {
  if (!is.null(bin.size) && is.numeric(bin.size)) {
    binned_df <- df %>%
      dplyr::mutate(
        bin = ggplot2::cut_width(start, width = bin.size, center = bin.size / 2, labels = FALSE) * bin.size
      ) %>%
      dplyr::group_by(.data$seqnames, .data$bin, .data$strand) %>%
      dplyr::summarize(score = mean(.data$score, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(start = .data$bin - (min(.data$bin) - 1), end = .data$bin, width = bin.size) %>%
      dplyr::select(dplyr::all_of(c("seqnames", "start", "end", "width", "strand", "score"))) %>%
      as.data.frame()
    return(binned_df)
  } else {
    return(df)
  }
}
