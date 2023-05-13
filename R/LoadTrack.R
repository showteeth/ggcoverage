#' Load Track File to Dataframe.
#'
#' @param track.file Track file, when \code{track.folder} is not NULL, determined by \code{track.folder}.
#' @param track.folder Track file folder. Default: NULL.
#' @param format Track file format, chosen from bam, wig, bw(bigwig), bedgraph(bedGraph) and txt.
#' @param region Region used to create coverage plot, eg: chr14:21,677,306-21,737,601 or chr14:21,677,306.
#' Default: "chr14:21,677,306-21,737,601"
#' @param extend Extend length of \code{region}. Default: 2000.
#' @param gtf.gr Granges object of GTF, created with \code{\link{import.gff}}. Default: NULL.
#' @param gene.name The name of gene. Default: HNRNPC.
#' @param gene.name.type Gene name type (filed of \code{gtf.gr}), chosen from gene_name and gene_id.
#' Default: gene_name.
#' @param meta.info Track file metadata. The columns should be: SampleName (\code{track.file} without suffix),
#' Type (sample with replicates information), Group (sample group). when \code{meta.file} is not NULL,
#' determined by \code{meta.file}.Default: NULL.
#' @param meta.file File contains track file metadata. Default: "".
#' @param bamcoverage.path The path to \code{bamCoverage}, used when \code{format} is bam. Default: NULL (auto-detect).
#' @param norm.method Methods to normalize the number of reads per bin, chosen from "RPKM", "CPM", "BPM", "RPGC", "None".
#' Default: RPKM.
#' @param single.nuc Logical value, whether to visualize at single nucleotide level. Default: FALSE.
#' @param single.nuc.region Region for \code{single.nuc}. Default: NULL
#' @param bin.size Size of the bins, in bases. Default: 50.
#' @param bc.extra.para Extra parameters for \code{bamCoverage}, eg: "--effectiveGenomeSize 2700000000 --ignoreForNormalization chrX"
#' @param n.cores The number of cores to be used for this job. Default:1.
#'
#' @return A dataframe.
#' @importFrom rtracklayer import
#' @importFrom Rsamtools indexBam
#' @importFrom utils read.csv
#' @importFrom GenomicAlignments alphabetFrequencyFromBam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom magrittr %>%
#' @importFrom dplyr select filter
#' @importFrom BiocParallel register MulticoreParam bplapply
#' @export
#'
#' @examples
#' library(ggcoverage)
#' library(utils)
#' meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
#' sample.meta <- utils::read.csv(meta.file)
#' # track folder
#' track.folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
#' # load bigwig file
#' track.df <- LoadTrackFile(
#'   track.folder = track.folder, format = "bw", region = "chr14:21,677,306-21,737,601",
#'   extend = 2000, meta.info = sample.meta
#' )
LoadTrackFile <- function(track.file, track.folder = NULL, format = c("bam", "wig", "bw", "bedgraph", "txt"),
                          region = "chr14:21,677,306-21,737,601", extend = 2000,
                          gtf.gr = NULL, gene.name = "HNRNPC", gene.name.type = c("gene_name", "gene_id"),
                          meta.info = NULL, meta.file = "",
                          bamcoverage.path = NULL, norm.method = c("RPKM", "CPM", "BPM", "RPGC", "None"),
                          single.nuc = FALSE, single.nuc.region = NULL, bin.size = 10, bc.extra.para = NULL, n.cores = 1) {
  # check parameters
  format <- match.arg(arg = format)
  gene.name.type <- match.arg(arg = gene.name.type)
  norm.method <- match.arg(arg = norm.method)

  # prepare track files
  if (!is.null(track.folder)) {
    track.file <- list.files(path = track.folder, full.names = TRUE, pattern = paste0(format, "$"))
  }

  # get track dataframe
  if (format %in% c("wig", "bw", "bedgraph")) {
    if (single.nuc) {
      stop("To visualize single nucleotide, please use bam file!")
    } else {
      # get used gr
      gr <- PrepareRegion(region = region, gtf.gr = gtf.gr, gene.name = gene.name, gene.name.type = gene.name.type, extend = extend)
      if (is.null(n.cores) || n.cores == 1) {
        # read track file
        track.list <- lapply(track.file, function(x) {
          # get basename
          track.file.base <- basename(x)
          # import wig, bigwig and bedgraph file
          single.track.df <- as.data.frame(rtracklayer::import(x, which = gr))
          single.track.df$TrackFile <- track.file.base
          return(single.track.df)
        })
      } else {
        # register
        BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
        # read track file
        track.list <- BiocParallel::bplapply(track.file, BPPARAM = BiocParallel::MulticoreParam(), FUN = function(x) {
          # get basename
          track.file.base <- basename(x)
          # import wig, bigwig and bedgraph file
          single.track.df <- as.data.frame(rtracklayer::import(x, which = gr))
          single.track.df$TrackFile <- track.file.base
          return(single.track.df)
        })
      }
    }
  } else if (format == "bam") {
    # create index
    if (is.null(n.cores) || n.cores == 1) {
      for (bam in track.file) {
        bam.index.file <- paste(bam, "bai", sep = ".")
        if (!file.exists(bam.index.file)) {
          message("Create index file for: ", basename(bam))
          Rsamtools::indexBam(bam)
        }
      }
    } else {
      # register
      BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
      index.flag <- BiocParallel::bplapply(track.file, BPPARAM = BiocParallel::MulticoreParam(), FUN = function(x) {
        bam.index.file <- paste(x, "bai", sep = ".")
        if (!file.exists(bam.index.file)) {
          message("Create index file for: ", basename(x))
          Rsamtools::indexBam(x)
        }
      })
    }
    if (single.nuc) {
      if (!is.null(single.nuc.region)) {
        single.nuc.region <- gsub(pattern = ",", replacement = "", x = single.nuc.region)
        single.nuc.region.chr <- unlist(strsplit(x = single.nuc.region, split = ":"))[1]
        single.nuc.region.se <- unlist(strsplit(x = single.nuc.region, split = ":"))[2]
        single.nuc.region.start <- unlist(strsplit(x = single.nuc.region.se, split = "-"))[1]
        single.nuc.region.end <- unlist(strsplit(x = single.nuc.region.se, split = "-"))[2]
        # load
        if (is.null(n.cores) || n.cores == 1) {
          track.list <- lapply(track.file, function(x) {
            single.track.df <- GenomicAlignments::alphabetFrequencyFromBam(x, param = single.nuc.region, baseOnly = TRUE) %>% as.data.frame()
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
          })
        } else {
          # register
          BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
          track.list <- BiocParallel::bplapply(track.file, BPPARAM = BiocParallel::MulticoreParam(), FUN = function(x) {
            single.track.df <- GenomicAlignments::alphabetFrequencyFromBam(x, param = single.nuc.region, baseOnly = TRUE) %>% as.data.frame()
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
          })
        }
      } else {
        stop("Please provide region for visualizing single nucleotide!")
      }
    } else {
      # require deeptools
      if (is.null(bamcoverage.path)) {
        bamcoverage.path <- Sys.which("bamCoverage")
        if (bamcoverage.path == "") {
          stop("Can not find bamCoverage automatically, please specify the path!")
        }
      } else {
        bamcoverage.path <- bamcoverage.path
      }

      # get used gr
      gr <- PrepareRegion(region = region, gtf.gr = gtf.gr, gene.name = gene.name, gene.name.type = gene.name.type, extend = extend)
      # read track file
      if (is.null(n.cores) || n.cores == 1) {
        track.list <- lapply(track.file, function(x) {
          # get basename
          track.file.base <- basename(x)
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
            stop("Run bamCoverage error!")
          }
          # import wig, bigwig and bedgraph file
          single.track.df <- as.data.frame(rtracklayer::import(out.bw.file, which = gr))
          single.track.df$TrackFile <- track.file.base
          return(single.track.df)
        })
      } else {
        # register
        BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
        track.list <- BiocParallel::bplapply(track.file, BPPARAM = BiocParallel::MulticoreParam(), FUN = function(x) {
          # get basename
          track.file.base <- basename(x)
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
            stop("Run bamCoverage error!")
          }
          # import wig, bigwig and bedgraph file
          single.track.df <- as.data.frame(rtracklayer::import(out.bw.file, which = gr))
          single.track.df$TrackFile <- track.file.base
          return(single.track.df)
        })
      }
    }
  } else if (format == "txt") {
    if (single.nuc) {
      stop("To visualize single nucleotide, please use bam file!")
    } else {
      # read track file
      if (is.null(n.cores) || n.cores == 1) {
        track.list <- lapply(track.file, function(x) {
          # get basename
          track.file.base <- basename(x)
          # import wig, bigwig and bedgraph file
          single.track.df <- utils::read.table(track.file, header = TRUE)
          single.track.df$TrackFile <- track.file.base
          return(single.track.df)
        })
      } else {
        # register
        BiocParallel::register(BiocParallel::MulticoreParam(workers = n.cores), default = TRUE)
        track.list <- BiocParallel::bplapply(track.file, BPPARAM = BiocParallel::MulticoreParam(), FUN = function(x) {
          # get basename
          track.file.base <- basename(x)
          # import wig, bigwig and bedgraph file
          single.track.df <- utils::read.table(track.file, header = TRUE)
          single.track.df$TrackFile <- track.file.base
          return(single.track.df)
        })
      }
    }
  }
  # get track dataframe
  track.df <- do.call(rbind, track.list)
  # remove width and strand
  if (all(c("width", "strand") %in% colnames(track.df))) {
    track.df <- track.df %>% dplyr::select(-c(width, strand))
  }

  # get metadata
  if (file.exists(meta.file)) {
    meta.info.used <- utils::read.csv(meta.file)
  } else if (!is.null(meta.info)) {
    meta.info.used <- meta.info
  } else {
    message("Sample without metadata!")
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

  return(track.df)
}
