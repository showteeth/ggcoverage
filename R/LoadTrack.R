#' Load Track File to Dataframe.
#'
#' @param track.file Track file, when \code{track.folder} is not NULL, determined by \code{track.folder}.
#' @param track.folder Track file folder. Default: NULL.
#' @param format Track file format, chosen from bam, wig, bw(bigwig), bedgraph(bedGraph).
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
#'
#' @return A dataframe.
#' @importFrom rtracklayer import
#' @importFrom Rsamtools indexBam
#' @importFrom utils read.csv
#' @importFrom GenomicAlignments alphabetFrequencyFromBam
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export
#'
#' @examples
#' library(ggcoverage)
#' sample.meta <- data.frame(
#'   SampleName = c("Chr18_MCF7_ER_1", "Chr18_MCF7_ER_2", "Chr18_MCF7_ER_3", "Chr18_MCF7_input"),
#'   Type = c("MCF7_ER_1", "MCF7_ER_2", "MCF7_ER_3", "MCF7_input"),
#'   Group = c("IP", "IP", "IP", "Input")
#' )
#' # track folder
#' track.folder <- system.file("extdata", "ChIP-seq", package = "ggcoverage")
#' # load bigwig file
#' track.df <- LoadTrackFile(
#'   track.folder = track.folder, format = "bw",
#'   meta.info = sample.meta
#' )
LoadTrackFile <- function(track.file, track.folder = NULL, format = c("bam", "wig", "bw", "bedgraph"), meta.info = NULL, meta.file = "",
                          bamcoverage.path = NULL, norm.method = c("RPKM", "CPM", "BPM", "RPGC", "None"),
                          single.nuc = FALSE, single.nuc.region = NULL, bin.size = 10, bc.extra.para = NULL) {
  # check parameters
  format <- match.arg(arg = format)
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
      # read track file
      track.list <- lapply(track.file, function(x) {
        # get basename
        track.file.base <- basename(x)
        # import wig, bigwig and bedgraph file
        single.track.df <- as.data.frame(rtracklayer::import(x))
        single.track.df$TrackFile <- track.file.base
        return(single.track.df)
      })
    }
  } else if (format == "bam") {
    # create index
    for (bam in track.file) {
      bam.index.file <- paste(bam, "bai", sep = ".")
      if (!file.exists(bam.index.file)) {
        message("Create index file for: ", basename(bam))
        Rsamtools::indexBam(bam)
      }
    }
    if (single.nuc) {
      if (!is.null(single.nuc.region)) {
        single.nuc.region <- gsub(pattern = ",", replacement = "", x = single.nuc.region)
        single.nuc.region.chr <- unlist(strsplit(x = single.nuc.region, split = ":"))[1]
        single.nuc.region.se <- unlist(strsplit(x = single.nuc.region, split = ":"))[2]
        single.nuc.region.start <- unlist(strsplit(x = single.nuc.region.se, split = "-"))[1]
        single.nuc.region.end <- unlist(strsplit(x = single.nuc.region.se, split = "-"))[2]
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

      # read track file
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
        single.track.df <- as.data.frame(rtracklayer::import(out.bw.file))
        single.track.df$TrackFile <- track.file.base
        return(single.track.df)
      })
    }
  }
  # get track dataframe
  track.df <- do.call(rbind, track.list)
  # remove width and strand
  track.df <- track.df %>% dplyr::select(-c(width, strand))

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

  return(track.df)
}
