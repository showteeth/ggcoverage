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
#' @param bin.size Size of the bins, in bases. Default: 50.
#' @param bc.extra.para Extra parameters for \code{bamCoverage}, eg: "--effectiveGenomeSize 2700000000 --ignoreForNormalization chrX"
#'
#' @return A dataframe.
#' @importFrom rtracklayer import
#' @importFrom Rsamtools indexBam
#' @importFrom utils read.csv
#' @export
#'
LoadTrackFile <- function(track.file, track.folder = NULL, format = c("bam", "wig", "bw", "bedgraph"), meta.info = NULL, meta.file = "",
                          bamcoverage.path = NULL, norm.method = c("RPKM", "CPM", "BPM", "RPGC", "None"),
                          bin.size = 10, bc.extra.para = NULL) {
  # check parameters
  format <- match.arg(arg = format)
  norm.method <- match.arg(arg = norm.method)

  # prepare track files
  if (!is.null(track.folder)) {
    track.file <- list.files(path = track.folder, full.names = TRUE, pattern = format)
  }

  # get track dataframe
  if (format %in% c("wig", "bw", "bedgraph")) {
    # read track file
    track.list <- lapply(track.file, function(x) {
      # get basename
      track.file.base <- basename(x)
      # import wig, bigwig and bedgraph file
      single.track.df <- as.data.frame(rtracklayer::import(x))
      single.track.df$TrackFile <- track.file.base
      return(single.track.df)
    })
  } else if (format == "bam") {
    # require deeptools
    if (is.null(bamcoverage.path)) {
      bamcoverage.path <- Sys.which("bamCoverage")
      if (bamcoverage.path == "") {
        stop("Can not find bamCoverage automatically, please specify the path!")
      }
    } else {
      bamcoverage.path <- bamcoverage.path
    }
    # create index
    for (bam in track.file) {
      bam.index.file <- paste(bam, "bai", sep = ".")
      if (!file.exists(bam.index.file)) {
        message("Create index file for: ", basename(bam))
        Rsamtools::indexBam(bam)
      }
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
