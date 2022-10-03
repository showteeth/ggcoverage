#' Get Consensus Peak from Replicates with MSPC.
#'
#' @param peak.file Peak files (two or more file: get consensus peak; one file: read directly)
#' obtained from peak caller, eg: MACS2 (without header).
#' @param peak.folder The folder contains peak files. Default: NULL.
#' @param mspc.path MSPC path. Default: NULL (conduct automatic detection).
#' @param rep.type Replicate type, chosen from bio (biological) and tec (technical).
#' @param stringency.threshold Set a threshold on p-values, where peaks with p-value lower than this threshold, are considered stringent. Default: 1e-8.
#' @param weak.threshold Set a threshold on p-values, such that peaks with p-value between this and stringency threshold, are considered weak peaks. Default: 1e-4.
#' @param gamma Set the combined stringency threshold. Peaks with combined p-value below this threshold are confirmed. Default: 1e-8.
#' @param alpha Set the threshold for Benjamini-Hochberg multiple testing correction. Default: 0.05.
#' @param min.overlap.num Set the minimum number of overlapping peaks required before MSPC combines their p-value. Default: 1.
#' @param multiple.intersections When multiple peaks from a sample overlap with a given peak,
#' this argument defines which of the peaks to be considered: the one with lowest p-value, or the one with highest p-value?
#' Chosen from Lowest and Highest. Default: Lowest.
#' @param parallelism.degree Set the number of parallel threads MSPC can utilize simultaneously when processing data. Default: 1.
#'
#' @return A dataframe contains all consensus peak.
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' # library(ggcoverage)
#' # peak.file <- system.file("extdata", "ChIP-seq", "consensus.peak", package = "ggcoverage")
#' # peak.df <- GetConsensusPeak(peak.file = peak.file)
GetConsensusPeak <- function(peak.file, peak.folder = NULL, mspc.path = NULL, rep.type = c("bio", "tec"), stringency.threshold = 1e-8,
                             weak.threshold = 1e-4, gamma = 1e-8, alpha = 0.05, min.overlap.num = 1,
                             multiple.intersections = c("Lowest", "Highest"), parallelism.degree = 1) {

  # check parameters
  rep.type <- match.arg(arg = rep.type)
  multiple.intersections <- match.arg(arg = multiple.intersections)

  # prepare input files
  if (!is.null(peak.folder)) {
    peak.file <- list.files(path = peak.folder, full.names = TRUE)
  }
  # check peak file number
  if (length(peak.file) < 1) {
    stop("Peak file number is less than or equal to one!")
  } else if (length(peak.file) == 1) {
    # read file directly, do not get consensus peaks
    consensus.peak.df <- read.table(file = peak.file, sep = "\t", header = FALSE)
    consensus.peak.df <- consensus.peak.df[, 1:5]
    colnames(consensus.peak.df) <- c("chr", "start", "stop", "name", "score")
  } else {
    # peak file para
    input.para <- paste0("-i ", paste(peak.file, collapse = " -i "))

    # get mspc path
    if (is.null(mspc.path)) {
      # detect MSPC path
      mspc.path <- Sys.which("mspc")
      if (mspc.path == "") {
        stop("Can not find MSPC automatically, please specify the path!")
      }
    } else {
      mspc.path <- mspc.path
    }

    # get tmp folder
    out.folder <- tempdir()

    # full command
    mspc.cmd <- paste(
      mspc.path, input.para, "-r", rep.type, "-s", stringency.threshold,
      "-w", weak.threshold, "-g", gamma, "-a", alpha, "-c", min.overlap.num, "-m", multiple.intersections,
      "-d", parallelism.degree, "-o", out.folder
    )
    # change language information
    full.mspc.cmd <- paste0("export LC_ALL=en_US.UTF-8;", mspc.cmd)
    # run command
    message(paste("Calling MSPC: ", mspc.cmd))
    mspc.status <- system(full.mspc.cmd, intern = TRUE)
    mspc.status.code <- attr(mspc.status, "status")
    if (!is.null(mspc.status.code)) {
      stop("Run MSPC error!")
    }
    # obtain results
    if (!file.exists(file.path(out.folder, "ConsensusPeaks.bed"))) {
      out.base <- basename(out.folder)
      all.tmp.dirs <- sort(dir(path = dirname(out.folder), pattern = out.base, full.names = TRUE))
      out.folder <- all.tmp.dirs[length(all.tmp.dirs)]
    }
    consensus.peak.file <- file.path(out.folder, "ConsensusPeaks.bed")
    consensus.peak.df <- read.table(file = consensus.peak.file, sep = "\t", header = TRUE)
  }
  return(consensus.peak.df)
}
