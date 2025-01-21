#' Create Mass Spectrometry Protein Coverage Plot.
#'
#' @param coverage.df Protein coverage, for example output from Proteome Discoverer.
#' @param fasta.file Input reference protein fasta file.
#' @param protein.id The protein ID of exported coverage file. This should be unique and in \code{fasta.file}.
#' @param XCorr.threshold The cross-correlation threshold. Default: 2.
#' @param confidence The confidence level. Default: High.
#' @param contaminant Whether to remove contaminant peptides. Default: NULL (not remove).
#' @param remove.na Logical value, whether to remove NA value in Abundance column. Default: TRUE.
#' @param color The fill color of coverage plot. Default: grey.
#' @param mark.bare Logical value, whether to mark region where Abundance is zero or NA. Default: TRUE.
#' @param mark.color The color used for the marked region. Default: red.
#' @param mark.alpha The transparency used for the marked region. Default: 0.5.
#' @param show.table Logical value, whether to show coverage summary table. Default: TRUE.
#' @param table.position The position of the coverage summary table, choose from right_top, left_top, left_bottom, right_bottom.
#' Default: right_top.
#' @param table.size The font size of coverage summary table. Default: 4.
#' @param table.color The font color of coverage summary table. Default: black.
#' @param range.size The label size of range text, used when \code{range.position} is in. Default: 3.
#' @param range.position The position of y axis range, chosen from in (move y axis in the plot) and
#' out (normal y axis). Default: in.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' \donttest{
#'   if (requireNamespace("openxlsx", quietly = TRUE)) {
#'     library(ggplot2)
#'     library(openxlsx)
#'
#'     # import coverage dataframe with function from openxlsx
#'     coverage.file <- system.file(
#'       "extdata", "Proteomics", "MS_BSA_coverage.xlsx",
#'       package = "ggcoverage"
#'     )
#'     coverage.df <- read.xlsx(coverage.file)
#'     head(coverage.df)
#'
#'     # get fasta file
#'     fasta.file <- system.file(
#'       "extdata", "Proteomics", "MS_BSA_coverage.fasta",
#'       package = "ggcoverage"
#'     )
#'
#'     protein.id <- "sp|P02769|ALBU_BOVIN"
#'     ggprotein(
#'       coverage.df = coverage.df,
#'       fasta.file = fasta.file,
#'       protein.id = protein.id
#'     )
#'   }
#' }
ggprotein <- function(coverage.df, fasta.file, protein.id, XCorr.threshold = 2,
                      confidence = "High", contaminant = NULL, remove.na = TRUE,
                      color = "grey", mark.bare = TRUE, mark.color = "red", mark.alpha = 0.5,
                      show.table = TRUE, table.position = c("top_right", "top_left", "bottom_left", "bottom_right"),
                      table.size = 10, table.color = "black", range.size = 3, range.position = c("in", "out")) {
  # check parameters
  table.position <- match.arg(arg = table.position)
  range.position <- match.arg(arg = range.position)

  # ms protein plot
  protein.plot <- ggplot() +
    geom_protein(
      coverage.df = coverage.df, fasta.file = fasta.file, protein.id = protein.id,
      XCorr.threshold = XCorr.threshold, confidence = confidence, contaminant = contaminant,
      remove.na = remove.na, color = color, mark.bare = mark.bare, mark.color = mark.color,
      mark.alpha = mark.alpha, show.table = show.table, table.position = table.position,
      table.size = table.size, table.color = table.color, range.size = range.size, range.position = range.position
    )

  # add theme
  if (range.position == "in") {
    protein.plot +
      theme_protein()
  } else if (range.position == "out") {
    protein.plot +
      theme_protein2()
  }
}
