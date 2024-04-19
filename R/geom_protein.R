#' Layer for Protein Coverage Plot.
#'
#' @param coverage.file Exported protein coverage file, should be in excel.
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
#' @param table.size The font size of coverage summary table. Default: 12.
#' @param table.color The font color of coverage summary table. Default: black.
#' @param range.size The label size of range text, used when \code{range.position} is in. Default: 3.
#' @param range.position The position of y axis range, chosen from in (move y axis in the plot) and
#' out (normal y axis). Default: in.
#'
#' @return A ggplot2 object.
#' @importFrom openxlsx read.xlsx
#' @importFrom magrittr %>%
#' @importFrom dplyr filter group_by summarise arrange
#' @importFrom rlang .data
#' @importFrom Biostrings readAAStringSet
#' @importFrom stringr str_locate
#' @importFrom GenomicRanges reduce GRanges setdiff
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_rect geom_text aes aes_string
#' @importFrom scales scientific
#' @importFrom gridExtra ttheme_default tableGrob
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(ggcoverage)
#' coverage.file <- system.file(
#'   "extdata", "Proteomics", "MS_BSA_coverage.xlsx", package = "ggcoverage"
#' )
#' fasta.file <- system.file(
#'   "extdata", "Proteomics", "MS_BSA_coverage.fasta", package = "ggcoverage"
#' )
#' protein.id = "sp|P02769|ALBU_BOVIN"
#' ggplot() +
#'   geom_protein(coverage.file = coverage.file, fasta.file = fasta.file, protein.id = protein.id)
geom_protein <- function(coverage.file, fasta.file, protein.id, XCorr.threshold = 2,
                         confidence = "High", contaminant = NULL, remove.na = TRUE,
                         color = "grey", mark.bare = TRUE, mark.color = "red", mark.alpha = 0.5,
                         show.table = TRUE, table.position = c("top_right", "top_left", "bottom_right", "bottom_left"),
                         table.size = 4, table.color = "black", range.size = 3, range.position = c("in", "out")) {
  # check parameters
  table.position <- match.arg(arg = table.position)
  range.position <- match.arg(arg = range.position)

  # load coverage dataframe
  coverage.df <- openxlsx::read.xlsx(coverage.file)
  # remove suffix and prefix string
  coverage.df$Annotated.Sequence <- gsub(pattern = ".*\\.(.*)\\..*", replacement = "\\1", x = coverage.df$Annotated.Sequence)
  # filter converge according to confidence
  if (!is.null(confidence)) {
    coverage.df <- coverage.df[coverage.df[, "Confidence"] == confidence, ]
  }
  # filter converge according to contaminant
  if (!is.null(contaminant)) {
    coverage.df <- coverage.df[coverage.df[, "Contaminant"] == contaminant, ]
  }
  # filter converge according to cross-correlation
  if (!is.null(XCorr.threshold)) {
    xcorr.index <- grep(pattern = "XCorr", x = colnames(coverage.df))
    coverage.df <- coverage.df[coverage.df[, xcorr.index] >= XCorr.threshold, ]
  }
  # get abundance cols
  abundance.col <- grep(pattern = "Abundance", x = colnames(coverage.df), value = TRUE)
  # remove na abundance
  if (remove.na) {
    coverage.df <- coverage.df %>% dplyr::filter(!is.na(.data[[abundance.col]]))
  }
  # sum abundance of duplicated Annotated.Sequence
  coverage.df <- coverage.df %>%
    dplyr::group_by(.data[["Annotated.Sequence"]]) %>%
    dplyr::summarise(Abundance = sum(.data[[abundance.col]])) %>%
    as.data.frame()
  colnames(coverage.df) <- c("peptide", "abundance")
  # check the coverage dataframe
  if (nrow(coverage.df) == 0) {
    stop("There is no valid peptide, please check!")
  }

  # load genome fasta
  aa.set <- Biostrings::readAAStringSet(fasta.file)
  protein.index <- which(names(aa.set) == protein.id)
  if (length(protein.index) == 1) {
    aa.set.used <- aa.set[protein.index]
    aa.seq.used <- paste(aa.set.used)
  } else if (length(protein.index) > 1) {
    stop("Please check the protein.id you provided, there is more than one in provided fasta file!")
  } else {
    stop("Please check the protein.id you provided, it can't be found in provided fasta file!")
  }

  # get the region
  aa.anno.region <- sapply(coverage.df$peptide, function(x) {
    stringr::str_locate(pattern = x, aa.seq.used)
  }) %>%
    t() %>%
    as.data.frame()
  colnames(aa.anno.region) <- c("start", "end")

  # merge
  coverage.final <- merge(coverage.df, aa.anno.region, by.x = "peptide", by.y = 0, all.x = TRUE)
  coverage.final <- coverage.final %>% dplyr::arrange(.data[["start"]], .data[["end"]])
  coverage.final$ProteinID <- protein.id

  # get coverage positions
  coverage.pos <-
    GenomicRanges::reduce(GenomicRanges::GRanges(protein.id, IRanges::IRanges(coverage.final$start, coverage.final$end))) %>%
    as.data.frame()
  coverage.pos$strand <- NULL
  colnames(coverage.pos) <- c("ProteinID", "start", "end", "width")
  coverage.pos$Type <- "covered"
  # get coverage rate
  coverage.rate <- round(sum(coverage.pos$width) * 100 / nchar(aa.seq.used), 2)
  # non-cover position
  non.coverage.pos <-
    GenomicRanges::setdiff(
      GenomicRanges::GRanges(protein.id, IRanges::IRanges(1, nchar(aa.seq.used))),
      GenomicRanges::GRanges(protein.id, IRanges::IRanges(coverage.final$start, coverage.final$end))
    ) %>%
    as.data.frame()
  non.coverage.pos$strand <- NULL
  colnames(non.coverage.pos) <- c("ProteinID", "start", "end", "width")
  non.coverage.pos$Type <- "bare"
  # coverage summary
  coverage.summary <- rbind(coverage.pos, non.coverage.pos) %>% as.data.frame()

  # create whole range
  non.coverage.final <- data.frame(
    ProteinID = protein.id, peptide = "Empty", abundance = 0,
    start = non.coverage.pos$start, end = non.coverage.pos$end
  )
  coverage.final <- coverage.final[c("ProteinID", "peptide", "abundance", "start", "end")]
  coverage.final <- rbind(coverage.final, non.coverage.final) %>%
    as.data.frame() %>%
    dplyr::arrange(.data[["start"]], .data[["end"]])

  # coverage rect
  coverage.rect <- geom_rect(
    data = coverage.final, mapping = aes_string(
      xmin = "start", xmax = "end",
      ymin = "0", ymax = "abundance"
    ),
    show.legend = FALSE, fill = color
  )
  plot.ele <- list(coverage.rect)
  # mark bare
  if (mark.bare) {
    bare.rect <- geom_rect(
      data = non.coverage.pos, mapping = aes_string(
        xmin = "start", xmax = "end",
        ymin = "0", ymax = "Inf"
      ),
      show.legend = F, fill = mark.color, alpha = mark.alpha
    )
    plot.ele <- append(plot.ele, bare.rect)
  }
  # summary table
  if (show.table) {
    # table position

    if (table.position == "top_left") {
      table_xmin <- 0
      table_xmax <- nchar(aa.seq.used) * 0.5
      table_ymin <- max(coverage.final[, "abundance"]) * 0.5
      table_ymax <- max(coverage.final[, "abundance"]) * 1.0
    } else if (table.position == "top_right") {
      table_xmin <- nchar(aa.seq.used) * 0.5
      table_xmax <- nchar(aa.seq.used) * 1.0
      table_ymin <- max(coverage.final[, "abundance"]) * 0.5
      table_ymax <- max(coverage.final[, "abundance"]) * 1.0
    } else if (table.position == "bottom_left") {
      table_xmin <- 0
      table_xmax <- nchar(aa.seq.used) * 0.5
      table_ymin <- 0
      table_ymax <- max(coverage.final[, "abundance"]) * 0.5
    } else if (table.position == "bottom_right") {
      table_xmin <- nchar(aa.seq.used) * 0.5
      table_xmax <- nchar(aa.seq.used) * 1.0
      table_ymin <- 0
      table_ymax <- max(coverage.final[, "abundance"]) * 0.5
    }
    table_theme <- gridExtra::ttheme_default(
      base_size = table.size, base_colour = table.color
    )
    summary.table <- ggplot2::annotation_custom(
      grob = gridExtra::tableGrob(
        d = coverage.summary,
        theme = table_theme),
      xmin = table_xmin, xmax = table_xmax,
      ymin = table_ymin, ymax = table_ymax
    )
    plot.ele <- append(plot.ele, summary.table)
  }
  # range position
  if (range.position == "in") {
    # prepare range
    max.abundance <- max(pretty(coverage.final$abundance))
    abundance.range <- data.frame(label = paste0("[0, ", scales::scientific(max.abundance, digits = 2), "]"))
    range.text <- geom_text(
      data = abundance.range,
      mapping = aes(x = -Inf, y = Inf, label = label),
      hjust = 0,
      vjust = 1.5,
      size = range.size
    )
    plot.ele <- append(plot.ele, range.text)
  }
  return(plot.ele)
}
