#' Add Base and Amino Acid Annotation to Coverage Plot.
#'
#' @param bam.file BAM file.
#' @param fa.file Genome fasta file. Default: NULL.
#' @param bs.fa.seq BSgenome for species. Default: NULL.
#' @param chr.split Split between chromosome name and description in \code{fa.file}. Default: "[[:space:]]".
#' @param nuc.offset Offset of nucleotide to frequency plot. Default: -0.1.
#' @param nuc.size The size of nucleotide text. Default: 4.
#' @param nuc.padding Background padding of nucleotide annotation. Default: 0.05.
#' @param nuc.padding.r Radius of background padding. Default: 0.
#' @param nuc.color Color scheme for nucleotides. Default: "A": "#ff2b08", "C": "#009aff",
#' "G": "#ffb507", "T": "#00bc0d".
#' @param guide.line Nucleotide frequency guide line. Default: NULL (0.5).
#' @param guide.line.color The color of guide line. Default: "red".
#' @param guide.line.type The line type of guide line. Default: "dashed".
#' @param mark.type The mark type to highlight SNV position, choose from twill (add twill to position with SNV),
#' star (add star mark to position with SNV), and highlight (position without SNV is grey). Default: twill.
#' @param star.size The size of star when \code{mark.type} is star. Default: 1.
#' @param show.aa Logical value, whether to show amino acid. Default: TRUE.
#' @param sens Sense to translate: F for forward sense and R for reverse sense.
#' Parameter of \code{\link[Biostrings]{translate}}. Default: F.
#' @param numcode The ncbi genetic code number for translation.
#' Parameter of \code{\link[Biostrings]{translate}}. By default the standard genetic code is used.
#' @param NAstring How to translate amino-acids when there are ambiguous bases in codons.
#' Parameter of \code{\link[Biostrings]{translate}}. Default: X.
#' @param ambiguous If TRUE, ambiguous bases are taken into account so that for instance GGN is
#' translated to Gly in the standard genetic code. Parameter of \code{\link[Biostrings]{translate}}. Default: FALSE.
#' @param aa.color Color scheme for amino acids.
#' @param aa.border.color The border color of amino acid. Default: white.
#' @param aa.size The size of amino acid text. Default: 4.
#' @param aa.margin Top and bottom margin of amino acids. Default: 2.
#' @param aa.height The relative height of amino acid to base frequency plot. Default: 0.4.
#' @param plot.space Top and bottom margin. Default: 2.5.
#' @param plot.height The relative height of base and amino acid annotation to coverage plot. Default: 0.5.
#'
#' @return Plot.
#' @importFrom stats reshape
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomicAlignments alphabetFrequencyFromBam
#' @importFrom dplyr %>%
#' @importFrom Biostrings readDNAStringSet getSeq translate
#' @importFrom ggplot2 ggplot_add ggplot geom_bar geom_label unit aes_string geom_hline labs
#' geom_tile geom_text theme_classic theme element_blank scale_fill_manual
#' element_text element_rect margin scale_x_continuous scale_y_continuous coord_cartesian
#' @importFrom ggpattern geom_col_pattern scale_pattern_fill_manual
#' @export
#'
#' @examples
#' library("BSgenome.Hsapiens.UCSC.hg19")
#'
#' # get sample metadata
#' sample.meta <- data.frame(
#'   SampleName = c("tumorA.chr4.selected"),
#'   Type = c("tumorA"),
#'   Group = c("tumorA")
#' )
#'
#' # get bam file
#' bam.file <-
#'   system.file("extdata", "DNA-seq", "tumorA.chr4.selected.bam", package = "ggcoverage")
#'
#' # load bam file
#' track.df <- LoadTrackFile(
#'   track.file = bam.file,
#'   meta.info = sample.meta,
#'   single.nuc = TRUE,
#'   single.nuc.region = "chr4:62474235-62474295"
#' )
#'
#' # plot
#' ggcoverage(
#'   data = track.df,
#'   color = "grey",
#'   range.position = "out",
#'   single.nuc = TRUE,
#'   rect.color = "white"
#' ) +
#'   geom_base(
#'     bam.file = bam.file,
#'     bs.fa.seq = BSgenome.Hsapiens.UCSC.hg19
#'   )
#'
geom_base <- function(bam.file, fa.file = NULL, bs.fa.seq = NULL, chr.split = "[[:space:]]",
                      nuc.offset = -0.1, nuc.size = 4, nuc.padding = 0.05, nuc.padding.r = 0,
                      nuc.color = c("A" = "#ff2b08", "C" = "#009aff", "G" = "#ffb507", "T" = "#00bc0d"),
                      guide.line = NULL, guide.line.color = "red", guide.line.type = "dashed",
                      mark.type = "twill", star.size = 1,
                      show.aa = TRUE, sens = "F", numcode = 1, NAstring = "X", ambiguous = FALSE,
                      aa.color = c(
                        "D" = "#FF0000", "S" = "#FF2400", "T" = "#E34234", "G" = "#FF8000", "P" = "#F28500",
                        "C" = "#FFFF00", "A" = "#FDFF00", "V" = "#E3FF00", "I" = "#C0FF00", "L" = "#89318C",
                        "M" = "#00FF00", "F" = "#50C878", "Y" = "#30D5C8", "W" = "#00FFFF", "H" = "#0F2CB3",
                        "R" = "#0000FF", "K" = "#4b0082", "N" = "#800080", "Q" = "#FF00FF", "E" = "#8F00FF",
                        "*" = "#FFC0CB"
                      ), aa.border.color = "white", aa.size = 4, aa.margin = 2, aa.height = 0.4,
                      plot.space = 2.5, plot.height = 0.5) {
  structure(
    list(
      bam.file = bam.file, fa.file = fa.file, bs.fa.seq = bs.fa.seq, chr.split = chr.split,
      nuc.offset = nuc.offset, nuc.size = nuc.size, nuc.padding = nuc.padding, nuc.padding.r = nuc.padding.r,
      nuc.color = nuc.color, guide.line = guide.line, guide.line.color = guide.line.color, guide.line.type = guide.line.type,
      mark.type = mark.type, star.size = star.size,
      show.aa = show.aa, sens = sens, numcode = numcode, NAstring = NAstring, ambiguous = ambiguous,
      aa.color = aa.color, aa.border.color = aa.border.color, aa.size = aa.size, aa.margin = aa.margin, aa.height = aa.height,
      plot.space = plot.space, plot.height = plot.height
    ),
    class = "base"
  )
}

#' @export
ggplot_add.base <- function(object, plot, object_name) {
  # get plot data, plot data should contain bins
  if ("patchwork" %in% class(plot)) {
    plot.data <- plot[[1]]$layers[[1]]$data
  } else {
    plot.data <- plot$layers[[1]]$data
  }
  # prepare plot range
  plot.chr <- as.character(plot.data[1, "seqnames"])
  # plot.region.start <- plot.data[1, "start"]
  plot.region.start <- min(plot.data[, "start"])
  # notice that this is start, because the use of geom_bar instead of geom_rect
  # plot.region.end <- plot.data[nrow(plot.data), "start"]
  plot.region.end <- max(plot.data[, "start"])
  region <-
    GenomicRanges::GRanges(
      plot.chr,
      IRanges::IRanges(plot.region.start, plot.region.end)
    )

  # get parameters
  bam.file <- object$bam.file
  fa.file <- object$fa.file
  bs.fa.seq <- object$bs.fa.seq
  chr.split <- object$chr.split
  nuc.offset <- object$nuc.offset
  nuc.size <- object$nuc.size
  nuc.padding <- object$nuc.padding
  nuc.padding.r <- object$nuc.padding.r
  nuc.color <- object$nuc.color
  guide.line <- object$guide.line
  guide.line.color <- object$guide.line.color
  guide.line.type <- object$guide.line.type
  mark.type <- object$mark.type
  star.size <- object$star.size
  show.aa <- object$show.aa
  sens <- object$sens
  numcode <- object$numcode
  NAstring <- object$NAstring
  ambiguous <- object$ambiguous
  aa.color <- object$aa.color
  aa.border.color <- object$aa.border.color
  aa.size <- object$aa.size
  aa.margin <- object$aa.margin
  aa.height <- object$aa.height
  plot.space <- object$plot.space
  plot.height <- object$plot.height

  # get position AGCT frequency
  pos.nuc.freq <-
    GenomicAlignments::alphabetFrequencyFromBam(
      bam.file,
      param = region,
      baseOnly = TRUE
    )
  # filter out others
  pos.nuc.freq <-
    pos.nuc.freq[, c("A", "G", "C", "T")] %>% as.data.frame()
  # add chromosome and position
  pos.nuc.freq$Chr <- plot.chr
  pos.nuc.freq$Pos <- plot.region.start:plot.region.end

  # get region sequence
  if (is.null(bs.fa.seq)) {
    if (is.null(fa.file)) {
      stop("Please provide either fa.seq or fa.file!")
    } else {
      fa.seq <- Biostrings::readDNAStringSet(fa.file)
      # change fasta name
      names(fa.seq) <-
        sapply(strsplit(x = names(fa.seq), split = chr.split), "[", 1)
    }
  } else {
    fa.seq <- bs.fa.seq
  }
  region.seq <- Biostrings::getSeq(fa.seq, region)

  # get reference nuc
  pos.nuc.freq$Ref <- as.character(region.seq) %>%
    strsplit("") %>%
    unlist()
  # add fill column
  pos.nuc.freq$Total <-
    rowSums(pos.nuc.freq[, c("A", "G", "C", "T")])
  pos.nuc.freq$Fill <- ifelse(pos.nuc.freq$Total == 0, 1, 0)
  pos.nuc.freq.long <- reshape(
    pos.nuc.freq,
    varying = c("A", "G", "C", "T"),
    times = c("A", "G", "C", "T"),
    timevar = "Base",
    idvar = "id",
    v.names = "Freq",
    direction = "long",
    sep = "",
    new.row.names = 1:(nrow(pos.nuc.freq) * 4)
  )
  # get position with alt
  alt.pos <- pos.nuc.freq.long %>%
    dplyr::filter(.data$Ref == .data$Base &
      .data$Total != .data$Freq) %>%
    dplyr::pull(.data$Pos) %>%
    unique()
  alt.pos.nuc.freq.long <-
    pos.nuc.freq.long %>% dplyr::filter(.data$Pos %in% c(alt.pos))
  # get position without alt
  ref.pos <- pos.nuc.freq.long %>%
    dplyr::filter(.data$Ref == .data$Base &
      .data$Total == .data$Freq) %>%
    dplyr::pull(.data$Pos) %>%
    unique()
  ref.pos.nuc.freq.long <-
    pos.nuc.freq.long %>% dplyr::filter(.data$Pos %in% c(ref.pos))
  # create label offset
  pos.nuc.freq$Offset <- nuc.offset
  # add guide line
  if (is.null(guide.line)) {
    # use mean as guide line
    guide.line <- 0.5
  }

  # add mark
  if (length(alt.pos) >= 1) {
    if (mark.type %in% c("twill", "star")) {
      # create plot
      base.plot <- ggplot() +
        geom_bar(
          data = pos.nuc.freq.long,
          aes_string(x = "Pos", y = "Freq", fill = "Base"),
          position = "fill",
          stat = "identity",
          color = "white"
        )
      if (mark.type == "twill") {
        base.plot <- base.plot +
          ggpattern::geom_col_pattern(
            data = alt.pos.nuc.freq.long,
            aes_string(
              x = "Pos",
              y = "Freq",
              pattern = "Base",
              fill = "Base",
              pattern_fill = "Base",
              pattern_angle = "Base"
            ),
            position = "fill",
            colour = "black",
            pattern_density = 0.35,
            pattern_key_scale_factor = 1.3
          ) +
          ggpattern::scale_pattern_fill_manual(values = c(nuc.color, "white"))
      } else if (mark.type == "star") {
        base.plot <- base.plot +
          geom_point(
            data = alt.pos.nuc.freq.long,
            aes_string(x = "Pos", y = "1.01"),
            shape = 8,
            show.legend = FALSE,
            size = star.size
          )
      }
    } else if (mark.type == "highlight") {
      base.plot <- ggplot() +
        geom_bar(
          data = ref.pos.nuc.freq.long,
          aes_string(x = "Pos", y = "Freq"),
          position = "fill",
          stat = "identity",
          color = "white",
          fill = "grey"
        ) +
        geom_bar(
          data = alt.pos.nuc.freq.long,
          aes_string(x = "Pos", y = "Freq", fill = "Base"),
          position = "fill",
          stat = "identity",
          color = "white"
        )
    } else {
      stop(
        "The mark.type you provided is not valid, please choose from twill, star, highlight."
      )
    }
  } else {
    message("No SNV detected, do not add mark!")
    # create plot
    base.plot <- ggplot() +
      geom_bar(
        data = pos.nuc.freq.long,
        aes_string(x = "Pos", y = "Freq", fill = "Base"),
        position = "fill",
        stat = "identity",
        color = "white"
      )
  }

  if (show.aa) {
    # translate the three different reading frames
    list_translated <- lapply(c(0, 1, 2), function(reading_frame) {
      region_remain <- (width(region.seq) - reading_frame) %% 3
      region_framed <- Biostrings::subseq(
        x = region.seq,
        start = 1 + reading_frame,
        end = width(region.seq) - region_remain
      )
      region_aa <-
        Biostrings::translate(x = region_framed, no.init.codon = TRUE)
      region_aa_df <- AAPadding(
        len = width(region.seq),
        offset = reading_frame,
        aa.seq = strsplit(as.character(region_aa), "")[[1]]
      )
      region_aa_df$Pos <- pos.nuc.freq$Pos
      region_aa_df$y <- -0.2 * (1 + reading_frame)
      return(region_aa_df)
    })

    # get base plot
    base.plot <- base.plot +
      theme_base2(fill.color = nuc.color) +
      geom_label(
        data = pos.nuc.freq,
        aes_string(
          x = "Pos",
          y = "Offset",
          label = "Ref",
          fill = "Ref"
        ),
        show.legend = FALSE,
        size = nuc.size,
        label.padding = unit(nuc.padding, "lines"),
        label.r = unit(nuc.padding.r, "lines")
      ) +
      labs(y = "Base") +
      geom_hline(
        yintercept = guide.line,
        color = guide.line.color,
        linetype = guide.line.type
      )

    # make final AA plot
    options(digits = nchar(max(list_translated[[1]]$Pos)) + 1)
    aa_plot <- ggplot()
    for (reading_frame in seq_along(list_translated)) {
      region_aa_rect <-
        PrepareRect(df = list_translated[[reading_frame]], y.center = -0.2 * reading_frame)
      aa_plot <- aa_plot +
        geom_rect(
          data = region_aa_rect,
          aes_string(
            xmin = "xmin",
            xmax = "xmax",
            ymin = "ymin",
            ymax = "ymax",
            fill = "aa"
          ),
          color = aa.border.color
        ) +
        geom_text(
          data = list_translated[[reading_frame]],
          aes_string(x = "Pos", y = "y", label = "anno")
        )
    }
    aa_plot <- aa_plot +
      labs(y = "AA") +
      theme_aa(margin.len = aa.margin, fill.color = aa.color)
    final.plot <-
      patchwork::wrap_plots(base.plot,
        aa_plot,
        ncol = 1,
        heights = c(1, aa.height)
      )
  } else {
    # create plot without amino acid
    final.plot <- base.plot +
      theme_base(margin.len = plot.space, fill.color = nuc.color) +
      geom_label(
        data = pos.nuc.freq,
        aes_string(
          x = "Pos",
          y = "Offset",
          label = "Ref",
          fill = "Ref"
        ),
        show.legend = FALSE,
        size = nuc.size,
        label.padding = unit(nuc.padding, "lines"),
        label.r = unit(nuc.padding.r, "lines")
      ) +
      labs(y = "Base") +
      geom_hline(
        yintercept = guide.line,
        color = guide.line.color,
        linetype = guide.line.type
      )
  }

  # assemble plot
  patchwork::wrap_plots(
    plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    final.plot,
    ncol = 1,
    heights = c(1, plot.height)
  )
}
