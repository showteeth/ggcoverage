#' Add Links to Coverage Plot.
#'
#' @param link.file File contains region link information.
#' @param file.type The type of \code{link.file}, choose from bedpe, pairs. Default: bedpe.
#' @param score.col Column index contains score information, used when \code{file.type} is bedpe. Default: NULL.
#' @param score.threshold The score threshold, used when \code{score.col} is not NULL. Default: NULL.
#' @param score.color The score color vector. The length should be three, the first represents the lowest score, the second represents the
#' middle score, the third represents the highest score. Default: c("blue", "grey", "red").
#' @param scale.range Scale the height of link according to width, should be greater than or equal to 1 (not scale). Default: 10.
#' @param plot.space Top and bottom margin. Default: 0.1.
#' @param plot.height The relative height of link to coverage plot. Default: 0.2.
#' @param show.rect Logical value, whether to add rect border to the plot. Default: FALSE.
#'
#' @return Plot.
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame start end
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom utils read.table
#' @importFrom scales rescale
#' @importFrom ggforce geom_bezier
#' @importFrom ggplot2 ggplot_add ggplot aes_string scale_color_gradient2 labs theme_classic theme element_blank element_rect
#' element_text margin scale_y_continuous scale_x_continuous expansion coord_cartesian
#' @importFrom patchwork wrap_plots
#' @references \url{https://stuartlab.org/signac/articles/cicero.html}
#' @export
#'
#' @examples
#' library(ggcoverage)
#' # create test dataframe
#' # (random, but use seed to obtain same result every time)
#' set.seed(123)
#' df <- data.frame(
#'   seqnames = "chr2L", start = seq(from = 8000000, to = 8300000, by = 1000),
#'   end = seq(from = 8001000, to = 8301000, by = 1000), score = sample(1:100, 301, replace = TRUE),
#'   Type = "Example", Group = "Example"
#' )
#' # get links
#' link.file = system.file("extdata", "HiC", "HiC_link.bedpe", package = "ggcoverage")
#'
#' # create plot
#' ggcoverage(
#'   data = df, color = "grey",
#'   mark.region = NULL, range.position = "out"
#' ) +
#'   geom_link(link.file = link.file, file.type = "bedpe", show.rect = TRUE)
#'
geom_link <- function(link.file, file.type = "bedpe", score.col = NULL, score.threshold = NULL,
                      score.color = c("blue", "grey", "red"), scale.range = 10,
                      plot.space = 0.1, plot.height = 0.2, show.rect = FALSE) {
  structure(list(
    link.file = link.file, file.type = file.type, score.col = score.col, score.threshold = score.threshold,
    score.color = score.color, scale.range = scale.range, plot.space = plot.space, plot.height = plot.height, show.rect = show.rect
  ),
  class = "link"
  )
}


#' @export
ggplot_add.link <- function(object, plot, object_name) {
  # get plot data
  # track.data <- plot$layers[[1]]$data
  # get plot data, plot data should contain bins
  if ("patchwork" %in% class(plot)) {
    track.data <- plot[[1]]$layers[[1]]$data
  } else {
    track.data <- plot$layers[[1]]$data
  }
  # prepare plot range
  # the plot region are not normal, so start is minimum value
  plot.range.chr <- track.data[1, "seqnames"]
  # plot.range.start <- track.data[1, "start"]
  plot.range.start <- min(track.data[, "start"])
  # plot.range.end <- track.data[nrow(track.data), "end"]
  plot.range.end <- max(track.data[, "end"])
  plot.range.gr <- GenomicRanges::GRanges(
    seqnames = plot.range.chr,
    ranges = IRanges::IRanges(plot.range.start, plot.range.end)
  )
  # get parameters
  link.file <- object$link.file
  file.type <- object$file.type
  score.col <- object$score.col
  score.threshold <- object$score.threshold
  score.color <- object$score.color
  scale.range <- object$scale.range
  plot.space <- object$plot.space
  plot.height <- object$plot.height
  show.rect <- object$show.rect

  # check parameter
  if (length(score.color) >= 3) {
    score.color <- score.color[1:3]
  } else {
    warning("The score.color you provided is smaller than 3, please check! Use default (blue, grey, red) first.")
    score.color <- c("blue", "grey", "red")
  }

  # prepare dataframe
  if (file.type == "bedpe") {
    # bedpe: https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    # read bedpe file
    link.df <- utils::read.table(file = link.file, sep = "\t")
    if (!is.null(score.col)) {
      if (score.col > ncol(link.df)) {
        stop("The score column index is bigger than whole dataframe columns, please provide valid score column!")
      } else {
        link.df <- link.df[, c(1, 2, 3, 4, 5, 6, score.col)]
        colnames(link.df) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "score")
      }
    } else {
      link.df <- link.df[, c(1, 2, 3, 4, 5, 6)]
      colnames(link.df) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
    }
    # filter link dataframe
    link.df <- link.df[link.df$chr1 == link.df$chr2, ]
    # detect chrs
    if (length(unique(link.df$chr1)) > 1) {
      warning("The bedpe file provided contains multi-chromosomes, use first!")
      link.df <- link.df[link.df$chr1 == link.df$chr1[1], ]
    }
    # filter threshold
    if ("score" %in% colnames(link.df)) {
      if (!is.null(score.threshold)) {
        link.df <- link.df[link.df$score > score.threshold, ]
      }
    }
    # calculate center
    r1.center <- (link.df$start1 + link.df$end1) / 2
    r2.center <- (link.df$start2 + link.df$end2) / 2
    # change position
    point.start.vec <- ifelse(r1.center < r2.center, r1.center, r2.center)
    point.end.vec <- ifelse(r1.center < r2.center, r2.center, r1.center)
    # create link point dataframe
    link.point.df <- data.frame(
      chr = unique(link.df$chr1),
      start = point.start.vec,
      end = point.end.vec
    )
    # add score
    if ("score" %in% colnames(link.df)) {
      link.point.df$score <- link.df$score
    }
  } else if (file.type == "pairs") {
    # pair format: https://pairtools.readthedocs.io/en/latest/formats.html
    # read pairs file
    link.df <- utils::read.table(file = link.file, sep = "\t")
    # select used dataframe
    link.df <- link.df[, c(2, 3, 4, 5)]
    colnames(link.df) <- c("chr1", "start1", "chr2", "start2")
    # filter link dataframe
    link.df <- link.df[link.df$chr1 == link.df$chr2, ]
    # detect chrs
    if (length(unique(link.df$chr1)) > 1) {
      warning("The bedpe file provided contains multi-chromosomes, use first!")
      link.df <- link.df[link.df$chr1 == link.df$chr1[1], ]
    }
    # change position
    point.start.vec <- ifelse(link.df$start1 < link.df$start2, link.df$start1, link.df$start2)
    point.end.vec <- ifelse(link.df$start1 < link.df$start2, link.df$start2, link.df$start1)
    # create link point dataframe
    link.point.df <- data.frame(
      chr = unique(link.df$chr1),
      start = point.start.vec,
      end = point.end.vec
    )
  }
  # convert link.point.df to genomic ranges
  link.point.gr <- GenomicRanges::makeGRangesFromDataFrame(df = link.point.df, keep.extra.columns = TRUE)
  # filter link gr
  link.point.df <- as.data.frame(IRanges::subsetByOverlaps(x = link.point.gr, ranges = plot.range.gr))
  # remove links outside region
  link.point.df <- link.point.df[link.point.df$start >= GenomicRanges::start(x = plot.range.gr) &
    link.point.df$end <= GenomicRanges::end(x = plot.range.gr), ]
  rownames(link.point.df) <- 1:nrow(link.point.df)
  # check dataframe
  if (nrow(link.point.df) < 1) {
    warning("There is no valid links in given region!")
    # create empty plot
    link.basic.plot <-
      ggplot(data = link.plot.df)
  } else {
    # prepare plot dataframe
    link.point.df$group <- seq_len(length.out = nrow(x = link.point.df))
    link.point.plot <- link.point.df
    link.point.plot$width <- link.point.df$end - link.point.df$start
    # scale width to range
    link.point.plot$rw <- scales::rescale(link.point.plot$width, to = c(1, scale.range))

    # prepare plot dataframe
    link.plot.df <- data.frame(
      x = c(
        link.point.plot$start,
        link.point.plot$start,
        link.point.plot$end,
        link.point.plot$end
      ),
      y = c(
        rep(x = 0, nrow(x = link.point.plot)),
        -link.point.plot$rw,
        -link.point.plot$rw,
        rep(x = 0, nrow(x = link.point.plot))
      ),
      group = rep(x = link.point.plot$group, 4)
    )
    # add score and create basic plot
    if ("score" %in% colnames(link.point.plot)) {
      # add score
      link.plot.df$score <- rep(link.point.plot$score, 4)
      # create plot
      link.basic.plot <-
        ggplot(data = link.plot.df) +
        ggforce::geom_bezier(
          mapping = aes_string(x = "x", y = "y", group = "group", color = "score")
        ) +
        scale_color_gradient2(
          low = score.color[1], mid = score.color[2], high = score.color[3],
          limits = c(min(link.plot.df$score), max(link.plot.df$score)),
          n.breaks = 3
        )
    } else {
      link.basic.plot <-
        ggplot(data = link.plot.df) +
        ggforce::geom_bezier(
          mapping = aes_string(x = "x", y = "y", group = "group")
        )
    }
  }
  # create plot
  link.plot <-
    link.basic.plot +
    labs(y = "Links") +
    theme_link(x.range = c(plot.range.start, plot.range.end), margin.len = plot.space, show.rect = show.rect)
  # assemble plot
  patchwork::wrap_plots(plot + theme(plot.margin = margin(t = plot.space, b = plot.space)),
    link.plot,
    ncol = 1, heights = c(1, plot.height)
  )
}
