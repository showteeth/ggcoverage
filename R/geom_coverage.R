#' Layer for Coverage Plot.
#'
#' @param data Track prepared by \code{\link{FormatTrack}}.
#' @param mapping Set of aesthetic mappings created by \code{aes} or \code{aes_}. Default: NULL.
#' @param color Track color. Default: NULL (select automatically).
#' @param rect.color The color of every bin. Default: NA.
#' @param facet.key Sample type key to create coverage plot. Default: Type.
#' @param facet.order The order of coverage plot. Default: NULL.
#' @param facet.color The color of sample text. Default: NULL (select automatically).
#' @param group.key Group of samples. Default: NULL.
#' @param range.size The label size of range text, used when \code{range.position} is in. Default: 3.
#' @param range.position The position of y axis range, chosen from in (move y axis in the plot) and
#' out (normal y axis). Default: in.
#' @param mark.region Mark region on the plot. Default: NULL.
#' @param mark.color The color of marked region. Default: "grey".
#' @param mark.alpha The alpha of marked region. Default: 0.5.
#' @param show.mark.label Logical value, whether to show mark label (use label column in \code{mark.region}). Default: TRUE.
#' @param mark.label.size The label size of mark label. Default: 4.
#'
#' @return Layers of ggplot2.
#' @importFrom ggplot2 aes_string scale_fill_manual geom_rect geom_text aes
#' @importFrom rlang .data
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang as_label
#' @importFrom stats as.formula
#' @importFrom ggh4x facet_wrap2 strip_themed elem_list_rect
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#' @examples
#' library(ggcoverage)
#' library(utils)
#' library(ggplot2)
#' meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
#' sample.meta <- utils::read.csv(meta.file)
#' # track folder
#' track.folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
#' # load bigwig file
#' track.df <- LoadTrackFile(
#'   track.folder = track.folder, format = "bw",
#'   meta.info = sample.meta
#' )
#' ggplot() +
#'   geom_coverage(data = track.df, color = "auto", mark.region = NULL)
geom_coverage <- function(data, mapping = NULL,
                          color = NULL, rect.color = NA, facet.key = "Type", facet.order = NULL, facet.color = NULL,
                          group.key = "Group", range.size = 3, range.position = c("in", "out"),
                          mark.region = NULL, mark.color = "grey", mark.alpha = 0.5,
                          show.mark.label = TRUE, mark.label.size = 4) {
  # check parameters
  range.position <- match.arg(arg = range.position)

  # get mapping and color
  if (is.null(mapping)) {
    if (!is.null(color)) {
      mapping <- aes_string(xmin = "start", xmax = "end", ymin = "0", ymax = "score", fill = "Type")
      if (length(color) != length(unique(data$Type))) {
        warning("The color you provided is not as long as Type column in data, select automatically!")
        # sample group with same color
        tmp.color <- AutoColor(data = data, n = 9, name = "Set1", key = group.key)
        # change Type color
        fill.color.df <- merge(unique(data[c("Type", group.key)]), data.frame(color = tmp.color), by.x = group.key, by.y = 0)
        fill.color <- fill.color.df$color
        names(fill.color) <- fill.color.df$Type
      } else {
        fill.color <- color
        if (is.null(names(fill.color))) {
          names(fill.color) <- unique(data$Type)
        }
      }
      sacle_fill_cols <- scale_fill_manual(values = fill.color)
    } else {
      mapping <- aes_string(xmin = "start", xmax = "end", ymin = "0", ymax = "score")
      sacle_fill_cols <- NULL
    }
  } else {
    if ("fill" %in% names(mapping)) {
      fill.str <- rlang::as_label(mapping$fill)
      fill.str.len <- length(unique(data[, fill.str]))
      if (is.null(color) | length(color) != fill.str.len) {
        # sample group with same color
        tmp.color <- AutoColor(data = data, n = 9, name = "Set1", key = group.key)
        # change Type color
        fill.color.df <- merge(unique(data[c(fill.str, group.key)]), data.frame(color = tmp.color), by.x = group.key, by.y = 0)
        fill.color <- fill.color.df$color
        names(fill.color) <- fill.color.df[, fill.str]
      } else {
        fill.color <- color
        if (is.null(names(fill.color))) {
          names(fill.color) <- unique(data[, fill.str])
        }
      }
      sacle_fill_cols <- scale_fill_manual(values = fill.color)
    } else {
      sacle_fill_cols <- NULL
    }
  }

  # facet order
  if (is.null(facet.order)) {
    facet.order <- unique(data[, facet.key])
  }
  data[, facet.key] <- factor(data[, facet.key], levels = facet.order)

  # facet color
  if (is.null(facet.color)) {
    facet.color <- AutoColor(data = data, n = 12, name = "Set3", key = facet.key)
  }

  # facet formula
  facet.formula <- as.formula(paste0("~ ", facet.key))
  region.rect <- geom_rect(data = data, mapping = mapping, show.legend = F, colour = rect.color)
  region.facet <- facet_wrap2(
    facets = facet.formula, ncol = 1, scales = "free_y", strip.position = "right",
    strip = strip_themed(background_y = elem_list_rect(
      fill = facet.color
    ))
  )
  plot.ele <- list(region.rect, region.facet)

  if (!is.null(sacle_fill_cols)) {
    plot.ele <- append(plot.ele, sacle_fill_cols)
  }

  # add range
  ymax.str <- rlang::as_label(mapping$ymax)
  if (range.position == "in") {
    data.range <- data %>%
      dplyr::group_by(.data[[facet.key]]) %>%
      dplyr::summarise(max_score = max(.data[[ymax.str]]))
    data.range$max_score <- sapply(data.range$max_score, CeilingNumber)
    data.range$label <- paste0("[0, ", data.range$max_score, "]")
    region.range <- geom_text(
      data = data.range,
      mapping = aes(x = -Inf, y = Inf, label = label),
      hjust = 0,
      vjust = 1.5,
      size = range.size
    )
    plot.ele <- append(plot.ele, region.range)
  }

  # add rect
  if (!is.null(mark.region)) {
    # get valid mark region
    region.start <- data[1, "start"]
    region.end <- data[nrow(data), "end"]
    valid.region.list <- list()
    for (r in 1:nrow(mark.region)) {
      if (mark.region[r, "start"] <= region.end & mark.region[r, "end"] >= region.start) {
        if (mark.region[r, "end"] >= region.end) {
          mark.region[r, "end"] <- region.end
        }
        if (mark.region[r, "start"] <= region.start) {
          mark.region[r, "start"] <- region.start
        }
        valid.region.list[[r]] <- mark.region[r, ]
      }
    }
    valid.region.df <- do.call(rbind, valid.region.list) %>% as.data.frame()
    colnames(valid.region.df) <- colnames(mark.region)

    region.mark <- geom_rect(
      data = valid.region.df,
      aes_string(xmin = "start", xmax = "end", ymin = "0", ymax = "Inf"),
      fill = mark.color, alpha = mark.alpha
    )
    plot.ele <- append(plot.ele, region.mark)
    # add rect label
    if (show.mark.label) {
      if ("label" %in% colnames(valid.region.df)) {
        # create mark region label
        region.label <- valid.region.df
        region.label[, facet.key] <- facet.order[1]
        region.mark.label <- geom_text_repel(
          data = region.label,
          mapping = aes(x = end, y = Inf, label = label),
          hjust = -0.3,
          vjust = 0.5,
          size = mark.label.size
        )
        plot.ele <- append(plot.ele, region.mark.label)
      } else {
        warning("label is not in provided mark.region!")
      }
    }
  }
  return(plot.ele)
}
