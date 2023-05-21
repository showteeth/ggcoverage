#' Layer for Coverage Plot.
#'
#' @param data Track prepared by \code{\link{FormatTrack}}.
#' @param mapping Set of aesthetic mappings created by \code{aes} or \code{aes_}. Default: NULL.
#' @param color Track color. Default: NULL (select automatically).
#' @param rect.color The color of every bin. Default: NA.
#' @param single.nuc Logical value, whether to visualize at single nucleotide level (use bar plot). Default: FALSE.
#' @param plot.type The type of the plot, choose from facet (separate plot for every sample) and
#' joint (combine all sample in a single plot). Default: facet.
#' @param facet.key Sample type key to create coverage plot. Used in both facet and joint plot. Default: Type.
#' @param joint.avg Logical value, whether to show average coverage across \code{group.key}. Default: FALSE.
#' @param facet.order The order of coverage plot. Default: NULL.
#' @param facet.color The color of sample text. Default: NULL (select automatically).
#' @param facet.y.scale The shared type of y-axis scales across facets, choose from free (facets have different y-axis scales),
#' fixed (facets have same y-axis scales). Default: free.
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
#' @importFrom ggplot2 aes_string scale_fill_manual geom_rect geom_text aes geom_step
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
#' # library(ggcoverage)
#' # library(utils)
#' # library(ggplot2)
#' # meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
#' # sample.meta <- utils::read.csv(meta.file)
#' # track folder
#' # track.folder <- system.file("extdata", "RNA-seq", package = "ggcoverage")
#' # load bigwig file
#' # track.df <- LoadTrackFile(
#' #   track.folder = track.folder, format = "bw",
#' #   meta.info = sample.meta
#' # )
#' # ggplot() +
#' #   geom_coverage(data = track.df, color = "auto", mark.region = NULL)
geom_coverage <- function(data, mapping = NULL, color = NULL, rect.color = NA,
                          single.nuc = FALSE, plot.type = c("facet", "joint"),
                          facet.key = "Type", joint.avg = FALSE, facet.order = NULL,
                          facet.color = NULL, facet.y.scale = c("free", "fixed"),
                          group.key = "Group", range.size = 3, range.position = c("in", "out"),
                          mark.region = NULL, mark.color = "grey", mark.alpha = 0.5,
                          show.mark.label = TRUE, mark.label.size = 4) {
  # check parameters
  plot.type <- match.arg(arg = plot.type)
  facet.y.scale <- match.arg(arg = facet.y.scale)
  range.position <- match.arg(arg = range.position)

  # get mapping and color
  if (is.null(mapping)) {
    if (plot.type == "facet") {
      if (!is.null(color)) {
        if (length(color) != length(unique(data[, facet.key]))) {
          warning("The color you provided is not as long as ", facet.key, " column in data, select automatically!")
          # sample group with same color
          tmp.color <- AutoColor(data = data, n = 9, name = "Set1", key = group.key)
          # change facet key color color
          if (facet.key == group.key) {
            fill.color.df <- merge(unique(data[c(facet.key)]), data.frame(color = tmp.color), by.x = group.key, by.y = 0)
          } else {
            fill.color.df <- merge(unique(data[c(facet.key, group.key)]), data.frame(color = tmp.color), by.x = group.key, by.y = 0)
          }
          fill.color <- fill.color.df$color
          names(fill.color) <- fill.color.df[, facet.key]
        } else {
          fill.color <- color
          if (is.null(names(fill.color))) {
            names(fill.color) <- unique(data[, facet.key])
          }
        }
        sacle_fill_cols <- scale_fill_manual(values = fill.color)
      } else {
        sacle_fill_cols <- NULL
      }
      if (!single.nuc) {
        mapping <- aes_string(xmin = "start", xmax = "end", ymin = "0", ymax = "score", fill = facet.key)
      } else {
        mapping <- aes_string(x = "start", y = "score", fill = facet.key)
      }
    } else if (plot.type == "joint") {
      if (!is.null(color)) {
        if (length(color) != length(unique(data[, group.key]))) {
          warning("The color you provided is not as long as ", group.key, " column in data, select automatically!")
          # sample group with same color
          tmp.color <- AutoColor(data = data, n = 9, name = "Set1", key = group.key)
          # change group key color
          color.color.df <- merge(unique(data[c(group.key)]), data.frame(color = tmp.color), by.x = group.key, by.y = 0)
          color.color <- color.color.df$color
          names(color.color) <- color.color.df[, group.key]
        } else {
          color.color <- color
          if (is.null(names(color.color))) {
            names(color.color) <- unique(data[, group.key])
          }
        }
        sacle_color_cols <- scale_color_manual(values = color.color)
      } else {
        sacle_color_cols <- NULL
      }
      if (single.nuc) {
        stop("The joint mode is not available for SNV!")
      } else {
        mapping <- aes_string(x = "start", y = "score", color = group.key)
      }
    }
  } else {
    if (plot.type == "facet") {
      message("For SNV, the mapping should contains start, score, fill. For others, the mapping should contains start, end, score, fill.")
      if ("fill" %in% names(mapping)) {
        fill.str <- rlang::as_label(mapping$fill)
        fill.str.len <- length(unique(data[, fill.str]))
        if (is.null(color) | length(color) != fill.str.len) {
          # sample group with same color
          tmp.color <- AutoColor(data = data, n = 9, name = "Set1", key = group.key)
          # change color
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
    } else if (plot.type == "joint") {
      message("For joint visualization, the mapping should contains start, score, color.")
      if ("color" %in% names(mapping)) {
        color.str <- rlang::as_label(mapping$color)
        color.str.len <- length(unique(data[, color.str]))
        if (is.null(color) | length(color) != color.str.len) {
          # sample group with same color
          tmp.color <- AutoColor(data = data, n = 9, name = "Set1", key = group.key)
          # change color
          if (color.str == group.key) {
            color.color.df <- merge(unique(data[c(color.str)]), data.frame(color = tmp.color), by.x = group.key, by.y = 0)
          } else {
            color.color.df <- merge(unique(data[c(color.str, group.key)]), data.frame(color = tmp.color), by.x = group.key, by.y = 0)
          }
          color.color <- color.color.df$color
          names(color.color) <- color.color.df[, color.str]
        } else {
          color.color <- color
          if (is.null(names(color.color))) {
            names(color.color) <- unique(data[, color.str])
          }
        }
        sacle_color_cols <- scale_color_manual(values = color.str)
      } else {
        sacle_color_cols <- NULL
      }
    }
  }

  # plot type
  if (plot.type == "facet") {
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
    if (!single.nuc) {
      region.rect <- geom_rect(data = data, mapping = mapping, show.legend = F, colour = rect.color)
    } else {
      region.rect <- geom_bar(
        data = data, mapping = mapping, show.legend = F, colour = rect.color,
        stat = "identity"
      )
    }
    # prepare facet scale
    if (facet.y.scale == "free") {
      facet.ys <- "free_y"
    } else if (facet.y.scale == "fixed") {
      facet.ys <- "fixed"
    }
    region.facet <- facet_wrap2(
      facets = facet.formula, ncol = 1, scales = facet.ys, strip.position = "right",
      strip = strip_themed(background_y = elem_list_rect(
        fill = facet.color
      ))
    )
    plot.ele <- list(region.rect, region.facet)

    # color the track
    if (!is.null(sacle_fill_cols)) {
      plot.ele <- append(plot.ele, sacle_fill_cols)
    }

    # add range text to track
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
  } else if (plot.type == "joint") {
    data.split <- split(x = data, f = data[, facet.key])
    if (joint.avg) {
      # process data
      data.list <- lapply(data.split, function(ds) {
        ds.list <- apply(ds, 1, function(x) {
          x.start <- seq(x["start"], x["end"])
          x.len <- length(x.start)
          x.df <- data.frame(
            seqnames = rep(x["seqnames"], x.len), start = x.start,
            score = rep(x["score"], x.len)
          )
          x.df[, facet.key] <- x[facet.key]
          x.df[, group.key] <- x[group.key]
          x.df
        })
        ds.final <- do.call(rbind, ds.list)
        rownames(ds.final) <- NULL
        return(ds.final)
      })
      data.final <- do.call(rbind, data.list)
      rownames(data.final) <- NULL
      # as numeric
      data.final$score <- as.numeric(data.final$score)
      data.final$start <- as.numeric(data.final$start)
      # get average
      data.final.mean <- data.final %>%
        dplyr::group_by(.data[[group.key]], .data[["start"]]) %>%
        dplyr::summarise(score = mean(.data[["score"]]))
      # create avg plot
      avg.plot <- geom_line(data = data.final.mean, mapping = mapping)
      plot.ele <- list(avg.plot)
    } else {
      # preocess data
      data.split.final <- lapply(data.split, function(ds) {
        ds <- ds[order(ds$start), ]
        ds.ll <- ds[nrow(ds), ]
        if (ds.ll$start < ds.ll$end) {
          ds.ll.new <- ds.ll
          ds.ll.new$start <- ds.ll.new$end
          ds.final <- rbind(ds, ds.ll.new) %>% as.data.frame()
        } else {
          ds.final <- ds
        }
        rownames(ds.final) <- NULL
        return(ds.final)
      })
      # get final data
      data.final <- do.call(rbind, data.split.final) %>% as.data.frame()
      rownames(data.final) <- NULL

      # create step plot
      region.step <- geom_step(data = data.final, mapping = mapping)
      plot.ele <- list(region.step)
    }
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
        if (plot.type == "facet") {
          region.label[, facet.key] <- facet.order[1]
        }
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
