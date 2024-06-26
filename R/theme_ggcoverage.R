# theme for ggcoverage: suitable for range position is in
#' Theme for geom_coverage.
#'
#' @param space The space between facets. Default: 0.2.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme unit element_blank annotate rel scale_y_continuous expansion
#' scale_x_continuous
#' @importFrom scales comma
#' @export
#'
theme_coverage <- function(space = 0.2) {
  list(
    theme_classic(),
    theme(
      panel.spacing.y = unit(x = space, units = "line"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank()
    ),
    ggplot2::annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = rel(1)),
    scale_y_continuous(
      limits = ~ c(pretty(.x)[1], tail(pretty(.x), 1)),
      expand = expansion(mult = c(0))
    ),
    scale_x_continuous(labels = scales::comma, expand = c(0, 0))
  )
}

# theme for ggcoverage: suitable for range position is out
#' Theme for geom_coverage.
#'
#' @param space The space between facets. Default: 0.2.
#'
#' @return List of layers.
#' @importFrom ggplot2 scale_y_continuous expansion theme_classic theme unit element_blank annotate rel
#' scale_x_continuous
#' @importFrom scales comma
#' @export
#'
theme_coverage2 <- function(space = 0.2) {
  list(
    scale_y_continuous(
      limits = ~ range(pretty(.x)),
      breaks = ~ .x[2],
      expand = expansion(mult = c(0))
    ),
    theme_classic(),
    theme(
      panel.spacing.y = unit(x = space, units = "line"),
      axis.title = element_blank()
    ),
    ggplot2::annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = rel(1)),
    scale_x_continuous(labels = scales::comma, expand = c(0, 0))
  )
}

# theme for geom_gc
#' Theme for geom_gc.
#'
#' @param x.range X axis ranges.
#' @param margin.len Top and bottom margin.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_y_continuous scale_x_continuous coord_cartesian
#' @export
#'
theme_gc <- function(x.range, margin.len) {
  list(
    theme_classic(),
    theme(
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      plot.margin = margin(t = margin.len, b = margin.len)
    ),
    scale_y_continuous(position = "right"),
    scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(xlim = x.range)
  )
}

#' Theme for geom_base.
#'
#' @param margin.len Top and bottom margin.
#' @param fill.color Fill color.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text margin
#' scale_y_continuous scale_x_continuous coord_cartesian scale_fill_manual
#' @export
theme_base <- function(margin.len, fill.color) {
  list(
    theme_classic(),
    theme(
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = margin.len, b = margin.len)
    ),
    scale_y_continuous(expand = c(0, 0), position = "right"),
    scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(clip = "off"),
    scale_fill_manual(values = fill.color)
  )
}

#' Theme for geom_base without margin.
#'
#' @param fill.color Fill color.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text
#' scale_y_continuous scale_x_continuous coord_cartesian scale_fill_manual
#' @export
theme_base2 <- function(fill.color) {
  list(
    theme_classic(),
    theme(
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "none"
    ),
    scale_y_continuous(expand = c(0, 0), position = "right"),
    scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(clip = "off"),
    scale_fill_manual(values = fill.color)
  )
}

#' Theme for geom_base with Amino Acid.
#'
#' @param margin.len Top and bottom margin.
#' @param fill.color Fill color.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text margin
#' scale_y_continuous scale_x_continuous coord_cartesian scale_fill_manual
#' @export
theme_aa <- function(margin.len, fill.color) {
  list(
    theme_classic(),
    theme(
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = margin.len, b = margin.len)
    ),
    scale_y_continuous(expand = c(0, 0), position = "right"),
    scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(clip = "off"),
    scale_fill_manual(values = fill.color)
  )
}

# theme for geom_gene
#' Theme for geom_gene.
#'
#' @param overlap.gene.gap The gap between gene groups.
#' @param group.num The number of groups.
#' @param fill.color Fill color.
#' @param x.range X axis ranges.
#' @param margin.len Top and bottom margin.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_y_continuous scale_color_manual scale_x_continuous coord_cartesian
#' @export
#'
theme_gene <- function(overlap.gene.gap, group.num, fill.color, x.range, margin.len) {
  list(
    theme_classic(),
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      plot.margin = margin(t = margin.len, b = margin.len)
    ),
    scale_y_continuous(
      limits = c(1 - overlap.gene.gap, 1 + overlap.gene.gap * (group.num)),
      expand = c(0, 0), position = "right"
    ),
    scale_color_manual(values = fill.color),
    scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(xlim = x.range)
  )
}

# theme for geom_transcript
#' Theme for geom_transcript.
#'
#' @param overlap.tx.gap The gap between transcript groups.
#' @param group.num The number of groups.
#' @param fill.color Fill color.
#' @param x.range X axis ranges.
#' @param margin.len Top and bottom margin.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_y_continuous scale_color_manual scale_x_continuous coord_cartesian
#' @export
#'
theme_transcript <- function(overlap.tx.gap, group.num, fill.color, x.range, margin.len) {
  list(
    theme_classic(),
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      plot.margin = margin(t = margin.len, b = margin.len)
    ),
    scale_y_continuous(
      limits = c(1 - overlap.tx.gap, 1 + overlap.tx.gap * (group.num)),
      expand = c(0, 0), position = "right"
    ),
    scale_color_manual(values = fill.color),
    scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(xlim = x.range)
  )
}


# theme for ideogram
#' Theme for geom_ideogram.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank scale_x_continuous scale_y_continuous
#' @export
#'
theme_ideogram <- function() {
  list(
    theme_classic(),
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    ),
    scale_x_continuous(expand = c(0, 0)),
    scale_y_continuous(expand = c(0, 0), limits = c(-0.2, 16))
  )
}

# theme for peak
#' Theme for geom_peak.
#'
#' @param margin.len Top and bottom margin.
#' @param x.range X axis ranges.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_x_continuous scale_y_continuous coord_cartesian
#' @export
#'
theme_peak <- function(margin.len, x.range) {
  list(
    theme_classic(),
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      plot.margin = margin(t = margin.len, b = margin.len)
    ),
    scale_y_continuous(
      limits = c(1 - 0.1, 1 + 0.1),
      expand = c(0, 0), position = "right"
    ),
    scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(xlim = x.range)
  )
}

#' Theme for geom_tad2.
#'
#' @param x.range X axis ranges.
#' @param margin.len Top and bottom margin.
#' @param show.rect Logical value, whether to add rect border to the plot.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_x_continuous scale_y_continuous coord_cartesian
#' @export
theme_tad2 <- function(x.range, margin.len, show.rect) {
  if (show.rect) {
    list(
      theme_classic(),
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.margin = margin(t = margin.len, b = margin.len)
      ),
      scale_y_continuous(expand = c(0, 0), position = "right"),
      scale_x_continuous(expand = c(0, 0)),
      coord_cartesian(xlim = x.range)
    )
  } else {
    list(
      theme_classic(),
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.margin = margin(t = margin.len, b = margin.len)
      ),
      scale_y_continuous(expand = c(0, 0), position = "right"),
      scale_x_continuous(expand = c(0, 0)),
      coord_cartesian(xlim = x.range)
    )
  }
}

#' Theme for geom_tad.
#'
#' @param margin.len Top and bottom margin.
#' @param show.rect Logical value, whether to add rect border to the plot.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_x_continuous scale_y_continuous
#' @export
theme_tad <- function(margin.len, show.rect) {
  if (show.rect) {
    list(
      theme_classic(),
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.margin = margin(t = margin.len, b = margin.len)
      ),
      scale_y_continuous(expand = c(0, 0), position = "right"),
      scale_x_continuous(expand = c(0, 0))
    )
  } else {
    list(
      theme_classic(),
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_blank(),
        plot.margin = margin(t = margin.len, b = margin.len)
      ),
      scale_y_continuous(expand = c(0, 0), position = "right"),
      scale_x_continuous(expand = c(0, 0))
    )
  }
}

#' Theme for geom_link.
#'
#' @param x.range X axis ranges.
#' @param margin.len Top and bottom margin.
#' @param show.rect Logical value, whether to add rect border to the plot.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_x_continuous scale_y_continuous coord_cartesian expansion
#' @export
theme_link <- function(x.range, margin.len, show.rect) {
  if (show.rect) {
    list(
      theme_classic(),
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.margin = margin(t = margin.len, b = margin.len)
      ),
      scale_y_continuous(expand = expansion(add = c(0.5, 0)), position = "right"),
      scale_x_continuous(expand = c(0, 0)),
      coord_cartesian(xlim = x.range)
    )
  } else {
    list(
      theme_classic(),
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        plot.margin = margin(t = margin.len, b = margin.len)
      ),
      scale_y_continuous(expand = expansion(add = c(0.5, 0))),
      scale_x_continuous(expand = c(0, 0)),
      coord_cartesian(xlim = x.range)
    )
  }
}

#' Theme for geom_cnv.
#'
#' @param x.range X axis ranges.
#' @param margin.len Top and bottom margin.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_x_continuous coord_cartesian
#' @export
theme_cnv <- function(x.range, margin.len) {
  list(
    theme_classic(),
    theme(
      axis.line.x = element_blank(),
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.title = element_blank(),
      plot.margin = margin(t = margin.len, b = margin.len)
    ),
    scale_x_continuous(expand = c(0, 0)),
    scale_y_continuous(position = "right"),
    coord_cartesian(xlim = x.range)
  )
}

# theme for ggprotein: suitable for range.position = "in"
#' Theme for geom_protein.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme unit element_blank annotate rel scale_y_continuous expansion
#' scale_x_continuous
#' @importFrom scales comma
#' @export
#'
theme_protein <- function() {
  list(
    theme_classic(),
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank()
    ),
    ggplot2::annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = rel(1)),
    scale_y_continuous(expand = expansion(mult = c(0))),
    scale_x_continuous(expand = c(0, 0))
  )
}

# theme for ggprotein: suitable for range.position = "out"
#' Theme for geom_protein.
#'
#' @return List of layers.
#' @importFrom ggplot2 scale_y_continuous expansion theme_classic theme unit element_blank annotate rel
#' scale_x_continuous
#' @importFrom scales comma scientific
#' @export
#'
theme_protein2 <- function() {
  list(
    scale_y_continuous(
      limits = ~ range(pretty(.x)),
      breaks = ~ .x[2],
      expand = expansion(mult = c(0)),
      labels = function(x) format(x, scientific = TRUE, digits = 2)
    ),
    theme_classic(),
    theme(
      axis.title = element_blank()
    ),
    ggplot2::annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = rel(1)),
    scale_x_continuous(expand = c(0, 0))
  )
}

#' Theme for geom_feature.
#'
#' @param margin.len Top and bottom margin.
#' @param x.range X axis ranges.
#' @param fill.color Fill color for different feature type.
#'
#' @return List of layers.
#' @importFrom ggplot2 theme_classic theme element_blank element_text element_rect margin
#' scale_x_continuous scale_y_continuous coord_cartesian scale_color_manual
#' @export
#'
theme_feature <- function(margin.len, x.range, fill.color) {
  list(
    theme_classic(),
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y.right = element_text(color = "black", angle = 90, vjust = 0.5),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      plot.margin = margin(t = margin.len, b = margin.len)
    ),
    scale_y_continuous(
      limits = c(1 - 0.1, 1 + 0.1),
      expand = c(0, 0), position = "right"
    ),
    scale_x_continuous(expand = c(0, 0)),
    coord_cartesian(xlim = x.range),
    scale_color_manual(values = fill.color)
  )
}
