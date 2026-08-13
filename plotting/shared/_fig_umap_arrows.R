# UMAP corner-arrow axis indicator (replaces full-length axes).
# Use:  p + umap_arrow_axes(data = df, x_col = "UMAP_1", y_col = "UMAP_2")
# Drops axis text/title/line; draws a small L of arrows at bottom-left.
umap_arrow_axes <- function(data, x_col, y_col,
                            frac = 0.16,
                            label_x = "UMAP 1", label_y = "UMAP 2",
                            text_size = 2.4, line_size = 0.45,
                            arrow_mm = 1.6,
                            inset_frac = 0.02) {
  xr <- range(data[[x_col]], na.rm = TRUE)
  yr <- range(data[[y_col]], na.rm = TRUE)
  x0 <- xr[1] + inset_frac * diff(xr)
  y0 <- yr[1] + inset_frac * diff(yr)
  x1 <- x0 + frac * diff(xr)
  y1 <- y0 + frac * diff(yr)
  arr <- grid::arrow(length = grid::unit(arrow_mm, "mm"),
                     ends = "last", type = "closed")
  fam <- if (exists("my_font")) my_font else "sans"
  list(
    annotate("segment", x = x0, xend = x1, y = y0, yend = y0,
             arrow = arr, linewidth = line_size, color = "black"),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y1,
             arrow = arr, linewidth = line_size, color = "black"),
    annotate("text", x = (x0 + x1) / 2, y = y0, label = label_x,
             vjust = 2.4, size = text_size, family = fam),
    annotate("text", x = x0, y = (y0 + y1) / 2, label = label_y,
             angle = 90, vjust = -1.4, size = text_size, family = fam),
    theme(axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.line  = element_blank(),
          panel.grid = element_blank())
  )
}
