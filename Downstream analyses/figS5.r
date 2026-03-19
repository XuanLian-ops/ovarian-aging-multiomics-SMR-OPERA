########################################
# 1. Load Required Package
########################################

library(circlize)
########################################
# 2. Input Data
########################################

# Input file (gene interaction)
input_file <- "tables1.csv"

# Output file
output_file <- "chord_plot.pdf"

# Read data
data <- read.csv(input_file, header = TRUE)

# Keep required columns: gene1, gene2, score
plot_data <- data[, c("node1", "node2", "combined_score")]

colnames(plot_data) <- c("from", "to", "value")


########################################
# 3. Define Colors
########################################

# Blue → White → Red gradient (same as original)
color_ramp <- colorRamp(c("#3581b7", "#FFFFFF", "#f36569"))
link_colors <- rgb(color_ramp(seq(0, 1, length.out = 100)), maxColorValue = 255)

# Map score to color
plot_data$color <- link_colors[
  as.numeric(cut(plot_data$value, breaks = 100))
]

# Sector (gene) colors
genes <- unique(c(plot_data$from, plot_data$to))

# Assign colors (blue for first half, red for second half)
sector_colors <- c(
  rep("#3581b7", length(genes) %/% 2),
  rep("#f36569", length(genes) - length(genes) %/% 2)
)

names(sector_colors) <- genes


########################################
# 4. Plot Chord Diagram
########################################

pdf(output_file, width = 10, height = 10)

circos.clear()
circos.par(start.degree = 90, gap.degree = 2)

chordDiagram(
  x = plot_data[, c("from", "to", "value")],
  grid.col = sector_colors,
  col = plot_data$color,
  transparency = 0.4,
  directional = 0,
  link.lwd = 1,
  link.lty = 0,
  link.border = NA,
  annotationTrack = "grid"
)


########################################
# 5. Add Labels
########################################

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  },
  bg.border = NA
)

dev.off()