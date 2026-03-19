########################################
# 1. Load Required Packages
########################################

# Install if needed
# install.packages("circlize")
# install.packages("RColorBrewer")
# install.packages("ggplot2")

library(circlize)
library(RColorBrewer)
library(ggplot2)


########################################
# 2. Set Working Directory
########################################

setwd("D:/keti/circos")


########################################
# 3. Input Files
########################################

files <- list(
  chr = "example_data_chromosome_general.csv",
  
  pqtl_points = "pqtl_points.csv",
  pqtl_labels = "pqtl_labels.csv",
  
  eqtl_points = "eqtl_points.csv",
  eqtl_labels = "eqtl_labels.csv",
  
  sqtl_points = "sqtl_points.csv",
  sqtl_labels = "sqtl_labels.csv",
  
  mqtl_points = "mqtl_points.csv",
  mqtl_labels = "mqtl_labels.csv"
)


########################################
# 4. Load Data
########################################

chr <- read.csv(files$chr)

pqtl <- read.csv(files$pqtl_points)
pqtllabel <- read.csv(files$pqtl_labels)

eqtl <- read.csv(files$eqtl_points)
eqtllabel <- read.csv(files$eqtl_labels)

sqtl <- read.csv(files$sqtl_points)
sqtllabel <- read.csv(files$sqtl_labels)

mqtl <- read.csv(files$mqtl_points)
mqtllabel <- read.csv(files$mqtl_labels)


########################################
# 5. Initialize Circos Plot
########################################

circos.clear()

# Chromosome colors (alternating)
chr_col <- rep(c("#DCDCDC", "#F5F5F5"), length.out = length(unique(chr$chr)))

# Gap between chromosomes
gap <- rep(1, length(chr_col))

circos.par(
  start.degree = 90,
  gap.degree = gap
)

# Initialize chromosome layout
circos.genomicInitialize(
  chr,
  tickLabelsStartFromZero = FALSE,
  major.by = 400000000,
  labels.cex = 0.6
)

# Draw chromosome track
circos.genomicTrackPlotRegion(
  chr,
  ylim = c(-2, 2),
  track.height = 0.05,
  bg.col = chr_col,
  bg.border = NA
)


########################################
# 6. Add QTL Tracks (Scatter Layers)
########################################

# PQTL
circos.track(
  factors = pqtl$chr,
  x = pqtl$start,
  y = pqtl$name1,
  track.height = 0.13,
  ylim = c(0, 10),
  bg.col = "#FFFAF0",
  bg.border = NA,
  panel.fun = function(x, y) {
    circos.points(x, y, pch = 16, col = "grey", cex = 0.4)
    circos.lines(CELL_META$xlim, c(4.56, 4.56), lty = 2, col = "red")
  }
)

# EQTL
circos.track(
  factors = eqtl$chr,
  x = eqtl$start,
  y = eqtl$name1,
  track.height = 0.18,
  ylim = c(0, 27),
  bg.col = "#FFFAFA",
  bg.border = NA,
  panel.fun = function(x, y) {
    circos.points(x, y, pch = 16, col = eqtl$name2, cex = 0.4)
    circos.lines(CELL_META$xlim, c(5.5, 5.5), lty = 2, col = "red")
    circos.lines(CELL_META$xlim, c(4.65, 4.65), lty = 2, col = "blue")
  }
)

# SQTL
circos.track(
  factors = sqtl$chr,
  x = sqtl$start,
  y = sqtl$value1,
  track.height = 0.18,
  ylim = c(0, 36),
  bg.col = "#FFF5EE",
  bg.border = NA,
  panel.fun = function(x, y) {
    circos.points(x, y, pch = 16, col = sqtl$color, cex = 0.4)
    circos.lines(CELL_META$xlim, c(6.07, 6.07), lty = 2, col = "red")
    circos.lines(CELL_META$xlim, c(4.95, 4.95), lty = 2, col = "blue")
  }
)

# MQTL
circos.track(
  factors = mqtl$chr,
  x = mqtl$start,
  y = mqtl$value1,
  track.height = 0.16,
  ylim = c(0, 10),
  bg.col = "#FDF5E6",
  bg.border = NA,
  panel.fun = function(x, y) {
    circos.points(x, y, pch = 16, col = "grey", cex = 0.4)
    circos.lines(CELL_META$xlim, c(4.31, 4.31), lty = 2, col = "red")
  }
)


########################################
# 7. Add Labels (General Function)
########################################

add_labels <- function(label_data, track_index) {
  for (i in 1:nrow(label_data)) {
    
    # Point
    circos.points(
      label_data[i, 2],
      label_data[i, 5],
      sector.index = label_data[i, 1],
      track.index = track_index,
      col = label_data[i, 6],
      pch = 16,
      cex = 0.4
    )
    
    # Text
    circos.text(
      label_data[i, 2],
      label_data[i, 5],
      label_data[i, 4],
      sector.index = label_data[i, 1],
      track.index = track_index,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.4,
      adj = c(-0.2, 0)
    )
  }
}

# Apply labels
add_labels(pqtllabel, 3)
add_labels(eqtllabel, 4)
add_labels(sqtllabel, 5)
add_labels(mqtllabel, 6)


########################################
# 8. Export Plot (Optional)
########################################

# Save as PDF
pdf("circos_plot.pdf", width = 8, height = 8)
circos.clear()

# (Re-run plotting code here if exporting)

dev.off()