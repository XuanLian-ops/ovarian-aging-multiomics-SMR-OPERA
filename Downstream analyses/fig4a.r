library(ggplot2)
library(readr)
library(dplyr)
data <- read_csv("D:/2026/genome biology/GO.csv")
data$Description <- factor(data$Description, levels = rev(data$Description))
data$LogP_abs <- abs(data$LogP)
ggplot(data, aes(x = InTerm_InList, y = Description, fill = LogP_abs)) +
  geom_col(width = 0.75) +
  scale_fill_gradient(low = "#C8DCE9", high = "#4480B0") +
  labs(
    x = "Number of Genes",
    y = "",
    fill = "-LogP"
  ) +
  theme_bw()
scale_fill_gradient(
  low = "#d0e1f2",
  high = "#08519c"
)
ggsave("D:/2026/genome biology/enrichment_barplot.pdf", width = 8, height = 6, dpi = 300)