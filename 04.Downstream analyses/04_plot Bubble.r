library(ggplot2)
library(readr)
library(dplyr)
data <- read_csv("D:/2026/genome biology/pathway.csv")
data$logP <- -log10(data$`Raw p`)
data$ratio <- data$hits / data$expected
data$pathway <- factor(data$pathway,
                       levels = rev(data$pathway))
ggplot(data, aes(x = logP,
                 y = pathway,
                 color = logP,
                 size = ratio)) +
  
  geom_point(alpha = 0.9) +
  
  scale_color_gradient(
    low = "#2b6cb0",
    high = "#d73027"
  ) +
  
  scale_size(range = c(4,12)) +
  
  labs(
    x = "-log10(P)",
    y = "Pathway",
    color = "-log10(P)",
    size = "Enrichment Ratio"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12)
  )
ggplot(data, aes(x = logP,
                 y = pathway,
                 color = logP,
                 size = ratio)) +
  
  geom_point(alpha = 0.9) +
  
  scale_color_gradient(
    low = "#2b6cb0",
    high = "#d73027"
  ) +
  
  scale_size(range = c(4,12)) +
  
  scale_x_continuous(limits = c(0, max(data$logP) * 1.1)) +
  
  labs(
    x = "-log10(P)",
    y = "Pathway",
    color = "-log10(P)",
    size = "Enrichment Ratio"
  ) +
  
  theme_bw()
ggsave("D:/2026/genome biology/pathway.pdf",
       width = 8,
       height = 8,
       dpi = 300)
