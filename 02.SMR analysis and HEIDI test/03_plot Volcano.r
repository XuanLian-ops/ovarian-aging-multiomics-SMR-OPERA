########################################
# 1. Load Packages
########################################

library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)

########################################
# 2. Input Files
########################################

files <- list(
  PQTL = "pqtl_discovery.txt",
  
  EQTL_dis = "eqtl_discovery.txt",
  EQTL_rep = "eqtl_ovary.txt",
  
  SQTL_dis = "sqtl_discovery.txt",
  SQTL_rep = "sqtl_ovary.txt",
  
  MQTL = "mqtl_discovery.txt"
)

output_file <- "multiomics_volcano_plot.pdf"

# threshold
threshold <- 0.0000275482093663911

########################################
# 3. Data Processing
########################################

process_data <- function(df) {
  
  df %>%
    mutate(
      logP = -log10(p_SMR),
      
      significant = p_SMR < threshold,
      
      color = case_when(
        significant & b_SMR > 0 ~ "#f36569",  # red
        significant & b_SMR < 0 ~ "#3581b7",  # blue
        TRUE ~ "grey"
      )
    )
}


########################################
# 4. Plot Functions
########################################

# ---- Single Omics (PQTL / MQTL) ----
plot_single <- function(file, title_text) {
  
  df <- read.table(file, header = TRUE)
  df <- process_data(df)
  
  p <- ggplot() +
    
    # background (non-significant)
    geom_point(
      data = df %>% filter(!significant),
      aes(x = b_SMR, y = logP),
      color = "grey",
      size = 1
    ) +
    
    # significant points
    geom_point(
      data = df %>% filter(significant),
      aes(x = b_SMR, y = logP, color = color),
      size = 2
    ) +
    
    # labels
    geom_text_repel(
      data = df %>% filter(significant),
      aes(x = b_SMR, y = logP, label = Gene),
      size = 3,
      max.overlaps = 100
    ) +
    
    scale_color_identity() +
    
    # reference lines
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed") +
    
    labs(
      title = title_text,
      x = "Effect Size (b_SMR)",
      y = expression(-log[10](p[SMR]))
    ) +
    
    theme_bw()
  
  return(p)
}


# ---- Dual Omics (EQTL / SQTL) ----
plot_dual <- function(file_dis, file_rep, title_text) {
  
  dis <- read.table(file_dis, header = TRUE)
  rep <- read.table(file_rep, header = TRUE)
  
  dis <- process_data(dis)
  rep <- process_data(rep)
  
  p <- ggplot() +
    
    # background (non-significant discovery)
    geom_point(
      data = dis %>% filter(!significant),
      aes(x = b_SMR, y = logP),
      color = "grey",
      size = 1
    ) +
    
    # replication (transparent layer)
    geom_point(
      data = rep,
      aes(x = b_SMR, y = logP, color = color),
      size = 1.5,
      alpha = 0.3
    ) +
    
    # discovery significant
    geom_point(
      data = dis %>% filter(significant),
      aes(x = b_SMR, y = logP, color = color),
      size = 2
    ) +
    
    # labels
    geom_text_repel(
      data = dis %>% filter(significant),
      aes(x = b_SMR, y = logP, label = Gene),
      size = 3,
      max.overlaps = 100
    ) +
    
    scale_color_identity() +
    
    # reference lines
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed") +
    
    labs(
      title = title_text,
      x = "Effect Size (b_SMR)",
      y = expression(-log[10](p[SMR]))
    ) +
    
    theme_bw()
  
  return(p)
}


########################################
# 5. Generate Plots
########################################

p_pqtl <- plot_single(files$PQTL, "PQTL")

p_eqtl <- plot_dual(
  files$EQTL_dis,
  files$EQTL_rep,
  "EQTL"
)

p_sqtl <- plot_dual(
  files$SQTL_dis,
  files$SQTL_rep,
  "SQTL"
)

p_mqtl <- plot_single(files$MQTL, "MQTL")


########################################
# 6. Combine Panels
########################################

final_plot <- (p_pqtl | p_eqtl) / (p_sqtl | p_mqtl)


########################################
# 7. Export Figure
########################################

ggsave(
  filename = output_file,
  plot = final_plot,
  width = 10,
  height = 8
)