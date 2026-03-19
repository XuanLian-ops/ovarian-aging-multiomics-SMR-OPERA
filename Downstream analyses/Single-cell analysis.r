library(Seurat)

# Load data: obtain Seurat objects for ovary from human, monkey, and mouse
load('oa_wang.Rdata')   # Human ovary Seurat object
load('sce_monkey.Rdata') # Monkey ovary Seurat object
load('scRNA_mouse_ovary.Rdata') # Mouse ovary Seurat object

# Define colors for cell types
cell_colors <- c(
  "GC" = "#c1c88a",     # yellow-green
  "Oo" = "#dee68c",     # light yellow-green
  "endo" = "#d1464a",   # red
  "mono" = "#f5c1be",   # pink
  "NK&T" = "#4f7f9d",   # blue
  "SMC" = "#c6b15e",    # light brown-yellow
  "T&S" = "#3d6c5c"     # dark green
)

# Assign human object to a common variable (adjust if needed)
scRNA <- oa_wang.Rdata   # Note: This line may need correction; original used load() but then assigned incorrectly

# Plot UMAP
DimPlot(scRNA, 
        reduction = "umap", 
        group.by = "celltype", 
        label = TRUE, 
        pt.size = 1.5) + 
  scale_color_manual(values = cell_colors) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# Load package for bubble plots
library(scRNAtoolVis)    

# Marker genes to explore (human genes)
x1 <- c('MSH6','MC1R','RAD52','RBBP8','XRCC5','BABAM2','SIRT1',
        'RAD18','FANCM','SPATA18')

# Marker genes for mouse (if applicable)
x2 <- c('Msh6','Mc1r','Rad52','Rbbp8','Xrcc5','Babam2','Sirt1',
        'Rad18','Fancm','Spata18')

# Bubble plot grouped by cell type
jjDotPlot(object = scRNA, id = 'celltype', 
          gene = x1,                       # change to x2 for mouse
          xtree = FALSE, ytree = FALSE,
          scale = TRUE, x.text.angle = 90,
          dot.col = c('#692F7C','#D96558','#F6C63C'),
          midpoint = 0, same.pos.label = TRUE,
          tile.geom = TRUE) + coord_flip() + coord_flip()

# Bubble plot grouped by condition (e.g., young vs old)
jjDotPlot(object = scRNA, id = 'group', 
          gene = x1,                       # change to x2 for mouse
          xtree = FALSE, ytree = FALSE,
          scale = TRUE, x.text.angle = 30,
          dot.col = c('#692F7C','#D96558','#F6C63C'),
          midpoint = 0, same.pos.label = TRUE,
          tile.geom = TRUE) + coord_flip() + coord_flip()

# Bubble plot for a specific cell type (e.g., GC) grouped by condition
jjDotPlot(object = subset(scRNA, celltype == 'GC'), id = 'group', 
          gene = x1,                       # change to x2 for mouse
          xtree = FALSE, ytree = FALSE,
          scale = TRUE, x.text.angle = 30,
          dot.col = c('#692F7C','#D96558','#F6C63C'),
          midpoint = 0, same.pos.label = TRUE,
          tile.geom = TRUE) + coord_flip() + coord_flip()

# Differential expression analysis between conditions (e.g., old vs young)
diff_gene <- FindMarkers(scRNA, group.by = 'group',
                         ident.1 = 'old', ident.2 = 'young')