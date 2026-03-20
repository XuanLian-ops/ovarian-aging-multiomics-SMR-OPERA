#### Data loading
library(Seurat)
load("oa_wang.Rdata")  ## Human ovary
load("sce_monkey.Rdata")  ## Monkey ovary
load("scRNA_mouse_ovary.Rdata") ## Mouse ovary
scRNA <- oa_wang
#scRNA <- sce_monkey
#scRNA <- scRNA_mouse_ovary


#### UMAP visualization
#Cell-type distributions were visualized using UMAP embedding.

cell_colors <- c(
  "GC" = "#c1c88a",
  "Oo" = "#dee68c",
  "endo" = "#d1464a",
  "mono" = "#f5c1be",
  "NK&T" = "#4f7f9d",
  "SMC" = "#c6b15e",
  "T&S" = "#3d6c5c"
)
#  UMAP
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

####  Gene expression visualization
#Expression patterns of prioritized genes were visualized using dot plots.

library(scRNAtoolVis)
genes <- c('MSH6','MC1R','RAD52','RBBP8','XRCC5',
           'BABAM2','SIRT1','RAD18','FANCM','SPATA18')

genes2=c('Msh6','Mc1r','Rad52','Rbbp8','Xrcc5','Babam2','Sirt1',
     'Rad18','Fancm','Spata18')


jjDotPlot(object=scRNA,id = 'celltype',
          gene=genes, ## mouse genes2
          xtree=F,ytree = F,
          scale = T,x.text.angle=90,
          dot.col =c('#692F7C','#D96558','#F6C63C'  ) 
          ,midpoint =0,same.pos.label=T,
          tile.geom=T) + coord_flip()+ coord_flip()

jjDotPlot(object=scRNA,id = 'group',# 
          gene=genes, ## mouse genes2
          xtree=F,ytree = F,
          scale = T,x.text.angle=30,
          dot.col =c('#692F7C','#D96558','#F6C63C'  ) 
          ,midpoint =0,same.pos.label=T,
          tile.geom=T) + coord_flip()+ coord_flip()


jjDotPlot(object=subset(scRNA,celltype=='GC'),id = 'group',
          gene=genes, ## mouse genes2
          xtree=F,ytree = F,
          scale = T,x.text.angle=30,
          dot.col =c('#692F7C','#D96558','#F6C63C'  ) 
          ,midpoint =0,same.pos.label=T,
          tile.geom=T) + coord_flip()+ coord_flip()


#### Differential expression analysis
#Differential expression analysis between age groups was performed using Seurat:
  

diff_gene=FindMarkers(scRNA,group.by = 'group',
                      ident.1 = 'old',ident.2 = 'young')


