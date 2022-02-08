source("packages.R")
sv4 = readRDS("/workdir/dropbox/seurat_object_v4.RDS")
sv4.meta = sv4@meta.data
Idents(sv4) <- sv4$SCT_snn_res.0.6
sv4$seurat_clusters = sv4$SCT_snn_res.0.6

Idents(sv4)
png(filename = "figures/sv4__res0.6__umap.png", width = 1920, height = 1080)
DimPlot(sv4, label = T, label.size = 9, raster = F) + ggtitle("SV4 res0.6 UMAP")
dev.off()

tcell.clusts = c("0","1","23","11","16","25","6","4","12","28","13","5","20","3","15")
# 
# lv.feats2 = c("CD3E","CD247","CD4","CD8A","TRBC1","TRGC1","TRDC","CCR7","CD28","SELL","CD69","CCL5","IFNG","GATA3","IRF4","SPI1","KLRB1","AHR","CTLA4","FOXP3","NCAM1","FCGR3A","KLRK1","NCR3","LILRB1","CD14","LGALS2","S100A8","PF4","MS4A1","CD1C","HLA-DRA","ITGAX","CD34")
# 
# DefaultAssay(sv4)
# 
# png(filename = "figures/sv4__dotPlot-gex-clusters.png", width = 1000, height = 1000)
# DotPlot(sv4, assay="RNA", features=lv.feats2, col.min=-1.5,
#         group.by = "SCT_snn_res.0.6",cols="Spectral", dot.scale=3) + RotatedAxis() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("sv4__dotPlot-gex-clusters") 
# dev.off()


# DotPlot(sv4, assay="RNA", features=lv.feats2, col.min=-1.5,
#         group.by = "SCT_snn_res.0.6",cols="Spectral", dot.scale=3, idents = "0", split.by = "orig.ident") + RotatedAxis() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("sv4__dotPlot-gex-clusters") 

# select t-cell barcodes
sv4.meta = sv4@meta.data
sv4.tcell.barcodes.orig.idents = sv4.meta %>% 
  rownames_to_column("BC") %>% 
  filter(SCT_snn_res.0.6 %in% tcell.clusts) %>% 
  select(BC, orig.ident, SCT_snn_res.0.6)

cells.use = unname(unlist(sv4.tcell.barcodes.orig.idents %>% select(BC)))

sv4.tcells = subset(sv4, cells = cells.use)

saveRDS(sv4.tcells, "data/sv4__TCELLS__.RDS")

