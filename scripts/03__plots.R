Idents(object = tcells.sv4) <- tcells.sv4$SCT_snn_res.0.4
res0.4 = DimPlot(tcells.sv4 , group.by = "SCT_snn_res.0.4", label = T, label.size = 8, raster = F) + ggtitle("res 0.4")

Idents(object = tcells.sv4) <- tcells.sv4$SCT_snn_res.0.5
res0.5 = DimPlot(tcells.sv4 , group.by = "SCT_snn_res.0.5", label = T, label.size = 8, raster = F) + ggtitle("res 0.5")

Idents(object = tcells.sv4) <- tcells.sv4$SCT_snn_res.0.6
res0.6 = DimPlot(tcells.sv4 , group.by = "SCT_snn_res.0.6", label = T, label.size = 8, raster = F) + ggtitle("res 0.6")

Idents(object = tcells.sv4) <- tcells.sv4$SCT_snn_res.0.8
res0.8 = DimPlot(tcells.sv4 , group.by = "SCT_snn_res.0.8", label = T, label.size = 8, raster = F) + ggtitle("res 0.8")

png(filename = "figures/TCELLS__sct__harmony__umap.png", width = 1929, height = 1080)
res0.4 + res0.5 + res0.6 + res0.8
dev.off()


Idents(object = tcells.sv4) <- tcells.sv4$SCT_snn_res.0.5
tcells.sv4$seurat_clusters = tcells.sv4$SCT_snn_res.0.5
tcells.sv4.meta = tcells.sv4@meta.data


lv.feats2 = c("CD3E","CD247","CD4","CD8A","TRBC1","TRGC1","TRDC",
              "CCR7","CD28","SELL","CD69","CCL5","IFNG","GATA3","IRF4",
              "SPI1","KLRB1","AHR","CTLA4","FOXP3","NCAM1","FCGR3A","KLRK1",
              "NCR3","LILRB1","CD14","LGALS2","S100A8","PF4","MS4A1","CD1C",
              "HLA-DRA","ITGAX","CD34")

DefaultAssay(tcells.sv4) <- "SCT"
png(filename = "figures/SCT__TCELLS__sv4__dotPlot-gex-clusters.png", width = 1000, height = 1000)
DotPlot(tcells.sv4, assay="SCT", features=lv.feats2, col.min=-1.5,
        group.by = "SCT_snn_res.0.5",cols="Spectral", dot.scale=3) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Tcells__sv4__dotPlot-gex-clusters")
dev.off()

DefaultAssay(tcells.sv4) <- "RNA"
png(filename = "figures/RNA__TCELLS__sv4__dotPlot-gex-clusters.png", width = 1000, height = 1000)
DotPlot(tcells.sv4, assay="RNA", features=lv.feats2, col.min=-1.5,
        group.by = "SCT_snn_res.0.5",cols="Spectral", dot.scale=3) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Tcells__sv4__dotPlot-gex-clusters")
dev.off()


