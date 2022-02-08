source("packages.R")
sv4 = readRDS("/workdir/dropbox/seurat_object_v4.RDS")
sv4.meta = sv4@meta.data
Idents(sv4) <- sv4$SCT_snn_res.0.6
sv4$seurat_clusters <- sv4$SCT_snn_res.0.6

Idents(sv4)
DimPlot(sv4, label = T, label.size = 9)


