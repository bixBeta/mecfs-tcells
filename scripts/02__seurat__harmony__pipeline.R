.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")

library("Seurat")
library("glmGamPoi")
library("harmony")
library("patchwork")
library("ggplot2")

# ran using command line

sobj = readRDS("data/sv4__TCELLS__.RDS")

sobj.list <- SplitObject(sobj, split.by="orig.ident")

sobj.list <- lapply(X = sobj.list, FUN = function(x){SCTransform(object = x, method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = FALSE)})

var.features <- SelectIntegrationFeatures(object.list = sobj.list, nfeatures = 3000)

sobj.sct <- merge(x = sobj.list[[1]], y = sobj.list[2:length(sobj.list)], merge.data=TRUE)

VariableFeatures(sobj.sct) <- var.features

# save.image("A__sct__harmony__split__by__orig.ident.Rdata")


sobj.sct <- RunPCA(sobj.sct, verbose = FALSE, dims = 1:50)

sobj.sct <- RunHarmony(sobj.sct, assay.use="SCT", group.by.vars = "orig.ident")

sobj.sct <- RunUMAP(sobj.sct, reduction = "harmony", dims = 1:50)

sobj.sct <- FindNeighbors(sobj.sct, reduction = "harmony", dims = 1:50) %>% FindClusters(resolution = c(0.4, 0.5, 0.6, 0.8))


saveRDS(sobj.sct, file = "data/TCELLS_sct__harmony__ndim50__res0.4_0.5_0.6_0.8.RDS")

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.4
res0.4 = DimPlot(sobj.sct , group.by = "SCT_snn_res.0.4", label = T, label.size = 5) + ggtitle("res 0.4")

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.5
res0.5 = DimPlot(sobj.sct , group.by = "SCT_snn_res.0.5", label = T, label.size = 5) + ggtitle("res 0.5")

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.6
res0.6 = DimPlot(sobj.sct , group.by = "SCT_snn_res.0.6", label = T, label.size = 5) + ggtitle("res 0.6")

Idents(object = sobj.sct) <- sobj.sct$SCT_snn_res.0.8
res0.8 = DimPlot(sobj.sct , group.by = "SCT_snn_res.0.8", label = T, label.size = 5) + ggtitle("res 0.8")

png(filename = "figures/TCELLS__sct__harmony__umap.png", width = 3840, height = 2160)
res0.4 + res0.5 + res0.6 + res0.8
dev.off()


# save.image("data/CMD__TCELLS__sct__harmony__split__by__orig.ident.Rdata")
DefaultAssay(tcells.sv4) <- "RNA"
Idents(tcells.sv4)
allmarkers.res0.5 = FindAllMarkers(tcells.sv4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(allmarkers.res0.5, file = "results/Tcells_seurat_v4_findAllMarkers_minpct.25_lfc.25.csv", quote = F)
saveRDS(allmarkers.res0.5, file = "results/Tcells_seurat_v4_findAllMarkers_minpct.25_lfc.25.RDS") 
