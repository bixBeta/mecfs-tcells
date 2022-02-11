# get sankey for sv4 tcell mega v tcell subset reclustered object
library(tibble)
library(dplyr)
library(networkD3)
library(tidyr)
library(reshape2)
library(viridis)

tcell.clusts = c("0","1","23","11","16","25","6","4","12","28","13","5","20","3","15")

sv4.tcell.barcodes.orig.idents = sv4.meta %>% 
  rownames_to_column("BC") %>% 
  filter(SCT_snn_res.0.6 %in% tcell.clusts) %>% 
  select(BC, SCT_snn_res.0.6) %>% column_to_rownames("BC") 

tcellSub.v.tcellMega = as.matrix(table(tcells.sv4$SCT_snn_res.0.5, sv4.tcell.barcodes.orig.idents$SCT_snn_res.0.6))

colnames(tcellSub.v.tcellMega) = paste0("MEGA_", colnames(tcellSub.v.tcellMega))
row.names(tcellSub.v.tcellMega) = paste0("SUB_", row.names(tcellSub.v.tcellMega))

getSankey <- function(matrix, floor){
  
  x = as.data.frame(matrix) 
  colnames(x) = c("source", "target", "value")
  
  x = x %>% filter(value > floor)
  
  nodes = as.data.frame(c(as.character(x$source), as.character(x$target)) %>% unique())
  
  colnames(nodes) = "name"
  
  x$IDsource=match(x$source, nodes$name)-1 
  x$IDtarget=match(x$target, nodes$name)-1
  
  y = sankeyNetwork(Links = x, Nodes = nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=12, nodePadding=20, iterations = 0)
  
  return(y)
  
}



getSankey(matrix = tcellSub.v.tcellMega, floor = 10)

pheatmap(tcellSub.v.tcellMega, scale = "none", display_numbers = T, cluster_rows = F, cluster_cols = F, angle_col = 45, 
         number_color = "#FF0001", number_format = "%.0f", color = RColorBrewer::brewer.pal(5,"Blues"))


