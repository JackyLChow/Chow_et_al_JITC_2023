# final, fully annotated sce, contains TSNE and UMAP calculated for all cells
sce_final <- SingleCellExperiment()
for(file_ in list.files("~/_Projects/Chow Nat Comm/data/sce_final", full.names = T)){
  if(nrow(sce_final) == 0){
    sce_final <- readRDS(file_)
  } else {
    sce_final <- cbind(sce_final, readRDS(file_))
  }
}