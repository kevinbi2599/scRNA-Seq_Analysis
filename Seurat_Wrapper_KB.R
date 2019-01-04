##Automated Seurat object construction -> filtering -> normalization -> PCA 
#mat    matrix
#allgenes   T or F, use all genes in PCA or subset of variable genes, recommend F
#delim      column name field delimiter (i.e. '_')
#field      number of delimited field to use for identity grouping (i.e. 2)
#doregress    T or F, regress out effects from nUMI and percent.mito, recommend F
#maxgenes     max number of genes a cell can have, recommend 10000
#save_seurobj_for_qcplot T or F, save and store object for output_qcplots, plotting various QC metrics

library(Seurat)
library(dplyr)
ribo <- readRDS("/Users/kevinbi/Documents/scRNASeq/ribosomal_genes_GO0005840.RDS")

seurat_10x_autofrommat_complex <- function(mat, allgenes, delim = NULL, field = NULL, doregress, vars_to_regress = NULL, maxgenes, save_seurobj_for_qcplot){
  pbmc.data <- mat
  genes <- row.names(pbmc.data)
  genes <- gsub("hg19_","", genes)
  row.names(pbmc.data) <- genes
  if(is.null(delim)){
    pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                               project = "10X_PBMC")
  } else {
    pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
                               project = "10X_PBMC", names.field = field, names.delim = delim)}
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
  rp.genes <- intersect(ribo, rownames(pbmc@data))
  percent.mito <- colSums(pbmc@raw.data[mito.genes, ])/colSums(pbmc@raw.data)
  percent.rp <- colSums(pbmc@raw.data[rp.genes, ])/colSums(pbmc@raw.data)
  log10_nUMI <- log(pbmc@meta.data$nUMI, 10)
  pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
  pbmc <- AddMetaData(object = pbmc, metadata = percent.rp, col.name = "percent.rp")
  pbmc <- AddMetaData(object = pbmc, metadata = log10_nUMI, col.name = "log10_nUMI")
  if(isTRUE(save_seurobj_for_qcplot)){
    SeurObj_ForQCPlot <<- pbmc}
  VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                      low.thresholds = c(200, -Inf), high.thresholds = c(maxgenes, 0.1))
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 4, y.cutoff = 0.5)
  length(x = pbmc@var.genes)
  if(isTRUE(doregress)){
    pbmc <- ScaleData(object = pbmc, vars.to.regress = vars_to_regress, model.use = 'negbinom')} 
  else {pbmc <- ScaleData(object = pbmc)}
  if(isTRUE(allgenes)){
    pbmc <- RunPCA(object = pbmc, pc.genes = row.names(as.matrix(pbmc@data)), do.print = TRUE, pcs.print = 1:5, 
                   genes.print = 5, pcs.compute = 60)
  } else {
    pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
                   genes.print = 5, pcs.compute = 60)}
  return(pbmc)
}


#Return percentage of cells in a given batch assigned to each cluster
batchcontribution_tocluster <- function(object){
  seur <- object
  dat <- table(seur@ident, seur@meta.data$orig.ident)
  numbatches <- length(colnames(dat))
  colsums <- colSums(dat)
  results <- sweep(dat, 2, colsums, '/')
  return(results)
}



##calculate percent expression of a single gene across dataset
#vec    expression vector
calc_percent_exp_onegene <- function(vec){
  percent <- length(which(vec > 0))/length(vec)
  return(percent)
}

##filter genes by percent expression across cells
#mat    matrix
#low    low cutoff, recommend 0.01
#high   high cutoff, recommend 0.95-0.99 to remove ribosomal genes
percent_filter <- function(mat, low, high){
  percentages <- apply(mat, 1, calc_percent_exp_onegene)
  goodgenes <- which(percentages > low & percentages < high)
  filtered_mat <- mat[goodgenes,]
  return(filtered_mat)
}