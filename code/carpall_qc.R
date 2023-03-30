#### Load Libraries ####
library(Seurat)
library(plyr);library(dplyr)
library(ggplot2)
library(reshape2)

#' Basic QC and filtering functions.  
#' Initial filtering and metadata
#'
#' Calculates metadata, drop cells based on QC metrics.  Assumes Seurat V4.
#'
#' @param data dgCMatrix (i.e. output from Read10X function). 
#' @param ID Identifier to save plots and objects. 
#' @param sPhaseGenes Genes used to identify S phase genes.  If NULL, will be loaded from Seurat.
#' @param g2mPhaseGenes Genes used to identify G2M phase genes.  If NULL, will be loaded from Seurat.
#' @param percentmt Maximum fraction of expression from Mitochondrial genes allow per cell.
#' @param nFeature_RNA_min Minimum number of genes for a cell to express and pass QC.
#' @param nCountRNA_min Minimum number of UMIs for a cell to have and pass QC.
#' @param npcs1 Number of PCs to use for the rough clustering. 
#' @param clustres Resolution for fine clustering.
#' @param finalNN This determines the number of neighboring points used in local approximations of manifold structure.
#' @param subCAR Subset the object for CAR T-cells only. If subCAR=TRUE, subset the object. 
#' @param saverds Do you want to save the final seurat object as an .rds file? If saverds=YES, save the file. 
#' @return A Seurat object with cells dropped that fail the QC filters.\

##Requirements##
#excludeGenes: Genes that will be dropped from the count matrix. Should be loaded prior to running the function. 
#cart_index: CAR T-cell metadata. Should be loaded prior to running the function. 
#barcode_metadata: Dataframe of barcodes that are CAR T-cells with Cell Type Labels. Should be loaded prior to running the function. 
excludeGenes <- as.character(read.table("excludeGenes.tsv", header = F, sep = "\t")$V1)
cart_index <- read.table(file="CARPALL_scRNAseq_CART_metadata.txt", header=T, sep="\t", stringsAsFactors = F)
barcode_metadata <- read.table(file="CARTcell_barcodes.txt", header=T, sep="\t", stringsAsFactors = F)

single_cell_processing <- function(data, ID, nFeature_RNA_min, percentmt, nCountRNA_min, npcs1, clustres, finalMinDist, finalNN, subCAR, saverds){
  
  #### Read Count Matrix ####
  data_seurat <- CreateSeuratObject(data, project = "SeuratProject", assay = "RNA", names.field = 1, names.delim = "_", meta.data = NULL)
  
  #Calculate % MT content 
  data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
  
  #Add CellID to Metadata 
  data_seurat@meta.data$CellID <- rownames(data_seurat@meta.data)
  RMS_meta <- data_seurat@meta.data
  
  #### Standard Pre-processing ####
  keeprows <- RMS_meta[RMS_meta$nCount_RNA > nCountRNA_min,] #Number of UMIs must be greater than nCountRNA_min
  keeprows <- keeprows[keeprows$nFeature_RNA > nFeature_RNA_min,] #Number of Features must be greater than nFeature_RNA_min
  keeprows <- keeprows[keeprows$percent.mt <= percentmt,] #Mitochondrial Content should not be greater than percent.mt
  cellids_keep <- rownames(keeprows)
  data_filtered <- subset(data_seurat, subset = CellID %in% cellids_keep)
  
  #Subset for CAR T-cells 
  if(subCAR == TRUE){
    data_filtered <- subset(data_filtered, subset = CellID %in% barcode_metadata$CellID)
  } else {}

  #Exclude Ribosomal, Mitochondrial and other unimportant Genes
  keep_features=rownames(data_seurat)[which(!rownames(data_seurat)%in%excludeGenes)]
  data_filtered <- subset(data_filtered, features=keep_features)
  data_filtered = NormalizeData(object = data_filtered, normalization.method = "LogNormalize", scale.factor = 1e4)
  data_filtered = FindVariableFeatures(data_filtered, selection.method = "vst", nfeatures = 2000)
  #RMS_scaled <- ScaleData(RMS_normalized, features = rownames(RMS_normalized), verbose=T)
  data_filtered = ScaleData(data_filtered, verbose = T) #Memory issues if we scale using all genes. 
  data_filtered = RunPCA(data_filtered, npcs = npcs1, features=VariableFeatures(object = data_filtered), verbose = F)
  
  data_filtered = FindNeighbors(data_filtered, dims = 1:npcs1, verbose = T)
  data_filtered = FindClusters(data_filtered, resolution = clustres, verbose = T)
  data_filtered = RunUMAP(data_filtered, dims = 1:npcs1, verbose = T, min.dist = finalMinDist, n.neighbors = finalNN)
  
  s.genes = cc.genes.updated.2019$s.genes
  g2m.genes = cc.genes.updated.2019$g2m.genes
  data_filtered = CellCycleScoring(data_filtered, s.features = s.genes, 
                              g2m.features = g2m.genes, set.ident = TRUE)
  
  ## Add Metadata ##
  colnames(cart_index)[2] <- "orig.ident"
  data_filtered_meta <- data_filtered@meta.data
  getindex <- join(data_filtered_meta, cart_index, by="orig.ident")
  data_filtered@meta.data$Patient <- getindex$Patient
  data_filtered@meta.data$TCR_Sample <- getindex$TCR_Sample
  data_filtered@meta.data$Sort <- getindex$Sort
  data_filtered@meta.data$Source <- getindex$Source
  data_filtered@meta.data$Timepoint <- getindex$Timepoint
  data_filtered@meta.data$Timepoint_bin <- getindex$Timepoint_bin

  data_filtered@meta.data$CellType1 <- sapply(1:nrow(data_filtered_meta), function(x) ifelse(nrow(barcode_metadata[barcode_metadata$CellID == data_filtered_meta$CellID[x],])>0,
                                                                                             barcode_metadata[barcode_metadata$CellID == data_filtered_meta$CellID[x],]$CellType1,
                                                                                             NA))
  
  data_filtered@meta.data$CellType2 <- sapply(1:nrow(data_filtered_meta), function(x) ifelse(nrow(barcode_metadata[barcode_metadata$CellID == data_filtered_meta$CellID[x],])>0,
                                                                                             barcode_metadata[barcode_metadata$CellID == data_filtered_meta$CellID[x],]$CellType2,
                                                                                             NA))
  
  data_filtered@meta.data$CARTcell <- sapply(1:nrow(data_filtered_meta), function(x) ifelse(nrow(barcode_metadata[barcode_metadata$CellID == data_filtered_meta$CellID[x],])>0,
                                                                                             barcode_metadata[barcode_metadata$CellID == data_filtered_meta$CellID[x],]$CARTcell,
                                                                                             NA))
                                                                                                  
  Idents(data_filtered) <- data_filtered@meta.data$seurat_clusters
  if(saverds == "YES"){
    saveRDS(data_filtered, file=paste(sep=".", ID, nFeature_RNA_min, percentmt, nCountRNA_min, npcs1,"rds"))
    return(data_filtered)
  } else {
    return(data_filtered)
  }
}

carpall_cart <- single_cell_processing(data, "<id>_<date>", 300, 10, 1000, 75, 1, 0.5, 50, TRUE, "YES")

### Look at CAR T-cells ### 
FeaturePlot(carpall_cart, features="CAT-scFv")


