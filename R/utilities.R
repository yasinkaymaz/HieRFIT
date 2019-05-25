


#' Homology mapping via orhologous genes between mouse and rat.
Gmor <- function(RatGenes){
  # This function retrieves mouse homolog associated gene names of Rat genes.
  #library(biomaRt)
  ensembl.rat <- biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  R2M.ort <- biomaRt::getBM(attributes = c("external_gene_name",
                                  "mmusculus_homolog_associated_gene_name",
                                  "mmusculus_homolog_orthology_confidence",
                                  "mmusculus_homolog_orthology_type",
                                  "ensembl_gene_id",
                                  "mmusculus_homolog_ensembl_gene"),
                   filters = 'external_gene_name',
                   values = RatGenes,
                   uniqueRows = T,
                   mart = ensembl.rat)
  R2M.ort <- R2M.ort[which(R2M.ort$mmusculus_homolog_orthology_confidence == "1"), ]
  R2M.ort <- R2M.ort[which(R2M.ort$mmusculus_homolog_orthology_type == "ortholog_one2one"), ]
  return(R2M.ort)
}


RAMmerger <- function(RatObj, MouseObj){
  ort <- Gmor(RatGenes = rownames(RatObj@data))
  mm.data <- as.matrix(MouseObj@raw.data[which(rownames(MouseObj@raw.data) %in% ort$mmusculus_homolog_associated_gene_name),])
  rownames(mm.data) <- ort[match(rownames(mm.data), ort$mmusculus_homolog_associated_gene_name),]$external_gene_name
  rownames(mm.data) <- make.names(rownames(mm.data), unique = T)
  Mouse <- SeuratWrapper(ExpData = mm.data, ProjectLabel = "Mouse_data", NewMeta = MouseObj@meta.data, Normalize = T, dump.files = F)

  Rat.data <- as.matrix(RatObj@raw.data[ort[which(ort$mmusculus_homolog_associated_gene_name %in% rownames(zeisub.data)),]$external_gene_name,])
  rownames(Rat.data) <- make.names(rownames(Rat.data), unique = T)
  Rat <- SeuratWrapper(ExpData = Rat.data, ProjectLabel = "Rat_data",  NewMeta = RatObj@meta.data, Normalize = T, dump.files = F)

  ccaMergedMouseRat <- SeuratCCAmerger(listofObjects = c(Mouse, Rat))

  return(ccaMergedMouseRat)
}


TipBias <- function(tree, confmat){
  #confmat is table() of Prior/Prediction comparison
  err.rate <- NULL
  for(i in 1:length(tree$tip.label)){
    nn <- length(GetAncestorsPath(tree=tree, i)[[2]])
    leafErr <- 1-confmat[tree$tip.label[i],tree$tip.label[i]]/sum(confmat[tree$tip.label[i],])
    print(paste(i,nn,leafErr,sep = "    "))
    err.rate <- c(err.rate, leafErr)
  }
  err.rate <- data.frame(err.rate)
  rownames(err.rate) <- tree$tip.label
  return(err.rate)
}


#' A wrapper function for quick data processing with Seurat functions
#' This function allows to determine highly variable genes and scale their expression,
#' run PCA, tSNE, and cluster detection.
#' @param SeuratObj a Seurat S4 object.
#' @param scale.only.var a boolean variable to determine whether to scale entire data or only highly variable genes. Default is True.
#' @param PCs number of PCs to be included in the downstream analysis
#' @param vars2reg variables to be regressed out when scaling the expression data.
#' @param perp perplexity parameter to be passed to RunTSNE function from Seurat.
#' @keywords seurat quick
#' @export
#' @usage pbmc <- QuickSeurat(pbmc, scale.only.var=F, PCs=5, perp=20)
QuickSeurat <- function(SeuratObj, scale.only.var=T, PCs=20, perp=30, vars2reg) {

  SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)
  hv.genes <- head(rownames(SeuratObj@hvg.info), 1000)
  if (scale.only.var == TRUE) {
    if (!missing(vars2reg)) {
      SeuratObj <- ScaleData(SeuratObj, vars.to.regress = vars2reg, genes.use = hv.genes, do.par=T, num.cores = 8)
    }else{
      SeuratObj <- ScaleData(SeuratObj, genes.use = hv.genes, do.par=T, num.cores = 8)
    }
  }else{
    if (!missing(vars2reg)) {
      SeuratObj <- ScaleData(SeuratObj, vars.to.regress = vars2reg, do.par=T, num.cores = 8)
    }else{
      SeuratObj <- ScaleData(SeuratObj, do.par=T, num.cores = 8)
    }
  }
  SeuratObj <- RunPCA(SeuratObj, pc.genes = hv.genes, do.print = FALSE, pcs.compute=PCs)
  SeuratObj <- FindClusters(SeuratObj, reduction.type = "pca", dims.use = 1:PCs, resolution = 1, print.output = FALSE, save.SNN = TRUE, force.recalc = T)
  SeuratObj <- RunTSNE(SeuratObj, dims.use = 1:PCs, do.fast = TRUE, check_duplicates = FALSE)

  return(SeuratObj)

}

#' A function to downsample Seurat object based on cell type identity
#' @param SeuratObj a Seurat S4 object.
#' @param IdentityCol the column 'number' in the metadata slot showing the cell type identities.
#' @usage pbmc1 <- DownSizeSeurat(SeuratObj = pbmc, IdentityCol = 7)
DownSizeSeurat <- function(SeuratObj, IdentityCol, min_n=NULL){
  cells <- NULL
  classes <- table(SeuratObj@meta.data[,IdentityCol])
  print(classes)
  if(is.null(min_n)){
    min_n <- min(classes)
    print(min_n)
  }
  for(type in names(classes)){
    if( classes[type] > min_n ){
      cells <- c(cells, sample(rownames(SeuratObj@meta.data[which(SeuratObj@meta.data[,IdentityCol] == type), ]), size = min_n, replace = F))
    }else{
      cells <- c(cells, sample(rownames(SeuratObj@meta.data[which(SeuratObj@meta.data[,IdentityCol] == type), ]), size = classes[type], replace = F))
    }
  }
  downSobj <- SubsetData(object = SeuratObj, cells.use = cells, do.clean=T)
  return(downSobj)
}

SeuratWrapper <- function(ExpData, ProjectLabel, NewMeta, Normalize=T, suppressLog=F, scale.only.var=T, PCs=20, perp=30, dump.files=F, min.cells=0, min.genes=0) {

  if (Normalize == TRUE) {print("Assuming the input is in count ...")
  }else{
    print("Assuming the input is in TPM ...")
    if (suppressLog == TRUE) {
      print("not taking log ...")
    }else{
      ExpData <- log1p(ExpData)
    }
  }

  SeuratObj <- CreateSeuratObject(raw.data = ExpData, project = ProjectLabel, min.cells=min.cells, min.genes = min.genes)

  if (Normalize == TRUE) {
    SeuratObj <- NormalizeData(object = SeuratObj)
  }else{
    print("Not normalizing the data since TPM is assumed ... ")
  }

  SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)
  hv.genes <- head(rownames(SeuratObj@hvg.info), 1000)

  if (scale.only.var == TRUE) {
    SeuratObj <- ScaleData(SeuratObj, genes.use = hv.genes, do.par=T, num.cores = 8)
  }else{
    SeuratObj <- ScaleData(SeuratObj, do.par=T, num.cores = 8)
  }

  if (!missing(NewMeta)) {
    SeuratObj <- AddMetaData(SeuratObj, NewMeta[rownames(SeuratObj@meta.data), ])
  }else{
    print("No new meta file is provided. Skipping...")
  }

  SeuratObj <- RunPCA(SeuratObj, pc.genes = hv.genes, do.print = FALSE, pcs.compute=PCs)
  SeuratObj <- FindClusters(SeuratObj, reduction.type = "pca", dims.use = 1:PCs, resolution = 1, print.output = FALSE, save.SNN = TRUE, force.recalc = T)
  SeuratObj <- RunTSNE(SeuratObj, dims.use = 1:PCs, do.fast = TRUE,check_duplicates = FALSE, perplexity=perp)

  if (dump.files == T) {
    #Export the tSNE coordinates along with the Cluster assignment IDs
    rownames_to_column(as.data.frame(SeuratObj@dr$tsne@cell.embeddings))  %>%
      as.tibble() %>%
      add_column(Clusters=SeuratObj@meta.data$res.1) %>%
      dplyr::rename(Cellname = rowname) %>%
      as_data_frame() %>% write_csv(paste(ProjectLabel, "_tSNECoordinates_Clusters.csv", sep=""))

    #Export Normalized and Scaled Expression matrix for cells and genes in the analysis
    rownames_to_column(as.data.frame(as.matrix(SeuratObj@data))) %>%
      dplyr::rename(GeneName = rowname) %>%
      as_data_frame() %>%
      write_delim(paste(ProjectLabel, "_Normalized_Expression_matrix.txt", sep=""))

    #Save seurat object:
    save(SeuratObj, file=paste(ProjectLabel, ".seurat.Robj", sep=""))
  }

  return(SeuratObj)
}

SeuratCCAmerger <- function(listofObjects) {
  # Determine genes to use for CCA, must be highly variable in at least 2 datasets
  #ob.list <- list(zeisel, romanov, tasic, marques)
  ob.list <- listofObjects
  genesuse <- c()
  ids=NULL
  for (i in 1:length(ob.list)) {
    genesuse <- c(genesuse, head(rownames(ob.list[[i]]@hvg.info), 1000))
    ob.list[[i]]@meta.data$dataSource <- paste("id", i, sep="")
    ids <- c(ids, paste("id", i, sep=""))
  }
  genesuse <- names(which(table(genesuse) > 1))
  for (i in 1:length(ob.list)) {
    genesuse <- genesuse[genesuse %in% rownames(ob.list[[i]]@scale.data)]
  }

  if (length(ob.list) > 2) {
    # Run multi-set CCA
    integrated <- RunMultiCCA(ob.list, genes.use = genesuse, num.ccs = 15, add.cell.ids = ids)
    # Run rare non-overlapping filtering
    integrated <- CalcVarExpRatio(object = integrated, reduction.type = "pca", dims.use = 1:10, grouping.var = "dataSource")
    integrated <- SubsetData(integrated, subset.name = "var.ratio.pca", accept.low = 0.5)
  }else{
    #integrated <- RunCCA(object = ob.list[[1]], object2 = ob.list[[2]], genes.use = genesuse, num.cc = 15, add.cell.id = ids)
    integrated <- RunCCA(object = ob.list[[1]], object2 = ob.list[[2]], genes.use = genesuse, num.cc = 15)
  }
  # Alignment
  integrated <- AlignSubspace(integrated, reduction.type = "cca", dims.align = 1:10, grouping.var = "dataSource")
  # t-SNE and Clustering
  integrated <- FindClusters(integrated, reduction.type = "cca.aligned", dims.use = 1:10, save.SNN = T, resolution = 0.4)
  integrated <- RunTSNE(integrated, reduction.use = "cca.aligned", dims.use = 1:10)
  save(integrated, file="integrated.Aligned.seurat.Robj")
  return(integrated)
}
