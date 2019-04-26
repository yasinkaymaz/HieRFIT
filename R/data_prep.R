
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
#' @examples pbmc <- QuickSeurat(pbmc, scale.only.var=F, PCs=5, perp=20)
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

#' A function used internally for selecting genes based on their weight in the top principle components.
#' @description p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom.
#' Weighted gene picking depending on PC number: Initial PCs give more genes.
#' For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
#' @param ExpData an data matrix storing gene expression as genes in rows and samples in columns.
#' @param PCs the number of PCs to be looked at when selecting genes.
#' @param num the number of genes in total to be included.
#' @param spitPlots boolean to generate a pdf for PCA plots. Default is False.
#' @param prefix a prefix to separate each run.
#' @keywords PCA loadings selection
#' @export
#' @examples genes <- as.character(SelectGenesBestLoadings(ExpData = train, prefix = "test.run", pcs = 10, num = 2000))
SelectGenesBestLoadings <- function(ExpData, PCs, num, spitPlots=F, prefix="gene.select") {

  print("Performing PCA...")
  #Do PCA here on input data

  ExpData <- ExpData[, apply(ExpData, 2, var) != 0]
  ExpData <- ExpData[, which(matrixStats::colVars(as.matrix(ExpData)) > 0.05)]

  if (file.exists(paste(prefix,".train.prcomp.Rdata",sep=""))) {
    warning("An existing PCA data is found in the directory! Used it instead...")
    pcatrain <- get(load(paste(prefix,".train.prcomp.Rdata",sep="")))
  }else{
    pcatrain <- prcomp(ExpData, center = TRUE, scale=TRUE, rank. = PCs)
    save(pcatrain,file=paste(prefix,".train.prcomp.Rdata",sep=""))
  }

  print("Selecting the genes as best features...")

  trainingData <- get(load(paste(prefix,".trainingData.tmp.Rdata",sep = "")))
  pcadata <- data.frame(pcatrain$x, CellType = trainingData$CellType)
  class_n=length(levels(pcadata$CellType))

  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  plots <- list()
  pb <- txtProgressBar(min = 0, max = PCs, style = 3)
  for(i in 1:(PCs-1)) {
    orderedpcai <- pcatrain$rotation[order(abs(pcatrain$rotation[,i]),decreasing = TRUE),i]
    Gn <- round((PCs-i+1)*(num*2)/(PCs*(PCs+1)))
    TotalGenes <- as.numeric(TotalGenes) + Gn
    top <- data.frame(genes=names(head(orderedpcai,Gn)),bestgenes=head(orderedpcai,Gn))
    load <- rbind(load, top)
    setTxtProgressBar(pb,i+1)
    pci <- paste("PC",i,sep="")
    pcj <- paste("PC",i+1,sep="")
    #pc.p <- ggplot(pcadata, aes_string(x=pci, y=pcj, color="CellType"))+geom_point()
    pc.p <- ggplot(pcadata, aes_string(x=pci, y="CellType", color="CellType"))+geom_point()
    plots[[i]] <- pc.p
    cat("\n",'Picking the best genes from first', i+1,'PCs is completed.',"\n")
  }
  if (spitPlots == TRUE) {
    pdf(paste(prefix,"_PCAplots.pdf",sep=""),width = 10,height = class_n*PCs/2)
    multiplot(plotlist = plots)
    dev.off()
  }
  bestgenes <- unique(load$genes)
  return(bestgenes)
}


#' A function used internally for selecting genes based on their weight in the top principle components.
#' @description p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom.
#' Weighted gene picking depending on PC number: Initial PCs give more genes.
#' For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
#' @param ExpData an data matrix storing gene expression as genes in rows and samples in columns.
#' @param PCs the number of PCs to be looked at when selecting genes.
#' @param num the number of genes in total to be included.
#' @param spitPlots boolean to generate a pdf for PCA plots. Default is False.
#' @param prefix a prefix to separate each run.
#' @keywords PCA loadings selection
#' @export
#' @examples genes <- as.character(SelectGenesBestLoadings(ExpData = train, prefix = "test.run", pcs = 10, num = 2000))
FeatureSelector <- function(ExpData, PCs, num, spitPlots=F, prefix="gene.select") {

  print("Performing PCA...")
  #Do PCA here on input data
  #Select genes for PCA
  ExpData <- ExpData[, apply(ExpData, 2, var) != 0]
  ExpData <- ExpData[, which(matrixStats::colVars(as.matrix(ExpData)) > 0.1)]

  if (file.exists(paste(prefix,".train.prcomp.Rdata",sep=""))) {
    warning("An existing PCA data is found in the directory! Used it instead...")
    pcatrain <- get(load(paste(prefix,".train.prcomp.Rdata",sep="")))
  }else{
    pcatrain <- prcomp(ExpData, center = TRUE, scale=TRUE, rank. = PCs)
    save(pcatrain,file=paste(prefix,".train.prcomp.Rdata",sep=""))
  }

  print("Selecting the genes as best features...")

  trainingData <- get(load(paste(prefix,".trainingData.tmp.Rdata",sep = "")))
  pcadata <- data.frame(pcatrain$x, CellType = trainingData$CellType)
  class_n=length(levels(pcadata$CellType))
  #Select only PCs which are significanly separating at least one of the class.
  is_significant <- function(x) any(as.numeric(x) <= 0.05)

  #Create a table for evaluation of all PCs with a significance test.
  PCs.sig.table <- pcadata %>% group_by(CellType) %>%
    select(paste("PC",1,sep = ""):paste("PC",PCs,sep = "")) %>%
    summarize_all(list(~t.test(x = .)$p.value)) %>%
    as.data.frame()
  rownames(PCs.sig.table) <- PCs.sig.table$CellType
  PCs.sig.table <- PCs.sig.table[,-c(1)]


  PCs.sig <- pcadata %>% group_by(CellType) %>%
    select(paste("PC",1,sep = ""):paste("PC",PCs,sep = "")) %>%
    summarize_all(list(~t.test(x = .)$p.value)) %>%
    select_if(., is_significant) %>% colnames()

  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  plots <- list()
  pb <- txtProgressBar(min = 0, max = PCs, style = 3)
  for (i in 1:length(PCs.sig)){
    orderedpcai <- pcatrain$rotation[order(abs(pcatrain$rotation[,PCs.sig[i]]),decreasing = TRUE),PCs.sig[i]]
    Gn <- round((PCs-i+1)*(num*2)/(PCs*(PCs+1)))
    print(Gn)
    TotalGenes <- as.numeric(TotalGenes) + Gn
    top <- data.frame(genes=names(head(orderedpcai,Gn)),bestgenes=head(orderedpcai,Gn))
    load <- rbind(load, top)
    setTxtProgressBar(pb,i+1)
    pci <- paste("PC",i,sep="")
    pcj <- paste("PC",i+1,sep="")
    pc.p <- ggplot(pcadata, aes_string(x=pci, y="CellType", color="CellType"))+geom_point()
    plots[[i]] <- pc.p
    cat("\n",'Picking the best genes from PC', PCs.sig[i],' is completed.',"\n")
  }

  if (spitPlots == TRUE) {
    pdf(paste(prefix,"_PCA_eval_p.pdf",sep=""),width = class_n/2, height = class_n/2)
    pheatmap::pheatmap(-1*log10(PCs.sig.table+.000001),cluster_cols = F,treeheight_row = 0,cellheight = 10,cellwidth = 10)
    dev.off()
    pdf(paste(prefix,"_PCA_class.pdf",sep=""),width = 10,height = class_n*PCs/2)
    multiplot(plotlist = plots)
    dev.off()
  }
  bestgenes <- unique(load$genes)
  return(bestgenes)
}


#' An internal function to prepare a training/test dataset for model generation.
#' @param ExpData a Normalized expression data matrix, genes in rows and samples in columns. Example: .
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param PCs the number of PCs to be looked at when selecting genes.
#' @param featureGeneSet an optional list of genes to be used in the training matrix rather than PC based selection. If this is provided, PC selection is skipped.
#' @param prefix a prefix to separate each run.
#' @keywords Training preparation
#' @export
#' @examples trainingData <- prepareDataset(ExpData = as.matrix(SeuratObject@data), ClassLabels = SeuratObject@meta.data$CellTypes, PCs = 10)
prepareDataset <- function(ExpData, ClassLabels, PCs, featureGeneSet, prefix="prep.run", ...) {
  #Transpose the matrix cols <--> rows t()
  #Keep the data in matrix form, otherwise randomForest will throw error: 'Error: protect(): protection stack overflow'
  trainingData <- as.data.frame(t(ExpData))
  #It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
  names(trainingData) <- make.names(names(trainingData))
  trainingData$CellType <- ClassLabels
  #Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
  trainingData$CellType <- factor(trainingData$CellType)

  train <- droplevels(trainingData[,!(names(trainingData) %in% c("CellType"))])
  indx <- sapply(train, is.factor)
  train[indx] <- lapply(train[indx], function(x) as.numeric(as.character(x)))

  save(trainingData, file=paste(prefix,".trainingData.tmp.Rdata",sep = ""))
  rm(trainingData)
  gc()

  if (missing(featureGeneSet)) {
    #Perform PCA on data (optional: only on training portion) and select genes for features
    #genes <- as.character(SelectGenesBestLoadings(ExpData = train, prefix = prefix, PCs = PCs, num = 2000, ...))
    genes <- as.character(FeatureSelector(ExpData = train, prefix = prefix, PCs = PCs, num = 1000, ...))
  }else{
    genes <- make.names(featureGeneSet)
  }#closes.if.missing.featureset

  trainingData <- get(load(paste(prefix,".trainingData.tmp.Rdata",sep = "")))

  trainingData.postPCA <- droplevels(trainingData[,c(genes,"CellType")])

  save(trainingData.postPCA, file=paste(prefix, ".trainingData.postPCA.data", sep = ""))
  file.remove(paste(prefix, ".trainingData.tmp.Rdata", sep = ""))

  return(trainingData.postPCA)
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

  pdf(paste(ProjectLabel,".plots.pdf", sep=""), width=8, height = 8)
  PCAPlot(SeuratObj, dim.1 = 1, dim.2 = 2)
  PCElbowPlot(SeuratObj, num.pc = PCs)
  TSNEPlot(SeuratObj, do.label = TRUE)
  dev.off()

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
  }

  save(SeuratObj, file=paste(ProjectLabel, ".seurat.Robj", sep=""))
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

