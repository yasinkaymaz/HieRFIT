

ParApp <- function() {
  '"
  Sets the parameters used throughout the run.
  "'
  PCs <- 20 # number of principle components to be analyzed.
  Pred_n <- 1000 # number of predictors aimed to be collected.
  prefix <- ""
}


NodesAcc <- function(HieRMod, ...){
  #Extract the node accuracy metrics:
  nodeStats <- NULL
  for(i in names(HieRMod@model)){
    mtry <- HieRMod@model[[as.character(i)]]$finalModel$mtry
    nodeStats <- rbind(nodeStats,
                       cbind(node=i,
                             HieRMod@model[[as.character(i)]]$results[which(HieRMod@model[[as.character(i)]]$results$mtry == mtry),],
                             NodeLabel=HieRMod@tree[[1]]$node.label[ as.numeric(i) - length(HieRMod@tree[[1]]$tip.label)],
                             classSize=length(HieRMod@model[[as.character(i)]]$levels)))
  }
  nodeStats <- nodeStats[which(nodeStats$NodeLabel %in% HieRMod@tree[[1]]$node.label), ]
  rownames(nodeStats) <- nodeStats$NodeLabel
  nodeAcc <- round(nodeStats$Accuracy*100, digits = 1)
  names(nodeAcc) <- nodeStats$NodeLabel

  return(nodeAcc)
}

#' This function calculates a scaled Kullback-Leibler divergence
#' @param probs a list of observed probability scores.
#' @return KLscaled = KLe/KLmax, where KLe is the empirical divergence given
#' the distributions while KLmax is the maximum KLe value that can be achieved
#' given the number of items.
KLeCalc <- function(probs){
  class_n <- length(probs)
  null <- c(1, rep(0, class_n-1))
  KLe <- entropy::KL.empirical(y1 = as.numeric(probs), y2 = rep(1/class_n, class_n))
  KLmax <- entropy::KL.empirical(y1 = as.numeric(null), y2 = rep(1/class_n, class_n))
  KLscaled <- KLe/KLmax
  return(KLscaled)
}


#' A function to calculate Asymmetric Entropy
#' @param p measured probability array
#' @param w empirically calculated W set
#' @return U set of certainty values for each probility outcome in p given w.
GetCertaintyArray <- function(p, w){
  U <- numeric()
  for(i in 1:length(p)){
    w[i] <- w[i]+1e-10#to prevent math err.
    if(p[i] > w[i]){
      lambda=1
    }else{
      lambda=-1
    }
    U <- c(U, (lambda*(p[i]-w[i])^2)/( ((1-2*w[i])*p[i])+w[i]^2 ) )
  }
  U <- as.numeric(U)
  return(U)
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



#' A function to generate a random tree using 'ape' package.
#' Returns a tree object.
#' @param LN The number of leafs to be in the tree. Default is 8.
RandTreeSim <- function(LN=8, furcation="binary"){
  if (furcation == "binary"){
    tree <- ape::rtree(n = LN, br = NULL)
    plot(tree, edge.width = 2)
    tree$edge
  } else if (furcation == "multi"){#Fix this.
    tiplabs <- paste("t", seq(1:LN), sep = "")
    while (length(tiplabs) > 0){
      sub <- sample(tiplabs, size = sample(seq(1:length(tiplabs)-1)), replace = F)
      print(sub)
      tiplabs <- tiplabs[which(!tiplabs %in% sub)]
    }
    tree <- ape::read.tree(text="(((L, K), E, F), (G, H));")
  }
  return(tree)
}

#' An internal function to prepare a training/test dataset for model generation.
#' @param Data a Normalized expression data matrix, genes in rows and samples in columns.
#' @param Predictors the predictor feature list selected by FeatureSelector.
#' @param ClassLabels [optional] A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param alpa variation cutoff for filtering data
#' @keywords data preparation
#' @export
#' @usage trainingData <- DataReshaper(Data = as.matrix(SeuratObject@data), Predictors = genes, ClassLabels = SeuratObject@meta.data$CellTypes)
DataReshaper <- function(Data, Predictors, ClassLabels, alpa=0.1, ...) {
  if(missing(Predictors)){#prepare the data for PCA
    TData <- as.data.frame(t(Data))
    indx <- sapply(TData, is.factor)
    TData[indx] <- lapply(TData[indx], function(x) as.numeric(as.character(x)))
    #Filter candidate predictors before PCA
    TData <- TData[, apply(TData, 2, var) != 0]
    TData <- droplevels(TData[, which(matrixStats::colVars(as.matrix(TData)) > alpa)])
    return(TData)

  } else {

    if (missing(ClassLabels)) {
      #Then the output is for Query
      QueData <- as.data.frame(t(Data), col.names=rownames(Data))
      colnames(QueData) <- make.names(colnames(QueData))
      QueData_sub <- droplevels(QueData[, which(colnames(QueData) %in% Predictors)])
      #missing Predictors
      mp <- Predictors[which(!Predictors %in% colnames(QueData))]
      #Add missing predictors into QueData by setting to 0.
      mp_df <- data.frame(matrix(0,
                                 ncol = length(mp),
                                 nrow = length(colnames(Data))))
      colnames(mp_df) <- mp
      QueData <- cbind(QueData_sub, mp_df)
      QueData <- QueData[, Predictors]
      return(QueData)

    } else {

      RefData <- as.data.frame(t(Data), col.names=rownames(Data))
      colnames(RefData) <- make.names(colnames(RefData))
      #convert factors to numeric
      indx <- sapply(RefData, is.factor)
      RefData[indx] <- lapply(RefData[indx], function(x) as.numeric(as.character(x)))
      RefData$ClassLabels <- factor(make.names(ClassLabels))
      RefData <- droplevels(RefData[, c(Predictors, "ClassLabels")])

      return(RefData)
    }
  }
}

#' A wrapper function for random forest.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(RefData).
#' @param mod.meth The model training method, "rf" for random forest.
RandForestWrap <- function(RefData, ClassLabels, prefix, mod.meth, train.control, thread=5, ...){
  #library(doParallel)
  library(caret)
  library(doMC)
  registerDoMC(cores = thread)
  #1. Select the predictors.
  P_dicts <- FeatureSelector(Data = RefData,
                             ClassLabels = ClassLabels,
                             PCs = 10,
                             num = 200,
                             prefix = prefix,
                             doPlots = F, ...)
  #2. Prepare the reference data.
  TrainData <- DataReshaper(Data = RefData,
                            Predictors = P_dicts,
                            ClassLabels = ClassLabels, ...)
  print(RefData[1:5, 1:5])
  #3. Train the model.
  #cl <- makePSOCKcluster(ncor)
  #registerDoParallel(cl)
  model <- caret::train(ClassLabels~.,
                        data = TrainData,
                        method = mod.meth,
                        norm.votes = TRUE,
                        importance = FALSE,
                        proximity = FALSE,
                        outscale = FALSE,
                        preProcess = c("center", "scale"),
                        ntree=50,
                        trControl = train.control)
  #stopCluster(cl)

  return(model)
}

#' A wrapper function for random forest.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param mod.meth The model training method, "svmLinear" for support vector machine.
SvmWrap <- function(RefData, ClassLabels, prefix, mod.meth, train.control, ncor=5, ...){
  library(doParallel)
  #1. Select the predictors.
  P_dicts <- FeatureSelector(Data = RefData,
                             ClassLabels = ClassLabels,
                             PCs = 10,
                             num = 200,
                             prefix = prefix,
                             doPlots = F, ...)
  #2. Prepare the reference data.
  TrainData <- DataReshaper(Data = RefData,
                            Predictors = P_dicts,
                            ClassLabels = ClassLabels, ...)

  #3. Train the model.
  cl <- makePSOCKcluster(ncor)
  registerDoParallel(cl)
  model <- caret::train(ClassLabels~.,
                        data = TrainData,
                        trControl = train.control,
                        method = mod.meth,
                        norm.votes = TRUE,
                        importance = TRUE,
                        proximity = TRUE,
                        preProcess = c("center", "scale"),
                        tuneLength = 10, ...)
  stopCluster(cl)

  return(model)
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
  if(class(SeuratObj)[1] == "seurat"){
    downSobj <- SubsetData(object = SeuratObj, cells.use = cells, do.clean=T)
  }else if(class(SeuratObj)[1] == "Seurat"){
    downSobj <- subset(x = SeuratObj, cells = cells)
  }
  return(downSobj)
}

#' A function to downsample a refdata table based on classLabels
#' @param RefData a data table with features as columns (last column being ClassLabels), instances in the rows.
#' @param IdentityCol the name of the column in the refdata storing class labels. default is "ClassLabels"
#' @param min_n min number of samples to downsample each class. default is the size of the minority class.
#' @usage RefData_d <- DownSampleRef(RefData = RefData)
DownSampleRef <- function(RefData, IdentityCol="ClassLabels", min_n=NULL){
  samples <- NULL
  classes <- table(RefData[, IdentityCol])
  print(classes)
  if(is.null(min_n)){
    min_n <- min(classes)
    print(min_n)
  }
  for(type in names(classes)){
    if( classes[type] > min_n ){
      samples <- c(samples, sample(rownames(RefData[which(RefData[, IdentityCol] == type), ]), size = min_n, replace = F))
    }else{
      samples <- c(samples, sample(rownames(RefData[which(RefData[, IdentityCol] == type), ]), size = min_n, replace = T))
    }
  }
  RefData_d <- RefData[samples, ]
  return(RefData_d)
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

FlatRF <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions") {

  library(caret)
  library(randomForest)
  library(tidyverse)

  if (!missing(SeuratObject)) {
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{
    print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
  }#Closes missing(SeuratObj)
  colnames(testExpSet) <- make.names(colnames(testExpSet))
  #Prepare Test Expression set
  testsub <- testExpSet[,which(colnames(testExpSet) %in% model$finalModel$xNames)]
  missingGenes <- model$finalModel$xNames[which(!model$finalModel$xNames %in% colnames(testExpSet))]
  print(model$finalModel$importance[missingGenes,])

  missingGenes.df <- data.frame(matrix(0, ncol = length(missingGenes), nrow = length(rownames(testExpSet))))
  colnames(missingGenes.df) <- missingGenes
  TestData <- cbind(testsub, missingGenes.df)
  TestData <- TestData[,model$finalModel$xNames]
  mmDGm <- mean(model$finalModel$importance[missingGenes,])
  mmDGf <- mean(model$finalModel$importance[which(!model$finalModel$xNames %in% missingGenes),])

  cat("Number of Features (genes) to be considered is", length(colnames(testsub)), '\n',
      "Number of missing Features set to zero is", length(missingGenes), '\n',
      "Mean MeanDecreaseGini of missing genes is", mmDGm, '\n',
      "Mean MeanDecreaseGini of Featured genes is", mmDGf, '\n',
      sep = ' ')

  if (!is.nan(mmDGm)) {
    if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
      warning("A significant portion of features are missing...")
    }
  }


  rm(testsub, missingGenes, missingGenes.df)
  gc()
  #Predict
  library(entropy)
  testPred <- as.data.frame(predict(model, TestData, type = "prob"))
  class.n <- length(model$finalModel$classes)
  testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
  testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class.n, class.n)) )
  testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
  testPred$Prediction <- predict(model, TestData, type="raw")
  #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class.n)
  testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", "Detected")) %>% as.data.frame()
  #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class.n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

  if (missing(priorLabels)) {
    print("Prior class labels are not provided!")

  }else{
    #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    testPred <- cbind(testPred, priorLabels)

    #Plot the crosscheck here:
    #Crosscheck Predictions
    library(tidyverse)
    library(alluvial)
    library(ggalluvial)
    crx <- testPred %>% group_by(Prior, Intermediate, Prediction) %>% tally() %>% as.data.frame()

    p5 <- ggplot(crx,aes(y = n, axis1 = Prior, axis2 = Intermediate, axis3 = Prediction )) +
      geom_alluvium(aes(fill = Prediction), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("Prior", "Clusters", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class.n, base_width = 1.2*class.n)
  }

  if (!missing(SeuratObject)) {
    #update predictions in the meta.data slot
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(testPred))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, testPred)

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function


#' An internal function to collect model training parameters and direct them to model creation.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param mod.meth The model training for hierarchical random forest. Default is "hrf"
#' @param cv.k Fold cross validation. Default is 5.
#' @param tree A tree storing relatinship between the class labels. Default is null. required when mod.meth "hrf" is chosen.
#' @param save.int.f Boolean to save model. Default is False.
#' @keywords
#' @export
#' @usage rf.model <- Modeller(RefData = as.matrix(pbmc1@data), ClassLabels = pbmc1@meta.data$ClusterNames_0.6)
Modeller <- function(RefData, ClassLabels=NULL, mod.meth="rf", thread=NULL, tree=NULL, save.int.f=FALSE, ...){
  if(mod.meth == "hrf"){
    try(if(missing(tree)|| missing(ClassLabels) || missing(RefData))
      stop("Please, provide the required inputs!"))
    model <- HieRandForest(RefData = RefData,
                           ClassLabels = ClassLabels,
                           tree, thread = thread)
  }

  return(model)
}
