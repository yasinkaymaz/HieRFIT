#CreateHieR functions:

#' Reference class
#' @slot model A list of models to be used.
#' @slot tree A hierarchical tree representing class relationships.
#' @slot modtype type of the model method used to create object.
RefMod <- setClass(Class = "RefMod",
                   slots = c(model = "list",
                             tree = "list",
                             modtype = "character"))

#' The main function for creating a reference model.
#' @param RefData Reference data from which class labels will be projected on Query data.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix (Ref). Same length as colnames(Ref).
#' @param method The model training method, "rf" for random forest, "svmLinear" for support vector machine, "hrf" for hierarchical random forest. Default is "hrf"
#' @param TreeFile An input file to create a hierarchical tree for relationship between class labels. Default is null but required if method 'hrf' is chosen.
#' @usage  pbmc.refmod <- CreateRef(RefData = as.matrix(pbmc@data), ClassLabels = pbmc@meta.data$ClusterNames_0.6, TreeFile = pbmc3k_tree)
CreateHieR <- function(RefData, ClassLabels, TreeTable=NULL, method="hrf", thread=NULL){

  if(method == "hrf"){
    tree <- CreateTree(treeTable = TreeTable)
  }else{ tree <- NULL}

  ClassLabels <- FixLab(xstring = ClassLabels)

  #Create predictive model structure. No need to create if exists.
  print("Training model... This may take some time... Please, be patient!")
  model <- Modeller(RefData = RefData,
                    ClassLabels = ClassLabels,
                    mod.meth = method,
                    thread = thread,
                    tree = tree)
  refObj <- new(Class = "RefMod",
                model = model,
                modtype = method)
  refObj@tree[[1]] <- tree

  #foreach::registerDoSEQ()
  #Return a Ref object rather than single model.
  return(refObj)
}


#' An internal function to collect model training parameters and direct them to model creation.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param mod.meth The model training method, "rf" for random forest, "svmLinear" for support vector machine, "hrf" for hierarchical random forest. Default is "rf"
#' @param cv.k Fold cross validation. Default is 5.
#' @param tree A tree storing relatinship between the class labels. Default is null. required when mod.meth "hrf" is chosen.
#' @param save.int.f Boolean to save model. Default is False.
#' @keywords
#' @export
#' @usage rf.model <- Modeller(RefData = as.matrix(pbmc1@data), ClassLabels = pbmc1@meta.data$ClusterNames_0.6)
Modeller <- function(RefData, ClassLabels=NULL, mod.meth="rf", cv.k=2, thread=NULL, tree=NULL, save.int.f=FALSE, ...){
  prefix=paste(cv.k, "K-fold", mod.meth, sep=".")
  #k-fold Cross Validation
  train.control <- caret::trainControl(method="oob",
                                       returnData = FALSE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim=TRUE,
                                       sampling = "down")
  #Choose whether flat or hierarchical classifier.
  if(mod.meth == "rf"){
    model <- RandForestWrap(RefData = RefData,
                            ClassLabels = ClassLabels,
                            prefix = prefix,
                            mod.meth = mod.meth,
                            train.control = train.control,
                            thread = thread,
                            ...)
    model <- list(model)
  } else if(mod.meth == "svmLinear"){
    model <- SvmWrap(RefData = RefData,
                     ClassLabels = ClassLabels,
                     prefix = prefix,
                     mod.meth = mod.meth,
                     train.control = train.control,
                     thread = thread,
                     ...)
    model <- list(model)
  } else if(mod.meth == "hrf"){
    try(if(missing(tree)|| missing(ClassLabels) || missing(RefData))
      stop("Please, provide the required inputs!"))
    model <- HieRandForest(RefData = RefData,
                           ClassLabels = ClassLabels,
                           tree, thread = thread) # For now.
  }
  #save(model, file=paste(prefix,"model.Robj",sep = "."))
  return(model)
}

#' An internal function to create hierarchical classification models (hiemods) with random forest.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param tree A tree storing relatinship between the class labels. Default is null. required when mod.meth "hrf" is chosen.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
#' @keywords
#' @export
#' @usage
HieRandForest <- function(RefData, ClassLabels, tree, thread=3){
  library(doParallel)
  node.list <- DigestTree(tree = tree)
  hiemods <- vector("list", length = max(node.list))

  Rdata <- DataReshaper(Data = RefData, Predictors = make.names(rownames(RefData)), ClassLabels = ClassLabels)
  #thread=NULL #TEMPORARY for Skipping foreach.
  if(is.null(thread)){
    #Build a local classifier for each node in the tree. Binary or multi-class mixed.
    for(i in node.list){
      hiemods[[i]] <- NodeTrainer(Tdata = Rdata, tree = tree, node = i)
    } #closes the for loop.

  }else{# thread is specified. For now, use this only when running on bigMem machines.
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/"))

    print(paste("registered cores is", getDoParWorkers(), sep = " "))

    out <- foreach(i=node.list, .packages = c('caret'), .inorder = TRUE, .export = ls(.GlobalEnv)) %dopar% {
      NodeTrainer(Tdata = Rdata, tree = tree, node = i)
    }
    stopCluster(cl)

    for(x in 1:length(out)){
      hiemods[node.list[x]] <- out[x]
    }
  }

  names(hiemods) <- seq_along(hiemods)
  hiemods[sapply(hiemods, is.null)] <- NULL

  return(hiemods)
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

#' A function used internally for selecting genes based on their weight in the top principle components.
#' @description p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom.
#' Weighted gene picking depending on PC number: Initial PCs give more genes.
#' For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
#' @param Data an data matrix storing gene expression as genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param PCs the number of PCs to be looked at when selecting genes. Default is 40.
#' @param num the number of Predictors (genes) in total to be included. Default is 2000.
#' @param doPlots boolean to generate a pdf for PCA plots. Default is False.
#' @param prefix a prefix to separate each run.
#' @keywords PCA loadings selection
#' @export
#' @usage Predictors <- FeatureSelector(Data = as.matrix(SeuratObject@data), PCs = 10, num = 2000)
FeatureSelector <- function(Data, ClassLabels, PCs=40, num=2000, doPlots=F, prefix="Feature.select") {
  library(dplyr)
  #Transpose the matrix cols <--> rows t()
  TData <- DataReshaper(Data = Data)
  #do pca
  pcatrain <- prcomp(TData, center = TRUE, scale=TRUE, rank. = PCs)
  pcadata <- data.frame(pcatrain$x, ClassLabels = ClassLabels)
  class_n=length(levels(pcadata$ClassLabels))

  print("Testing the PCs for classification...")
  #Select only PCs which are significanly separating at least one of the class.
  is_significant <- function(x) any(as.numeric(x) <= 0.05)

  #Create a table for evaluation of all PCs with a significance test.
  PCs.sig.table <- pcadata %>% group_by(ClassLabels) %>%
    dplyr::select(paste("PC", 1, sep = ""):paste("PC", PCs, sep = "")) %>%
    summarize_all(list(~t.test(x = .)$p.value)) %>%
    as.data.frame()
  rownames(PCs.sig.table) <- PCs.sig.table$ClassLabels
  PCs.sig.table <- PCs.sig.table[, -c(1)]

  PCs.sig <- pcadata %>% group_by(ClassLabels) %>%
    dplyr::select(paste("PC", 1, sep = ""):paste("PC", PCs, sep = "")) %>%
    summarize_all(list(~t.test(x = .)$p.value)) %>%
    select_if(., is_significant) %>% colnames()

  print(PCs.sig)
  if(length(PCs.sig) == 0){PCs.sig <- c("PC1", "PC2", "PC3")}
  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  plots <- list()
  pb <- txtProgressBar(min = 0, max = length(PCs.sig), style = 3)
  for (i in 1:length(PCs.sig)){
    orderedpcai <- pcatrain$rotation[order(abs(pcatrain$rotation[, PCs.sig[i]]), decreasing = TRUE), PCs.sig[i]]
    Gn <- round((length(PCs.sig)-i+1)*(num*2)/(length(PCs.sig)*(length(PCs.sig)+1)))
    print(Gn)
    TotalGenes <- as.numeric(TotalGenes) + Gn
    top <- data.frame(genes=names(head(orderedpcai, Gn)), bestgenes=head(orderedpcai, Gn))
    load <- rbind(load, top)
    print(top)
    setTxtProgressBar(pb, i+1)
    pc.p <- ggplot(pcadata, aes_string(x=PCs.sig[i], y="ClassLabels", color="ClassLabels"))+geom_point()
    plots[[i]] <- pc.p
    cat("\n", 'Picking the best genes from', PCs.sig[i], ' is completed.', "\n")
  }

  if (doPlots == TRUE) {
    pdf(paste(prefix,"_PCA_eval_p.pdf",sep=""),width = PCs+5, height = class_n/1.2)
    pheatmap::pheatmap(-1*log10(PCs.sig.table+.000001), cluster_cols = F, treeheight_row = 0, cellheight = 10, cellwidth = 10)
    dev.off()
    pdf(paste(prefix,"_PCA_class.pdf",sep=""), width = 10, height = class_n*PCs/2)
    multiplot(plotlist = plots)
    dev.off()
  }

  print(paste("memory used:", pryr::mem_used()))
  bestgenes <- make.names(as.character(unique(load$genes)))
  print(bestgenes)
  return(bestgenes)
}

#' A function to create a tree in phylo format.
#' @param treeTable a data table/data.frame storing class relationships with intermediate class labels. Can be read from a tab separated file.
#' @usage treeTable <- read.delim("TaxonomyRank.tree.txt", header = F)
#' @usage tree <- CreateTree(treeTable)
CreateTree <- function(treeTable){
  #Create a tree:
  library(ape)
  library(data.tree)# '::' doesn't work!

  treeTable <- data.frame(lapply(treeTable, function(x) {gsub("\\+|-", ".", x)}))

  treeTable$pathString <- apply(cbind("TaxaRoot", treeTable), 1, paste0, collapse="/")
  tree <- as.phylo(as.Node(treeTable))
  return(tree)
}

NodeTrainer <- function(Tdata, tree, node){
  print(paste("Node id: ", node))
  node.Data <- SubsetTData(Tdata = Tdata, tree = tree, node = node)
  node.ClassLabels <- node.Data["ClassLabels", ]
  node.Data <- node.Data[which(!rownames(node.Data) %in% c("ClassLabels")), ]

  node.mod <- Modeller(RefData = node.Data,
                       ClassLabels = node.ClassLabels,
                       mod.meth="rf",
                       cv.k=2,
                       save.int.f=FALSE)
  return(node.mod)
}

