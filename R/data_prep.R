
ParApp <- function() {
  '"
  Sets the parameters used throughout the run.
  "'
  PCs <- 20 # number of principle components to be analyzed.
  Pred_n <- 1000 # number of predictors aimed to be collected.
  prefix <- ""
}

#' The CellRpred Class
#' @slot Prior
#' @slot Predictors
#' @slot Projection
#' @slot mod.name
CellRpred <- setClass(Class = "CellRpred",
                      slots = c(Prior = "character",
                               ClassProbilities = "data.frame",
                               Projection = "character",
                               mod.name = "character"))


#' The main function to project reference data on query in order to identify class labels.
#' @param Ref Reference data from which class labels will be projected on Query data.
#' @param Query Query data whose components will be labeled with reference data.
#' @param method The model training method, "rf" for random forest, "svmLinear" for support vector machine, "hrf" for hierarchical random forest. Default is "hrf"
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix (Ref). Same length as colnames(Ref).
#' @param TreeFile An input file to create a hierarchical tree for relationship between class labels. Default is null but required if method 'hrf' is chosen.
#' @param model optional input if model exist already. Default is null and generated from scratch. Input can be an caret model object or a .Rdata file.
#' @keywords
#' @export
#' @examples cpo <- CellProjector(Ref = as.matrix(pbmc@data), ClassLabels = pbmc@meta.data$ClusterNames_0.6, Query = as.matrix(pbmc1@data))
CellProjector <- function(Ref, ClassLabels, Query, TreeFile=NULL, model=NULL, method="hrf"){
  #Tree file. Required if method 'hrf' is chosen.
  tree <- CreateTree(TreeFile)
  #Determine the input (Ref/Query) properties.
  #Create predictive model structure. No need to create if exists.
  #Skip this step if already ran.
  if (!is.null(model)) {
    if(file.exists(model)){
      print("An existing file for a model is found in the directory! Using it...")
      model <- get(load(model))
    }else{
      model <- model
    }
  }else{
    print("Training model... This may take some time... Please, be patient!")
    model <- Modeller(ExpData = Ref, ClassLabels = ClassLabels, mod.meth = method, tree = tree)
  }

  if(method == "hrf"){
    nodes_P_all <- CTTraverser(Query = Query, tree = tree, hiemods = model)
    P_path_prod <- ClassProbCalculator(tree = tree, nodes_P_all = nodes_P_all)
    #Run uncertainty function
    #exclude first column with query ids.
    Prediction <- colnames(P_path_prod)[apply(P_path_prod, 1, which.max)]

  }else{
    P_path_prod <- Predictor(model = model, Query = Query)
    #Run uncertainty function
    Prediction <- colnames(P_path_prod)[apply(P_path_prod, 1, which.max)]
  }

  object <- new(Class = 'CellRpred',
                ClassProbilities = P_path_prod,
                Projection = Prediction,
                mod.name = method
                )

  return(object)
}


#' An internal function to collect model training parameters and direct them to model creation.
#' @param ExpData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param mod.meth The model training method, "rf" for random forest, "svmLinear" for support vector machine, "hrf" for hierarchical random forest. Default is "rf"
#' @param cv.k Fold cross validation. Default is 5.
#' @param tree A tree storing relatinship between the class labels. Default is null. required when mod.meth "hrf" is chosen.
#' @param save.int.f Boolean to save model. Default is False.
#' @keywords
#' @export
#' @examples rf.model <- Modeller(ExpData = as.matrix(pbmc1@data), ClassLabels = pbmc1@meta.data$ClusterNames_0.6)
Modeller <- function(ExpData, ClassLabels=NULL, mod.meth="rf", cv.k=5, thread=NULL, tree=NULL, save.int.f=FALSE, ...){
  prefix=paste(cv.k, "K-fold", mod.meth, sep=".")
  #k-fold Cross Validation
  train.control <- caret::trainControl(method="cv",
                                       number=cv.k,
                                       savePredictions = TRUE,
                                       classProbs =  TRUE)
  #Choose whether flat or hierarchical classifier.
  if(mod.meth == "rf"){
    model <- RandForestWrap(ExpData = ExpData,
                            ClassLabels = ClassLabels,
                            prefix = prefix,
                            mod.meth = mod.meth,
                            train.control = train.control,
                            ...)
  } else if(mod.meth == "svmLinear"){
    model <- SvmWrap(ExpData = ExpData,
                     ClassLabels = ClassLabels,
                     prefix = prefix,
                     mod.meth = mod.meth,
                     train.control = train.control,
                     ...)
  } else if(mod.meth == "hrf"){
    try(if(missing(tree)|| missing(ClassLabels) || missing(ExpData))
      stop("Please, provide the required inputs!"))
    model <- HieRandForest(ExpData = ExpData,
                           ClassLabels = ClassLabels,
                           tree, thread = thread) # For now.
  }
  save(model, file=paste(prefix,"model.Robj",sep = "."))
  return(model)
}

#' A wrapper function for random forest.
#' @param ExpData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param mod.meth The model training method, "rf" for random forest.
RandForestWrap <- function(ExpData=ExpData, ClassLabels=ClassLabels, prefix, mod.meth, train.control, ncor=5, ...){
  library(doParallel)
  #1. Select the predictors.
  P_dicts <- FeatureSelector(ExpData = ExpData,
                             ClassLabels = ClassLabels,
                             PCs = 10,
                             num = 200,
                             prefix = prefix,
                             spitPlots = F, ...)
  #2. Prepare the reference data.
  RefData <- DataReshaper(ExpData = ExpData,
                          Predictors = P_dicts,
                          ClassLabels = ClassLabels, ...)

  #3. Train the model.
  cl <- makePSOCKcluster(ncor)
  registerDoParallel(cl)
  model <- caret::train(ClassLabels~., data = RefData,
                        trControl = train.control,
                        method = mod.meth,
                        norm.votes = TRUE,
                        importance = TRUE,
                        proximity = TRUE,
                        #       preProcess = c("center", "scale"),
                        ntree=50, ...)
  stopCluster(cl)

  return(model)
}

#' A wrapper function for random forest.
#' @param ExpData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param mod.meth The model training method, "svmLinear" for support vector machine.
SvmWrap <- function(ExpData=ExpData, ClassLabels=ClassLabels, prefix, mod.meth, train.control, ncor=5, ...){
  library(doParallel)
  #1. Select the predictors.
  P_dicts <- FeatureSelector(ExpData = ExpData,
                             ClassLabels = ClassLabels,
                             PCs = 10,
                             num = 200,
                             prefix = prefix,
                             spitPlots = F, ...)
  #2. Prepare the reference data.
  RefData <- DataReshaper(ExpData = ExpData,
                          Predictors = P_dicts,
                          ClassLabels = ClassLabels, ...)

  #3. Train the model.
  cl <- makePSOCKcluster(ncor)
  registerDoParallel(cl)
  model <- caret::train(ClassLabels~., data = RefData,
                        trControl = train.control,
                        method = mod.meth,
                        norm.votes = TRUE,
                        importance = TRUE,
                        proximity = TRUE,
                        #       preProcess = c("center", "scale"),
                        tuneLength = 10, ...)
  stopCluster(cl)

  return(model)
}

#' An internal function to create hierarchical classification models (hiemods) with random forest.
#' @param ExpData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param tree A tree storing relatinship between the class labels. Default is null. required when mod.meth "hrf" is chosen.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
#' @keywords
#' @export
#' @examples
HieRandForest <- function(ExpData=ExpData, ClassLabels=ClassLabels, tree, thread=3){
  library(doParallel)
  node.list <- DigestTree(tree = tree)
  hiemods <- vector("list", length = max(node.list))

  Tdata <- DataReshaper(ExpData = ExpData, Predictors = make.names(rownames(ExpData)), ClassLabels = ClassLabels)

  if(is.null(thread)){
    #Build a local classifier for each node in the tree. Binary or multi-class mixed.
    for(i in node.list){
      hiemods[[i]] <- NodeTrainer(Tdata = Tdata, tree = tree, node = i)
    } #closes the for loop.

  }else{# thread is specified. For now, use this only when running on bigMem machines.
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    print(paste("registered cores is", getDoParWorkers(), sep = " "))

    out <- foreach(i=node.list, .packages = c('ggplot2'), .inorder = TRUE) %dopar% {
      NodeTrainer(Tdata = Tdata, tree = tree, node = i)
    }
    stopCluster(cl)

    for(x in 1:length(out)){
      hiemods[node.list[x]] <- out[x]
    }
  }

  names(hiemods) <- seq_along(hiemods)
  hiemods[sapply(hiemods, is.null)] <- NULL
  class(hiemods) <- "hie.rf"

  return(hiemods)
}

#' A function used internally for selecting genes based on their weight in the top principle components.
#' @description p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom.
#' Weighted gene picking depending on PC number: Initial PCs give more genes.
#' For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
#' @param ExpData an data matrix storing gene expression as genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param PCs the number of PCs to be looked at when selecting genes. Default is 40.
#' @param num the number of Predictors (genes) in total to be included. Default is 2000.
#' @param spitPlots boolean to generate a pdf for PCA plots. Default is False.
#' @param prefix a prefix to separate each run.
#' @keywords PCA loadings selection
#' @export
#' @examples Predictors <- FeatureSelector(ExpData = as.matrix(SeuratObject@data), PCs = 10, num = 2000)
FeatureSelector <- function(ExpData, ClassLabels, PCs=40, num=2000, spitPlots=F, prefix="Feature.select") {
  library(dplyr)
  #Transpose the matrix cols <--> rows t()
  TData <- DataReshaper(ExpData = ExpData)

  if (file.exists(paste(prefix,"train.prcomp.Rdata",sep="."))) {
    warning("An existing PCA data is found in the directory! Used it instead...")
#    pcatrain <- get(load(paste(prefix,"train.prcomp.Rdata",sep=".")))
    pcatrain <- prcomp(TData, center = TRUE, scale=TRUE, rank. = PCs)
  }else{
    print("Performing PCA...")
    pcatrain <- prcomp(TData, center = TRUE, scale=TRUE, rank. = PCs)
#    save(pcatrain,file=paste(prefix,"train.prcomp.Rdata",sep="."))
  }

  pcadata <- data.frame(pcatrain$x, ClassLabels = ClassLabels)
  class_n=length(levels(pcadata$ClassLabels))

  print("Testing the PCs for classification...")
  save(pcadata,file=paste(prefix,"train.pcadata.Rdata",sep="."))
  #Select only PCs which are significanly separating at least one of the class.
  is_significant <- function(x) any(as.numeric(x) <= 0.05)

  #Create a table for evaluation of all PCs with a significance test.
  PCs.sig.table <- pcadata %>% group_by(ClassLabels) %>%
    select(paste("PC", 1, sep = ""):paste("PC", PCs, sep = "")) %>%
    summarize_all(list(~t.test(x = .)$p.value)) %>%
    as.data.frame()
  rownames(PCs.sig.table) <- PCs.sig.table$ClassLabels
  PCs.sig.table <- PCs.sig.table[, -c(1)]

  PCs.sig <- pcadata %>% group_by(ClassLabels) %>%
    select(paste("PC", 1, sep = ""):paste("PC", PCs, sep = "")) %>%
    summarize_all(list(~t.test(x = .)$p.value)) %>%
    select_if(., is_significant) %>% colnames()

  print(PCs.sig)

  load <- NULL
  topbotnames <- NULL
  TotalGenes <- 0
  plots <- list()
  pb <- txtProgressBar(min = 0, max = length(PCs.sig), style = 3)
  for (i in 1:length(PCs.sig)){
    orderedpcai <- pcatrain$rotation[order(abs(pcatrain$rotation[, PCs.sig[i]]), decreasing = TRUE), PCs.sig[i]]
    Gn <- round((PCs-i+1)*(num*2)/(PCs*(PCs+1)))
    print(Gn)
    TotalGenes <- as.numeric(TotalGenes) + Gn
    top <- data.frame(genes=names(head(orderedpcai, Gn)), bestgenes=head(orderedpcai, Gn))
    load <- rbind(load, top)
    setTxtProgressBar(pb, i+1)
    pc.p <- ggplot(pcadata, aes_string(x=PCs.sig[i], y="ClassLabels", color="ClassLabels"))+geom_point()
    plots[[i]] <- pc.p
    cat("\n", 'Picking the best genes from', PCs.sig[i], ' is completed.', "\n")
  }

  if (spitPlots == TRUE) {
    pdf(paste(prefix,"_PCA_eval_p.pdf",sep=""),width = PCs+5, height = class_n/1.2)
    pheatmap::pheatmap(-1*log10(PCs.sig.table+.000001), cluster_cols = F, treeheight_row = 0, cellheight = 10, cellwidth = 10)
    dev.off()
    pdf(paste(prefix,"_PCA_class.pdf",sep=""), width = 10, height = class_n*PCs/2)
    multiplot(plotlist = plots)
    dev.off()
  }

  print(paste("memory used:", pryr::mem_used()))
  bestgenes <- make.names(as.character(unique(load$genes)))
  return(bestgenes)
}

#' An internal function to prepare a training/test dataset for model generation.
#' @param ExpData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param Predictors the predictor feature list selected by FeatureSelector.
#' @param ClassLabels [optional] A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param alpa variation cutoff for filtering data
#' @keywords data preparation
#' @export
#' @examples trainingData <- DataReshaper(ExpData = as.matrix(SeuratObject@data), Predictors = genes, ClassLabels = SeuratObject@meta.data$CellTypes)
DataReshaper <- function(ExpData, Predictors, ClassLabels, alpa=0.1, ...) {
  if(missing(Predictors)){#prepare the data for PCA
    TData <- as.data.frame(t(ExpData))
    indx <- sapply(TData, is.factor)
    TData[indx] <- lapply(TData[indx], function(x) as.numeric(as.character(x)))
    #Filter candidate predictors before PCA
    TData <- TData[, apply(TData, 2, var) != 0]
    TData <- droplevels(TData[, which(matrixStats::colVars(as.matrix(TData)) > alpa)])
    return(TData)

  } else {

    if (missing(ClassLabels)) {
      #Then the output is for Query
      QueData <- as.data.frame(t(ExpData), col.names=rownames(ExpData))
      colnames(QueData) <- make.names(colnames(QueData))
      QueData_sub <- droplevels(QueData[, which(colnames(QueData) %in% Predictors)])
      #missing Predictors
      mp <- Predictors[which(!Predictors %in% colnames(QueData))]
      #Add missing predictors into QueData by setting to 0.
      mp_df <- data.frame(matrix(0,
                               ncol = length(mp),
                               nrow = length(colnames(ExpData))))
      colnames(mp_df) <- mp
      QueData <- cbind(QueData_sub, mp_df)
      QueData <- QueData[, Predictors]
      return(QueData)

  } else {

    RefData <- as.data.frame(t(ExpData), col.names=rownames(ExpData))
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

#' A hierarchical Classifier Tree Traverser function.
#' @param Query is the input query data. rows are genes and columns are cells.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
#' @param hiemods models list from HieRandForest function.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
CTTraverser <- function(Query, tree, hiemods, thread=NULL){
  library(doParallel)
  node.list <- DigestTree(tree = tree)
  #Create a table for storing node probabilities.
  ProbTab <- data.frame(Queries = colnames(Query))

  if(is.null(thread)){
    #Build a local classifier for each node in the tree. Binary or multi-class mixed.
    for(i in node.list){
      nodeProb <- Predictor(model = hiemods[[as.character(i)]],
                            format = "prob",
                            Query = Query,
                            node = i)
      ProbTab <- cbind(ProbTab, nodeProb)
    } #closes the for loop.

  }else{# thread is specified. For now, use this only when running on bigMem machines.
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    print(paste("registered cores is", getDoParWorkers(), sep = " "))

    nodeProb <- foreach(i=node.list, .inorder = TRUE, .combine=cbind) %dopar% {
      Predictor(model = hiemods[[as.character(i)]],
                format = "prob",
                Query = Query,
                node = i)
    }
    stopCluster(cl)
    ProbTab <- cbind(ProbTab, nodeProb)
  }

  return(ProbTab)
}

#'
#' @param model
#' @param Query
#' @param format type of prediction output, "prob" or "resp".
Predictor <- function(model, Query, format="prob", node=NULL){
  P_dicts <- colnames(model$trainingData)
  P_dicts <- P_dicts[P_dicts != ".outcome"]
  QueData <- DataReshaper(ExpData = Query, Predictors = P_dicts)

  if(is.null(node)){
    QuePred <- as.data.frame(predict(model, QueData, type = "prob"))
  } else{
    if(format == "prob"){
      c_f <- length(model$levels) #correction factor; class size
      QuePred <- as.data.frame(predict(model, QueData, type = "prob"))
      QuePred <- QuePred*c_f
      colnames(QuePred) <- paste(node, colnames(QuePred), sep = "")
    } else {
      QuePred <- as.data.frame(predict(model, QueData, type = "raw"))
      colnames(QuePred) <- as.character(node)
    }
  }
  return(QuePred)
}


#' Function to generate a random tree using 'ape' package.
#' Returns a tree object.
#' @param LN The target number of leafs in the tree. Default is 8.
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

CreateTree <- function(x){
  #Create a tree:
  final.text <- "((CD8.T.cells, CD4.T.cells), (CD14..Monocytes, FCGR3A..Monocytes), Dendritic.cells, NK.cells, B.cells, Megakaryocytes);"
  tree <- ape::read.tree(text=final.text)
  L <- list(T_cells = "T.cells", Monocytes = "Monocytes", Pbmc = c("T.cells", "Monocytes"))
  tree <- ape::makeNodeLabel(tree, "u", nodeList = L)

  return(tree)
} # Fix dummy input

CheckTreeIntegrity <- function(){}

DigestTree <- function(tree) {
  all.nodes <- unique(tree$edge[,1])
  return(all.nodes)
}

#' A function to retrieve corresponding leafs of the child nodes of a given node.
#' @param tree A tree storing relatinship between the class labels.
#' @param node a particular non-terminal node in the tree.
GetChildNodeLeafs <- function(tree, node){
  #Get the children nodes.
  children <- tree$edge[which(tree$edge[, 1] == node), 2]
  #store ordered tips
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  #Then extract all leafs below these children
  c.tips <- list()
  for (c in children){
    if(c > length(tree$tip.label)){
      print(paste("extracting tips for node", c, sep=" "))
      c.tips[[c]] <- ape::extract.clade(tree, c)$tip.label
    }else{
      c.tips[[c]] <- tree$tip.label[ordered_tips][match(c, ordered_tips)]
    }
  }
  return(c.tips)
}

SubsetTData <- function(Tdata, tree, node){
  #This function subsets a Tdata with class labels (trainingData - output of DataReshaper)
  #And updates the class labels.
  # 1. Extract the data under the node. Subsample if necessary.
  childNodes <- GetChildNodeLeafs(tree = tree, node = node)

  SubTdata <- NULL
  #loop through child nodes that are not null in the list.
  for (i in which(lapply(childNodes, length) > 0)){

    Subdata <- droplevels(Tdata[which(Tdata$ClassLabels %in% childNodes[i][[1]]), ])

    if (i > length(tree$tip.label)){# if the node is not a leaf node, then
      if(!is.null(tree$node.label)){# if labels for subnodes exist
        labels <- c(tree$tip.label, tree$node.label)
        #Replace labels with subnode labels.
        Subdata$ClassLabels <- as.factor(labels[i])
        print(paste("For Node", i, "the label is:", labels[i], sep = " "))
      } else {
        #if subnodes don't exist, replace class tip labels with childnode label.
        Subdata$ClassLabels <- as.factor(i)
      }
    } else {#if the node is a terminal leaf
      #Replace class tip labels with Leaf labels.
      Subdata$ClassLabels <- as.factor(childNodes[i][[1]])
    }
    #Combine with all other child node data
    SubTdata <- rbind(SubTdata, Subdata)
  }
  return(t(SubTdata))
}

NodeTrainer <- function(Tdata, tree, node){
  print(paste("Node id: ", node))
  node.ExpData <- SubsetTData(Tdata = Tdata, tree = tree, node = node)
  node.ClassLabels <- node.ExpData["ClassLabels", ]
  node.ExpData <- node.ExpData[which(!rownames(node.ExpData) %in% c("ClassLabels")), ]

  node.mod <- Modeller(ExpData = node.ExpData,
                         ClassLabels = node.ClassLabels,
                         mod.meth="rf",
                         cv.k=2,
                         save.int.f=FALSE)
  return(node.mod)
}

GetAncestPath <- function(tree, class){
  path <- c()
  labs_l <- c(tree$tip.label, tree$node.label)
  Node <- match(class, labs_l)
  parent <- tree$edge[which(x = tree$edge[, 2] == Node), ][1]
  while(!is.na(parent)){
    path <- c(path, paste(parent, class, sep = ""))
    class <- labs_l[parent]
    parent <- tree$edge[which(x = tree$edge[, 2] == parent), ][1]
  }
  return(path)
}

#' A function to calculate class scores for all internal and tip node classes.
#' @param tree
#' @param nodes_P_all a table output from CTTraverser() which contains all class probabilities from every node.
#' @return P_path_prod a table for products of ancestor node scores.
ClassProbCalculator <- function(tree, nodes_P_all){
  #CTip_table <- data.frame(matrix(ncol = 0, nrow = length(rownames(Htable))))
  P_path_prod <- data.frame(row.names = rownames(nodes_P_all))
  clabs <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]

  for(cl in clabs){
    map <- GetAncestPath(tree = tree, class = cl)
    nt_prob <- data.frame(matrixStats::rowProds(as.matrix(nodes_P_all[, map])))
    colnames(nt_prob) <- cl
    P_path_prod <- cbind(P_path_prod, nt_prob)
  }
  return(P_path_prod)
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

#' A function to downsample Seurat object based on cell type identity
#' @param SeuratObj a Seurat S4 object.
#' @param IdentityCol the column 'number' in the metadata slot showing the cell type identities.
#' @example pbmc1 <- DownSizeSeurat(SeuratObj = pbmc, IdentityCol = 7)
DownSizeSeurat <- function(SeuratObj, IdentityCol, min_n=NULL){
  cells <- NULL
  classes <- table(SeuratObj@meta.data[,IdentityCol])
  print(classes)
  if(is.null(min_n)){
    min_n <- min(classes)
    print(min_n)
  }
  for(type in names(classes)){
    cells <- c(cells, sample(rownames(SeuratObj@meta.data[which(SeuratObj@meta.data[,IdentityCol] == type), ]), size = min_n, replace = F))
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

CellTyper <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions"){

  if(!missing(SeuratObject)){
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{
    print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
  }#Closes missing(SeuratObj)
  colnames(testExpSet) <- make.names(colnames(testExpSet))
  #Prepare Test Expression set
  testsub <- testExpSet[,which(colnames(testExpSet) %in% attributes(model$terms)$term.labels)]
  missingGenes <- attributes(model$terms)$term.labels[which(!attributes(model$terms)$term.labels %in% colnames(testExpSet))]
  print(missingGenes)
  missingGenes.df <- data.frame(matrix(0, ncol = length(missingGenes), nrow = length(rownames(testExpSet))))
  colnames(missingGenes.df) <- missingGenes
  TestData <- cbind(testsub, missingGenes.df)
  TestData <- TestData[,attributes(model$terms)$term.labels]
  cat("Number of Features (genes) to be considered is", length(colnames(testsub)), '\n', "Number of missing Features set to zero is", length(missingGenes), '\n', sep = ' ')

  rm(testsub, missingGenes, missingGenes.df)
  gc()
  #Predict
  # library(entropy)
  testPred <- as.data.frame(predict(model, TestData, type = "prob"))
  class_n <- length(model$classes)
  testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
  testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class_n, class_n)) )
  testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
  testPred$Prediction <- predict(model, TestData, type="raw")

  #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class_n)
  testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class_n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class_n), "Undetermined", "Detected")) %>% as.data.frame()
  #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class_n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

  if(missing(priorLabels)){
    print("Prior class labels are not provided!")

  }else{
    #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    testPred <- cbind(testPred, priorLabels)

    #Plot the crosscheck here:
    #Crosscheck Predictions
    #  library(tidyverse)
    #  library(alluvial)
    #  library(ggalluvial)
    crx <- testPred %>% group_by(Prior, Intermediate, Prediction) %>% tally() %>% as.data.frame()

    p5 <- ggplot(crx,aes(y = n, axis1 = Prior, axis2 = Intermediate, axis3 = Prediction )) +
      geom_alluvium(aes(fill = Prediction), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("Prior", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    cowplot::save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 16, base_width = 20)

  }

  if(!missing(SeuratObject)){

    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(testPred))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, testPred)

    PlotPredictions(SeuratObject = SeuratObject, model = model, outputFilename = outputFilename)

    return(SeuratObject)
  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function

