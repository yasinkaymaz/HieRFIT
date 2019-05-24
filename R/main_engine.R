
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
                               Evaluation = "data.frame"))

#' Reference class
#' @slot model A list of models to be used.
#' @slot tree A hierarchical tree representing class relationships.
#' @slot modtype type of the model method used to create object.
RefMod <- setClass(Class = "RefMod",
                   slots = c(model = "list",
                             tree = "list",
                             modtype = "character"))

#'
#' @param Ref Reference data from which class labels will be projected on Query data.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix (Ref). Same length as colnames(Ref).
#' @param method The model training method, "rf" for random forest, "svmLinear" for support vector machine, "hrf" for hierarchical random forest. Default is "hrf"
#' @param TreeFile An input file to create a hierarchical tree for relationship between class labels. Default is null but required if method 'hrf' is chosen.
#' @examples pbmc.refmod <- CreateRef(Ref = as.matrix(pbmc@data), ClassLabels = pbmc@meta.data$ClusterNames_0.6, TreeFile = pbmc3k_tree)
#' @examples
CreateRef <- function(Ref, ClassLabels, TreeTable=NULL, method="hrf"){

  if(method == "hrf"){
    tree <- CreateTree(treeTable = TreeTable)
  }

  #Replace white space with '_'
  ClassLabels <- gsub(ClassLabels, pattern = " ", replacement = "_")
  print(ClassLabels)
  #Create predictive model structure. No need to create if exists.
  print("Training model... This may take some time... Please, be patient!")
  model <- Modeller(ExpData = Ref,
                    ClassLabels = ClassLabels,
                    mod.meth = method,
                    tree = tree)
  refObj <- new(Class = "RefMod",
                model = model,
                modtype = method)
  refObj@tree[[1]] <- tree
  #Return a Ref object rather than single model.
  return(refObj)
}

#' The main function to project reference data on query in order to identify class labels.
#' @param Query Query data whose components will be labeled with reference data.
#' @param refMod optional input if model exist already. Default is null and generated from scratch. Input can be an caret model object or a .Rdata file.
#' @param xSpecies optional argument to specify cross species information transfer. Default is null. Possible options are 'rat2mouse', 'mouse2rat', 'mouse2human', human2mouse. With respect to model data.
#' @keywords
#' @export
#' @examples expRefObj <- get(load("data/exp_refObj.Rdata"))
#' @examples cpo <- HieRFIT(Query = as.matrix(pbmc1@data), refMod = expRefobj)
HieRFIT <- function(Query, refMod, xSpecies=NULL){

  if( !is.null(xSpecies)) {
    if(xSpecies == "mouse2rat"){ ## pay attention! swapped logic.
    print("Rat to mouse gene id conversion...")
    ort <- Gmor(RatGenes = rownames(Query))
    Query <- Query[which(rownames(Query) %in% ort$external_gene_name), ]
    rownames(Query) <- ort[match(rownames(Query), ort$external_gene_name),]$mmusculus_homolog_associated_gene_name
    }
  }

  if(refMod@modtype == "hrf"){
    nodes_P_all <- CTTraverser(Query = Query, tree = refMod@tree[[1]], hiemods = refMod@model)
    P_path_prod <- ClassProbCalculator(tree = refMod@tree[[1]], nodes_P_all = nodes_P_all)
    #Run uncertainty function
    #exclude first column with query ids.
    Prediction <- colnames(P_path_prod)[apply(P_path_prod, 1, which.max)]
    #ScoreEvals
    ScoreEvals <- ScoreEval(P_path_prod = P_path_prod)
  }else{
    P_path_prod <- Predictor(model = refMod@model[[1]], Query = Query)
    #Run uncertainty function
    Prediction <- colnames(P_path_prod)[apply(P_path_prod, 1, which.max)]
    #ScoreEvals
    ScoreEvals <- ScoreEval(P_path_prod = P_path_prod)
  }
  object <- new(Class = "CellRpred",
                ClassProbilities = P_path_prod,
                Projection = Prediction,
                Evaluation = ScoreEvals
                )

  return(object)
}

ScoreEval <- function(P_path_prod){
  #library(entropy)

  class_n <- length(colnames(P_path_prod))
  P_path_prod$Diff <- apply(P_path_prod, 1, function(x) max(x) - sort(x, partial=length(x) - 1)[length(x) - 1])
  P_path_prod$KLe <- apply(P_path_prod[, which(!names(P_path_prod) %in% c("Diff"))], 1, function(x) entropy::KL.empirical(y1 = as.numeric(x), y2 = rep(1/class_n, class_n)))
  P_path_prod$BestScore <- apply(P_path_prod[, which(!names(P_path_prod) %in% c("Diff","KLe"))], 1, function(x) max(x))
  PredEval <- P_path_prod
 return(PredEval)
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
                            model <- list(model)
  } else if(mod.meth == "svmLinear"){
    model <- SvmWrap(ExpData = ExpData,
                     ClassLabels = ClassLabels,
                     prefix = prefix,
                     mod.meth = mod.meth,
                     train.control = train.control,
                     ...)
                     model <- list(model)
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
                             doPlots = F, ...)
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
                             doPlots = F, ...)
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

#' A function used internally for selecting genes based on their weight in the top principle components.
#' @description p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom.
#' Weighted gene picking depending on PC number: Initial PCs give more genes.
#' For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
#' @param ExpData an data matrix storing gene expression as genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param PCs the number of PCs to be looked at when selecting genes. Default is 40.
#' @param num the number of Predictors (genes) in total to be included. Default is 2000.
#' @param doPlots boolean to generate a pdf for PCA plots. Default is False.
#' @param prefix a prefix to separate each run.
#' @keywords PCA loadings selection
#' @export
#' @examples Predictors <- FeatureSelector(ExpData = as.matrix(SeuratObject@data), PCs = 10, num = 2000)
FeatureSelector <- function(ExpData, ClassLabels, PCs=40, num=2000, doPlots=F, prefix="Feature.select") {
  library(dplyr)
  #Transpose the matrix cols <--> rows t()
  TData <- DataReshaper(ExpData = ExpData)
  #do pca
  pcatrain <- prcomp(TData, center = TRUE, scale=TRUE, rank. = PCs)
  pcadata <- data.frame(pcatrain$x, ClassLabels = ClassLabels)
  class_n=length(levels(pcadata$ClassLabels))

  print("Testing the PCs for classification...")
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
