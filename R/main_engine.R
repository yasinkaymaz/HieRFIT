
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
HieRFIT <- setClass(Class = "HieRFIT",
                      slots = c(Prior = "character",
                               ClassProbs = "data.frame",
                               ClassWeights = "data.frame",
                               CertaintyValues = "data.frame",
                               ClassScores = "data.frame",
                               Evaluation = "data.frame"))

#' Reference class
#' @slot model A list of models to be used.
#' @slot tree A hierarchical tree representing class relationships.
#' @slot modtype type of the model method used to create object.
RefMod <- setClass(Class = "RefMod",
                   slots = c(model = "list",
                             tree = "list",
                             modtype = "character"))


HieMetrics <- setClass(Class = "HieMetrics",
                       slots = c(Pvotes = "data.frame",
                                 QueWs = "data.frame",
                                 QueCers = "data.frame",
                                 Scores = "data.frame")
                       )

#' The main function for creating a reference model.
#' @param Ref Reference data from which class labels will be projected on Query data.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix (Ref). Same length as colnames(Ref).
#' @param method The model training method, "rf" for random forest, "svmLinear" for support vector machine, "hrf" for hierarchical random forest. Default is "hrf"
#' @param TreeFile An input file to create a hierarchical tree for relationship between class labels. Default is null but required if method 'hrf' is chosen.
#' @usage  pbmc.refmod <- CreateRef(Ref = as.matrix(pbmc@data), ClassLabels = pbmc@meta.data$ClusterNames_0.6, TreeFile = pbmc3k_tree)
CreateHieR <- function(Ref, ClassLabels, TreeTable=NULL, method="hrf", thread=NULL){

  if(method == "hrf"){
    tree <- CreateTree(treeTable = TreeTable)
  }else{ tree <- NULL}

  ClassLabels <- FixLab(x = ClassLabels)

  #Create predictive model structure. No need to create if exists.
  print("Training model... This may take some time... Please, be patient!")
  model <- Modeller(ExpData = Ref,
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

#' The main function to project reference data on query in order to identify class labels.
#' @param Query Query data whose components will be labeled with reference data.
#' @param refMod optional input if model exist already. Default is null and generated from scratch. Input can be an caret model object or a .Rdata file.
#' @param Prior prior class labels if exist. For cross comparison. Should correspond row order of the Query.
#' @param xSpecies optional argument to specify cross species information transfer. Default is null. Possible options are 'rat2mouse', 'mouse2rat', 'mouse2human', human2mouse. With respect to model data.
#' @keywords
#' @export
#' @usage expRefObj <- get(load("data/exp_refObj.Rdata"))
#' @usage cpo <- HieRFIT(Query = as.matrix(pbmc1@data), refMod = expRefobj)
HieRFIT <- function(Query, refMod, Prior=NULL, xSpecies=NULL, alpha=.9){

  if (class(Query) == "seurat" | class(Query) == "Seurat" ){
    Query_d <- as.matrix(Query@data)}else{Query_d <- Query}

  if( !is.null(xSpecies)) {
    if(xSpecies == "mouse2rat"){ ## pay attention! swapped logic.
    print("Rat to mouse gene id conversion...")
    ort <- Gmor(RatGenes = rownames(Query_d))
    Query_d <- Query_d[which(rownames(Query_d) %in% ort$external_gene_name), ]
    rownames(Query_d) <- ort[match(rownames(Query_d), ort$external_gene_name),]$mmusculus_homolog_associated_gene_name
    }
  }

  #Query backround
  #Query_bg <- RandomizeR(Query_d)

  if(refMod@modtype == "hrf"){
    # To Do: Implement parallel here as well!

    HieMetObj <- CTTraverser(Query = Query_d, tree = refMod@tree[[1]], hiemods = refMod@model)
    #P_path_prod <- ClassProbCalculator(tree = refMod@tree[[1]], nodes_P_all = HieMetObj@Pvotes)
    #Run uncertainty function
    #nodes_P_all_bg <- CTTraverser(Query = Query_bg, tree = refMod@tree[[1]], hiemods = refMod@model)
    #P_path_prod_bg <- ClassProbCalculator(tree = refMod@tree[[1]], nodes_P_all = nodes_P_all_bg)

  }else{#FIX THIS
    P_path_prod <- Predictor(model = refMod@model[[1]], Query = Query_d)
    #Run uncertainty function
    #P_path_prod_bg <- Predictor(model = refMod@model[[1]], Query = Query_bg)
  }

  #Evaluate scores and run uncertainty function, then, project the class labels.
  ScoreEvals <- ScoreEval(ScoreObs = HieMetObj@Scores, ProbCert = HieMetObj@QueCers, alpha=alpha)

  if (class(Query) == "seurat" | class(Query) == "Seurat" ){
    Query@meta.data <- cbind(Query@meta.data[, which(!colnames(Query@meta.data) %in% colnames(ScoreEvals))],
                             ScoreEvals)
    object <- Query
  }else{
    object <- new(Class = "HieRFIT",
                  Prior = as.character(Prior),
                  ClassProbs = HieMetObj@Pvotes,
                  ClassWeights = HieMetObj@QueWs,
                  CertaintyValues = HieMetObj@QueCers,
                  ClassScores = HieMetObj@Scores,
                  Evaluation = ScoreEvals)
  }

  return(object)
}



#' A function for evalating the uncertainty.
#' @param ScoreObs P_path_prod for observed scores
#' @param ProbCert Certainty scores.
ScoreEval <- function(ScoreObs, ProbCert, alpha=.9){

  df <- data.frame(row.names = rownames(ScoreObs))
  for(i in 1:length(ProbCert[,1])){
    candits <- colnames(ProbCert)[ProbCert[i,] > alpha]
    if(length(candits) == 0){
      classL <- "Undetermined"
      classU <- max(ProbCert[i,])
      classS <- max(ScoreObs[i,])
    }else if(length(candits) == 1){
      classL <- candits
      classU <- ProbCert[i,classL]
      classS <- ScoreObs[i,classL]
    }else{
      classL <- colnames(ScoreObs[i,candits])[apply(ScoreObs[i,candits], 1, which.max)]
      classU <- ProbCert[i,classL]
      classS <- ScoreObs[i,classL]
    }

    df <- rbind(df, data.frame(row.names = rownames(ProbCert[i,]),
                               Score = classS,
                               Certainty = classU,
                               Projection = classL))

  }


 #  library(dplyr)
 #  colMax <- function(x) sapply(x, max, na.rm = TRUE)
 #  #Substract the Background_max.i from Scores
 #  bg_max <- colMax(ScoreBg)
 #  ScoreObs <- ScoreObs - t(as.data.frame(bg_max))[rep(1, each=nrow(ScoreObs)), ]
 #
 #  ScoreEval <- ScoreObs
 #  ScoreEval$BestScore <- apply(ScoreObs, 1, function(x) max(x))
 #  ScoreEval$Projection <- colnames(ScoreObs)[apply(ScoreObs, 1, which.max)]
 #  ScoreEval <- ScoreEval %>%
 #    mutate(Projection = if_else( (BestScore <= 0),
 #                                 "Undetermined",
 #                                 as.character(Projection) )) %>% as.data.frame()
  return(df)
}


#' A function to tally votes from the model locally trained.
#' @param model
#' @param QueData is the prepared data matrix ready to be used in predict
#' @param format type of prediction output, "prob" or "resp".
PvoteR <- function(model, QueData, format="prob", node=NULL){
  if(is.null(node)){
    QuePvotes <- as.data.frame(predict(model, QueData, type = "prob", scale=T, center=T))
  } else{
    if(format == "prob"){
      QuePvotes <- as.data.frame(predict(model, QueData, type = "prob", scale=T, center=T))
      colnames(QuePvotes) <- paste(node, colnames(QuePvotes), sep = "")
    } else {
      QuePvotes <- as.data.frame(predict(model, QueData, type = "raw", scale=T, center=T))
      colnames(QuePvotes) <- as.character(node)
    }
  }
  return(QuePvotes)
}

#' A function to calculate gravitity center (weight center) of the probability distributions
#' among classes of the node. By permutation.
#' @param model of the node.
#' @param QueData same matrix used in PvoteR().
#' @return QueWeights a set of probability weightes per class to be used in asymetrix entropy estimations.
graWeighteR <- function(model, QueData){
  #Randomizing only feature space
  QueData_R <- NULL
  for(i in 1:1){
    QueData_R <- rbind(QueData_R, RandomizeR(df = QueData))
  }
  Ws <- colMeans(PvoteR(model = model, QueData = QueData_R))
  QueWeights <- t(as.data.frame(Ws))[rep(1, each=nrow(QueData)), ]
  #rownames(QueWeights) <- rownames(QueData)
  QueWeights <- as.data.frame(QueWeights)
  return(QueWeights)
}

#' Certainty Estimation function. Calculates certainty values for each class probability.
#' @param qP Observed probability for each class
#' @param qW Probability weights estimated by graWeighteR()
#' @return QueCers
ceR <- function(qP, qW){
  QueCers <- NULL
  for(i in 1:length(qP[,1])){
    QueCers <- rbind(QueCers, GetCertaintyArray(p = qP[i,], w = qW[i,]))
  }
  #rownames(QueCers) <- rownames(qP)
  colnames(QueCers) <- colnames(qP)
  QueCers <- as.data.frame(QueCers)
  return(QueCers)
}


#' The predictor function.
#' @param model
#' @param QueData
#' @param format type of prediction output, "prob" or "resp".
scoR <- function(model, QueData, format="prob", node=NULL){
  # P_dicts <- colnames(model$trainingData)
  # P_dicts <- P_dicts[P_dicts != ".outcome"]
  # QueData <- DataReshaper(ExpData = Query, Predictors = P_dicts)

  if(is.null(node)){
    QuePred <- as.data.frame(predict(model, QueData, type = "prob", scale=T, center=T))
  } else{
    if(format == "prob"){
      c_f <- length(model$levels) #correction factor; class size
      QuePred <- as.data.frame(predict(model, QueData, type = "prob", scale=T, center=T))
      #      QuePred <- as.data.frame(t(apply(QuePred, 1, function(x) KLeCalc(x)*x)))
      QuePred <- QuePred*c_f
      colnames(QuePred) <- paste(node, colnames(QuePred), sep = "")
    } else {
      QuePred <- as.data.frame(predict(model, QueData, type = "raw", scale=T, center=T))
      colnames(QuePred) <- as.character(node)
    }
  }
  return(QuePred)
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
#' @usage rf.model <- Modeller(ExpData = as.matrix(pbmc1@data), ClassLabels = pbmc1@meta.data$ClusterNames_0.6)
Modeller <- function(ExpData, ClassLabels=NULL, mod.meth="rf", cv.k=2, thread=NULL, tree=NULL, save.int.f=FALSE, ...){
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
  #save(model, file=paste(prefix,"model.Robj",sep = "."))
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
  print(RefData[1:5, 1:5])
  #3. Train the model.
  cl <- makePSOCKcluster(ncor)
  registerDoParallel(cl)
  model <- caret::train(ClassLabels~., data = RefData,
                        method = mod.meth,
                        norm.votes = TRUE,
                        importance = TRUE,
                        proximity = TRUE,
                        preProcess = c("center", "scale"),
                        ntree=50,
                        #sampsize=rep(1,length(unique(ClassLabels))),
                        #strata= RefData$ClassLabels,
                        #classwt=table(ClassLabels)/sum(table(ClassLabels)),
                        #trControl = train.control,
                        ...)
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
                        preProcess = c("center", "scale"),
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
#' @usage Predictors <- FeatureSelector(ExpData = as.matrix(SeuratObject@data), PCs = 10, num = 2000)
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

#' An internal function to prepare a training/test dataset for model generation.
#' @param ExpData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param Predictors the predictor feature list selected by FeatureSelector.
#' @param ClassLabels [optional] A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param alpa variation cutoff for filtering data
#' @keywords data preparation
#' @export
#' @usage trainingData <- DataReshaper(ExpData = as.matrix(SeuratObject@data), Predictors = genes, ClassLabels = SeuratObject@meta.data$CellTypes)
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

