#HieRFIT functions:

#' The HieRFIT Class
#' @slot Prior
#' @slot ClassProbs
#' @slot ClassWeights
#' @slot CertaintyValues
#' @slot ClassScores
#' @slot Evaluation
HieRFIT <- setClass(Class = "HieRFIT",
                    slots = c(Prior = "character",
                              ClassProbs = "data.frame",
                              ClassWeights = "data.frame",
                              CertaintyValues = "data.frame",
                              ClassScores = "data.frame",
                              Evaluation = "data.frame"))

#' HieMetrics class
#' @slot Pvotes
#' @slot QueWs
#' @slot QueCers
#' @slot Scores
HieMetrics <- setClass(Class = "HieMetrics",
                       slots = c(Pvotes = "data.frame",
                                 QueWs = "data.frame",
                                 QueCers = "data.frame",
                                 Scores = "data.frame")
)

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

#' A hierarchical Classifier Tree Traverser function.
#' @param Query is the input query data. rows are genes and columns are cells.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
#' @param hiemods models list from HieRandForest function.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
CTTraverser <- function(Query, tree, hiemods, thread=NULL){
  library(doParallel)
  node.list <- DigestTree(tree = tree)
  #Create a table for storing node probabilities.
  Scores <- data.frame(row.names = colnames(Query))
  Pvotes <- data.frame(row.names = colnames(Query))
  QueWs <- data.frame(row.names = colnames(Query))
  QueCers <- data.frame(row.names = colnames(Query))

  if(is.null(thread)){
    #
    for(i in node.list){
      nodeModel <- hiemods[[as.character(i)]][[1]]
      #Create QueData:
      #P_dicts <- colnames(nodeModel$trainingData)
      #P_dicts <- P_dicts[P_dicts != ".outcome"]
      P_dicts <- nodeModel$finalModel$xNames
      nodeQueData <- DataReshaper(Data = Query, Predictors = P_dicts)
      #Calculate node Scores:
      nodeScores <- scoR(model = nodeModel,
                         format = "prob",
                         QueData = nodeQueData,
                         node = i)

      #Tally votes for class from the local model:
      nodePvotes <- PvoteR(model = nodeModel, QueData = nodeQueData)
      #Calculate the probability weights of each class by random permutation:
      nodeQueWs <- graWeighteR(model = nodeModel, QueData = nodeQueData )
      #Estimate Certainty of prediction probabilities per class:
      nodeQueCers <- ceR(qP = nodePvotes, qW = nodeQueWs)

      Scores <- cbind(Scores, nodeScores)
      Pvotes <- cbind(Pvotes, nodePvotes)
      QueWs <- cbind(QueWs, nodeQueWs)
      QueCers <- cbind(QueCers, nodeQueCers)
    } #closes the for loop.
    S_path_prod <- ClassProbCalculator(tree = tree, nodes_P_all = Scores)
    HieMetrxObj <- new(Class = "HieMetrics",
                       Pvotes = Pvotes,
                       QueWs = QueWs,
                       QueCers = QueCers,
                       Scores = S_path_prod)

  }else{# FIX THIS PART! thread is specified. For now, use this only when running on bigMem machines.
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    print(paste("registered cores is", getDoParWorkers(), sep = " "))

    nodeProb <- foreach(i=node.list, .inorder = TRUE, .combine=cbind) %dopar% {
      scoR(model = hiemods[[as.character(i)]][[1]],
           format = "prob",
           QueData = Query,
           node = i)
    }
    stopCluster(cl)
    ProbTab <- cbind(ProbTab, nodeProb)
  }

  return(HieMetrxObj)
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
  # QueData <- DataReshaper(Data = Query, Predictors = P_dicts)

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
