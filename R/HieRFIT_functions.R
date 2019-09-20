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
#' @param ortoDict optional argument to specify an existing orthology gene table for inter-species projection. Can be used to avoid re-running Biomart query which takes time.
#' @keywords
#' @export
#' @usage expRefObj <- get(load("data/exp_refObj.Rdata"))
#' @usage cpo <- HieRFIT(Query = as.matrix(pbmc1@data), refMod = expRefobj)
HieRFIT <- function(Query, refMod, Prior=NULL, xSpecies=NULL, ortoDict=NULL, alpha=.9){

  if (class(Query) == "seurat" | class(Query) == "Seurat" ){
    Query_d <- as.matrix(Query@data)}else{Query_d <- Query}

  rownames(Query_d) <- FixLab(xstring = rownames(Query_d))

  if( !is.null(xSpecies)) {
    if(xSpecies == "mouse2rat"){ ## pay attention! flipped logic.
      if(is.null(ortoDict)){
        print("Rat to mouse gene id conversion...")
        ortoDict <- Gmor(RatGenes = rownames(Query_d))
        assign("ortoDict", ortoDict, .GlobalEnv)
      }
      ortoDict$external_gene_name <- FixLab(ortoDict$external_gene_name)
      ortoDict$mmusculus_homolog_associated_gene_name <- FixLab(ortoDict$mmusculus_homolog_associated_gene_name)
      Query_d <- Query_d[which(rownames(Query_d) %in% ortoDict$external_gene_name), ]
      rownames(Query_d) <- ortoDict[match(rownames(Query_d), ortoDict$external_gene_name),]$mmusculus_homolog_associated_gene_name
    }
  }
  if(refMod@modtype == "hrf"){
    HieMetObj <- CTTraverser(Query = Query_d, tree = refMod@tree[[1]], hiemods = refMod@model)
  }else{#FIX THIS
    P_path_prod <- Predictor(model = refMod@model[[1]], Query = Query_d)
  }
  #Evaluate scores and run uncertainty function, then, project the class labels.
  ScoreEvals <- ScoreEval(ScoreObs = HieMetObj@Scores, tree=refMod@tree[[1]], ProbCert = HieMetObj@QueCers, alpha=alpha)

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

#' A hierarchical Classifier Tree Traverser function.
#' @param Query is the input query data. rows are genes and columns are cells.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
#' @param hiemods models list from HieRandForest function.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
CTTraverser <- function(Query, tree, hiemods, thread=NULL){

  node.list <- DigestTree(tree = tree)
  #Create a table for storing node probabilities.
  Scores <- data.frame(row.names = colnames(Query))
  Pvotes <- data.frame(row.names = colnames(Query))
  QueWs <- data.frame(row.names = colnames(Query))
  QueCers <- data.frame(row.names = colnames(Query))
  fi=node.list[length(node.list)]+1
  node.list <- c(node.list, fi)
  if(is.null(thread)){
    for(i in node.list){
      #nodeModel <- hiemods[[as.character(i)]][[1]]
      nodeModel <- hiemods[[as.character(i)]]
      #Create QueData:
      P_dicts <- FixLab(xstring = nodeModel$finalModel$xNames)
      nodeQueData <- Query[which(rownames(Query) %in% P_dicts), ]
      nodeQueData <- t(nodeQueData)
      #Add the missing features matrix
      mp <- P_dicts[which(!P_dicts %in% colnames(nodeQueData))]
      mp_df <- data.frame(matrix(0,
                                 ncol = length(mp),
                                 nrow = length(colnames(Query))))
      colnames(mp_df) <- mp
      nodeQueData <- cbind(nodeQueData, mp_df)

      #Tally votes for class from the local model:
      nodePvotes <- PvoteR(model = nodeModel, QueData = nodeQueData)
      #Calculate the probability weights of each class by random permutation:
      nodeQueWs <- graWeighteR(model = nodeModel, QueData = nodeQueData )
      #Estimate Certainty of prediction probabilities per class:
      nodeQueCers <- ceR(qP = nodePvotes, qW = nodeQueWs)
      #Calculate node Scores:
      #nodeScores <- scoR(model = nodeModel,
      #                   format = "prob",
      #                   QueData = nodeQueData,
      #                   node = i)
      nodeScores <- nodePvotes/nodeQueWs
      colnames(nodeScores) <- paste(i, colnames(nodeScores), sep = "")
      Scores <- cbind(Scores, nodeScores)
      if(i != fi){
        Pvotes <- cbind(Pvotes, nodePvotes)
        QueWs <- cbind(QueWs, nodeQueWs)
        QueCers <- cbind(QueCers, nodeQueCers)
      }
    } #closes the for loop.
    #S_path_prod <- ClassProbCalculator(tree = tree, nodes_P_all = Scores)
    S_path_prod <- ClassProbCalculator2(tree = tree, nodes_P_all = Scores)
    HieMetrxObj <- new(Class = "HieMetrics",
                       Pvotes = Pvotes,
                       QueWs = QueWs,
                       QueCers = QueCers,
                       Scores = S_path_prod)

  }else{# FIX THIS PART! thread is specified. For now, use this only when running on bigMem machines.
    library(doParallel)
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    print(paste("registered cores is", getDoParWorkers(), sep = " "))

    nodeProb <- foreach(i=node.list, .inorder = TRUE, .combine=cbind) %dopar% {
      scoR(model = hiemods[[as.character(i)]],
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
  QueData_R <- RandomizeR(df = QueData, n = 20)
  pvts_R <- PvoteR(model = model, QueData = QueData_R)
  Ws <- apply(pvts_R, 2, mean) + apply(pvts_R, 2, sd)
  #Ws <- apply(pvts_R, 2, mean) + apply(pvts_R, 2, sd)/sqrt(dim(pvts_R)[1])
  #Ws <- apply(pvts_R, 2, mean)
  #Ws <- colMeans(PvoteR(model = model, QueData = QueData_R))
  QueWeights <- t(as.data.frame(Ws))[rep(1, each=nrow(QueData)), ]
  QueWeights <- as.data.frame(QueWeights)
  return(QueWeights)
}

#' Certainty Estimation function. Calculates certainty values for each class probability.
#' @param qP Observed probability for each class
#' @param qW Probability weights estimated by graWeighteR()
#' @return QueCers
ceR <- function(qP, qW){
  QueCers <- NULL
  qL <- qP > qW
  qL <- ifelse(qL == TRUE, 1, -1)
  rownames(qW) <- rownames(qP)#Otherwise, duplicate rownames of qW gives error.
  QueCers <- (qL*(qP-qW)^2)/( ((1-2*qW)*qP)+qW^2 )
  colnames(QueCers) <- colnames(qP)
  QueCers <- as.data.frame(QueCers)
  return(QueCers)
}

#' The predictor function. REDUNDANT!!!
#' @param model
#' @param QueData
#' @param format type of prediction output, "prob" or "resp".
scoR <- function(model, QueData, format="prob", node=NULL){

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

#' A function to find the full path of all true ancestors:
#' @param PCertVector Probability certainty table.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
CandidateDetector <- function(PCertVector, tree){

  labs_l <- c(tree$tip.label, tree$node.label)#The order is important! tips first. Don't change!#New
  labs_l <- labs_l[!labs_l %in% "TaxaRoot"] #Double-check#New
  Path_nodes_of_candits <- NULL
  CandidNodes <- NULL
  for(node.lab in labs_l){
    AncPath <- GetAncestPath(tree = tree, class = node.lab, labels = T)
    if(all(as.logical(PCertVector[AncPath])) & !(AncPath[1] %in% Path_nodes_of_candits)){
      CandidNodes <- c(CandidNodes, AncPath[1])
      Path_nodes_of_candits <- unique(c(Path_nodes_of_candits, AncPath[2:length(AncPath)]))
    }
  }
  return(CandidNodes[!CandidNodes %in% Path_nodes_of_candits])
}

#' A function to find the full path of all true ancestors:
#' @param PCertVector Probability certainty table.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
CandidateDetector2 <- function(PCertVector, tree, alpha=0){
  labs_l <- c(tree$tip.label, tree$node.label)#The order is important! tips first. Don't change!#New
  labs_l <- labs_l[!labs_l %in% "TaxaRoot"] #Double-check#New
  Path_nodes_of_candits <- NULL
  CandidNodes <- NULL
  for(node.lab in labs_l){
    AncPath <- GetAncestPath(tree = tree, class = node.lab, labels = T)
    if( (mean(as.numeric(PCertVector[AncPath])) > alpha) &
        !(AncPath[1] %in% Path_nodes_of_candits) &
        (as.numeric(PCertVector[AncPath[1]]) > alpha)
        ){
      CandidNodes <- c(CandidNodes, AncPath[1])
      Path_nodes_of_candits <- unique(c(Path_nodes_of_candits, AncPath[2:length(AncPath)]))
    }
  }
  return(CandidNodes)
}

#' A function for evalating the uncertainty.
#' @param ScoreObs P_path_prod for observed scores
#' @param ProbCert Certainty scores.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
ScoreEval <- function(ScoreObs, ProbCert, tree, alpha=0){

  #ProbCert.logic <- ProbCert > alpha#New

  df <- data.frame(row.names = rownames(ScoreObs))
  for(i in 1:length(ProbCert[,1])){
    #candits <- colnames(ProbCert)[ProbCert[i,] > alpha]
    #candits <- CandidateDetector(PCertVector = ProbCert.logic[i,], tree=tree)
    candits <- CandidateDetector2(PCertVector = ProbCert[i,], tree=tree, alpha=alpha)
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
      #classL <- colnames(ProbCert[i,candits])[apply(ProbCert[i,candits], 1, which.max)]
      classU <- ProbCert[i,classL]
      classS <- ScoreObs[i,classL]
    }

    df <- rbind(df, data.frame(row.names = rownames(ProbCert[i,]),
                               Score = classS,
                               Certainty = classU,
                               Projection = classL))
  }
  return(df)
}

#' A function to calculate class scores for all internal and tip node classes.
#' @param tree
#' @param nodes_P_all a table output from CTTraverser() which contains all class probabilities from every node.
#' @return P_path_prod a table for products of ancestor node scores.
ClassProbCalculator <- function(tree, nodes_P_all){
  P_path_prod <- data.frame(row.names = rownames(nodes_P_all))
  clabs <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]
  for(cl in clabs){
    map <- GetAncestPath(tree = tree, class = cl)
    if(length(map) > 1){
      nt_prob <- data.frame(apply(nodes_P_all[, map], 1, prod))
    }else{
      nt_prob <- data.frame(nodes_P_all[,map])
    }
    colnames(nt_prob) <- cl
    P_path_prod <- cbind(P_path_prod, nt_prob)
  }
  return(P_path_prod)
}

ClassProbCalculator2 <- function(tree, nodes_P_all){
  P_path_prod <- data.frame(row.names = rownames(nodes_P_all))
  clabs <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]
  fi=length(c(tree$tip.label, tree$node.label))+1
  for(cl in clabs){
    map <- GetAncestPath(tree = tree, class = cl)
    if(cl %in% tree$tip.label){
      map <- c(map, paste(fi, cl, sep = ""))
    }
    if(length(map) > 1){
      nt_prob <- data.frame(apply(nodes_P_all[, map], 1, prod))
    }else{
      nt_prob <- data.frame(nodes_P_all[, map])
    }
    colnames(nt_prob) <- cl
    P_path_prod <- cbind(P_path_prod, nt_prob)
  }
  return(P_path_prod)
}
