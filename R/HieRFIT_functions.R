#HieRFIT functions:

#' The HieRFIT Class
#' @slot Prior a list of prior class labels if provided by the user.
#' @slot ClassProbs raw class probabilities from each local classifier.
#' @slot ClassWeights probability centers when the null is true for each class.
#' @slot CertaintyValues certainty scores calculated by certainty function.
#' @slot Evaluation table of scores and class projections for each sample.
#' @slot tree classification tree copied from the hiermod.
HieRFIT <- setClass(Class = "HieRFIT",
                    slots = c(Prior = "character",
                              ClassProbs = "data.frame",
                              ClassWeights = "data.frame",
                              ClassCertainties = "data.frame",
                              PathScores = "data.frame",
                              Evaluation = "data.frame",
                              tree = "list"))

#' HieMetrics class
#' @slot Pvotes raw class probabilities from each local classifier.
#' @slot QueWs probability centers when the null is true for each class.
#' @slot QueCers certainty scores calculated by certainty function.
HieMetrics <- setClass(Class = "HieMetrics",
                       slots = c(Pvotes = "data.frame",
                                 QueWs = "data.frame",
                                 QueCers = "data.frame",
                                 PathProbScores = "data.frame"))

#' The main function to project reference data on query in order to identify class labels.
#' @param Query Query data whose components will be labeled with reference data.
#' @param refMod optional input if model exist already. Default is null and generated from scratch. Input can be an caret model object or a .Rdata file.
#' @param Prior prior class labels if exist. For cross comparison. Should correspond row order of the Query.
#' @param Qspecies optional argument to specify cross species information transfer. Default is null. Possible options are 'rat2mouse', 'mouse2rat', 'mouse2human', human2mouse. With respect to model data.
#' @param ortoDict optional argument to specify an existing orthology gene table for inter-species projection. Can be used to avoid re-running Biomart query which takes time.
#' @keywords
#' @export
#' @usage refmod <- readRDS("data/exp_refObj.RDS")
#' @usage hierObj <- HieRFIT(Query = pbmc[["RNA"]]@data, refMod = refmod)
HieRFIT <- function(Query, refMod, Prior=NULL, Qspecies=NULL, ortoDict=NULL, binarize=FALSE, alpha=NULL, ...){
  options(warn=-1)
  if ( class(Query) == "Seurat" ){
    Query_d <- as.matrix(Query[["RNA"]]@data)
    }else{
      Query_d <- as.matrix(Query)
      }

  cat("Preparing the query data...\n")
  rownames(Query_d) <- FixLab(xstring = rownames(Query_d))

  Query_d <- Query_d[apply(Query_d, 1, var) != 0, ]

  if(binarize){
    Query_d <- as.matrix((Query_d > 0) + 0)
  }

  if(is.null(Qspecies)){
    Qspecies <- refMod@species
    cat(paste("Assuming that the query data is same species as refmod:", Qspecies, "\n"))
    }

  if(Qspecies != refMod@species) {
    if(is.null(ortoDict)){
      print(paste("Retrieving the orthologous genes between", Qspecies, "and", refMod@species, "..."))
      ref_f <- ExtractHierModfeatures(refMod = refMod)
      print(length(ref_f))
      ortoDict <- GetOrtologs(Genes_r = ref_f, species_r = refMod@species, species_q = Qspecies)
      assign("ortoDict", ortoDict, .GlobalEnv)
    }
    Query_d <- Query_d[which(rownames(Query_d) %in% ortoDict[, paste(Qspecies, "homolog_associated_gene_name", sep = "_")]), ]
    rownames(Query_d) <- ortoDict[match(rownames(Query_d), ortoDict[, paste(Qspecies, "homolog_associated_gene_name", sep = "_")]), 'external_gene_name']
    print(dim(Query_d))
  }

  HieMetObj <- CTTraverser(Query = Query_d, refMod = refMod, ...)

  #Evaluate scores and run uncertainty function, then, project the class labels.
  if(is.null(alpha)){alpha <- refMod@alphas}
  ScoreEvals <- ScoreEvaluate(ProbCert = HieMetObj@QueCers, ProbScores = HieMetObj@PathProbScores, tree=refMod@tree[[1]], alphaList = alpha, ...)
  if ( class(Query) == "Seurat" ){
    Query@meta.data <- cbind(Query@meta.data[, which(!colnames(Query@meta.data) %in% colnames(ScoreEvals))], ScoreEvals)
    object <- Query
  }else{
    object <- new(Class = "HieRFIT",
                  Prior = as.character(Prior),
                  ClassProbs = HieMetObj@Pvotes,
                  ClassWeights = HieMetObj@QueWs,
                  ClassCertainties = HieMetObj@QueCers,
                  PathScores = ScoreEvals$ScoresArray,
                  Evaluation = ScoreEvals$Evals,
                  tree = refMod@tree)
  }
  cat("Completed!\n")
  return(object)
}

#' A hierarchical Classifier Tree Traverser function.
#' @param Query is the input query data. rows are genes and columns are cells.
#' @param refMod model from HieRandForest function.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
CTTraverser <- function(Query, refMod, thread=NULL, ...){
  cat("Traversing the tree...\n")
  cat("--- Gathering the scores...\n")
  node.list <- DigestTree(tree = refMod@tree[[1]])
  #Create a table for storing node probabilities.
  Pvotes <- data.frame(row.names = colnames(Query))
  QueWs <- data.frame(row.names = colnames(Query))
  QueCers <- data.frame(row.names = colnames(Query))

  Query_R <- Shuffler(df = Query)

  if(is.null(thread)){
    for(i in node.list){
      nodeModel <- refMod@model[[as.character(i)]]
      nodeMlr <- refMod@mlr[[as.character(i)]]
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

      #Use randomized Query data
      nodeQueData_R <- Query_R[which(rownames(Query_R) %in% P_dicts), ]
      nodeQueData_R <- t(nodeQueData_R)
      #Add the missing features matrix
      nodeQueData_R <- cbind(nodeQueData_R, mp_df)

      #Calculate the probability weights of each class by random permutation:
      #nodeQueWs <- CalibratedgraWeighteR(model = nodeModel, mlr = nodeMlr, QueData = nodeQueData, ...)
      nodeQueWs <- CalibratedgraWeighteR(model = nodeModel, mlr = nodeMlr, QueData = nodeQueData_R, ...)
      #Tally votes for class from the local model:
      nodePvotes <- CalibratedPvoteR(model = nodeModel, mlr = nodeMlr, QueData = nodeQueData, ...)

      #Estimate Certainty of prediction probabilities per class:
      nodeQueCers <- ceR(qP = nodePvotes, qW = nodeQueWs)

      Pvotes <- cbind(Pvotes, nodePvotes)
      QueWs <- cbind(QueWs, nodeQueWs)
      QueCers <- cbind(QueCers, nodeQueCers)

    } #closes the for loop.
    QueScores <- Pvotes/QueWs
    PathProbScores <- ClassProbCalculator(tree = refMod@tree[[1]], nodes_P_all = QueScores)

    HieMetrxObj <- new(Class = "HieMetrics",
                       Pvotes = Pvotes,
                       QueWs = QueWs,
                       QueCers = QueCers,
                       PathProbScores = PathProbScores)

  }else{# FIX THIS PART! thread is specified. For now, use this only when running on bigMem machines.
    library(doParallel)
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    print(paste("registered cores is", getDoParWorkers(), sep = " "))

    nodeProb <- foreach(i=node.list, .inorder = TRUE, .combine=cbind) %dopar% {
      scoR(model = hiemods[[as.character(i)]],#Deprecated function
           format = "prob",
           QueData = Query,
           node = i)
    }
    stopCluster(cl)
    ProbTab <- cbind(ProbTab, nodeProb)
  }

  return(HieMetrxObj)
}

#' Tree traverser for collecting class probabilities.
#' @param Query is the input query data. rows are genes and columns are cells.
#' @param refMod model from HieRandForest function.
ProbTraverser <- function(Query, refMod, ...){

  rownames(Query) <- FixLab(xstring = rownames(Query))
  node.list <- DigestTree(tree = refMod@tree[[1]])
  #Create a table for storing node probabilities.
  Pvotes <- data.frame(row.names = colnames(Query))
  for(i in node.list){
    nodeModel <- refMod@model[[as.character(i)]]
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
    Pvotes <- cbind(Pvotes, nodePvotes)
  }

  return(Pvotes)
}

#' A function to tally votes from the model locally trained.
#' @param model local classifier.
#' @param QueData is the prepared data matrix ready to be used in predict
#' @param format type of prediction output, "prob" or "resp".
PvoteR <- function(model, QueData, format="prob", node=NULL){
  if(is.null(node)){
    QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "prob", scale=T, center=T))
  } else{
    if(format == "prob"){
      QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "prob", scale=T, center=T))
      colnames(QuePvotes) <- paste(node, colnames(QuePvotes), sep = "")
    } else {
      QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "raw", scale=T, center=T))
      colnames(QuePvotes) <- as.character(node)
    }
  }
  return(QuePvotes)
}

#' A function to tally votes from the model locally trained.
#' @param model local classifier.
#' @param QueData is the prepared data matrix ready to be used in predict
#' @param format type of prediction output, "prob" or "resp".
CalibratedPvoteR <- function(model, mlr, QueData, format="prob", package='caret', node=NULL, ...){

  if(package == 'scikit'){
    if(is.null(node)){
      QuePvotes <- as.data.frame(ski_predict(fit=model, newdata=QueData, type = "prob"))
    }else{
      if(format == "prob"){
        QuePvotes <- as.data.frame(ski_predict(fit=model, newdata=QueData, type = "prob"))
        colnames(QuePvotes) <- paste(node, colnames(QuePvotes), sep = "")
      }else{
        QuePvotes <- as.data.frame(ski_predict(fit=model, newdata=QueData, type = "raw"))
        colnames(QuePvotes) <- as.character(node)
      }
    }
  }else{
    if(is.null(node)){
      QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "prob", scale=T, center=T))
    }else{
      if(format == "prob"){
        QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "prob", scale=T, center=T))
        colnames(QuePvotes) <- paste(node, colnames(QuePvotes), sep = "")
      }else{
        QuePvotes <- as.data.frame(predict(model, newdata=QueData, type = "raw", scale=T, center=T))
        colnames(QuePvotes) <- as.character(node)
      }
    }
  }
  #Calibrate the probabilities:
  QuePvotes.cal <- predict(mlr, newdata=QuePvotes, type="prob")

  if(is.null(dim(QuePvotes.cal))){
    QuePvotes.cal <- data.frame(1-QuePvotes.cal, QuePvotes.cal)
    colnames(QuePvotes.cal) <- colnames(QuePvotes)#Double-check the colnames.
  }

  return(QuePvotes.cal)
}

#' A function to calculate gravitity center (weight center) of the probability distributions
#' among classes of the node. By permutation.
#' @param model of the node.
#' @param mlr multinominal logistic regression model for calibrating the probabilities with Platt scaling, aka sigmoid calibration.
#' @param QueData same matrix used in PvoteR().
#' @return QueWeights a set of probability weightes per class to be used in asymetrix entropy estimations.
CalibratedgraWeighteR <- function(model, mlr, QueData, ...){
  #Randomizing only feature space
  QueData_R <- Shuffler(df = QueData)
  pvts_R <- CalibratedPvoteR(model = model, mlr = mlr, QueData = QueData_R, ...)
  Ws <- apply(pvts_R, 2, mean)
  QueWeights <- t(as.data.frame(Ws))[rep(1, each=nrow(QueData)), ]
  QueWeights <- as.data.frame(QueWeights)
  return(QueWeights)
}

#' Certainty Estimation function. Calculates certainty values for each class probability.
#' @param qP Observed probability for each class
#' @param qW Probability weights estimated by graWeighteR()
#' @return QueCers
ceR <- function(qP, qW){
  qW <- qW + 1e-05#to prevent math error.
  QueCers <- NULL
  qL <- qP > qW
  qL <- ifelse(qL == TRUE, 1, -1)
  rownames(qW) <- rownames(qP)#Otherwise, duplicate rownames of qW gives error.
  QueCers <- (qL*(qP-qW)^2)/( ((1-2*qW)*qP)+qW^2 )
  colnames(QueCers) <- colnames(qP)
  QueCers <- as.data.frame(QueCers)
  return(QueCers)
}

#' A function for evalating the uncertainty.
#' @param ProbCert Certainty scores.
#' @param tree a tree topology with which hrf model was trained on. 'phylo' format.
#' @param full.tbl if true, only the entire score table is returned.
#' @param alphaList list of alpha threshold of each class.
ScoreEvaluate <- function(ProbCert, ProbScores, tree, full.tbl=FALSE, alphaList=NULL, ...){
  cat("--- Evaluating the scores...\n")
  labs_l <- c(tree$tip.label, tree$node.label)#The order is important! tips first. Don't change!#New
  labs_l <- labs_l[!labs_l %in% "TaxaRoot"]

  df <- data.frame(row.names = rownames(ProbCert))
  path_sums <- data.frame(row.names = rownames(ProbCert))
  paths_outs <- data.frame(row.names = rownames(ProbCert))
  all_sibs <- data.frame(row.names = rownames(ProbCert))
  for(cl in labs_l){
    AncPath <- GetAncestPath(tree = tree, class = cl, labels = T)

    clSibs <- c()
    for(c in AncPath){
      sib <- GetSiblings(tree = tree, class = c)
      clSibs <- append(clSibs, sib)
    }

    if(length(clSibs) > 1){
      sibs_out <- data.frame(apply(ProbCert[, clSibs], 1, sum))
    }else{
      sibs_out <- data.frame(ProbCert[, clSibs])
    }

    PathOuts <- c("TaxaRoot_OutGroup")#NEW
    if(length(AncPath) > 1){
      cpaths_sum <- data.frame(apply( ProbCert[, AncPath], 1, sum))

      for(n in AncPath[-1]){
        outg <- grep(n, grep("OutGroup", names(ProbCert), value = T), value = T)
        PathOuts <- append(PathOuts, outg)
      }

      if(length(PathOuts) > 1){
        cpaths_out <- data.frame(apply(ProbCert[, PathOuts], 1, sum))
      }else{
        cpaths_out <- data.frame(ProbCert[, PathOuts])
      }

    }else{
      cpaths_sum <- data.frame( ProbCert[, AncPath])
      cpaths_out <- data.frame(rep(0, length(ProbCert[,1])), row.names = rownames(ProbCert))
    }
    colnames(cpaths_sum) <- cl
    colnames(cpaths_out) <- cl
    colnames(sibs_out) <- cl

    path_sums <- cbind(path_sums, cpaths_sum)
    paths_outs <- cbind(paths_outs, cpaths_out)
    all_sibs <- cbind(all_sibs, sibs_out)
  }

  CertScore = path_sums - paths_outs - all_sibs### ****

  if(full.tbl){
    return(CertScore)
  }else{
    if(!is.null(alphaList)){
      cat("--- Computing the path scores...\n")
      #cl.CertScore.cmax <- apply(CertScore, 1, GetCls5, alphaList=alphaList, tree=tree, class=TRUE)
      #CertScoremax <- apply(CertScore, 1, GetCls5, alphaList=alphaList, tree=tree, class=FALSE)
      cl.CertScore.cmax <- sapply(seq_len(nrow(CertScore)), GetCls_beta, Upath=CertScore, Ppath=ProbScores, alphaList=alphaList, tree=tree, class=TRUE)
      CertScoremax <- sapply(seq_len(nrow(CertScore)), GetCls_beta, Upath=CertScore, Ppath=ProbScores, alphaList=alphaList, tree=tree, class=FALSE)
    }else{
      cl.CertScore.cmax <- colnames(CertScore)[apply(CertScore, 1, which.max)]
      CertScoremax <- apply(CertScore, 1, function(x) max(x))
    }
    df <- data.frame(Score = CertScoremax,
                     Projection = cl.CertScore.cmax)

    return(list(Evals = df, ScoresArray=CertScore))
  }
}


GetCls_beta <- function(i, Upath, Ppath, alphaList, tree, class, ...){
  #Select the candidates:
  log.list <- vector()
  for(c in names(alphaList)){
    x1=as.numeric(Upath[i,])
    x2=as.numeric(alphaList[[c]]["50%",])
    ed <- dist(rbind(x1,x2))
    logic.c <- ed/length(x1) < 1.00
    names(logic.c) <- c
    log.list <- c(log.list, logic.c)
  }

  if(any(log.list)){
    candits <- names(which(log.list))
    #Pick the best class among the candidates
    if(length(candits) == 1){
      classL <- candits
      }else{
        classL <- colnames(Ppath[i, candits])[apply(Ppath[i, candits], 1, which.max)]
      }
    classS <- Ppath[i, classL]
  }else{
    classL <- "Undetermined"
    classS <- max(Ppath[i,])
  }

  if(class){
    return(classL)
  }else{
    return(classS)
  }
}

#' A function to calculate class scores for all internal and tip node classes.
#' @param tree
#' @param nodes_P_all a table output from CTTraverser() which contains all class probabilities from every node.
#' @return P_path_prod a table for products of ancestor node scores.
ClassProbCalculator <- function(tree, nodes_P_all){
  P_path_prod <- data.frame(row.names = rownames(nodes_P_all))
  clabs <- c(tree$tip.label, tree$node.label)[tree$edge[, 2]]
  for(cl in clabs){
    map <- GetAncestPath(tree = tree, class = cl,labels = T)

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

