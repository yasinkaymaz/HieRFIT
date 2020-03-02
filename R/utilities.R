#Utility functions:

NoiseInject <- function(RefData, ClassLabels, refmod){
  NoisedRef <- list()
  Caldata <- RefData
  rownames(Caldata) <- FixLab(rownames(Caldata))#Not necessary
  features <- NULL
  for( i in names(refmod@model)){ features <- append(features, refmod@model[[i]]$finalModel$xNames)}
  features <- unique(features)
  Caldata <- Caldata[features, ]
  CalY <- FixLab(ClassLabels)
  Caldata.noised <- Caldata

  # max(table(CalY))  #For balancing the class data
  # Caldata[, CalY == unique(CalY)[1]]

  rands <- seq(1:20)
  CaldataFull <- Caldata
  for(r in rands){
    for(x in 1:length(Caldata[1,])){
      nz_idx <- which(!Caldata[,x] %in% c(0))
      if(r > length(nz_idx)){
        r = length(nz_idx)-1
      }
      idx.to.set.zero <- sample(nz_idx, r)
      Caldata.noised[idx.to.set.zero, x] <- 0
    }
    CaldataFull <- cbind(CaldataFull, Caldata.noised)
  }

  TaxaOutdata <- Caldata[, sample(colnames(Caldata), size = 500)]#Replace 500 with an adaptive.
  TaxaOutdata <- Shuffler(df = TaxaOutdata)
  CaldataFull <- cbind(CaldataFull, TaxaOutdata)

  colnames(CaldataFull) <- make.names(colnames(CaldataFull), unique = TRUE)
  CalYFull <- c(rep(CalY, length(rands)+1), rep("TaxaOut", 500))

  NoisedRef[["data"]] <- CaldataFull
  NoisedRef[["Y"]] <- CalYFull

  return(NoisedRef)
}


ProbTraverser <- function(Query, refmod){

  rownames(Query) <- FixLab(xstring = rownames(Query))
  node.list <- DigestTree(tree = refmod@tree[[1]])
  #Create a table for storing node probabilities.
  Pvotes <- data.frame(row.names = colnames(Query))
    for(i in node.list){
      nodeModel <- refmod@model[[as.character(i)]]
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


GetSigMods2 <- function(RefData, ClassLabels, refmod){

  Cal.data <- NoiseInject(RefData = RefData, ClassLabels = ClassLabels, refmod = refmod)
  Pvotes <- ProbTraverser(Query = Cal.data[["data"]], refmod = refmod)
  labs_l <- c(refmod@tree[[1]]$tip.label, refmod@tree[[1]]$node.label)
  votes <- data.frame(Prior=Cal.data[["Y"]], Pvotes)
  node.list <- DigestTree(tree = refmod@tree[[1]])
  for(i in node.list){
    nodeModel <- refmod@model[[as.character(i)]]
    outGname <- grep("OutGroup", nodeModel$finalModel$classes, value=T)
    Leafs <- NULL
    classes <- nodeModel$finalModel$classes[! nodeModel$finalModel$classes %in% outGname]
    node.dict <- list()
    Labs <- votes$Prior
    for(l in classes){
      if(l %in% refmod@tree[[1]]$tip.label){
        node.dict[[l]] <- l
      }else{
        leaf <- GetChildNodeLeafs(tree = refmod@tree[[1]], node = match(l, labs_l))
        leaf <- leaf[which(lapply(leaf, length)>0)]
        for(j in 1:length(leaf)){
          Leafs <- append(Leafs, leaf[[j]])
        }
        node.dict[[l]] <- Leafs
        Leafs <- NULL
      }
      Labs <- ifelse(Labs %in% node.dict[[l]], l, as.character(Labs))
    }
    Labs <- ifelse(!Labs %in% c(names(node.dict), node.dict), outGname, as.character(Labs))
    print(table(Labs))
    dfr.full <- data.frame(Prior=Labs, votes[, nodeModel$finalModel$classes])
    mlrmod <- nnet::multinom(Prior ~ ., data = dfr.full, model=FALSE)
    #refmod@model[[as.character(i)]]$mlr <- stripGlmModel(mlrmod)
    refmod@model[[as.character(i)]]$mlr <- mlrmod
  }
  return(refmod)
}


GetSigMods <- function(RefData, ClassLabels, refmod){

  Cal.data <- NoiseInject(RefData = RefData, ClassLabels = ClassLabels, refmod = refmod)
  Pvotes <- ProbTraverser(Query = Cal.data[["data"]], refmod = refmod)

  votes <- data.frame(Prior=Cal.data[["Y"]], Pvotes)
  classes <- unique(votes$Prior)
  labs_l <- c(refmod@tree[[1]]$tip.label, refmod@tree[[1]]$node.label)
  sigmoidMods <- list()

  for(n in 1:length(labs_l)){
    if(labs_l[n] == "TaxaRoot"){next}
    Leafs <- NULL
    if(labs_l[n] %in% refmod@tree[[1]]$tip.label){
      Leafs <- labs_l[n]
    }else{
      leaf <- GetChildNodeLeafs(tree = refmod@tree[[1]], node = n)
      leaf <- leaf[which(lapply(leaf, length)>0)]
      for(j in 1:length(leaf)){
        Leafs <- append(Leafs, leaf[[j]])
      }
    }
    votes <- votes[order(votes[[labs_l[n]]]),]
    Labs <- ifelse(votes$Prior %in% Leafs, 1, 0)
    dfr <- data.frame(votes[[labs_l[n]]], Labs)
    colnames(dfr) <- c("x","y")
    # training a logistic regression model on the cross validation dataset
    glmmod <- glm(y~x, data = dfr, family = binomial, maxit = 100, y=FALSE, model=FALSE)
    sigmoidMods[[labs_l[n]]] <- stripGlmModel(glmmod)
  }
  outs <- grep("OutGroup", names(votes), value = T)
  for(x in 1:length(outs)){
    int.node <- gsub("_OutGroup", "", outs[x])
    n <- match(int.node, labs_l)
    Leafs <- NULL
    leaf <- GetChildNodeLeafs(tree = refmod@tree[[1]], node = n)
    leaf <- leaf[which(lapply(leaf, length)>0)]
    for(j in 1:length(leaf)){
      Leafs <- append(Leafs, leaf[[j]])
    }
    votes <- votes[order(votes[[outs[x]]]),]
    Labs <- ifelse(votes$Prior %in% Leafs, 0, 1)#Labels (1,0) reversed since outgroup
    dfr <- data.frame(votes[[outs[x]]], Labs)
    colnames(dfr) <- c("x","y")
    # training a logistic regression model on the cross validation dataset
    glmmod <- glm(y~x, data = dfr, family = binomial, maxit = 100, y=FALSE, model=FALSE)
    sigmoidMods[[outs[x]]] <- stripGlmModel(glmmod)
  }

  return(sigmoidMods)
}



stripGlmModel = function(cm) {
  cm$y = c()
  cm$model = c()
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  return(cm)
}


DigestTree <- function(tree) {
  all.nodes <- unique(tree$edge[,1])
  return(all.nodes)
}

NodePredictorImportance <- function(treeTable, RefMod){

  tree <- CreateTree(treeTable = treeTable)
  node.list <- DigestTree(tree = tree)
  pdf("out.pdf")
  for(i in node.list){
    varImpPlot(RefMod@model[[as.character(i)]][[1]]$finalModel,
               n.var = 10,
               main = paste("Important predictors of node", tree$node.label[i-length(tree$tip.label)], sep="\n"))
    par()
  }
  dev.off()

}


ExtractHierModfeatures <- function(RefMod){
  node.list <- DigestTree(tree = RefMod@tree[[1]])
  Tufs <- NULL
  for(i in node.list){
    Tufs <- c(Tufs, RefMod@model[[as.character(i)]]$finalModel$xNames)
  }
  Tufs <- unique(Tufs)
  return(Tufs)
}

SaveHieRMod <- function(refMod, fileName="Reference.HierMod"){
  lapply(refMod@model,
    function(x) rm(list=ls(envir = attr(x$terms, ".Environment")),
    envir = attr(x$terms, ".Environment")))
  lapply(refMod@model,
    function(x) environment(x$terms) <- NULL)
  saveRDS(refMod, file = paste(fileName,".RDS", sep = ""))
}

LoadHieRMod <- function(fileName){
  refMod <- readRDS(file = fileName)
  lapply(refMod@model,
    function(x) attr(x$terms, ".Environment") <- globalenv() )
  return(refMod)
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
      #print(paste("extracting tips for node", c, sep=" "))
      c.tips[[c]] <- ape::extract.clade(tree, c)$tip.label
    }else{
      c.tips[[c]] <- tree$tip.label[ordered_tips][match(c, ordered_tips)]
    }
  }
  return(c.tips)
}

GetAncestPath <- function(tree, class, labels=FALSE){
  path <- c()
  labs_l <- c(tree$tip.label, tree$node.label)
  Node <- match(class, labs_l)
  parent <- tree$edge[which(x = tree$edge[, 2] == Node), ][1]
  while(!is.na(parent)){
    if(labels){
    path <- c(path, class)
    }else{
      path <- c(path, paste(parent, class, sep = ""))
    }
    class <- labs_l[parent]
    parent <- tree$edge[which(x = tree$edge[, 2] == parent), ][1]
  }
  return(path)
}


#' A function to fix class labels.
#' @param xstring is a list of class labels in character.
FixLab <- function(xstring){
  #Replace white space with '_'
  xstring <- gsub(xstring, pattern = " ", replacement = "_")
  xstring <- gsub(xstring, pattern = "\\+|-|/", replacement = ".")
  xstring <- gsub(xstring, pattern = "`|,", replacement = "")
  return(xstring)
}

CVRunner <- function(Ref, ClassLabels, TreeTable=NULL, cv_k=5, method="hrf"){

  # Create K-fold data split:
  CMtab <- vector("list", length = cv_k)
  flds <- caret::createFolds(ClassLabels, k = cv_k, list = TRUE, returnTrain = FALSE)
  hPRF_cv <- do.call(rbind.data.frame,
                     lapply(1:cv_k,
                            function(i){
                              #Create hieR
                              trainClassLabels <- ClassLabels[-flds[[i]]]
                              trainRef <- Ref[, -flds[[i]]]

                              refmod <- CreateHieR(RefData = trainRef,
                                                   ClassLabels = trainClassLabels,
                                                   TreeTable = TreeTable)
                              #Hierfit
                              testClassLabels <- ClassLabels[flds[[i]]]
                              testRef <- Ref[, flds[[i]]]


                              testObj <- HieRFIT(Query = testRef, refMod = refmod, Prior = testClassLabels)
                              PriorPostTable <- data.frame(Prior=testObj@Prior, Projection = testObj@Evaluation$Projection)
                              hPRF.out <- hPRF(tpT = PriorPostTable, tree = refmod@tree[[1]])
                              hPRF.out <- c(Tool="HieRFIT",
                                            CV_k=i,
                                            hPRF.out)

                              print(hPRF.out)
                              cc <- t(hPRF.out)
                              return(cc)
                            }
                     )
  )
  return(hPRF_cv)
}

#' A function to evaluate performance of HieRFIT with various Certainty thresholds.
#' @param RefSeuObj reference seurat object.
#' @param IdentityCol the column 'number' in the metadata slot showing the cell type identities.
#' @param RefMod
#' @param Uinter number of intervals to assess certainty threshold.
#' @param perm_n # of permutations.
#' @param samp_n # samples to draw at each permutation.
EvaluateCertainty <- function(RefSeuObj, IdentityCol, RefMod, Uinter=20, perm_n=10, samp_n=100){
  hPRF_table <- numeric()
  for(p in 1:perm_n){
    testcells <- sample(rownames(RefSeuObj@meta.data), samp_n)
    for( ai in seq(1:Uinter)/Uinter){
      testObj <- HieRFIT(Query = as.matrix(RefSeuObj@data)[, testcells],
                         refMod = RefMod,
                         Prior = RefSeuObj@meta.data[testcells, IdentityCol],
                         alpha = ai)

      PriorPostTable <- data.frame(Prior = testObj@Prior,
                                   Projection = testObj@Evaluation$Projection)
      hPRF.out <- hPRF(tpT = PriorPostTable, tree = RefMod@tree[[1]])
      hPRF.out <- c(Perm=p,
                    U_threshold=ai,
                    hPRF.out,
                    Undetermined=nrow(PriorPostTable[PriorPostTable$Projection == "Undetermined",])/length(PriorPostTable[,1]))
      print(hPRF.out)
      hPRF_table <- rbind(hPRF_table, hPRF.out)
    }
  }
  return(hPRF_table)
}

Rander <- function(df){
  if(dim(df)[1] > 1000){Sn <- 1000}else{Sn <- dim(df)[1]}
  dfr <- as.matrix(df[sample(row.names(df), size = Sn),])
  for(i in 1:length(dfr[1,])){#for each col:
    for(j in 1:length(dfr[,1])){#for each row:
      dfr[,i] <- sample(dfr[,i])
      dfr[j,] <- sample(dfr[j,])
    }
  }
  return(as.data.frame(dfr))
}

RandomizeR <- function(df, n=10){
  #set.seed(192939)
  dfRand <- NULL
  for (i in 1:n){
    dfR <- df[sample(nrow(df)), sample(ncol(df))]
    rownames(dfR) <- rownames(df)
    colnames(dfR) <- colnames(df)
    dfRand <- rbind(dfRand, dfR)
  }
  dfRand <- dfRand[sample(rownames(dfRand),size = nrow(df)),]
  return(dfRand)
}

FakeRandomizeR <- function(df, seed.num=192939){
  set.seed(seed.num)
  dfRand <- df
  rownames(dfRand) <- sample(rownames(df))
  colnames(dfRand) <- sample(colnames(df))
  return(dfRand)
}

Shuffler <- function(df){
  dfr <- t(apply(df, 1, sample))
  dfr <- apply(dfr, 2, sample)
  colnames(dfr) <- colnames(df)
  rownames(dfr)<-rownames(df)
  return(as.data.frame(dfr))
}

#' Homology mapping via orthologous genes between mouse and rat.
Gmor <- function(RatGenes){
  # This function retrieves mouse homolog associated gene names of Rat genes.
  #library(biomaRt)
  ensembl.rat <- biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl",
                                  host = "www.ensembl.org",
                                  ensemblRedirect = FALSE)
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



GetSiblings <- function(tree, class, labels=FALSE){

  labs_l <- c(tree$tip.label, tree$node.label)
  Node <- match(class, labs_l)
  parent <- tree$edge[which(x = tree$edge[, 2] == Node), ][1]

  #labs_s <- c(tree$tip.label, tree$node.label)
  labs_s <- labs_l[!labs_l %in% c(class, "TaxaRoot")]

  siblings <- c()
  for(cl in labs_s){
    Node.sib <- match(cl, labs_l)
    par.sib <- tree$edge[which(x = tree$edge[, 2] == Node.sib), ][1]

    if( par.sib == parent){
      #print(cl)
      siblings <- append(siblings, cl)
    }
  }

  return(siblings)
}


#' A function to determine the size of intersection between ancestors of True class and Predicted class
#' @param t true class
#' @param p predicted class
EvalPred <- function(t, p, tree){
  #Fix the label characters:
  t <- FixLab(t)
  p <- FixLab(p)
  #make a list of node indexes of the entire tree:
  labs_l <- c(tree$tip.label, tree$node.label)
  #look up the index of prior(t) and predicted(p) labels:
  Node.t <- match(t, labs_l)
  Node.p <- match(p, labs_l)
  #look up the parent node indexes of prior(t) and predicted(p) label indexes:
  parent.t <- tree$edge[which(x = tree$edge[, 2] == Node.t), ][1]
  parent.p <- tree$edge[which(x = tree$edge[, 2] == Node.p), ][1]
  #extract the list of children node indexes: #can be multiple children.
  children <- tree$edge[which(x = tree$edge[, 1] == Node.t), 2]

  if(t %in% tree$node.label){#if the prior node label is an internal node not a leaf.
  #extract the grandChildren node labels if exist.
  grandChilds <- c(ape::extract.clade(tree, t)$node.label, ape::extract.clade(tree, t)$tip.label)
  #Exclude children and self node labels.
  grandChilds <- grandChilds[!grandChilds %in% c(t, labs_l[children])]
  }else{
    grandChilds <- NULL
  }
  #Look up entire path for ancestors. This returns node index and node labels concatinated: e.g. "8B" "7A"
  Ancestors <- GetAncestPath(tree = tree, class = t)
  if(is.na(Node.p) || is.na(Node.t)){
    out <- "NotDefined"
  }else if(any(grep(p, Ancestors))){
    if(t == p){
      out <- "Correct_node"
    }else if(labs_l[parent.t] == p){
      out <- "Correct_parent_node"
    }else{
      out <- "Correct_ancestral_node"
    }
  }else{
    if(parent.t == parent.p){
      out <- "Incorrect_node_sibling"
    }else if(Node.p %in% children){
      out <- "Correct_children_node"
    }else if(p %in% grandChilds){
      out <- "Correct_grandchildren_node"
    }else{
    out <- "Incorrect_clade"
    }
  }
  return(out)
}


#' A function to determine the size of intersection between ancestors of True class and Predicted class
#' @param t true class
#' @param p predicted class
IntSectSize <- function(t, p, tree){
  Ti <- GetAncestPath(tree = tree, class = FixLab(t))
  Pi <- GetAncestPath(tree = tree, class = FixLab(p))
  intL <- length(intersect(Ti, Pi))
  return(intL)
}

#' A function for Hierarchical Precision, Recall, and F-measure.
#' @param tpT PriorPostTable: a table with two columns of which first is Prior and second is Post-prediction.
#' @param tree tree topology in phylo format.
#' @param BetaSq Beta coefficient
#' @param ND_term The label used for undetermined class types.
hPRF <- function(tpT, tree, BetaSq=1, ND_term="Undetermined"){
  # To Do: Consider Undetermined class!
  tpT$Int <- apply(tpT, 1, function(x) IntSectSize(t = x[1], p = x[2], tree = tree))
  tpT$PiL <- apply(tpT, 1, function(x) length(GetAncestPath(tree = tree, class = FixLab(x[2]) )))
  tpT$TiL <- apply(tpT, 1, function(x) length(GetAncestPath(tree = tree, class = FixLab(x[1]) )))

  hP <- sum(tpT$Int)/sum(tpT$PiL)
  hR <- sum(tpT$Int)/sum(tpT$TiL)
  hF <- (BetaSq+1)*hP*hR/(BetaSq*hP+hR)
  #Calculate ND rate
  ND.rate <- dim(tpT[tpT[, 2] == ND_term, ])[1]/dim(tpT)[1]
  #Filter out ND predictions
  tpT.size <- dim(tpT)[1]
  tpT <- tpT[tpT[, 2] != ND_term, ]
  #Calculate correctness
  metrics <- c("Correct_node", "Correct_parent_node", "Correct_ancestral_node",
  "Incorrect_node_sibling", "Correct_children_node", "Correct_grandchildren_node",
  "Incorrect_clade", "NotDefined")
  tpT$Eval <- apply(tpT, 1, function(x) EvalPred(t = x[1], p = x[2], tree = tree) )
  evals <- table(tpT$Eval)/tpT.size
  mm <- metrics[!metrics %in% names(evals)]
  mm.x <- rep(0, length(mm))
  names(mm.x) <- mm
  evals <- c(evals, mm.x)
  evals <- evals[order(names(evals))]
  return(c(Precision=hP, Recall=hR, Fmeasure=hF, evals, UndetectedRate=ND.rate))
}


evaluate <- function(TrueLabels, PredLabels, Indices = NULL, HierModPath=NULL){
  "
  This function was taken from https://github.com/tabdelaal/scRNAseq_Benchmark
  "
  #true_lab <- unlist(read.csv(TrueLabelsPath))
  #pred_lab <- unlist(read.csv(PredLabelsPath))
  true_lab <- TrueLabels
  pred_lab <- PredLabels

  true_lab <- FixLab(xstring = true_lab)
  pred_lab <- FixLab(xstring = pred_lab)

  if (! is.null(Indices)){
    true_lab <- true_lab[Indices]
    pred_lab <- pred_lab[Indices]
  }

  if(!is.null(HierModPath)){

    if(class(HierModPath)[1] == "RefMod"){
      refmod <- HierModPath
    }else{
      suppressPackageStartupMessages(library(HieRFIT))
      refmod <- LoadHieRMod(fileName=HierModPath)
    }

    hPRFtab <- hPRF(tpT = as.data.frame(cbind(true_lab, pred_lab)), tree = refmod@tree[[1]])

    }else{hPRFtab <- NULL}

  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))

  unique_all <- unique(c(unique_true,unique_pred))
  conf <- table(true_lab,pred_lab)
  pop_size <- rowSums(conf)

  pred_lab = gsub('Node..','Node',pred_lab)

  conf_F1 <- table(true_lab,pred_lab,exclude = c('Undetermined',
                                                 'unassigned',
                                                 'Unassigned',
                                                 'Unknown',
                                                 'rand',
                                                 'Node',
                                                 'Int.Node',
                                                 'ambiguous',
                                                 'unknown'))

  F1 <- vector()
  sum_acc <- 0

  for (i in c(1:length(unique_true))){
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(findLabel)){
      prec <- conf_F1[i,findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i,findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] <- (2*prec*rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i,findLabel]
    } else {
      F1[i] = 0
    }
  }

  pop_size <- pop_size[pop_size > 0]

  names(F1) <- names(pop_size)

  med_F1 <- median(F1)
  mean_F1 <- mean(F1)

  total <- length(pred_lab)
  num_unlab <- sum(pred_lab == 'Undetermined') +
    sum(pred_lab == 'unassigned') +
    sum(pred_lab == 'Unassigned') +
    sum(pred_lab == 'rand') +
    sum(pred_lab == 'Unknown') +
    sum(pred_lab == 'unknown') +
    sum(pred_lab == 'Node') +
    sum(pred_lab == 'Int.Node') +
    sum(pred_lab == 'ambiguous')
  per_unlab <- num_unlab / total

  num_Interlab <- sum(pred_lab == 'Node') +
    sum(pred_lab == 'Int.Node')
  per_Interlab <- num_Interlab / total

  acc <- sum_acc/sum(conf_F1)

  result <- list(Conf = conf, MeanF1=mean_F1, MedF1 = med_F1, F1 = F1, Acc = acc, PercInter= per_Interlab, PercUnl = per_unlab, PopSize = pop_size, hPRF=hPRFtab)

  return(result)
}
