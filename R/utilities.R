#Utility functions:

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
                              trainClassLabels <- ClassLabels[-flds[[1]]]
                              trainRef <- Ref[, -flds[[1]]]

                              refmod <- CreateHieR(RefData = trainRef,
                                                   ClassLabels = trainClassLabels,
                                                   TreeTable = TreeTable)
                              #Hierfit
                              testClassLables <- ClassLabels[flds[[1]]]
                              testRef <- Ref[, flds[[1]]]


                              testObj <- HieRFIT(Query = testRef, refMod = refmod, Prior = testClassLables)
                              PriorPostTable <- data.frame(Prior=testObj@Prior, Projection = testObj@Evaluation$Projection)
                              hPRF.out <- hPRF(tpT = PriorPostTable, tree = refmod@tree[[1]])
                              hPRF.out <- c(Tool="HieRFIT",
                                            CV_k=i,
                                            hPRF.out)

                              CMtab[[i]] <- PriorPostTable
                              print(hPRF.out)
                              cc <- t(hPRF.out)
                              return(cc)
                            }
                     )
  )
  return(list(metrics=hPRF_cv, Confusion=CMtab))
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

RandomizeR <- function(df, n=10){
  set.seed(192939)
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
      out <- "Correct_sibling_node"
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
  "Correct_sibling_node", "Correct_children_node", "Correct_grandchildren_node",
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
