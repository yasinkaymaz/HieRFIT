#Utility functions:

ParApp <- function() {
  '"
  Sets the parameters used throughout the run.
  "'
  PCs <- 20 # number of principle components to be analyzed.
  Pred_n <- 1000 # number of predictors aimed to be collected.
  prefix <- ""
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

#' A function to fix class labels.
#' @param xstring is a list of class labels in character.
FixLab <- function(xstring){
  #Replace white space with '_'
  xstring <- gsub(xstring, pattern = " ", replacement = "_")
  xstring <- gsub(xstring, pattern = "\\+|-", replacement = ".")
  return(xstring)
}

SubsetTData <- function(Tdata, tree, node, transpose=FALSE){
  #This function subsets a Tdata with class labels (trainingData - output of DataReshaper)
  #And updates the class labels.
  # 1. Extract the data under the node. Subsample if necessary.
  childNodes <- GetChildNodeLeafs(tree = tree, node = node)

  SubTdata <- NULL
  #loop through child nodes that are not null in the list.
  for (i in which(lapply(childNodes, length) > 0)){

    Subdata <- droplevels(Tdata[which(Tdata$ClassLabels %in% childNodes[i][[1]]), ])
    print(childNodes[i][[1]])
    print(class(childNodes[i][[1]]))
    print(Subdata[1:4,1:4])
    if (i > length(tree$tip.label)){# if the node is not a leaf node, then
      if(!is.null(tree$node.label)){# if labels for subnodes exist
        labels <- c(tree$tip.label, tree$node.label)
        #Replace labels with subnode labels.
        print(i)
        print(labels[i])
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
  #Here balance the class sizes with median:
  #SubTdata <- DownSampleRef(RefData = SubTdata, min_n = round(median(table(SubTdata$ClassLabels))))
  if(transpose){
    return(t(SubTdata))
  }else{
    return(SubTdata)
  }

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

CVRunner <- function(Ref, ClassLabels, TreeTable=NULL, cv_k=5, method="hrf"){

  # Create K-fold data split:
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
                              hPRF.out <- c(CV_k=i,
                                            hPRF.out,
                                            Undetermined=nrow(PriorPostTable[PriorPostTable$Projection == "Undetermined",])/length(PriorPostTable[,1]))
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

RandomizeR <- function(df){
  #set.seed(192939)
  dfR <- df[sample(nrow(df)), sample(ncol(df))]
  rownames(dfR) <- rownames(df)
  colnames(dfR) <- colnames(df)
  return(dfR)
}

#' Homology mapping via orhologous genes between mouse and rat.
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
IntSectSize <- function(t, p, tree){
  Ti <- GetAncestPath(tree = tree, class = FixLab(t))
  Pi <- GetAncestPath(tree = tree, class = FixLab(p))
  intL <- length(intersect(Ti, Pi))
  return(intL)
}

#' A function for Hierarchical Precision, Recall, and F-measure.
#' @param tpT PriorPostTable: a table with two columns of which first is Prior and second is Post-prediction.
#' @param tree tree topology in phylo format.
hPRF <- function(tpT, tree, BetaSq=1){
  # To Do: Consider Undetermined class!
  tpT$Int <- apply(tpT, 1, function(x) IntSectSize(t = x[1], p = x[2], tree = tree))
  tpT$PiL <- apply(tpT, 1, function(x) length(GetAncestPath(tree = tree, class = FixLab(x[2]) )))
  tpT$TiL <- apply(tpT, 1, function(x) length(GetAncestPath(tree = tree, class = FixLab(x[1]) )))

  hP <- sum(tpT$Int)/sum(tpT$PiL)
  hR <- sum(tpT$Int)/sum(tpT$TiL)
  hF <- (BetaSq+1)*hP*hR/(BetaSq*hP+hR)

  return(c(Precision=hP, Recall=hR, Fmeasure=hF))
}

