#CreateHieR functions:

#' Reference class
#' @slot model A list of models to be used.
#' @slot tree A hierarchical tree representing class relationships.
#' @slot modtype type of the model method used to create object.
RefMod <- setClass(Class = "RefMod",
                   slots = c(model = "list",
                             tree = "list",
                             modtype = "character"))#Add treeTable

#' The main function for creating a reference model.
#' @param RefData Reference data from which class labels will be projected on Query data.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix (Ref). Same length as colnames(Ref).
#' @param method The model training method, "rf" for random forest, "svmLinear" for support vector machine, "hrf" for hierarchical random forest. Default is "hrf"
#' @param TreeFile An input file to create a hierarchical tree for relationship between class labels. Default is null but required if method 'hrf' is chosen.
#' @usage  pbmc.refmod <- CreateRef(RefData = as.matrix(pbmc@data), ClassLabels = pbmc@meta.data$ClusterNames_0.6, TreeFile = pbmc3k_tree)
CreateHieR <- function(RefData, ClassLabels, TreeTable=NULL, method="hrf", thread=NULL, ...){

  if(method == "hrf"){
    tree <- CreateTree(treeTable = TreeTable)
  }else{ tree <- NULL}

  ClassLabels <- FixLab(xstring = ClassLabels)

  try(if(missing(tree)|| missing(ClassLabels) || missing(RefData))
    stop("Please, provide the required inputs!"))

    print("Training model... This may take some time... Please, be patient!")

  try(if(!identical(sort(tree$tip.label),
                   sort(unique(ClassLabels))))
    stop("Please, check class labels in the reference and tree file."))

  model <- HieRandForest(RefData = RefData,
                         ClassLabels = ClassLabels,
                         tree, thread = thread, ...)
  refObj <- new(Class = "RefMod",
                model = model,
                modtype = method)
  refObj@tree[[1]] <- tree

  #foreach::registerDoSEQ()
  return(refObj)
}

#' A function to create a tree in phylo format.
#' @param treeTable a data table/data.frame storing class relationships with intermediate class labels. Can be read from a tab separated file.
#' @usage treeTable <- read.delim("TaxonomyRank.tree.txt", header = F)
#' @usage tree <- CreateTree(treeTable)
CreateTree <- function(treeTable){
  #Create a tree:
  library(ape)#TRY REPLACING
  library(data.tree)# '::' doesn't work!

  treeTable <- data.frame(lapply(treeTable, function(x) {gsub("\\+|-|/", ".", x)}))

  treeTable$pathString <- apply(cbind("TaxaRoot", treeTable), 1, paste0, collapse="/")
  tree <- as.phylo(as.Node(treeTable))
  return(tree)
}

#' An internal function to create hierarchical classification models (hiemods) with random forest.
#' @param RefData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param tree A tree storing relatinship between the class labels. Default is null. required when mod.meth "hrf" is chosen.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
#' @keywords
#' @export
#' @usage
HieRandForest <- function(RefData, ClassLabels, tree, thread=3, ...){

  node.list <- DigestTree(tree = tree)
  hiemods <- vector("list", length = max(node.list))
  var_cutoff=0.1

  Rdata <- RefData[which(apply(RefData, 1, var) > var_cutoff),]
  Rdata <- as.data.frame(t(Rdata))
  colnames(Rdata) <- FixLab(xstring = colnames(Rdata))
  Rdata$ClassLabels <- factor(make.names(ClassLabels))

  pb <- txtProgressBar(min = 0, max = length(node.list), style = 3)
  p=1
  if(is.null(thread)){
    #Build a local classifier for each node in the tree. Binary or multi-class mixed.
    for(i in node.list){

      print(paste("Training local node", tree$node.label[i-length(tree$tip.label)], sep=" "))
      setTxtProgressBar(pb, p)

      #hiemods[[i]] <- NodeTrainer(Rdata = Rdata, tree = tree, node = i, ...)
      hiemods[[i]] <- NodeTrainer2(Rdata = Rdata, tree = tree, node = i, ...)
      p=p+1
    } #closes the for loop.

  }else{# thread is specified. For now, use this only when running on bigMem machines.
    library(doParallel)
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths("/n/home13/yasinkaymaz/biotools/Rlibs/"))
    print(paste("registered cores is", getDoParWorkers(), sep = " "))
    out <- foreach(i=node.list, .packages = c('caret'), .inorder = TRUE, .export = ls(.GlobalEnv)) %dopar% {
      NodeTrainer(Rdata = Rdata, tree = tree, node = i, ...)
    }
    stopCluster(cl)

    for(x in 1:length(out)){
      hiemods[node.list[x]] <- out[x]
    }
  }
  FlatTreeTable <- data.frame(V1=names(table(ClassLabels)))
  tree <- CreateTree(treeTable = FlatTreeTable)
  fi <- DigestTree(tree = tree)
  #i is the last node in the node.list above:
  hiemods[[i+1]] <- NodeTrainer(Rdata = Rdata, tree = tree, node = fi)
  names(hiemods) <- seq_along(hiemods)
  hiemods[sapply(hiemods, is.null)] <- NULL

  return(hiemods)
}

#' A function to create a local classifier for a given node.
#' @param Rdata
#' @param tree
#' @param node
#' @param f_n number of features to be included in local classifier.
NodeTrainer <- function(Rdata, tree, node, f_n=200, tree_n=500, ...){

  node.Data <- SubsetTData2(Tdata = Rdata, tree = tree, node = node)
  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))
  node.Data <- node.Data[, apply(node.Data, 2, var) != 0]
  #First select the highly variable genes that correlate with the PCs
  P_dict <- FeatureSelector(Data = node.Data,
                             ClassLabels = node.ClassLabels,
                             num = 2000,
                             ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  #Then, select the genes as predictors if statistically DE between the classes.
  P_dict <- FeatureSelector2(Data = node.Data,
                            ClassLabels = node.ClassLabels,
                            num = f_n,
                            ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  node.Data$ClassLabels <- node.ClassLabels

  #Generate outgroup data for the node:
  #1. check the node: is.internal? is.not.Root?
  labs_l <- c(tree$tip.label, tree$node.label)

  if(labs_l[node] != "TaxaRoot"){
  parent.t <- tree$edge[which(x = tree$edge[, 2] == node), ][1]
  for(n in 1:length(labs_l)){
    if(labs_l[n] != "TaxaRoot"){
      parent.p <- tree$edge[which(x = tree$edge[, 2] == n), ][1]
      if(parent.t == parent.p & node != n){
        print(paste("Training for", labs_l[node], "and its sibling is:", labs_l[n]))
        #SubsetTData
        node.outData <- SubsetTData2(Tdata = Rdata, tree = tree, node = n)
        node.outData <- droplevels(subset(node.outData, select=c(P_dict)))
        node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
        if(dim(node.outData)[1] > 500){
          node.outData <- node.outData[sample(rownames(node.outData), size = 500), ]
        }
        node.Data <- rbind(node.Data, node.outData)
       }
     }
   }
  }
  #--#
  train.control <- caret::trainControl(method="oob",
                                       returnData = FALSE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim = TRUE)
  node.mod <- caret::train(ClassLabels~.,
                           data = node.Data,
                           method = "rf",
                           norm.votes = TRUE,
                           importance = FALSE,
                           proximity = FALSE,
                           outscale = FALSE,
                           preProcess = c("center", "scale"),
                           ntree = tree_n,
                           trControl = train.control, ...)

  return(node.mod)
}

NodeTrainer2 <- function(Rdata, tree, node, f_n=200, tree_n=500, ...){

  node.Data <- SubsetTData(Tdata = Rdata, tree = tree, node = node)
  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))
  node.Data <- node.Data[, apply(node.Data, 2, var) != 0]
  #First select the highly variable genes that correlate with the PCs
  P_dict <- FeatureSelector(Data = node.Data,
                             ClassLabels = node.ClassLabels,
                             num = 2000,
                             ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  #Then, select the genes as predictors if statistically DE between the classes.
  P_dict <- FeatureSelector2(Data = node.Data,
                            ClassLabels = node.ClassLabels,
                            num = f_n,
                            ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  node.Data$ClassLabels <- node.ClassLabels

  #Generate outgroup data for the node:
  #1. check the node: is.not.Root?
  labs_l <- c(tree$tip.label, tree$node.label)
  if(labs_l[node] != "TaxaRoot"){
    childNodes <- GetChildNodeLeafs(tree = tree, node = node)
    childLeafs <- NULL
    for (i in which(lapply(childNodes, length) > 0)){
      childLeafs <- c(childLeafs, childNodes[i][[1]])
    }
    outGroupLeafs <- tree$tip.label[!tree$tip.label %in% childLeafs]
    node.outData <- droplevels(Rdata[which(Rdata$ClassLabels %in% outGroupLeafs), ])
    if(dim(node.outData)[1] > 500){
      node.outData <- node.outData[sample(rownames(node.outData), size = 500), ]#500 can be replaced with a better #
    }
    node.outData <- droplevels(subset(node.outData, select=c(P_dict)))
    node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
    node.Data <- rbind(node.Data, node.outData)
  }
  #--#
  train.control <- caret::trainControl(method="oob",
                                       returnData = FALSE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim = TRUE)
  node.mod <- caret::train(ClassLabels~.,
                           data = node.Data,
                           method = "rf",
                           norm.votes = TRUE,
                           importance = FALSE,
                           proximity = FALSE,
                           outscale = FALSE,
                           preProcess = c("center", "scale"),
                           ntree = tree_n,
                           trControl = train.control, ...)

  return(node.mod)
}

#' A function used internally for selecting genes based on their weight in the top principle components.
#' @description p: pca object, pcs: number of PCs to include, num: number of genes from top and bottom.
#' Weighted gene picking depending on PC number: Initial PCs give more genes.
#' For example; the number of genes are taken from PC i is calculated by (pcs-i+1)*(num*2)/(pcs*(pcs+1))
#' @param Data an data matrix storing gene expression as genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param PC_n the number of PCs to be looked at when selecting genes. Default is 40.
#' @param num the number of Predictors (genes) in total to be included. Default is 2000.
#' @param ... parameters to be passed down to subfunctions such as f_n, tree_n, and PC_n,
#' for "number of features per local classifier", "number of trees per local classsifier",
#' and "number of PC space for feature search", respectively.
#' @keywords PCA loadings selection
#' @export
#' @usage Predictors <- FeatureSelector(Data = as.matrix(SeuratObject@data), PCs = 10, num = 2000)
FeatureSelector <- function(Data, ClassLabels, PC_n = 40, num = 200, ...) {
  #do pca
  pcatrain <- prcomp(Data, center = TRUE, scale=TRUE, rank. = PC_n)
  pcadata <- data.frame(pcatrain$x, ClassLabels = ClassLabels)
  cls <- levels(ClassLabels)

  ptab <- NULL
  for(i in 1:PC_n){
    PC.stats <- NULL
    for(c in cls){
      p.cl <- t.test(pcadata[pcadata$ClassLabels == c, i],
                     pcadata[pcadata$ClassLabels != c, i])$p.value
      PC.stats <- c(PC.stats, p.cl)
    }
    names(PC.stats) <- cls
    pc.col <- paste("PC", i, sep = "")
    ptab <- cbind(ptab, pc.col=PC.stats)
  }

  ptab <- as.data.frame(ptab)
  colnames(ptab) <- colnames(pcadata)[-length(pcadata[1,])]
  ptab <- ptab*length(cls)*PC_n#Correct for the multiple test. Bonferroni.
  #Select only PCs which are significanly separating at least one of the class.
  PCs.sig <- colnames(ptab[, apply(ptab < 0.05, 2 ,any)])

  if(length(PCs.sig) == 0){PCs.sig <- paste(rep("PC", 3), 1:3, sep="")}

  #Calculate the variance proportions explained by each selected PC.
  vars <- apply(pcatrain$x, 2, var)
  props <- vars[PCs.sig] / sum(vars[PCs.sig])
  props <- round(props*num)

  pick.rot <- NULL
  for(i in PCs.sig){ #For each principle component
    #sort variables based on their absolute component rotations values
    pc.rot <- pcatrain$rotation[order(abs(pcatrain$rotation[, i]), decreasing = TRUE), i]
    #exclude already selected variables.
    pc.rot <- pc.rot[!names(pc.rot) %in% pick.rot]
    #select top N variables. N is proportional to variance explained by overall PC i.
    pick.rot <- c(pick.rot, names(head(pc.rot, props[i])))
  }
  return(pick.rot)
}


FeatureSelector2 <- function(Data, ClassLabels, num = 200, ...) {
  cls <- levels(ClassLabels)
  #Data.cl <- data.frame(Data, ClassLabels = ClassLabels)# If you do this way, it cannot find genes, gives error. Weird!
  Data.cl <- Data
  Data.cl$ClassLabels <- ClassLabels
  ptab <- NULL
  for(i in names(Data)){
    fea.stats <- NULL
    if(length(cls) > 0){# Do this if only more than 2 classes exist. Redundant
      for(c in cls){
        p.cl <- wilcox.test(Data.cl[Data.cl$ClassLabels == c, i],
                            Data.cl[Data.cl$ClassLabels != c, i])$p.value
        fea.stats <- c(fea.stats, p.cl)
      }
      names(fea.stats) <- cls
    }else{
      c <- cls[1]
      p.cl <- wilcox.test(Data.cl[Data.cl$ClassLabels == c, i],
                          Data.cl[Data.cl$ClassLabels != c, i])$p.value
      fea.stats <- c(fea.stats, p.cl)
      names(fea.stats) <- c
    }
    fea.stats <- as.data.frame(t(fea.stats))
    rownames(fea.stats) <- i
    ptab <- rbind(ptab, fea.stats)
  }

  pick.fea <- NULL
  for(c in names(ptab)){
    fea.p <- ptab[order(abs(ptab[, c]), decreasing = FALSE), ]
    #exclude already selected variables.
    fea.p <- fea.p[!rownames(fea.p) %in% pick.fea,]
    #select top N variables. N is proportional to each class.
    pick.fea <- c(pick.fea, rownames(head(fea.p, round(num/length(cls)))))
  }
  return(pick.fea)
}


#' A function to slide data according to the class hierarchies. And updates the class labels.
#' @param Tdata
#' @param tree
#' @param node
#' @param transpose
SubsetTData <- function(Tdata, tree, node){
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
  return(SubTdata)
}

SubsetTData2 <- function(Tdata, tree, node){

  labs_l <- c(tree$tip.label, tree$node.label)

  if(labs_l[node] %in% tree$tip.label){
    childNodes <- list(labs_l[node])
    SubTdata <- droplevels(Tdata[which(Tdata$ClassLabels %in% childNodes[[1]]), ])

  }else{

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
 }
  return(SubTdata)
}
