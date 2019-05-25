
#' An internal function to create hierarchical classification models (hiemods) with random forest.
#' @param ExpData a Normalized expression data matrix, genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the ExpData matrix. Same length as colnames(ExpData).
#' @param tree A tree storing relatinship between the class labels. Default is null. required when mod.meth "hrf" is chosen.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
#' @keywords
#' @export
#' @usage
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

  return(hiemods)
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

DigestTree <- function(tree) {
  all.nodes <- unique(tree$edge[,1])
  return(all.nodes)
}

#' A function to create a tree in phylo format.
#' @param treeTable a data table/data.frame storing class relationships with intermediate class labels. Can be read from a tab separated file.
#' @usage treeTable <- read.delim("TaxonomyRank.tree.txt", header = F)
#' @usage tree <- CreateTree(treeTable)
CreateTree <- function(treeTable){
  #Create a tree:
  library(ape)
  library(data.tree)# '::' doesn't work!

  treeTable <- data.frame(lapply(treeTable, function(x) {gsub("\\+|-", ".", x)}))

  treeTable$pathString <- apply(cbind("TaxaRoot", treeTable), 1, paste0, collapse="/")
  tree <- as.phylo(as.Node(treeTable))
  return(tree)
}

CheckTreeIntegrity <- function(){}

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
