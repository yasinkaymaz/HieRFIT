#CreateHieR functions:

#' Reference class
#' @slot model A list of models to be used.
#' @slot tree A hierarchical tree representing class relationships.
#' @slot mlr multinominal logistic regression model for probability calibration.
#' @slot alhas list of alpha threshold of each class.
#' @slot species stores the species of the data, default is "hsapiens".
RefMod <- setClass(Class = "RefMod",
                   slots = c(model = "list",
                             tree = "list",
                             mlr = "list",
                             alphas = "list",
                             species = "character"))#Add treeTable

#' The main function for creating a reference model.
#' @param RefData Reference data from which class labels will be projected on Query data.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix (Ref). Same length as colnames(Ref).
#' @param Tree An input file to create a hierarchical tree for relationship between class labels. If stays null, a de novo tree is created.
#' @param species The organism from which the reference data generated. Allowed species are scientific species names, e.g. 'Homo sapiens'.
#' @param thread number of workers to be used for parallel processing. Default is Null, so serial processing.
#' @usage  pbmc.refmod <- CreateRef(RefData = as.matrix(pbmc@data), ClassLabels = pbmc@meta.data$ClusterNames_0.6, TreeFile = pbmc3k_tree)
CreateHieR <- function(RefData, ClassLabels, Tree=NULL, species = "hsapiens", thread=NULL, injectNoise=FALSE, binarize=FALSE, alpha=0.05, ...){
  options(warn=-1)
  if(missing(ClassLabels) || missing(RefData)){
    stop("Please, provide the required inputs! ClassLabels or RefData is missing.")
  }
  cat(paste("Species is", species, "\n"))
  cat("Preparing the reference data...\n")
  RefData <- as.matrix(RefData)
  #Downsample:
  samp.idx <- DownSampleIdx(RefData = RefData, ClassLabels = ClassLabels, ...)
  ClassLabels <- ClassLabels[samp.idx]
  RefData <- RefData[, samp.idx]

  RefData <- RefData[apply(RefData, 1, var) != 0, ]
  rownames(RefData) <- FixLab(xstring = rownames(RefData))
  ClassLabels <- FixLab(xstring = ClassLabels)

  if(binarize){
    RefData <- as.matrix((RefData > 0) + 0)
  }

  if(is.null(Tree)){
    tree <- CreateDeNovoTree(RefData, ClassLabels)
  }else{
    if(class(Tree) == "phylo"){tree <- Tree}else{tree <- CreateTree(treeTable = Tree)}
  }

  try(if(!identical(sort(tree$tip.label),
                   sort(unique(ClassLabels))))
    stop("Please, check class labels in the reference and tree file."))

  if(injectNoise == TRUE){
    noised.data <- NoiseInjecter(RefData = RefData, ClassLabels = ClassLabels, ...)
    RefData <- noised.data[["data"]]
    ClassLabels <- noised.data[["Y"]]
  }

  cat("Training model... This may take some time... Please, be patient!\n")
  model <- HieRandForest(RefData = RefData,
                         ClassLabels = ClassLabels,
                         tree = tree, thread = thread, ...)

  refObj <- new(Class = "RefMod",
                model = model,
                species = species)

  refObj@tree[[1]] <- tree
  #Get sigmoid functions for probability calibrations:
  refObj <- GetSigMods(RefData = RefData, ClassLabels = ClassLabels, refMod = refObj, ...)
  #Determine the alphas
  #refObj@alphas <- DetermineAlphas(refMod = refObj, RefData = RefData, Prior = ClassLabels, ...)
  refObj@alphas <- DetermineAlphaArray(refMod = refObj, RefData = RefData, Prior = ClassLabels, ...)
  #refObj@alphas <- alpha
  cat("Successfull!\n")
  #foreach::registerDoSEQ()
  return(refObj)
}


#' A function to determine alpha thresholds of each class in the tree.
#' @param refMod reference model. This will be used for self-projection.
#' @param RefData Reference data from which class labels will be projected on Query data (in this case self-projection).
#' @param Prior prior class labels if exist. For cross comparison. Should correspond row order of the Query.
DetermineAlphaArray <- function(refMod, RefData, Prior, ...){
  # cat("temp processing...")
  # RefData <- as.matrix(RefData)
  # samp.idx <- DownSampleIdx(RefData = RefData, ClassLabels = Prior, ...)
  # Prior <- Prior[samp.idx]
  # RefData <- RefData[, samp.idx]
  # RefData <- RefData[apply(RefData, 1, var) != 0, ]
  # rownames(RefData) <- FixLab(xstring = rownames(RefData))
  # Prior <- FixLab(xstring = Prior)
  cat("Now, determining the alpha tresholds...\n")
  Caldata <- NoiseInjecter(RefData = RefData,
                           ClassLabels = Prior,
                           refMod = refMod, ...)

  rownames(Caldata$data) <- FixLab(xstring = rownames(Caldata$data))
  #Fix this step when cross-species is the case
  HieMetObj <- CTTraverser(Query = Caldata$data, refMod = refMod)
  df <- ScoreEvaluate(ProbCert = HieMetObj@QueCers,
                      tree = refMod@tree[[1]],
                      full.tbl = TRUE)
  #dfs <- data.frame(Prior=Caldata$Y, df)
  alpha.list <- list()
  labs_l <- c(refMod@tree[[1]]$tip.label, refMod@tree[[1]]$node.label)
  for(class in labs_l[!labs_l %in% "TaxaRoot"]){
    if(!class %in% refMod@tree[[1]]$tip.label){
      Node <- match(class, labs_l)
      leaf <- GetChildNodeLeafs(tree = refMod@tree[[1]], node = Node)
      leaf <- unlist(leaf[which(lapply(leaf, length)>0)])
    }else{
      leaf <- class
    }
    alpha.list[[class]] <- colQuarts(df[which(Caldata$Y %in% leaf), ])
  }
  # pdf("CertScores_test.pdf",width = 12, height = 6)
  # for (class in unique(Prior)){
  #   df.sub <- as.data.frame(dfs)%>% filter(Prior==class) %>% select(-Prior)
  #   alpha.list[[class]] <- colQuarts(df.sub)
  #
  #   p <- as.data.frame(dfs)%>% filter(Prior==class) %>% select(-Prior) %>%
  #     reshape2::melt() %>% ggplot(aes(variable, value))+geom_boxplot()+
  #     theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  #     ggtitle(class)+
  #     geom_hline(yintercept=0.05, linetype="dashed", color = "red")
  #   print(p)
  #
  # }
  # dev.off()

  return(alpha.list)
}

#' A function to determine alpha thresholds of each class in the tree.
#' @param refMod reference model. This will be used for self-projection.
#' @param RefData Reference data from which class labels will be projected on Query data (in this case self-projection).
#' @param Prior prior class labels if exist. For cross comparison. Should correspond row order of the Query.
DetermineAlphas <- function(refMod, RefData, Prior, ...){
  cat("Now, determining the alpha tresholds...\n")
  Caldata <- NoiseInjecter(RefData = RefData,
                           ClassLabels = Prior,
                           refMod = refMod, ...)

  rownames(Caldata$data) <- FixLab(xstring = rownames(Caldata$data))
  #Fix this step when cross-species is the case
  HieMetObj <- CTTraverser(Query = Caldata$data, refMod = refMod)
  df <- ScoreEvaluate(ProbCert = HieMetObj@QueCers,
                      tree = refMod@tree[[1]],
                      full.tbl = TRUE)

  labs_l <- c(refMod@tree[[1]]$tip.label, refMod@tree[[1]]$node.label)
  alpha.list <- list()
  for(class in labs_l[!labs_l %in% "TaxaRoot"]){
    if(!class %in% refMod@tree[[1]]$tip.label){
      Node <- match(class, labs_l)
      leaf <- GetChildNodeLeafs(tree = refMod@tree[[1]], node = Node)
      leaf <- unlist(leaf[which(lapply(leaf, length)>0)])
    }else{
      leaf <- class
    }
    #cer.mean <- mean(df[which(Caldata$Y %in% leaf), class])
    #cer.sd <- sd(df[which(Caldata$Y %in% leaf), class])
    #cer.min <- min(df[which(Caldata$Y %in% leaf), class])
    cer.max <- max(df[which(Caldata$Y %in% leaf), class])

    data <- df[which(Caldata$Y %in% leaf), class]
    #data
    #lowerq = quantile(data)[2]
    #upperq = quantile(data)[4]
    #iqr = upperq - lowerq #inter quartile range

    #extreme.threshold.upper = (iqr * 3) + upperq
    #extreme.threshold.lower = lowerq - (iqr * 3)

    #data <- df[which(Caldata$Y %in% leaf), class]
    #lowerq = quantile(data)[2]
    #upperq = quantile(data)[4]

    data.not <- df[which(!Caldata$Y %in% leaf), class]
    #lowerq.not = quantile(data.not)[2]
    #upperq.not = quantile(data.not)[4]

    #iqr = upperq - lowerq #inter quartile range

    #extreme.threshold.upper = (iqr * 3) + upperq
    #extreme.threshold.lower = lowerq - (iqr * 3)

    df_10 <- rbind(data.frame(Class="1", Cert=data), data.frame(Class="0", Cert=data.not))
    df_10 <- df_10[order(df_10$Cert), ]

    C1 <- 0; C0 <- 0; Diff <- c()
    for(i in 1:length(df_10[,1]) ){
      if(df_10[i,]$Class == 1){
        C1 <- C1 +1
      }else{
        C0 <- C0 +1
      }
      Diff <- c(Diff, abs(C1 - C0))
    }
    df2 <- cbind(df_10, Diff=Diff)

    SmartCutoff <- as.numeric(head(df2[order(df2$Diff, decreasing = T),],1)['Cert'])

    #alpha.list[[class]]['Up.ext'] <- extreme.threshold.upper
    #alpha.list[[class]]['Low.ext'] <- extreme.threshold.lower
    alpha.list[[class]]['Low.ext'] <- SmartCutoff
    alpha.list[[class]]['Up.ext'] <- cer.max
    #alpha.list[[class]]['Low.ext'] <- cer.min
    }

  return(alpha.list)
}


#' A function to create a tree in phylo format.
#' @param treeTable a data table/data.frame storing class relationships with intermediate class labels. Can be read from a tab separated file.
#' @usage treeTable <- read.delim("TaxonomyRank.tree.txt", header = F)
#' @usage tree <- CreateTree(treeTable)
CreateTree <- function(treeTable){
  #Create a tree:
  library(ape)#TRY REPLACING
  library(data.tree)# '::' doesn't work!
  cat("Creating a tree with user table input...\n")
  treeTable <- data.frame(lapply(treeTable, function(x) {gsub("\\+|-|/", ".", x)}))

  treeTable$pathString <- apply(cbind("TaxaRoot", treeTable), 1, paste0, collapse="/")
  tree <- as.phylo(as.Node(treeTable))
  cat("Done!\n")

  return(tree)
}

#' A function to create a tree from the reference data based on euclidean distances between class types.
#' @param RefData Reference data from which class labels will be projected on Query data.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix (Ref). Same length as colnames(Ref).
CreateDeNovoTree <- function(RefData, ClassLabels){
  cat("Creating a de novo tree...\n")
  refData_mean <- NULL
  for(ct in unique(ClassLabels)){
    dt <- RefData[, ClassLabels == ct]
    ctmean <- rowMeans(as.matrix(dt))
    refData_mean <- cbind(refData_mean, ctmean)
  }
  colnames(refData_mean) <- unique(ClassLabels)

  cl_dists <- dist(t(refData_mean))
  clusts <- hclust(cl_dists)
  tree <- ape::as.phylo(clusts)
  tree$node.label <- c("TaxaRoot", paste("Int.Node", seq(tree$Nnode-1), sep = "."))
  cat("Done!\n")
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
HieRandForest <- function(RefData, ClassLabels, tree, thread, RPATH=NULL, ...){

  node.list <- DigestTree(tree = tree)
  hiemods <- vector("list", length = max(node.list))

  Rdata <- as.data.frame(t(RefData))
  colnames(Rdata) <- FixLab(xstring = colnames(Rdata))
  Rdata$ClassLabels <- factor(ClassLabels)

  pb <- txtProgressBar(min = 0, max = length(node.list), style = 3)
  p=1
  if(is.null(thread)){
    #Build a local classifier for each node in the tree. Binary or multi-class mixed.
    for(i in node.list){
      cat("\n")
      #print(paste("Training local node", tree$node.label[i - length(tree$tip.label)], sep=" "))
      #cat(paste("Training local node", tree$node.label[i - length(tree$tip.label)], sep=" "))
      setTxtProgressBar(pb, p)
      hiemods[[i]] <- NodeTrainer(Rdata = Rdata, tree = tree, node = i, ...)
      p=p+1
    } #closes the for loop.
    cat("\nDone!\n")

  }else{# thread is specified. For now, use this only when running on bigMem machines.
    library(doParallel)
    cl <- makePSOCKcluster(thread)
    registerDoParallel(cl)
    print(paste("registered cores is", getDoParWorkers(), sep = " "))
    out <- foreach(i=node.list, .packages = c('caret'), .inorder = TRUE, .export = ls(.GlobalEnv)) %dopar% {
      NodeTrainer(Rdata = Rdata, tree = tree, node = i, ...)
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

#' A function to train local nodes.
#' @param Rdata the transposed data for the node.
#' @param tree the tree input.
#' @param f_n the total number of features to be selected.
#' @param tree_n the total number of tree in each random forest.
#' @param min_n the cealing threshold for each class sample when training.
NodeTrainer <- function(Rdata, tree, node, f_n=200, tree_n=500, min_n = 500, ...){

  node.Data <- SubsetTData(Tdata = Rdata, tree = tree, node = node)
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
    if(dim(node.outData)[1] > min_n){
      node.outData <- node.outData[sample(rownames(node.outData), size = min_n), ]
    }
    node.outData <- droplevels(subset(node.outData, select=c(colnames(node.Data))))
  }else{
    #Create OurGroup for the TaxaRoot
    if(dim(Rdata)[1] > min_n){s_n <- min_n}else{s_n <- dim(Rdata)[1]}
    node.outData <- droplevels(Rdata[sample(rownames(Rdata), size = s_n), ])
    node.outData <- droplevels(subset(node.outData, select=c(colnames(node.Data))))
    node.outData <- droplevels(subset(node.outData, select=-c(ClassLabels)))
    node.outData <- Shuffler(df = node.outData)
  }
  node.outData$ClassLabels <- paste(labs_l[node], "OutGroup", sep="_")
  node.Data <- rbind(node.Data, node.outData)

  node.ClassLabels <- node.Data$ClassLabels
  node.Data <- droplevels(subset(node.Data, select=-c(ClassLabels)))

  #Downsize samples:
  include_idx <- NULL
  classes <- table(node.ClassLabels)

  for(type in names(classes)){
    type_idx <- which(node.ClassLabels == type)
    if(length(type_idx) > min_n){
      include_idx <- c(include_idx, sample(type_idx, min_n, replace = F))
    }else{
      include_idx <- c(include_idx, type_idx)
    }
  }
  node.ClassLabels <- node.ClassLabels[include_idx]
  node.Data <- node.Data[include_idx, ]

  node.Data <- node.Data[, apply(node.Data, 2, var) != 0]


  #First select the highly variable genes that correlate with the PCs
  P_dict <- FeaSelectPCA(Data = node.Data,
                         ClassLabels = node.ClassLabels,
                         num = 2000,
                         ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))
  #Then, select the genes as predictors if statistically DE between the classes.
  P_dict <- FeaSelectWilcox(Data = node.Data,
                            ClassLabels = node.ClassLabels,
                            num = f_n,
                            ...)
  node.Data <- droplevels(subset(node.Data, select=c(P_dict)))

  train.control <- caret::trainControl(method="oob",
                                       returnData = FALSE,
                                       savePredictions = "none",
                                       returnResamp = "none",
                                       allowParallel = TRUE,
                                       classProbs =  TRUE,
                                       trim = TRUE)

  node.mod <- caret::train(x = node.Data,
                           y = node.ClassLabels,
                           method = "rf",
                           norm.votes = TRUE,
                           importance = TRUE,
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
FeaSelectPCA <- function(Data, ClassLabels, PC_n=40, num=200, ...) {
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

#' A function for applying Wilcox test for candidate features.
#' @param Data an data matrix storing gene expression as genes in rows and samples in columns.
#' @param ClassLabels A list of class labels for cells/samples in the Data matrix. Same length as colnames(Data).
#' @param num the number of Predictors (genes) in total to be included. Default is 2000.
FeaSelectWilcox <- function(Data, ClassLabels, num = 200, ...) {
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
    fea.p <- fea.p[!rownames(fea.p) %in% pick.fea, ]
    #select top N variables. N is proportional to each class.
    pick.fea <- c(pick.fea, rownames(head(fea.p, round(num/length(cls)))))
  }
  return(pick.fea)
}

#' A function to slide data according to the class hierarchies. And updates the class labels.
#' @param Tdata transposed data.
#' @param tree tree input.
#' @param node a particular node of the tree.
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

#' A function for generating the sigmoid functions.
#' @param RefData input reference data.
#' @param ClassLabels input class labels.
#' @param refMod reference model.
GetSigMods <- function(RefData, ClassLabels, refMod, ...){
  cat("Now calibrating the nodes...\n")
  #Cal.data <- NoiseInject(RefData = RefData, ClassLabels = ClassLabels, refMod = refMod)
  Cal.data <- NoiseInjecter(RefData = RefData, ClassLabels = ClassLabels, refMod = refMod)
  Pvotes <- ProbTraverser(Query = Cal.data[["data"]], refMod = refMod, ...)
  labs_l <- c(refMod@tree[[1]]$tip.label, refMod@tree[[1]]$node.label)
  votes <- data.frame(Prior=Cal.data[["Y"]], Pvotes)
  node.list <- DigestTree(tree = refMod@tree[[1]])
  for(i in node.list){
    nodeModel <- refMod@model[[as.character(i)]]
    outGname <- grep("OutGroup", nodeModel$finalModel$classes, value=T)
    Leafs <- NULL
    classes <- nodeModel$finalModel$classes[! nodeModel$finalModel$classes %in% outGname]
    node.dict <- list()
    Labs <- votes$Prior
    for(l in classes){
      if(l %in% refMod@tree[[1]]$tip.label){
        node.dict[[l]] <- l
      }else{
        leaf <- GetChildNodeLeafs(tree = refMod@tree[[1]], node = match(l, labs_l))
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
    dfr.full <- data.frame(Prior=Labs, votes[, nodeModel$finalModel$classes])
    mlrmod <- nnet::multinom(Prior ~ ., data = dfr.full, model=FALSE, trace=FALSE)
    #refMod@model[[as.character(i)]]$mlr <- stripGlmModel(mlrmod)
    #refMod@model[[as.character(i)]]$mlr <- mlrmod
    #refMod@mlr[[as.character(i)]] <- stripGlmModel(mlrmod)
    refMod@mlr[[as.character(i)]] <- mlrmod
  }
  cat("Done!\n")
  return(refMod)
}


#' A function for generating the sigmoid functions.
#' @param RefData input reference data.
#' @param ClassLabels input class labels.
#' @param NoiseRate the rate at which selected features are noised by setting their values to zero.
#' @param refMod reference model to inject noise to only selected features. If Null, the entire data is used to inject noise.
NoiseInjecter <- function(RefData, ClassLabels, NoiseRate=0.75, refMod=NULL, ...){
  NoisedRefData <- list()

  features <- NULL
  if(is.null(refMod)){
    features <- rownames(RefData)
  }else{
    for(i in names(refMod@model)){ features <- append(features, refMod@model[[i]]$finalModel$xNames)}
  }
  features <- unique(features)
  RefData <- RefData[features, ]

  #by default, rate is 75%
  #Indices of non-zero counts:
  NonZeroInds <- which(RefData!=0, arr.ind = T)
  #Total number of non-zeros
  nZ_c <- dim(NonZeroInds)[1]
  #Indices of values in the count matrix to be set zero
  noise_idx <- NonZeroInds[sample(1:nZ_c, nZ_c*NoiseRate),]
  #set the selected counts to zero
  X <- RefData
  X[noise_idx] <- 0
  #extract the noised samples
  X <- X[, unique(noise_idx[, 'col'])]
  Y <- ClassLabels[unique(noise_idx[, 'col'])]
  #Create the taxaout data
  min_n <- 500
  if(dim(RefData)[2] > min_n){s_n <- min_n}else{s_n <- dim(RefData)[2]}
  taxa.outData <- RefData[, sample(colnames(RefData), size = s_n)]
  taxa.outData <- Shuffler(df = taxa.outData)
  #combine them all
  NoisedRefData[["data"]] <- cbind(RefData, X, taxa.outData)
  NoisedRefData[["Y"]] <- c(ClassLabels, Y, rep("TaxaRoot_OutGroup", s_n))
  colnames(NoisedRefData[["data"]]) <- make.names(colnames(NoisedRefData[["data"]]), unique = TRUE)

  return(NoisedRefData)
}
