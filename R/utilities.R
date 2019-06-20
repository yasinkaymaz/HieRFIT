

#' A function to calculate Asymmetric Entropy
#' @param p measured probability array
#' @param w empirically calculated W set
#' @return U set of certainty values for each probility outcome in p given w.
GetCertaintyArray <- function(p, w){
  U <- numeric()
  for(i in 1:length(p)){
    w[i] <- w[i]+1e-10#to prevent math err.
    if(p[i] > w[i]){
      Beta=1
    }else{
        Beta=-1
        }
    #H <- c(H, ((p[i]*(1-p[i])) / (p[i]-2*w[i]*p[i]+w[i]^2)))
    #U <- c(U, (Beta*(p[i]-w[i])^2)/( ((1-2*w[i])*p[i])+w[i]^2 ) )
    U <- c(U, Beta*( 1 - ( (p[i]*(1-p[i])) / (p[i]-2*w[i]*p[i]+w[i]^2)  ) )  )
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

                              refmod <- CreateHieR(Ref = trainRef,
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
      print(RefSeuObj@meta.data[testcells, ]$ColClassLabs)
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

#' A wrapper function for quick data processing with Seurat functions
#' This function allows to determine highly variable genes and scale their expression,
#' run PCA, tSNE, and cluster detection.
#' @param SeuratObj a Seurat S4 object.
#' @param scale.only.var a boolean variable to determine whether to scale entire data or only highly variable genes. Default is True.
#' @param PCs number of PCs to be included in the downstream analysis
#' @param vars2reg variables to be regressed out when scaling the expression data.
#' @param perp perplexity parameter to be passed to RunTSNE function from Seurat.
#' @keywords seurat quick
#' @export
#' @usage pbmc <- QuickSeurat(pbmc, scale.only.var=F, PCs=5, perp=20)
QuickSeurat <- function(SeuratObj, scale.only.var=T, PCs=20, perp=30, vars2reg) {

  SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)
  hv.genes <- head(rownames(SeuratObj@hvg.info), 1000)
  if (scale.only.var == TRUE) {
    if (!missing(vars2reg)) {
      SeuratObj <- ScaleData(SeuratObj, vars.to.regress = vars2reg, genes.use = hv.genes, do.par=T, num.cores = 8)
    }else{
      SeuratObj <- ScaleData(SeuratObj, genes.use = hv.genes, do.par=T, num.cores = 8)
    }
  }else{
    if (!missing(vars2reg)) {
      SeuratObj <- ScaleData(SeuratObj, vars.to.regress = vars2reg, do.par=T, num.cores = 8)
    }else{
      SeuratObj <- ScaleData(SeuratObj, do.par=T, num.cores = 8)
    }
  }
  SeuratObj <- RunPCA(SeuratObj, pc.genes = hv.genes, do.print = FALSE, pcs.compute=PCs)
  SeuratObj <- FindClusters(SeuratObj, reduction.type = "pca", dims.use = 1:PCs, resolution = 1, print.output = FALSE, save.SNN = TRUE, force.recalc = T)
  SeuratObj <- RunTSNE(SeuratObj, dims.use = 1:PCs, do.fast = TRUE, check_duplicates = FALSE)

  return(SeuratObj)

}

#' A function to downsample Seurat object based on cell type identity
#' @param SeuratObj a Seurat S4 object.
#' @param IdentityCol the column 'number' in the metadata slot showing the cell type identities.
#' @usage pbmc1 <- DownSizeSeurat(SeuratObj = pbmc, IdentityCol = 7)
DownSizeSeurat <- function(SeuratObj, IdentityCol, min_n=NULL){
  cells <- NULL
  classes <- table(SeuratObj@meta.data[,IdentityCol])
  print(classes)
  if(is.null(min_n)){
    min_n <- min(classes)
    print(min_n)
  }
  for(type in names(classes)){
    if( classes[type] > min_n ){
      cells <- c(cells, sample(rownames(SeuratObj@meta.data[which(SeuratObj@meta.data[,IdentityCol] == type), ]), size = min_n, replace = F))
    }else{
      cells <- c(cells, sample(rownames(SeuratObj@meta.data[which(SeuratObj@meta.data[,IdentityCol] == type), ]), size = classes[type], replace = F))
    }
  }
  downSobj <- SubsetData(object = SeuratObj, cells.use = cells, do.clean=T)
  return(downSobj)
}

#' A function to downsample a refdata table based on classLabels
#' @param RefData a data table with features as columns (last column being ClassLabels), instances in the rows.
#' @param IdentityCol the name of the column in the refdata storing class labels. default is "ClassLabels"
#' @param min_n min number of samples to downsample each class. default is the size of the minority class.
#' @usage RefData_d <- DownSampleRef(RefData = RefData)
DownSampleRef <- function(RefData, IdentityCol="ClassLabels", min_n=NULL){
  samples <- NULL
  classes <- table(RefData[, IdentityCol])
  print(classes)
  if(is.null(min_n)){
    min_n <- min(classes)
    print(min_n)
  }
  for(type in names(classes)){
    if( classes[type] > min_n ){
      samples <- c(samples, sample(rownames(RefData[which(RefData[, IdentityCol] == type), ]), size = min_n, replace = F))
    }else{
      samples <- c(samples, sample(rownames(RefData[which(RefData[, IdentityCol] == type), ]), size = min_n, replace = T))
    }
  }
  RefData_d <- RefData[samples, ]
  return(RefData_d)
}
SeuratWrapper <- function(ExpData, ProjectLabel, NewMeta, Normalize=T, suppressLog=F, scale.only.var=T, PCs=20, perp=30, dump.files=F, min.cells=0, min.genes=0) {

  if (Normalize == TRUE) {print("Assuming the input is in count ...")
  }else{
    print("Assuming the input is in TPM ...")
    if (suppressLog == TRUE) {
      print("not taking log ...")
    }else{
      ExpData <- log1p(ExpData)
    }
  }

  SeuratObj <- CreateSeuratObject(raw.data = ExpData, project = ProjectLabel, min.cells=min.cells, min.genes = min.genes)

  if (Normalize == TRUE) {
    SeuratObj <- NormalizeData(object = SeuratObj)
  }else{
    print("Not normalizing the data since TPM is assumed ... ")
  }

  SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)
  hv.genes <- head(rownames(SeuratObj@hvg.info), 1000)

  if (scale.only.var == TRUE) {
    SeuratObj <- ScaleData(SeuratObj, genes.use = hv.genes, do.par=T, num.cores = 8)
  }else{
    SeuratObj <- ScaleData(SeuratObj, do.par=T, num.cores = 8)
  }

  if (!missing(NewMeta)) {
    SeuratObj <- AddMetaData(SeuratObj, NewMeta[rownames(SeuratObj@meta.data), ])
  }else{
    print("No new meta file is provided. Skipping...")
  }

  SeuratObj <- RunPCA(SeuratObj, pc.genes = hv.genes, do.print = FALSE, pcs.compute=PCs)
  SeuratObj <- FindClusters(SeuratObj, reduction.type = "pca", dims.use = 1:PCs, resolution = 1, print.output = FALSE, save.SNN = TRUE, force.recalc = T)
  SeuratObj <- RunTSNE(SeuratObj, dims.use = 1:PCs, do.fast = TRUE,check_duplicates = FALSE, perplexity=perp)

  if (dump.files == T) {
    #Export the tSNE coordinates along with the Cluster assignment IDs
    rownames_to_column(as.data.frame(SeuratObj@dr$tsne@cell.embeddings))  %>%
      as.tibble() %>%
      add_column(Clusters=SeuratObj@meta.data$res.1) %>%
      dplyr::rename(Cellname = rowname) %>%
      as_data_frame() %>% write_csv(paste(ProjectLabel, "_tSNECoordinates_Clusters.csv", sep=""))

    #Export Normalized and Scaled Expression matrix for cells and genes in the analysis
    rownames_to_column(as.data.frame(as.matrix(SeuratObj@data))) %>%
      dplyr::rename(GeneName = rowname) %>%
      as_data_frame() %>%
      write_delim(paste(ProjectLabel, "_Normalized_Expression_matrix.txt", sep=""))

    #Save seurat object:
    save(SeuratObj, file=paste(ProjectLabel, ".seurat.Robj", sep=""))
  }

  return(SeuratObj)
}

SeuratCCAmerger <- function(listofObjects) {
  # Determine genes to use for CCA, must be highly variable in at least 2 datasets
  #ob.list <- list(zeisel, romanov, tasic, marques)
  ob.list <- listofObjects
  genesuse <- c()
  ids=NULL
  for (i in 1:length(ob.list)) {
    genesuse <- c(genesuse, head(rownames(ob.list[[i]]@hvg.info), 1000))
    ob.list[[i]]@meta.data$dataSource <- paste("id", i, sep="")
    ids <- c(ids, paste("id", i, sep=""))
  }
  genesuse <- names(which(table(genesuse) > 1))
  for (i in 1:length(ob.list)) {
    genesuse <- genesuse[genesuse %in% rownames(ob.list[[i]]@scale.data)]
  }

  if (length(ob.list) > 2) {
    # Run multi-set CCA
    integrated <- RunMultiCCA(ob.list, genes.use = genesuse, num.ccs = 15, add.cell.ids = ids)
    # Run rare non-overlapping filtering
    integrated <- CalcVarExpRatio(object = integrated, reduction.type = "pca", dims.use = 1:10, grouping.var = "dataSource")
    integrated <- SubsetData(integrated, subset.name = "var.ratio.pca", accept.low = 0.5)
  }else{
    #integrated <- RunCCA(object = ob.list[[1]], object2 = ob.list[[2]], genes.use = genesuse, num.cc = 15, add.cell.id = ids)
    integrated <- RunCCA(object = ob.list[[1]], object2 = ob.list[[2]], genes.use = genesuse, num.cc = 15)
  }
  # Alignment
  integrated <- AlignSubspace(integrated, reduction.type = "cca", dims.align = 1:10, grouping.var = "dataSource")
  # t-SNE and Clustering
  integrated <- FindClusters(integrated, reduction.type = "cca.aligned", dims.use = 1:10, save.SNN = T, resolution = 0.4)
  integrated <- RunTSNE(integrated, reduction.use = "cca.aligned", dims.use = 1:10)
  save(integrated, file="integrated.Aligned.seurat.Robj")
  return(integrated)
}

FlatRF <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions") {

  library(caret)
  library(randomForest)
  library(tidyverse)

  if (!missing(SeuratObject)) {
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{
    print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
  }#Closes missing(SeuratObj)
  colnames(testExpSet) <- make.names(colnames(testExpSet))
  #Prepare Test Expression set
  testsub <- testExpSet[,which(colnames(testExpSet) %in% model$finalModel$xNames)]
  missingGenes <- model$finalModel$xNames[which(!model$finalModel$xNames %in% colnames(testExpSet))]
  print(model$finalModel$importance[missingGenes,])

  missingGenes.df <- data.frame(matrix(0, ncol = length(missingGenes), nrow = length(rownames(testExpSet))))
  colnames(missingGenes.df) <- missingGenes
  TestData <- cbind(testsub, missingGenes.df)
  TestData <- TestData[,model$finalModel$xNames]
  mmDGm <- mean(model$finalModel$importance[missingGenes,])
  mmDGf <- mean(model$finalModel$importance[which(!model$finalModel$xNames %in% missingGenes),])

  cat("Number of Features (genes) to be considered is", length(colnames(testsub)), '\n',
      "Number of missing Features set to zero is", length(missingGenes), '\n',
      "Mean MeanDecreaseGini of missing genes is", mmDGm, '\n',
      "Mean MeanDecreaseGini of Featured genes is", mmDGf, '\n',
      sep = ' ')

  if (!is.nan(mmDGm)) {
    if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
      warning("A significant portion of features are missing...")
    }
  }


  rm(testsub, missingGenes, missingGenes.df)
  gc()
  #Predict
  library(entropy)
  testPred <- as.data.frame(predict(model, TestData, type = "prob"))
  class.n <- length(model$finalModel$classes)
  testPred$Diff <- apply(testPred, 1, function(x) max(x)-sort(x,partial=length(x)-1)[length(x)-1])
  testPred$KLe <- apply(testPred[,which(!names(testPred) %in% c("Diff"))], 1, function(x) KL.empirical(y1 = as.numeric(x), y2 = rep(1/class.n, class.n)) )
  testPred$BestVotesPercent <- apply(testPred[,which(!names(testPred) %in% c("Diff","KLe"))],1, function(x) max(x)  )
  testPred$Prediction <- predict(model, TestData, type="raw")
  #Flag cell type prediction if Kullback-Leibler divergence value is higher than 0.5 OR the difference between the highest and the second highest percent vote (Diff) is higher than two time of random vote rate (2/class.n)
  testPred <- testPred %>% as.tibble() %>% mutate(Intermediate = Prediction ) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(Prediction = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", as.character(Prediction) )) %>% as.data.frame()
  testPred <- testPred %>% as.tibble() %>% mutate(PredictionStatus = if_else( (KLe <= 0.25) | (Diff <= 2/class.n), "Undetermined", "Detected")) %>% as.data.frame()
  #testPred <- testPred %>% as.tibble() %>% mutate(Prediction = ifelse( (KLe <= 0.5) | (Diff <= 2/class.n), "Unclassified", as.character(Prediction) )) %>% as.data.frame()

  if (missing(priorLabels)) {
    print("Prior class labels are not provided!")

  }else{
    #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    testPred <- cbind(testPred, priorLabels)

    #Plot the crosscheck here:
    #Crosscheck Predictions
    library(tidyverse)
    library(alluvial)
    library(ggalluvial)
    crx <- testPred %>% group_by(Prior, Intermediate, Prediction) %>% tally() %>% as.data.frame()

    p5 <- ggplot(crx,aes(y = n, axis1 = Prior, axis2 = Intermediate, axis3 = Prediction )) +
      geom_alluvium(aes(fill = Prediction), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "grey") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("Prior", "Clusters", "Int-Prediction", "Final-Prediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class.n, base_width = 1.2*class.n)
  }

  if (!missing(SeuratObject)) {
    #update predictions in the meta.data slot
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(testPred))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, testPred)

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function


