
CellTyperTrainer2 <- function(ExpData, ClassLabels, model.method="rf", run.name, PCs, cv.k=5) {
  library(randomForest)
  #library(rfUtilities)
  library(tidyverse)
  library(caret)

  if (missing(PCs)) {
    PCs <- length(unique(ClassLabels))
  }else{
    PCs <- PCs
  }

  ExpData <- as.matrix(ExpData)

  if (file.exists(paste(run.name,".trainingData.postPCA.data",sep = ""))) {
    print("Training data already exists...")
    trainingData <- get(load(paste(run.name,".trainingData.postPCA.data",sep = "")))
  }else{
    print("creating the training data...")
    trainingData <- prepareDataset(ExpData = ExpData, ClassLabels = ClassLabels, PCs = PCs, prefix = run.name)
  }

  #k-fold Cross Validation
  train.control <- trainControl(method="cv", number=cv.k, savePredictions = TRUE)
  model <- train(CellType~., data=trainingData, trControl=train.control, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500)
  save(model, file=paste(run.name,".RF_model.Robj",sep = ""))
  return(model)
}

CellTyper2 <- function(SeuratObject, testExpSet, model, priorLabels, outputFilename="plotpredictions") {

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

    PlotPredictions2(SeuratObject = SeuratObject, model = model, outputFilename = outputFilename)

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function

VarAsString <- function(x) {deparse(substitute(x))}

HTyper2 <- function(SeuratObject, testExpSet, models, priorLabels, outputFilename="plotpredictions") {

  #models is a list of of rf models
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

  Htable <- data.frame(cells = rownames(testExpSet))

  for (model in models) {
    modelname <- VarAsString(model)
    print(paste("Predicting with model",modelname,"...",sep = " "))

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

    #if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
    #  warning("A significant portion of features are missing...")
    #}

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

    Htable <- data.frame(Htable, modelname = testPred$Prediction)

  }#closes models for loop

  if (missing(priorLabels)) {
    print("Prior class labels are not provided!")

  }else{
    #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    Htable <- cbind(priorLabels,Htable)
    print(head(Htable))
    #Plot the crosscheck here:
    #Crosscheck Predictions
    library(tidyverse)
    library(alluvial)
    library(ggalluvial)
    Htable <- Htable %>% select(-cells) %>% droplevels()
    print(head(Htable))
    print(names(Htable))
    crx <- Htable %>% group_by_at(vars(one_of(names(Htable)))) %>% tally() %>% as.data.frame()
    print(head(crx))

    p5 <- ggplot(crx,aes_string(y = "n", axis1 = names(crx)[1], axis2 = names(crx)[2], axis3 = names(crx)[3], axis4 = names(crx)[4] )) +
      geom_alluvium(aes_string(fill = names(crx)[4]), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "red") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = names(crx), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    save_plot(filename = paste(outputFilename,".prediction-crosscheck.pdf",sep=""),plot = p5, base_height = 1.2*class.n, base_width = 1.2*class.n)
  }#closes missing PriorLabels

  if (!missing(SeuratObject)) {
    #update predictions in the meta.data slot
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(testPred))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, testPred)

    PlotPredictions2(SeuratObject = SeuratObject, model = model, outputFilename = outputFilename)

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(testPred) <- rownames(testExpSet)
    return(testPred)
  }#Closes missing(SeuratObj)
}#closes the function


HTyper22 <- function(SeuratObject, testExpSet, taxTable, models, priorLabels, outputFilename="plotpredictions") {

  #models is a list of of rf models
  library(caret)
  library(randomForest)
  library(tidyverse)
  if (missing(taxTable)) {stop("Please provide a proper taxanomy table as a dataframe with 'taxTable' ... exiting!")}else{taxtable <- taxTable}

  if (!missing(SeuratObject)) {
    testExpSet <- t(as.matrix(SeuratObject@data))
  }else{print("Expression matrix is provided...")
    testExpSet <- t(as.matrix(testExpSet))
  }#Closes missing(SeuratObj)

  colnames(testExpSet) <- make.names(colnames(testExpSet))

  Htable <- data.frame(cells = rownames(testExpSet))
  Bigtable <- data.frame(cells = rownames(testExpSet))
  i <- 1
  for (model in models) {
    modelname <- VarAsString(model)
    print(paste("Predicting with model",modelname,"...",sep = " "))

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

    #if ((mmDGm > mmDGf)|(length(colnames(testsub)) < length(missingGenes) )) {
    #  warning("A significant portion of features are missing...")
    #}

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

    colnames(testPred) <- paste(modelname,i, names(testPred),sep = ".")
    Bigtable <- cbind(Bigtable, testPred)
    Htable <- data.frame(Htable, modelname = testPred[,  paste(modelname,i,"Prediction",sep = ".")])

    i=i+1
  }#closes models for loop


  ConditionalProbTable <- data.frame(matrix(ncol = 0, nrow = length(rownames(Bigtable))))
  leafNames <- NULL
  #For each leaf:
  for(j in 1:dim(taxtable)[1]) {
    leafName <- paste(taxtable[j,dim(taxtable)[2]],sep = "")
    print(leafName)
    #Calculate the Conditional Probabilities
    nodeNames <- NULL
    for(i in 1:length(models.list)) {
      print(paste("model",i, taxtable[j,i], sep="."))
      nodeNames <- c(nodeNames, paste("model",i, taxtable[j,i], sep="."))
    }
    ConditionalProbTable <- cbind(ConditionalProbTable, matrixStats::rowProds(as.matrix(Bigtable[,nodeNames])))
    leafNames <- c(leafNames, leafName)
  }
  colnames(ConditionalProbTable) <- leafNames
  ConditionalProbTable$FinalBestProb <- apply(ConditionalProbTable, 1, function(x) max(x) )
  ConditionalProbTable$FinalPrediction <- colnames(ConditionalProbTable[,which(!colnames(ConditionalProbTable) %in% c("FinalBestProb"))])[apply(ConditionalProbTable[,which(!colnames(ConditionalProbTable) %in% c("FinalBestProb"))],1,which.max)]

  Htable$FinalPrediction <- ConditionalProbTable$FinalPrediction


  if (missing(priorLabels)) {print("Prior class labels are not provided!")}else{
    #Provided prior class labels (priorLabels) has to be a dataframe with same rownames as input testExpSet with one column storing labels.
    priorLabels <- as.data.frame(priorLabels)
    colnames(priorLabels) <- c("Prior")
    Htable <- cbind(priorLabels,Htable)
    print(head(Htable))
    save(Htable, file=paste(outputFilename,".prediction-crosscheck.Full_htable.Rdata",sep=""))

    #Plot the crosscheck here: For interactive Sankey diagram: https://www.r-graph-gallery.com/sankey-diagram/
    #Crosscheck Predictions
    library(tidyverse)
    library(alluvial)
    library(ggalluvial)
    Htable <- Htable %>% select(-cells) %>% droplevels()
    print(head(Htable))
    print(names(Htable))
    crx <- Htable %>% group_by_at(vars(one_of(names(Htable)))) %>% tally() %>% as.data.frame()
    print(head(crx))
    crx.f <- crx %>% mutate(freq = n*100 / sum(n)) %>% filter(freq > 1)

    p5 <- ggplot(crx,aes_string(y = "n", axis1 = names(crx)[1], axis2 = names(crx)[2], axis3 = names(crx)[3], axis4 = names(crx)[4], axis5 = names(crx)[5], axis6 = names(crx)[6] )) +
      geom_alluvium(aes_string(fill = names(crx)[6]), width = 1/12) +
      geom_stratum(width = 1/12, fill = "black", color = "red") +
      geom_label(stat = "stratum", label.strata = TRUE) +
      scale_x_discrete(limits = c("PriorLabels","Zeisel.Tax.Rank1","Zeisel.Tax.Rank2","Zeisel.Tax.Rank3","Zeisel.Tax.Rank4","FinalPrediction"), expand = c(.05, .05)) +
      ggtitle("Predictions Cross-Check")

    p6 <- ggplot(crx.f,aes_string(y = "n", axis1 = names(crx.f)[1], axis2 = names(crx.f)[2], axis3 = names(crx.f)[3], axis4 = names(crx.f)[4], axis5 = names(crx.f)[5], axis6 = names(crx.f)[6])) +
      geom_alluvium(aes_string(fill = names(crx.f)[6]), width = 0, knot.pos = 1/4) +
      guides(fill = FALSE)+
      geom_stratum(width = 1/12, fill = "grey", color = "red") +
      geom_label(stat = "stratum", label.strata = TRUE ) +
      ylab("Frequency")+
      scale_x_discrete(limits = c("PriorLabels","Zeisel.Tax.Rank1","Zeisel.Tax.Rank2","Zeisel.Tax.Rank3","Zeisel.Tax.Rank4","FinalPrediction"), expand = c(.05, .05)) +
      ggtitle(paste("Predictions","Cross-Check",sep = " "))

    pdf(paste(outputFilename,".prediction-crosscheck.pdf",sep=""),width = 40,height = 30)
    print(p5)
    print(p6)
    dev.off()
    save(crx, file=paste(outputFilename,".prediction-crosscheck.htable.Rdata",sep=""))
  }#closes missing PriorLabels

  if (!missing(SeuratObject)) {
    #update predictions in the meta.data slot
    SeuratObject@meta.data <- SeuratObject@meta.data[,which(!colnames(SeuratObject@meta.data) %in% colnames(ConditionalProbTable))]
    SeuratObject@meta.data <- cbind(SeuratObject@meta.data, ConditionalProbTable)

    PlotPredictions22(SeuratObject = SeuratObject, outputFilename = outputFilename)

    return(SeuratObject)

  }else{
    print("Prediction output is being exported ...")
    rownames(Bigtable) <- rownames(testExpSet)
    return(Bigtable)
  }#Closes missing(SeuratObj)
}#closes the function

TrainPrep <- function(model, modeSeurat, RankLabellist ) {
  imp <- as.data.frame(model$finalModel$importance)
  features <- rownames(head(imp[order(imp$MeanDecreaseGini, decreasing=T),],200))
  cells <- rownames(model$finalModel$votes)
  RLR <- modeSeurat@meta.data[,RankLabellist]
  colns <- NULL
  for(i in 1:(length(RankLabellist)-1 ) ) {colns <- c(colns,c(paste("R",i,sep="")))}
  print(colns)
  colnames(RLR) <- c(colns,"CellType")
  trainingData <- as.data.frame(t(as.matrix(modeSeurat@data)))
  #It is important to replace '-' with '.' in the names, otherwise, th rf function will throw error: 'Error in eval(predvars, data, env) : object 'RP11-206L10.2' not found'
  names(trainingData) <- make.names(names(trainingData))
  trainingData <- cbind(trainingData[,features], RLR)
  #Convert categorical variables to factors, otherwise, it will throw error: 'Error in y - ymean : non-numeric argument to binary operator'
  trainingData$CellType <- factor(trainingData$CellType)
  gc()
  return(trainingData)
}

RecursiveTrainer <- function(trainingData, model.method="rf", run.name, cv.k=5) {
  library(randomForest)
  library(rfUtilities)
  library(tidyverse)
  library(caret)

  #k-fold Cross Validation
  #train.control <- trainControl(method="cv", number=cv.k, savePredictions = TRUE)
  model <- train(CellType~., data=trainingData, method=model.method, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500)

  save(model, file=paste(run.name,".RFrec_model.Robj", sep = ""))

  return(model)
}



#Test functions

CellTyperTrainer <- function(ExpData, ClassLabels, run.name, PCs){
  library(randomForest)
  #library(rfUtilities)
  library(tidyverse)

  if(missing(PCs)){
    PCs <- length(unique(ClassLabels))
  }else{
    PCs <- PCs
  }

  ExpData <- as.matrix(ExpData)

  if(file.exists(paste(run.name,".trainingData.postPCA.data",sep = ""))){
    print("Training data already exists...")
    trainingData <- get(load(paste(run.name,".trainingData.postPCA.data",sep = "")))
  }else{
    print("creating the training data...")
    trainingData <- prepareDataset(ExpData = ExpData, ClassLabels = ClassLabels, PCs = PCs, run.name = run.name)
  }

  #Added: "sampsize=c(table(trainingData$CellType))". Revisit this later to make sure it is working as expected...
  rf <- randomForest(CellType~., data = trainingData, norm.votes = TRUE, importance=TRUE, proximity = TRUE, ntree=500, sampsize=c(table(trainingData$CellType)))
  print(rf)

  save(rf, file=paste(run.name,".RF_model.Robj",sep = ""))
  return(rf)
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
