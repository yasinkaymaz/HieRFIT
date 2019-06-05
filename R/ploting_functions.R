
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  #This function has been copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' This function uses 'ape' package.
#' @param tree a tree object in phylo format.
#' @usage PlotTopo(RandTreeSim(LN = sample(5:20, 1)))
PlotTopo <- function(tree){### Update with data.tree plotting
  ape::plot.phylo(x = tree, direction = "downwards", show.node.label = TRUE)
  #ape::nodelabels()
}

PlotPredictions <- function(SeuratObject, model, save.pdf=T, outputFilename="plotpredictions") {
  #Evaluate model prediction accuracy:
  conf.mat <- model$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>%
    mutate(freq = 100*value/sum(value))

  class.n = length(model$classes)

  pdf(paste(outputFilename,".pdf",sep=""),width= 1.5*class.n, height = 1.5*class.n)


  require(gridExtra)

  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) +
    geom_tile(color = "white")+
    coord_equal()+
    scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_discrete(name ="Predicted Cell Types")+
    scale_x_discrete(name ="Cell Types Classes")

  errorSize <- as.data.frame(cbind(model$confusion[,"class.error"],
                                   head(colSums(model$confusion),-1)))
  colnames(errorSize) <- c("ClassError","ClassTrainingSize")
  errorSize$CellTypeClass <- rownames(errorSize)
  acc <- 100-100*colSums(model$confusion)["class.error"]/length(head(colSums(model$confusion),-1))
  names(acc) <- "accuracy"

  p2 <- ggplot(errorSize)  +
    geom_bar(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassTrainingSize),stat="identity", fill="tan1", colour="sienna3")+
    geom_point(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassError*max(errorSize$ClassTrainingSize)),stat="identity",size=3)+
    scale_y_continuous(sec.axis = sec_axis(~./max(errorSize$ClassTrainingSize),name="Class % Error rate (Dots)"))+
    labs(y="Class Size in Training (Bars)",title=paste("Model Prediction Accuracy is ",round(acc,digits = 2),"%", sep=""))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))+
    scale_x_discrete(name ="Cell Type Classes")

  grid.arrange(p1, p2, nrow=2)

  #Prediction outputs
  FeaturePlot(object = SeuratObject,
              features.plot = model$classes,
              cols.use = c("grey", "blue"),
              reduction.use = "tsne")

  TSNEPlot(SeuratObject, group.by="Prediction",do.label=T)

  FeaturePlot(SeuratObject, features.plot = "BestVotesPercent", no.legend = F, cols.use = c("gray","red"))

  FeaturePlot(SeuratObject, features.plot = "KLe", no.legend = F, cols.use = c("gray","purple"))

  p3 <- ggplot(data=SeuratObject@meta.data)+
    geom_histogram(aes(x=Prediction,fill=Prediction),stat = "count")+
    geom_violin(aes(x=Prediction,y=BestVotesPercent*max(table(SeuratObject@meta.data$Prediction)),fill=Prediction))+
    scale_y_continuous(sec.axis = sec_axis(~./max(table(SeuratObject@meta.data$Prediction)),name="Probability Scores (violin)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Number of Cells (bars)",title=paste("Prediction outcome", sep=""))

  p4 <- ggplot(SeuratObject@meta.data, aes(KLe, Diff, color= PredictionStatus))+ geom_point(size=.6)

  grid.arrange(p3, p4, nrow=2)

  dev.off()
}

PlotPredictions2 <- function(SeuratObject, model, priorLabels, outputFilename="plotpredictions") {
  #Evaluate model prediction accuracy:
  conf.mat <- model$finalModel$confusion %>% as.data.frame() %>% select(-class.error)
  conf.mat <- reshape2::melt(as.matrix(conf.mat)) %>% as.tibble() %>% group_by(Var1) %>% mutate(freq = 100*value/sum(value))
  class.n = length(model$finalModel$classes)
  errorSize <- as.data.frame(cbind(model$finalModel$confusion[,"class.error"], head(colSums(model$finalModel$confusion),-1)))
  colnames(errorSize) <- c("ClassError","ClassTrainingSize")
  errorSize$CellTypeClass <- rownames(errorSize)
  acc <- getTrainPerf(model)["TrainAccuracy"]*100
  di <- round(sqrt(class.n),digits = 0)+1

  p1 <- ggplot(conf.mat, aes(Var1, Var2, fill=freq)) +
    geom_tile(color = "white")+
    coord_equal()+
    scale_fill_gradient2(low = "white", high = "red", name="% Predictions")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_y_discrete(name ="Predicted Cell Types")+
    scale_x_discrete(name ="Cell Types Classes")

  p2 <- ggplot(errorSize)  +
    geom_bar(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassTrainingSize),stat="identity", fill="tan1", colour="sienna3")+
    geom_point(aes(x=reorder(CellTypeClass,-ClassTrainingSize), y=ClassError*max(errorSize$ClassTrainingSize)),stat="identity",size=3)+
    scale_y_continuous(sec.axis = sec_axis(~./max(errorSize$ClassTrainingSize),name="Class % Error rate (Dots)"))+
    labs(y="Class Size in Training (Bars)",title=paste("Model Prediction Accuracy is ",round(acc,digits = 2),"%", sep=""))+
    theme(axis.text.x = element_text(angle = 90,hjust = 1))+
    scale_x_discrete(name ="Cell Type Classes")

  #Prediction Performance
  p3 <- ggplot(data=SeuratObject@meta.data)+
    geom_histogram(aes(x=Prediction,fill=Prediction),stat = "count")+
    geom_violin(aes(x=Prediction,y=BestVotesPercent*max(table(SeuratObject@meta.data$Prediction)),fill=Prediction))+
    scale_y_continuous(sec.axis = sec_axis(~./max(table(SeuratObject@meta.data$Prediction)),name="Probability Scores (violin)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Number of Cells (bars)",title=paste("Prediction outcome", sep=""))

  p4 <- TSNEPlot(SeuratObject, group.by="Prediction",do.label=T, do.return = T)


  ff <- FeaturePlot(object = SeuratObject,
                    features.plot = model$finalModel$classes,
                    cols.use = c("grey", "blue"),
                    reduction.use = "tsne",do.return = T)

  q1 <- FeaturePlot(SeuratObject, features.plot = "BestVotesPercent", no.legend = F, cols.use = c("gray","red"),do.return = T)
  q2 <- FeaturePlot(SeuratObject, features.plot = "KLe", no.legend = F, cols.use = c("gray","purple"),do.return = T)
  q3 <- ggplot(SeuratObject@meta.data, aes(KLe, Diff, color= PredictionStatus))+ geom_point(size=.6)

  p1p2 <- cowplot::plot_grid(p1,NULL,NULL, p2, ncol = 2, nrow=2)
  save_plot(filename = paste(outputFilename,".mode.-stats.pdf",sep=""),plot = p1p2,base_height = 12, base_width = 12)

  p3p4 <- cowplot::plot_grid(p3,NULL,NULL,p4, ncol = 2, nrow=2)
  save_plot(filename = paste(outputFilename,".prediction-stats.pdf",sep=""),plot = p3p4,base_height = 12, base_width = 12)

  ff <- cowplot::plot_grid(plotlist=ff, nrow = di, ncol = di-1)
  save_plot(filename = paste(outputFilename,".prediction-probability-projection.pdf",sep=""),plot = ff,base_height = 12, base_width = 12)

  q1q2q3 <- cowplot::plot_grid(plotlist = list(q1$BestVotesPercent, q2$KLe, q3 ), ncol = 2, nrow = 2)
  save_plot(filename = paste(outputFilename,".prediction-quality.pdf",sep=""),plot = q1q2q3,base_height = 12, base_width = 12)

}

PlotPredictions22 <- function(SeuratObject, outputFilename="plotpredictions") {

  #Prediction Performance
  p3 <- ggplot(data=SeuratObject@meta.data)+
    geom_histogram(aes(x=FinalPrediction,fill=FinalPrediction),stat = "count")+
    geom_violin(aes(x=FinalPrediction,y=FinalBestProb*max(table(SeuratObject@meta.data$FinalPrediction)),fill=FinalPrediction))+
    scale_y_continuous(sec.axis = sec_axis(~./max(table(SeuratObject@meta.data$FinalPrediction)),name="Probability Scores (violin)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="right")+
    labs(y="Number of Cells (bars)",title=paste("Prediction outcome", sep=""))

  p4 <- TSNEPlot(SeuratObject, group.by="FinalPrediction",do.label=T, do.return = T)

  p3p4 <- cowplot::plot_grid(p3,p4, ncol = 1, nrow=2)
  save_plot(filename = paste(outputFilename,".prediction-stats.pdf",sep=""),plot = p3p4,base_height = 20, base_width = 24)

}

CrossCheck <- function(PriorPostTable, outputprefix=""){
  #This function takes a table with two columns of which first is for prior cell labels and second is for predicted cell class.
  library(tidyverse)
  library(alluvial)
  library(ggalluvial)

  crx <- PriorPostTable %>% group_by_at(vars(one_of(names(PriorPostTable))))%>% tally() %>% arrange(desc(n))

  p5 <- ggplot(crx,aes_string(y = "n", axis1 = names(crx)[1], axis2 = names(crx)[2] )) +
    geom_alluvium(aes_string(fill = names(crx)[2]), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "red") +
    geom_label(stat = "stratum", label.strata = TRUE) +
    scale_x_discrete(limits = c("PriorLabels","FinalPrediction"), expand = c(.05, .05)) +
    ggtitle("Predictions Cross-Check")
  #pdf(paste(outputprefix,".prediction-crosscheck.pdf",sep=""),width = 20,height = 15)
  #print(p5)
  #dev.off()
  return(p5)
}

PlotClusterTree <- function(object, ...) {
  if (length(x = object@cluster.tree) == 0) {
    stop("Phylogenetic tree does not exist, build using BuildClusterTree")
  }
  data.tree <- object@cluster.tree[[1]]
  plot.phylo(x = data.tree, direction = "downwards", ...)
  nodelabels()
}

PlotCellPredictionTree <- function(tree){

  ape::plot.phylo(tree,direction = "downwards")
  dd <- data.frame("L"=c(.8,.9,.5),"R"=c(.2,.1,.5))
  ape::nodelabels(thermo = dd, horiz = TRUE,cex=1.2)
}
