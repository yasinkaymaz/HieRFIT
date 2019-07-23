# HieRFIT
## Hierarchical Random Forest for Information Transfer

![](data/extra/HieRFIT_banner.jpg)

Increasing demand for data integration and cross-comparison in bioinformatics drive the need for new methods. One of the applications of this R package is determining major cell types of samples in the single cell RNAseq (scRNAseq) datasets. Common methods for deciding type of the cells often involve manually checking known tissue/cell specific gene expression levels. Sensitive methods for automatically determining cell types is currently lacking. Therefore, our effort is to develop a machine learning approach to predict cell labels using gene expression levels from scRNAseq datasets. This will allow researchers to find out existing cell types in their experimental outcomes and use these predictions to further fine tune their data for downstream analysis, such as marker gene selection, differential expression test, etc.

## How to install

```
install.packages("devtools")

library(devtools)

install_github("yasinkaymaz/HieRFIT")

```


### Quick tutorial

The goal is to project cell type information from an existing experiment to our new experiment. We wonder whether our experiment went well and recovered a diverse cell population or not. Let's say our data is a new PBMC scRNAseq data.

We pretend that our data is this another PBMC run from 10X Genomics:
[Test set data is from 10X Genomics](http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz)
Note: Please, unzip all files and rename features.tsv as genes.tsv!


To use as a reference, get the data [PBMC single cell RNAseq data with 2,700 cell](https://www.dropbox.com/s/kwd3kcxkmpzqg6w/pbmc3k_final.rds?dl=0). The final data is in Seurat object format and processed by following the [tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html)



```{r}
library(Seurat)
# Process the data:
new.pbmc.data <- Read10X("pbmc_10k_v3/filtered_feature_bc_matrix/")

#The function SeuratWrapper is a wrapper function for quick processing.
newPBMC <- SeuratWrapper(ExpData = new.pbmc.data, ProjectLabel = "newPBMC", Normalize = T, scale.only.var = T, PCs = 10, dump.files = F )

```


Load the reference data
```{r}
pbmc <- get(load("pbmc3k_final.Rda"))

# A tree file for cell type hierarchy:
treeTable <- read.delim("pbmc3k_taxa.txt", header = F)

```

Then, create the reference model using the cell class labels. Here, we can use a topology tree for defining relationship between the cell groups ("pbmc3k_tree").

```{r}
library(HieRFIT)

refmod <- CreateHieR(Ref = as.matrix(pbmc@data),
                          ClassLabels = pbmc@meta.data$ClusterNames_0.6,
                          TreeTable = treeTable)

#Project the cell class labels on the new dataset:
ProObj <- HieRFIT(Query = as.matrix(newPBMC@data), refMod = refmod)

newPBMC@meta.data$ProjectedCellTypes <- ProObj@Projection

```
Alternatively, you can directly update the Seurat object:

```{r}
#Project the cell class labels on the new dataset:
newPBMC <- HieRFIT(Query = newPBMC, refMod = refmod)

head(newPBMC@meta.data)

```
meta.data slot of the object will carry all class probailities as well as the final cell type prediction column, "Projection".

### How it works:

1. Feature/predictor selection from significant principal components (PCA),
2. Build a hierarchical decision tree,
3. Create local classifiers at each node,
3. Combine local classifiers into a reference model,
4. HieRFIT new data using reference.


## Projecting cell types inter-species

HieRFIT can take an argument "xSpecies" for inter-species cross projection of cell types. In this example below, we demonstrate how to use a HieR model, "refmodZr4" build on mouse brain single cell atlas to project cell type labels on to a rat single cell data. Using biomaRt package, we find the orthologous genes between two species.

```{r}
hippo006.HierObj <- HieRFIT(Query = as.matrix(hippo006@data), refMod = refmodZr4, xSpecies = "mouse2rat")
```
This command above also creates a dataframe called "ortoDict" and saves it in the global environment. This table can be fed into other HieRFIT runs using "ortoDict" argument, if same xSpecies is the case, to save time as follows:

```{r}
hippo006.HierObj <- HieRFIT(Query = as.matrix(hippo006@data),
                            refMod = refmodZr4,
                            xSpecies = "mouse2rat",
                            ortoDict = ortoDict)
```
