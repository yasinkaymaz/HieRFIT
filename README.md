# HieRFIT
## Hierarchical Random Forest for Information Transfer

![](extra/HieRFIT_banner.jpg)

One of the applications can be determining major cell types of samples in the single cell RNAseq (scRNAseq) datasets. Common methods for deciding type of the cells often involve manually checking known tissue/cell specific gene expression levels. Sensitive methods for automatically determining cell types is currently lacking. Therefore, our effort is to develop a machine learning approach to predict cell labels using gene expression levels from scRNAseq datasets. This will allow researchers to find out existing cell types in their experimental outcomes and use these predictions to further fine tune their data for downstream analysis, such as marker gene selection, differential expression test, etc.

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
new.pbmc.data <- Read10X("~/data/pbmc_10k_v3/filtered_feature_bc_matrix/")

#The function SeuratWrapper is a wrapper function for quick processing.
newPBMC <- SeuratWrapper(ExpData = new.pbmc.data, ProjectLabel = "newPBMC", Normalize = T, scale.only.var = T, PCs = 10, dump.files = F )

dim(newPBMC@meta.data)
```
[1] 2638    7


Load the reference data
```{r}
pbmc <- get(load("~/Documents/RFTyper/pbmc3k_final.Rda"))

head(pbmc@meta.data)

```
nGene nUMI orig.ident percent.mito res.0.6 ClusterNames_0.6 res.0.8
AAACATACAACCAC   781 2421   10X_PBMC  0.030177759       0      CD4 T cells       1
AAACATTGAGCTAC  1352 4903   10X_PBMC  0.037935958       2          B cells       3
AAACATTGATCAGC  1131 3149   10X_PBMC  0.008897363       0      CD4 T cells       1
AAACCGTGCTTCCG   960 2639   10X_PBMC  0.017430845       1  CD14+ Monocytes       5
AAACCGTGTATGCG   522  981   10X_PBMC  0.012244898       5         NK cells       6
AAACGCACTGGTAC   782 2164   10X_PBMC  0.016643551       0      CD4 T cells       1


Then, create the reference model using the cell class labels. Here, we can use a topology tree for defining relationship between the cell groups ("pbmc3k_tree").

```{r}
library(HieRFIT)

refmod <- CreateRef(Ref = as.matrix(pbmc@data),
                          ClassLabels = pbmc@meta.data$ClusterNames_0.6,
                          TreeFile = pbmc3k_tree
                          )

#Project the cell class labels on the new dataset:
ProObj <- HieRFIT(Query = as.matrix(newPBMC@data), refMod = refmod)

```
