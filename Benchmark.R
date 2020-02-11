setwd(".")



# Data Preprocessing


library(Seurat)
library(dplyr)
# Load dataset
sc10x.data <- Read10X(data.dir = "/nethome/jzhou417/run_count_3kpbmc/outs/filtered_feature_bc_matrix")
# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = sc10x.data))
dense.size
sparse.size <- object.size(x = sc10x.data)
sparse.size



sc10x <- CreateSeuratObject(counts = sc10x.data, project = "sc10x", min.cells = 3, min.features = 200)
sc10x
## sc10x <- subset(sc10x, features = sample(rownames(GetAssayData(sc10x)), size = 13711*0.1))


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
sc10x[["percent.mt"]] <- PercentageFeatureSet(sc10x, pattern = "^MT-")
head(sc10x@meta.data, 5)



# Visualize QC metrics as a violin plot
VlnPlot(sc10x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



#sc10x <- subset(sc10x, subset = nFeature_RNA > 50 & nFeature_RNA < 200 & percent.mt < 10)
sc10x <- subset(sc10x, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 10)


sc10x <- NormalizeData(sc10x, normalization.method = "LogNormalize", scale.factor = 10000)



sc10x <- FindVariableFeatures(sc10x, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(sc10x), 10)
# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(sc10x)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))



all.genes <- rownames(sc10x)
sc10x <- ScaleData(sc10x, features = all.genes)



saveRDS(sc10x, file = "./sc10x-3c.rds")


# Create single cell experiment object


library(SingleCellExperiment)
sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(sc10x@assays$RNA@counts),
        logcounts = log2(as.matrix(sc10x@assays$RNA@counts) + 1)
    )
)
saveRDS(sce, file = "./sc10x-3c-sce.rds")



library(Matrix)
# save sparse matrix
sparse.sce <- Matrix(sc10x@assays$RNA@counts , sparse = T )
writeMM(obj = sparse.sce, file="./matrix.mtx")
# save genes and cells names
write(x = rownames(sce), file = "./genes.tsv")
write(x = colnames(sce), file = "./barcodes.tsv")


# Start Running


sc10x <- readRDS("./sc10x-3c.rds")
sce <- readRDS("./sc10x-3c-sce.rds")


## 1. Seurat


run_seurat = function(sc10x){
  sc10x <- RunPCA(sc10x, features = VariableFeatures(object = sc10x))
  sc10x <- FindNeighbors(sc10x, dims = 1:5)
  sc10x <- FindClusters(sc10x, resolution = 1)
  sc10x <- RunUMAP(sc10x, dims = 1:5)
  cluster1.markers <- FindMarkers(sc10x, ident.1 = 1, min.pct = 0.25)
  write.csv(sc10x[["seurat_clusters"]], ('seurat-Rrs2.csv'))
}

## 2. sc3


library(SingleCellExperiment)
library(SC3)
library(scater)
run_sc3 = function(sce){
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  # define spike-ins
  isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)
  sce_after_sc3 <- sc3(sce, ks = 10, biology = TRUE)
  # sc3_interactive(sce_after_sc3)
  sc3_export_results_xls(sce_after_sc3, filename = "sc3_results.xls")
}


## 3. RaceID3


library(RaceID)
run_raceid3 = function(sce){
  sc <- SCseq(assay(sce))
  sc <- filterdata(sc, mintotal=1, minexpr = 1, minnumber = 1,
                   LBatch = NULL, knn = 10, CGenes = NULL, FGenes = NULL, ccor = 0.4,
                   bmode = "RaceID")
  sc@ndata = sc@expdata
  sc@genes = rownames(sc@ndata)
  sc@counts = rep(1,ncol(sce))
  names(sc@counts) = colnames(sc@ndata)
  sc@cluster$features = sc@genes
  sc <- compdist(sc, metric="pearson", FSelect = FALSE, knn = NULL)
  sc <- clustexp(sc, sat = TRUE, samp = NULL, cln = NULL, clustnr = 30,
                 bootnr = 50, rseed = 17000, FUNcluster = "kmedoids")
  sc <- findoutliers(sc, probthr = 0.001, outminc = 5, outlg = 2,
                     outdistquant = 0.95)
  colData(sce)$clustering_res = as.factor(sc@cpart)
  write.csv(as.matrix(sc@cpart), 'RaceID3.csv')
}

run_seurat(sc10x)
run_sc3(sce)
run_raceid3(sce)
