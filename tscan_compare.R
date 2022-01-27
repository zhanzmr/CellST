library(CancerInSilico)
library(gplots)
library(SingleCellExperiment)
library(scater)
library(reticulate)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(gss)
library(TSCAN)
library(igraph)

set.seed(34895)
load("~/Documents/big_data/single_cell/new_dy_network_07082021/sim2/sim2_for_tscan.RData")

## time label
count = data.frame(t(ge))
name = rownames(count)
time_label = c(rep(c(1:13),each = 80),rep(c(1:13),each = 80))
time = data.frame(name,time_label)
rownames(time) = time$name


procdata <- preprocess(ge, minexpr_percent = 0.1)
lpsmclust <- exprmclust(procdata)
plotmclust(lpsmclust, show_cell_names = F, cell_name_size = 2)

#making not good graph
pca_coord = as.data.frame(lpsmclust$pcareduceres)
pca_coord$time = time_label
cluster_center = as.data.frame(lpsmclust$clucenter)

gra_tscan = ggplot(pca_coord, aes(x=PC1, y=PC2, color=factor(time))) +
  geom_point() + annotate(geom="text", x=cluster_center[1,1], y=cluster_center[1,2], size = 7, label="S1", color="black") + 
  annotate(geom="text", x=cluster_center[2,1], y=cluster_center[2,2], size = 7, label="S3", color="black") + 
  annotate(geom="text", x=cluster_center[3,1], y=cluster_center[3,2], size = 7, label="S2", color="black") +
  annotate(geom="text", x=cluster_center[4,1], y=cluster_center[4,2], size = 7, label="S4", color="black") + 
  geom_segment(aes(x = cluster_center[1,1], y = cluster_center[1,2], xend = cluster_center[3,1], yend = cluster_center[3,2]), color = "black") +
  geom_segment(aes(x = cluster_center[3,1], y = cluster_center[3,2], xend = cluster_center[2,1], yend = cluster_center[2,2]), color = "black") +
  geom_segment(aes(x = cluster_center[2,1], y = cluster_center[2,2], xend = cluster_center[4,1], yend = cluster_center[4,2]), color = "black")

gra_tscan[["labels"]][["colour"]] = "Experimental Time"
gra_tscan = gra_tscan + theme_bw() + theme(legend.position="none")
gra_tscan

#geom_text(x=cluster_center[1,1], y=cluster_center[1,2], label="S1", label.size = 4)

?annotate
#+ geom_text(x=3, y=30, label="Scatter plot")

### Create not good graph
mclustobj = lpsmclust
x = 1
y = 2 
MSTorder = NULL
show_tree = T
show_cell_names = F
cell_name_size = 3 
markerexpr = NULL

color_by = "stage"
lib_info_with_pseudo <- data.frame(State = mclustobj$clusterid, 
                                   sample_name = names(mclustobj$clusterid))
lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
S_matrix <- mclustobj$pcareduceres
pca_space_df <- data.frame(S_matrix[, c(x, y)])
colnames(pca_space_df) <- c("pca_dim_1", "pca_dim_2")
pca_space_df$sample_name <- row.names(pca_space_df)
edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", 
                 by.y = "sample_name")
edge_df$markerexpr <- markerexpr[edge_df$sample_name]
#edge_df = edge_df[edge_df$sample_name %in% nn,]

edge_df$stage = 0
for (i in 1:nrow(edge_df)){
  temp = edge_df[i,]
  name22 = temp$sample_name
  temp2 = time[time$name == name22,]
  name33 = temp2$time_label
  if(identical(name33, character(0)) == F){
    edge_df$stage[i] = name33
  }else{
    edge_df$stage[i] = NA
  }
  print(i)
}
edge_df$stage = paste0("t",edge_df$stage)
edge_df = edge_df[complete.cases(edge_df), ]

g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2))
g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE, 
                    size = 3)
clucenter <- mclustobj$clucenter[, c(x, y)]
clulines <- NULL
if (is.null(MSTorder)) {
  allsp <- shortest.paths(mclustobj$MSTtree)
  longestsp <- which(allsp == max(allsp), arr.ind = T)
  MSTorder <- get.shortest.paths(mclustobj$MSTtree, 
                                 longestsp[1, 1], longestsp[1, 2])$vpath[[1]]
}
for (i in 1:(length(MSTorder) - 1)) {
  clulines <- rbind(clulines, c(clucenter[MSTorder[i], 
  ], clucenter[MSTorder[i + 1], ]))
}
clulines <- data.frame(x = clulines[, 1], xend = clulines[, 
                                                          3], y = clulines[, 2], yend = clulines[, 4])
g <- g + geom_segment(aes_string(x = "x", xend = "xend", 
                                 y = "y", yend = "yend", size = NULL), data = clulines, 
                      size = 1)
clucenter <- data.frame(x = clucenter[, 1], y = clucenter[, 
                                                          2], id = 1:nrow(clucenter))
g <- g + geom_text(aes_string(label = "id", x = "x", 
                              y = "y", size = NULL), data = clucenter, size = 10)
# g <- g + guides(colour = guide_legend(override.aes = list(size = 5))) + 
#   xlab(paste0("PCA_dimension_", x)) + ylab(paste0("PCA_dimension_", 
#                                                   y)) + theme(panel.border = element_blank(), axis.line = element_line()) + 
#   theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
#   theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
#   theme(legend.position = "top", legend.key.size = unit(0.3, 
#                                                         "in"), legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + 
#   theme(legend.key = element_blank()) + theme(panel.background = element_rect(fill = "white")) + 
#   theme(axis.text.x = element_text(size = 17, color = "darkred"), 
#         axis.text.y = element_text(size = 17, color = "black"), 
#         axis.title.x = element_text(size = 20, vjust = -1), 
#         axis.title.y = element_text(size = 20, vjust = 1), 
#         plot.margin = unit(c(1, 1, 1, 1), "cm"))

g <- g + guides(colour = guide_legend(override.aes = list(size = 5))) + 
  xlab(paste0("PCA", x)) + ylab(paste0("PCA", 
                                                  y))
g[["labels"]][["colour"]] = "Time Points"

g + theme_bw()



##################################
########## Tradeseq $$$$$$$$$$$$$
##################################
suppressPackageStartupMessages({
  library(slingshot); library(SingleCellExperiment)
  library(RColorBrewer); library(scales)
  library(viridis); library(UpSetR)
  library(pheatmap); library(msigdbr)
  library(fgsea); library(knitr)
  library(ggplot2); library(gridExtra)
  library(tradeSeq)
})

data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(crv, package = "tradeSeq")
data(celltype, package = "tradeSeq")


data("sce", package = "bioc2020trajectories")

shuffle <- sample(ncol(sce))
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce)$UMAP[shuffle, ],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:2)[factor(colData(sce)$pheno$treatment_id)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 1:2, bty = "n", 
       legend = levels(factor(colData(sce)$pheno$treatment_id)))
plot(reducedDims(sce)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(3, 4)[factor(colData(sce)$pheno$spatial_id)][shuffle], alpha = .5))
legend("topright", pch = 16, col = 3:4, bty = "n", legend = levels(factor(colData(sce)$pheno$spatial_id)))

layout(1)
par(mar = c(5, 4, 4, 2) + .1)

scores <- bioc2020trajectories::imbalance_score(
  rd = reducedDims(sce)$UMAP, 
  cl = colData(sce)$pheno$treatment_id,
  k = 20, smooth = 40)

grad <- viridis::plasma(10, begin = 0, end = 1)
names(grad) <- levels(cut(scores$scaled_scores, breaks = 10))
plot(reducedDims(sce)$UMAP, col = grad[cut(scores$scaled_scores, breaks = 10)],
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", cex = .8)
legend("topleft", legend = names(grad), col = grad, pch = 16, bty = "n", cex = 2 / 3)

library(slingshot)
sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$pheno$spatial_id,
                 start.clus = 'inner', approx_points = 150)



##### Test for simulation with 2 pathway
load("~/Documents/big_data/single-cell/analysis03092021/time_sim2/ge_SingleCellExperiment.RData")


sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$pheno$spatial_id,
                 start.clus = 'inner', approx_points = 150)




#### Monocle3 ####
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(monocle)

# For reproducibility
data(countMatrix, package = "tradeSeq")
counts <- as.matrix(countMatrix)
rm(countMatrix)
data(celltype, package = "tradeSeq")



###### Monocle plus tradeseq test ##########

library(Seurat)
library(SeuratObject)

object <- CreateSeuratObject(counts = ge, project = "simulation", 
                             min.cells = 3, min.features = 10)
object

object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 100)
top10 <- head(VariableFeatures(object), 10)
object <- ScaleData(object)

object <- RunPCA(object, features = VariableFeatures(object = object))

object <- FindNeighbors(object)
object <- FindClusters(object)
object <- RunUMAP(object, dims = 1:50)

DimPlot(object, reduction = "umap")


countPCA = data.frame(object@reductions$pca@cell.embeddings)
countUMAP = data.frame(object@reductions$umap@cell.embeddings)
cluster_label = data.frame(object@meta.data)

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(monocle)

celltype = cluster_label$seurat_clusters

pd <- data.frame(cells = colnames(ge), cellType = celltype)
rownames(pd) <- pd$cells
fd <- data.frame(gene_short_name = rownames(ge))
rownames(fd) <- fd$gene_short_name
cds <- newCellDataSet(ge, phenoData = new("AnnotatedDataFrame", data = pd),
                      featureData = new("AnnotatedDataFrame", data = fd))
cds <- estimateSizeFactors(cds)
cds <- reduceDimension(cds, max_components = 2)
cds <- orderCells(cds)
plot_cell_trajectory(cds)

sce <- fitGAM(cds, verbose = TRUE)

assoRes <- associationTest(sce)
head(assoRes)
startRes <- startVsEndTest(sce)

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, ge, gene = sigGeneStart)

graphic = plotSmoothers(sce, ge, gene = sigGeneStart)
graphic + theme_bw()
graphic$layers
graphic$layers[[1]] <- NULL
graphic + theme_bw()

endRes <- diffEndTest(sce)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
plotSmoothers(sce, ge, sigGene)

### plot true
fda_data = ge[rownames(ge) %in% sigGene]
fda = data.frame(name = colnames(ge), value = fda_data)
fda$time = c(rep(c(1:13),each = 80),rep(c(1:13),each = 80))
use = matrix(nrow = 160)
for (i in 1:13){
  temp = fda[fda$time == i,]
  use = cbind(use, temp$value)
}

use = use[,-1]
use = data.frame(use)
colnames(use) = paste0("time_",seq(1:13))
group_label = rep(c(1:2), each = 80)


library(ggplot2)
library(fdANOVA)
library(fda)

use = t(use)

gra = plotFANOVA(x = use, group.label = as.character(group_label),int = c(0.5,13.5), means = TRUE)

gra[["labels"]][["colour"]] = "Cells Pathway"
gra[["labels"]][["x"]] = "Time"
gra[["labels"]][["y"]] = "Expression"
gra + theme_bw() + ggtitle("Simulated m_67 Gene")

plotFANOVA


##### seurat ####

library(Seurat)
library(SeuratObject)

object <- CreateSeuratObject(counts = ge, project = "simulation", 
                             min.cells = 3, min.features = 10)
object

object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 100)
top10 <- head(VariableFeatures(object), 10)
object <- ScaleData(object)

object <- RunPCA(object, features = VariableFeatures(object = object))

object <- FindNeighbors(object)
object <- FindClusters(object, resolution = 0.6)
object <- RunUMAP(object, dims = 1:50)

DimPlot(object, reduction = "umap")


countPCA = data.frame(object@reductions$pca@cell.embeddings)
countUMAP = data.frame(object@reductions$umap@cell.embeddings)
cluster_label = data.frame(object@meta.data)

count = data.frame(t(ge))
name = rownames(count)
time_label = c(rep(c(1:13),each = 80),rep(c(1:13),each = 80))
time = data.frame(name,time_label)
rownames(time) = time$name


library(ggplot2)
time_label = time["time_label"]
plotdata_umap = merge(time_label, countUMAP, by="row.names")  # merge by row names (by=0 or by="row.names")
#stage = unique(plotdata_umap$Stage) 

gra = ggplot(plotdata_umap, aes(x=UMAP_1, y=UMAP_2)) + geom_point(color=factor(plotdata_umap$time_label))
gra

sp<-ggplot(plotdata_umap, aes(x=UMAP_1, y=UMAP_2, color=factor(plotdata_umap$time_label))) + geom_point(size = 0.3)
sp[["labels"]][["colour"]] = "Time Points"
sp + theme_bw()


count = data.frame(t(ge))
name = rownames(count)
time_label = c(rep(c(1:13),each = 80),rep(c(1:13),each = 80))
time = data.frame(name,time_label)
rownames(time) = time$name

# np = import("numpy")
# ot = import("ot")
# 
# ## cell model
# 
# simple_mod <- suppressMessages(inSilicoCellModel(initialNum=80, runTime=96,
#                                                  density=0.1, outputIncrement=3, randSeed=267))
# 
# ######## first pathway data ###############
# mitosisGeneNames <- paste("m_", seq(1:100), sep="")
# mitosisExpression <- function(model, cell, time)
# {
#   ifelse(getCellPhase(model, time, cell) == "M", 1, 0)
# }
# 
# pwyMitosis <- new("Pathway", genes=mitosisGeneNames,
#                   expressionScale=mitosisExpression)
# 
# allGenes <- c(mitosisGeneNames)
# geneMeans <- 2 + rexp(length(allGenes), 1/20)
# data <- t(pmax(sapply(geneMeans, rnorm, n=25, sd=2), 0))
# rownames(data) <- allGenes
# 
# # calibrate pathways
# pwyMitosis <- calibratePathway(pwyMitosis, data)
# 
# ### Simulate RNA-seq Data
# params <- new("GeneExpressionParams")
# params@nCells <- 80 # sample 30 cells at each time point to measure activity
# params@sampleFreq <- 3 # measure activity every 6 hours
# params@RNAseq <- TRUE # 
# params@singleCell <- FALSE # generate bulk data
# params@perError <- 0 # parameter for simulated noise
# 
# pwyMitosis@expressionScale = function(model, cell, time)
# {
#   window <- c(max(time - 2, 0), min(time + 2, model@runTime))
#   a1 <- getAxisLength(model, window[1], cell)
#   a2 <- getAxisLength(model, window[2], cell)
#   if (is.na(a1)) a1 <- 0 # in case cell was just born
#   return(ifelse(a2 < a1, 1, 0))
# }
# 
# pwys <- c(pwyMitosis)
# ge_rna <- inSilicoGeneExpression(simple_mod, pwys, params)$expression
# 
# ############# second pathway #####################
# contactInhibitionGeneNames <- paste("ci_", seq(1:100), sep="")
# contactInhibitionExpression <- function(model, cell, time)
# {
#   getLocalDensity(model, time, cell, 3.3)
# }
# pwyContactInhibition <- new("Pathway", genes=contactInhibitionGeneNames,
#                             expressionScale=contactInhibitionExpression)
# 
# allGenes <- c(contactInhibitionGeneNames)
# geneMeans <- 2 + rexp(length(allGenes), 1/20)
# data <- t(pmax(sapply(geneMeans, rnorm, n=25, sd=2), 0))
# rownames(data) <- allGenes
# 
# # calibrate pathways
# pwyContactInhibition <- calibratePathway(pwyContactInhibition, data)
# 
# ### Simulate RNA-seq Data
# params <- new("GeneExpressionParams")
# params@nCells <- 80 # sample 30 cells at each time point to measure activity
# params@sampleFreq <- 3 # measure activity every 6 hours
# params@RNAseq <- TRUE # 
# params@singleCell <- FALSE # generate bulk data
# params@perError <- 0 # parameter for simulated noise
# 
# pwys <- c(pwyContactInhibition)
# ge_rna2 <- inSilicoGeneExpression(simple_mod, pwys, params)$expression
# 
# #########################################
# ######### single cell ###################
# #########################################
# 
# #### pathway 1 ###############
# params@randSeed <- 288 
# params@RNAseq <- TRUE
# params@singleCell <- TRUE
# params@dropoutPresent <- F
# pwys <- c(pwyMitosis)
# ge1 <- inSilicoGeneExpression(simple_mod, pwys, params)$expression
# 
# 
# #### pathway 2 ###############
# params@randSeed <- 288 
# params@RNAseq <- TRUE
# params@singleCell <- TRUE
# params@dropoutPresent <- F
# pwys <- c(pwyContactInhibition)
# ge2 <- inSilicoGeneExpression(simple_mod, pwys, params)$expression
# 
# col = colnames(ge2)
# col2 <- paste(col, "_2", sep="")
# colnames(ge2) = col2
# 
# ge = cbind(ge1,ge2)
# save(ge, file = "ge_33points.RData")







