# ##### Test for simulation with 2 pathway
# load("ge_SingleCellExperiment.RData")
# 
# sce <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = colData(sce)$pheno$spatial_id,
#                  start.clus = 'inner', approx_points = 150)
# 
# 
# 
# 
# #### Monocle3 ####
# library(tradeSeq)
# library(RColorBrewer)
# library(SingleCellExperiment)
# library(monocle)
# 
# # For reproducibility
# data(countMatrix, package = "tradeSeq")
# counts <- as.matrix(countMatrix)
# rm(countMatrix)
# data(celltype, package = "tradeSeq")

#load("~/Documents/big_data/single_cell/new_dy_network_07082021/sim2/sim2_for_tscan.RData")
#load("~/Documents/big_data/single_cell/new_dy_network_07082021/sim2/ge_SingleCellExperiment.RData")

library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
gene_name = paste0("g_", seq(1:100))
colnames(count_contactInhibition) = gene_name
colnames(count_mitosisGene) = gene_name
com_cell = rbind(count_contactInhibition,count_mitosisGene)
com_cell = com_cell[,-c(102,101)]
tt = t(com_cell)

sce <- SingleCellExperiment(list(counts=t(com_cell)))
sce <- logNormCounts(sce)

assays(sce)$norm <- FQnorm(assays(sce)$counts)

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')
plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 10)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'kmeans', reducedDim = 'PCA')

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

library(tradeSeq)
sce <- fitGAM(sce)
ATres <- associationTest(sce)
startRes <- startVsEndTest(sce)

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[9]]
plotSmoothers(sce, ge, gene = sigGeneStart)

names(sce)

graphic = plotSmoothers(sce, ge, gene = sigGeneStart)
graphic + theme_bw()
graphic$layers
graphic$layers[[1]] <- NULL
new_gra = graphic + theme_bw()
new_gra
data = new_gra[["data"]]

# endRes <- diffEndTest(sce)
# o <- order(endRes$waldStat, decreasing = TRUE)
# sigGene <- names(sce)[o[1]]
# plotSmoothers(sce, ge, sigGene)

aa = new_gra[["layers"]][[1]][["data"]]
aa$time = seq(1:100)
plot(aa$time, aa$gene_count, type = "l")
aa2 = new_gra[["layers"]][[3]][["data"]]
aa2$time = seq(1:100)
plot(aa2$time, aa2$gene_count, type = "l")

sigGene = "g_50"
ge = t(mucl)
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

gra[["labels"]][["colour"]] = "Cells Pathways"
gra[["labels"]][["x"]] = "Time"
gra[["labels"]][["y"]] = "Expression"
gra + theme_bw() + ggtitle("Simulated m_60 Gene")

ssano = data.frame(t(use))
ssano$label = group_label

ss1 = ssano[ssano$label == 1,]
ss2 = ssano[ssano$label == 2,]
ss1 = ss1[,-14]
ss2 = ss2[,-14]

# ssanova
ss1$id <- 1:nrow(ss1)
plot_data <- melt(ss1,id.var="id")
yy <- plot_data$value
cc <- plot_data$id
fit.random <- ssanova(yy~cc, alpha=0.03)
x = data.frame(cc = seq(1, 13, length = 100))
pred1 <- predict(fit.random, newdata = x, se.fit = T )
plot(x$cc,pred2$fit, type = "l", xlab = "Time_points", ylab = "Gene Set Expression", main = "Mitosis Mean Curve", ylim=c(0.5,2))

ss2$id <- 1:nrow(ss2)
plot_data <- melt(ss2,id.var="id")
yy <- plot_data$value
cc <- plot_data$id
fit.random <- ssanova(yy~cc, alpha=0.05)
x = data.frame(cc = seq(1, 13, length = 100))
pred2 <- predict(fit.random, newdata = x, se.fit = T )
plot(x$cc,pred2$fit, type = "l", xlab = "Time_points", ylab = "Gene Set Expression", main = "Mitosis Mean Curve", ylim=c(0.5,2))


##### tradeseq data
# tra1 = data[data$lineage == 1,]
# tra1$time = seq(1:1121)
# tra2 = data[data$lineage == 2,]
# tra2$time = seq(1:959)
# 
# plot_plot = rbind(tra1, tra2)
# 
# all_gra = ggplot() + 
#   geom_line(data = plot_plot, aes(x = time,y = gene_count, group = lineage), size=0.8, color = factor(temp$group.label)) + 
#   ylab("Expression") + xlab("Time Points")
# 
# all_gra = all_gra + theme_bw()
# all_gra

plot1 = data.frame(time = seq(1:100), value = pred1$fit, group = 1)
plot2 = data.frame(time = seq(1:100), value = pred2$fit, group = 2)
plot_plot = rbind(plot1, plot2)

plot_data2 = rbind(aa,aa2)

all_gra = ggplot() +
  geom_line(data = plot_plot, aes(x = time,y = value, group = group), size=0.8, color = factor(plot_plot$group)) +
  geom_line(data = plot_data2, aes(x = time,y = gene_count, group = lineage), size=0.8, color = factor(plot_data2$lineage), linetype = "dashed") + ylab("Expression") + xlab("Time Points")

all_gra = all_gra + theme_bw()
all_gra



#### Plot together
# temp = gra[["data"]]
# 
# yhat = graphic[["plot_env"]][["yhat"]]
# temp2= data.frame(x = seq(1,13,length.out = 100), y = yhat) 
# 
# 
# all_gra = ggplot() + 
#   geom_line(data = temp, aes(x = t,y = log(values1), group = group.label), size=0.8, color = factor(temp$group.label)) + 
#   geom_line(data = temp2, aes(x = x,y = log(y-1)), size=0.8, color = "orange", linetype = "dashed") + ylab("Expression") + xlab("Time Points")
# 
# all_gra = all_gra + theme_bw()
# all_gra


# geom_line(data = plot2,aes(x = x,y = y), size=1, color = "blue") + 
#   geom_line(data = plot3,aes(x = x,y = y), size=1, color = "black") + 
#   geom_line(data = plot4,aes(x = id,y = value,group = variable), size=0.8, color = "red") + 
#   geom_line(data = plot5,aes(x = x,y = y), size=1, color = "blue") + 
#   geom_line(data = plot6,aes(x = x,y = y), size=1, color = "black") + xlab("Time Points") + 
#   ylab("Expression") + ggtitle("Individual Cell Trajectories")
# 
# 
# 
# 
# 
# data(crv, package="tradeSeq")
# data(countMatrix, package="tradeSeq")
# gamList <- fitGAM(counts = as.matrix(countMatrix),
#                   sds = crv,
#                   nknots = 5)

set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE)

gamList <- fitGAM(counts = as.matrix(countMatrix),
                  sds = crv,
                  nknots = 5)

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

# # Run PCA then UMAP on the data
# cds <- preprocess_cds(cds, method = "PCA")
# cds <- reduce_dimension(cds, preprocess_method = "PCA",
#                         reduction_method = "UMAP")
# 
# # First display, coloring by the cell types from Paul et al
# plot_cells(cds, label_groups_by_cluster = FALSE, cell_size = 1,
#            color_cells_by = "cellType")
# 
# # Running the clustering method. This is necessary to the construct the graph
# cds <- cluster_cells(cds, reduction_method = "UMAP")
# # Visualize the clusters
# plot_cells(cds, color_cells_by = "cluster", cell_size = 1)
# 
# # Construct the graph
# # Note that, for the rest of the code to run, the graph should be fully connected
# cds <- learn_graph(cds, use_partition = FALSE)
# 
# # We find all the cells that are close to the starting point
# cell_ids <- colnames(cds)[pd$cellType ==  "Multipotent progenitors"]
# closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
# closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
# closest_vertex <- closest_vertex[cell_ids, ]
# closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
# mst <- principal_graph(cds)$UMAP
# root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
# 
# # We compute the trajectory
# cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
# plot_cells(cds, color_cells_by = "pseudotime")
# 
# 
# library(magrittr)
# # Get the closest vertice for every cell
# y_to_cells <-  principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
#   as.data.frame()
# y_to_cells$cells <- rownames(y_to_cells)
# y_to_cells$Y <- y_to_cells$V1
# 
# # Get the root vertices
# # It is the same node as above
# root <- cds@principal_graph_aux$UMAP$root_pr_nodes
# 
# # Get the other endpoints
# endpoints <- names(which(igraph::degree(mst) == 1))
# endpoints <- endpoints[!endpoints %in% root]
# 
# # For each endpoint
# cellWeights <- lapply(endpoints, function(endpoint) {
#   # We find the path between the endpoint and the root
#   path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
#   path <- as.character(path)
#   # We find the cells that map along that path
#   df <- y_to_cells[y_to_cells$Y %in% path, ]
#   df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
#   colnames(df) <- endpoint
#   return(df)
# }) %>% do.call(what = 'cbind', args = .) %>%
#   as.matrix()
# rownames(cellWeights) <- colnames(cds)
# pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights),
#                      nrow = ncol(cds), byrow = FALSE)
# 
# 
# sce <- fitGAM(counts = cds@assays@data@listData[["counts"]],
#               pseudotime = pseudotime,
#               cellWeights = cellWeights)

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

pseudotime <- matrix(pseudotime(cds), ncol = ncol(cellWeights), nrow = ncol(cds), byrow = FALSE)

sce <- fitGAM(counts = cds@assays@data@listData[["counts"]],
              pseudotime = pseudotime,
              cellWeights = cellWeights)

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


