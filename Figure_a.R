library(CancerInSilico)
library(gplots)
library(SingleCellExperiment)
library(scater)
library(reticulate)
library(reshape2)
library(ggplot2)
library(gss)
library(pheatmap)
conda_list()
use_condaenv(condaenv = 'ot', required = TRUE)
# import python module
np = import("numpy")
ot = import("ot")

#load("~/Documents/big_data/single_cell/new_dy_network_07082021/sim2/sim2_for_tscan.RData")
#load("~/Documents/big_data/single_cell/new_dy_network_07082021/sim2/ge_SingleCellExperiment.RData")

#choose only one pathway for plot
sub = t(ge)
sub = sub[1:1040,]
sub = t(sub)

# umap
library(Seurat)
library(SeuratObject)

object <- CreateSeuratObject(counts = sub, project = "simulation", 
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

mucl = data.frame(t(ge))
cell_time = rep(1:13, each = 80)
pathway = rep(1:2, each = 1040)
plot_data = data.frame(cell = rownames(mucl), expression = rowMeans(mucl), time = cell_time, pathway = pathway)


gra_b = ggplot(plot_data, aes(x=time, y=expression)) + 
  geom_point(data = plot_data, aes(x=time, y=expression, color = factor(pathway)), size = 2)

color = factor(time)

gra_b[["labels"]][["colour"]] = "Cell Type"
gra_b[["labels"]][["y"]] = "UMAP2"
gra_b[["labels"]][["x"]] = "UMAP1"
gra_b + theme_bw()

### plot line
rowname1 = rownames(count_contactInhibition)
rowname2 = paste0(rowname1,"_2")
rownames(count_mitosisGene) = rowname2
gene_name = paste0("g_", seq(1:100))
x = count_contactInhibition[1:100]
y = count_mitosisGene[1:100]
colnames(x) = gene_name
colnames(y) = gene_name

mucl = rbind(x,y)

cell_time = rep(1:13, each = 80)
pathway = c("Contact Inhibition", "Mitosis")
pathway = rep(pathway, each = 1040)
plot_data = data.frame(cell = rownames(mucl), expression = rowMeans(mucl), time = cell_time, pathway = pathway)

traj1 = traj_contactInhibition
#traj2 = paste0(as.matrix(traj_mitosisGene), "_2")
traj2 <- lapply(traj_mitosisGene, function(x) paste0(x, "_2"))
traj2 <- data.frame(matrix(unlist(traj2), nrow=length(traj2), byrow=T))
traj2 = data.frame(t(traj2))

time_col = paste0("time", seq(1:13))
colnames(traj1) = time_col
colnames(traj2) = time_col

traj_total = rbind(traj1,traj2)

library(purrr)
nn = as_vector(traj_total)

traject = seq(1:160)
traject_umap <- paste("traject ", traject, sep="")
umap_group = rep(traject_umap,13)
umap_plot = data.frame(cell = nn, group = umap_group)

umap = merge(plot_data, umap_plot, by="cell")

index = sample(160, 40)
index = paste0("traject ", index)

subset_umap = umap[umap$group %in% index,]

##add line
line_connect1 <- subset_umap[subset_umap$time == 1 | subset_umap$time == 2, ]
line_connect2 <- subset_umap[subset_umap$time == 2 | subset_umap$time == 3, ]
line_connect3 <- subset_umap[subset_umap$time == 3 | subset_umap$time == 4, ]
line_connect4 <- subset_umap[subset_umap$time == 4 | subset_umap$time == 5, ]
line_connect5 <- subset_umap[subset_umap$time == 5 | subset_umap$time == 6, ]
line_connect6 <- subset_umap[subset_umap$time == 6 | subset_umap$time == 7, ]
line_connect7 <- subset_umap[subset_umap$time == 7 | subset_umap$time == 8, ]
line_connect8 <- subset_umap[subset_umap$time == 8 | subset_umap$time == 9, ]
line_connect9 <- subset_umap[subset_umap$time == 9 | subset_umap$time == 10, ]
line_connect10 <- subset_umap[subset_umap$time == 10 | subset_umap$time == 11, ]
line_connect11 <- subset_umap[subset_umap$time == 11 | subset_umap$time == 12, ]
line_connect12 <- subset_umap[subset_umap$time == 12 | subset_umap$time == 13, ]

gra2 = ggplot(subset_umap, aes(x = time, y = expression))+
  geom_point(data = subset_umap, aes(x = time, y = expression, color = factor(pathway))) + 
  geom_line(data = line_connect1, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect2, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect3, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect4, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect5, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect6, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect7, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect8, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect9, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect10, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect11, aes(x = time, y = expression, group = group), size = 0.1) +
  geom_line(data = line_connect12, aes(x = time, y = expression, group = group), size = 0.1) +
  xlab("Time Points") + ylab("Exrpession")
gra2[["labels"]][["colour"]] = "Pathway Type"
gra2 = gra2 + theme_bw() + theme(legend.position="none")


library(ggpubr)
ggarrange(gra2, gra_dd, labels = c("a", "b") ,ncol = 2, nrow = 1)
ggarrange(all_gra, gra_tscan, labels = c("c", "d") ,ncol = 2, nrow = 1)
