library(ggplot2)
library(fdANOVA)
library(fda)
library(dplyr)
source("~/Documents/big_data/single-cell/new_dy_network_07082021/fun_code/fun.R")
setwd("~/Documents/big_data/single-cell/analysis03092021/zebrafish_data")
load("trajectory_cellnames_zebrafish.RData")
load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/count_200_gene.RData")
load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/time_label.RData")
load("cluster_metadata_zebrafish.RData")
load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/time_label.RData")



cell_id = cell_name_traj[,12]
sub_metadata <- subset(metadata, rownames(metadata) %in% cell_id)
sub_cluster <- data.frame(rownames(sub_metadata),sub_metadata[,5])
colnames(sub_cluster) <- c("last_time","cluster")
trajectory_cluster = cell_name_traj
trajectory_cluster$cluster = 0
for (ii in 1:311){
  name = trajectory_cluster[ii,12]
  id = as.numeric(as.character(sub_cluster[sub_cluster$last_time==name,][,2]))
  trajectory_cluster$cluster[ii] = id
  print(ii)
}

trajectory_cluster = trajectory_cluster[-which(trajectory_cluster$cluster == 19),]
cell_name_traj = trajectory_cluster[,-ncol(trajectory_cluster)]

#### get gene expression for each trajectory

## log count ##
data_count = count
gene_name = colnames(count)
data_count$index = rownames(data_count)
time_label$index = rownames(time_label)
time_label$time = as.numeric(factor(time_label$Stage))
time = time_label[c("index", "time")]
rownames(time) <- NULL
data_count = merge(data_count,time, by = "index")
data_count <- data_count %>% relocate(index, .before = time)

nTime = max(data_count$time)
nGenes = ncol(count)

trajectory_name = cell_name_traj
colnames(trajectory_name) = paste0("time_", seq(1:nTime))

trajectory_expression = get_traj_exp(cell_name_traj = trajectory_name, all = data_count, nTime = nTime, nGenes = nGenes, gene_name = gene_name)

#### run functional anova test ####
# gene_id = colnames(count)
# gene_id2 = gene_id
# 
# library(foreach)
# library(doParallel)
# 
# cl <- makeCluster(6)
# registerDoParallel(cl)
# 
# group_label <- trajectory_cluster$cluster
# 
# r <- foreach(gene = gene_id2, .combine=cbind) %dopar% {
#   library(ggplot2)
#   library(fdANOVA)
#   library(fda)
#   
#   test_data = matrix(nrow = length(trajectory_expression), ncol = nTime)
#   for (i in 1:length(trajectory_expression)){
#     temp = trajectory_expression[[i]]
#     temp2 = temp[rownames(temp)==gene,]
#     test_data[i,] = temp2
#   }
#   
#   colnames(test_data) = paste("time",seq(1:12))
#   rownames(test_data) = paste("trajectory",seq(1:length(trajectory_expression)))
#   use = t(test_data)
#   fanova <- fanova.tests(x = use, group.label = group_label, test = "FP", parallel = T)
#   pp = fanova[["FP"]][["pvalueFP"]]
#   print(gene)
#   
# }


###### for loop ############


gene_id = colnames(count)
gene_id2 = "DBX1A"

group_label <- trajectory_cluster$cluster
p_values = c()

for (gene in gene_id2){
  test_data = matrix(nrow = length(trajectory_expression), ncol = nTime)
  for (i in 1:length(trajectory_expression)){
    temp = trajectory_expression[[i]]
    temp2 = temp[rownames(temp)==gene,]
    test_data[i,] = temp2
  }
  
  colnames(test_data) = paste("time",seq(1:12))
  rownames(test_data) = paste("trajectory",seq(1:length(trajectory_expression)))
  use = t(test_data)
  fanova <- fanova.tests(x = use, group.label = group_label, test = "FP", parallel = T)
  pp = fanova[["FP"]][["pvalueFP"]]
  p_values = c(p_values,pp)
  print(gene)
}

first_200 = p_values

##### genes analysis ######

p_values = t(first_200)
#mucl = data.frame(gene = gene_id, p_value = t(r))
p_adjust = p.adjust(p_values, method = "bonferroni")

mucl = data.frame(gene = gene_id, p_value = first_200, p_adjust = p_adjust)

sign = mucl[mucl$p_adjust < 0.01,]
#save(sign,mucl, file = "fanova_genes_no19.RData")

# write.table(sign$gene, file = "gene_list.txt", append = FALSE, sep = " ", dec = ".",
#             row.names = F, col.names = F, quote = F)
####### Divided Gene ontology #######

functional_annotation <- read.delim(
  "~/Documents/big_data/single-cell/analysis03092021/zebrafish_data/functional_annotation.txt")

library(xtable)
gene_function = functional_annotation[c("Term","Count","Bonferroni")]
xtable(head(gene_function))


###### Plot significant gene
sign_gene = sign$gene
gene = "ALPI.1"

test_data = matrix(nrow = 309, ncol = 12)
for (i in 1:length(trajectory_expression)){
  temp = trajectory_expression[[i]]
  temp2 = temp[rownames(temp)==gene,]
  test_data[i,] = temp2
}

colnames(test_data) = paste("time",seq(1:12))
rownames(test_data) = paste("trajectory",seq(1:309))

use = t(test_data)
group_label <- trajectory_cluster$cluster

# plotFANOVA(x = use, int = c(0.5,12.5))
# plotFANOVA(x = use, group.label = as.character(group_label),int = c(0.5,12.5))
# plotFANOVA(x = use, group.label = as.character(group_label),int = c(0.5,12.5), separately = TRUE)
gra = plotFANOVA(x = use, group.label = as.character(group_label),int = c(0.5,12.5), means = TRUE)

gra[["labels"]][["colour"]] = "Trajectory Group"
gra[["labels"]][["x"]] = "Time"
gra[["labels"]][["y"]] = "Expression"

gra + theme_bw() + ggtitle("ALPI.1 Gene")

write.csv(sign, file = "function_annotion.csv")
write.csv(cluster1, file = "cluster1.csv")
write.csv(cluster_2, file = "cluster2.csv")


# 
#   
#   test_data = matrix(nrow = 311, ncol = 12)
#   for (i in 1:length(expression)){
#     temp = expression[[i]]
#     temp2 = temp[rownames(temp)==gene,]
#     test_data[i,] = temp2
#   }
#   
#   colnames(test_data) = paste("time",seq(1:12))
#   rownames(test_data) = paste("trajectory",seq(1:311))
#   use = t(test_data)
#   fanova <- fanova.tests(x = use, group.label = group_label, test = "FP", parallel = T)
#   fanova[["FP"]][["pvalueFP"]]
#   
#   print(gene)
#   
# }
# 


head(seurat_df, 5)


##### prepare dataset
gene_id = colnames(seurat_df[,-ncol(seurat_df)])
gene = gene_id[1]
p_values = c()

test_data = matrix(nrow = 311, ncol = 12)
for (i in 1:length(expression)){
  temp = expression[[i]]
  temp2 = temp[rownames(temp)==gene,]
  test_data[i,] = temp2
  print(paste0("add row ",i))
}

colnames(test_data) = paste("time",seq(1:12))
rownames(test_data) = paste("trajectory",seq(1:311))

use = t(test_data)
group_label <- trajectory_cluster$cluster

plotFANOVA(x = use, int = c(0.5,12.5))
plotFANOVA(x = use, group.label = as.character(group_label),int = c(0.5,12.5))
plotFANOVA(x = use, group.label = as.character(group_label),int = c(0.5,12.5), separately = TRUE)
plotFANOVA(x = use, group.label = as.character(group_label),int = c(0.5,12.5), means = TRUE)

fanova <- fanova.tests(x = use, group.label = group_label)
fanova



# seurat_df = data.frame(t(seurat_df))
# seurat_df$index = rownames(seurat_df)
# unique(trajectory_cluster$cluster)
# 
# expression = list()
# for (row in 1:309){
# 
#   name_set = c()
#   row_num = row
#   name = as.character(cell_name_traj[row_num,1])
#   value = as.numeric(t(seurat_df[seurat_df$index == name,]))
#   value = value[c(1:(length(value)-1))]
#   gene_value = as.data.frame(value)
#   name_set = c(name_set,name)
# 
#   for (nn in 2:12){
#     name = as.character(cell_name_traj[row_num,nn])
#     value = as.numeric(t(seurat_df[seurat_df$index == name,]))
#     value = value[c(1:(length(value)-1))]
#     value = as.data.frame(value)
#     gene_value = cbind(gene_value,value)
#     name_set = c(name_set,name)
#   }
# 
#   gga = as.matrix(gene_value)
#   rownames(gga) = colnames(seurat_df)[1:2000]
#   colnames(gga) = name_set
#   expression[[row]] = gga
#   print(row)
# }





