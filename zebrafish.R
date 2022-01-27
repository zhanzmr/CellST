suppressMessages(library(scater))
suppressMessages(library(CombMSC))
suppressMessages(library(fda))
suppressMessages(library(gss))
suppressMessages(library(CombMSC))
suppressMessages(library(dplyr))
suppressMessages(library(devtools))
suppressMessages(library(ESCO))

source("~/Documents/big_data/single-cell/new_dy_network_07082021/fun_code/fun.R")
source("~/Documents/big_data/single-cell/new_dy_network_07082021/fun_code/ssFCL.R")

load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/count_200_gene.RData")
load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/trajectory_cellnames_zebrafish.RData")
load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/time_label.RData")
load("~/Documents/big_data/single-cell/analysis03092021/zebrafish_data/cluster_metadata_zebrafish.RData")
load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/biotime.RData")


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

trajectory_expression = get_traj_exp(cell_name_traj = trajectory_name, all = data_count,nTime = nTime, nGenes = nGenes, gene_name = gene_name)

##### Split cluster

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

### make dynamic network for each of the clusters
cluster_label = trajectory_cluster$cluster
label_dyn = unique(cluster_label)

current_cluster = label_dyn[7]

temp_trajectory = trajectory_cluster[trajectory_cluster$cluster == current_cluster, ]
index = as.numeric(rownames(temp_trajectory))

temp_trajec_expression = trajectory_expression[index]
gene_expression = get_gene_exp(expression = temp_trajec_expression, nGenes = nGenes, gene_name = gene_name)
names(gene_expression) = gene_name

temp_biotime = biotime[,index]

res = subsets(200,2)
res = res[-9367,]
mse_all = c()
#name = gene_name
fit_fit = list()
names_name = c()
somePDFPath = "zebrafish_dy_7.pdf"
pdf(file=somePDFPath)  
for (jj in 1:nrow(res)){
  gene_gene = res[jj,]
  name = c(gene_name[gene_gene[1]], gene_name[gene_gene[2]])
  name1 = name[1]
  name2 = name[2]
  # ###### true ##########
  # gene_gene = res[jj,]
  # gene_name = c(name[gene_gene[1]], name[gene_gene[2]])
  # rho = c()
  # for (i in 1:nTime){
  #   temp = time[[i]]
  #   name1 = gene_name[1]
  #   name2 = gene_name[2]
  #   sub = subset(temp, rownames(temp) %in% name1)
  #   value = sub[, name2]
  #   rho = c(rho,value)
  # }
  # ss = data.frame(x = nTime_seq, y = rho)
  # fitss = ssanova(y~x, data = ss,alpha = 0.4)
  # hh = data.frame(x = seq(1,nTime,length = 101))
  # est_true = predict(fitss, newdata = hh)
  ###### Estimate
  gene1 = gene_expression[[gene_gene[1]]]
  gene2 = gene_expression[[gene_gene[2]]]
  
  K = nTime #num of time points
  xdomain = c(0,50)
  
  #sigma = 0.5
  pt = c()
  wt = c()
  for (i in 1:ncol(temp_biotime)){
    t.i = temp_biotime[,i]
    w.i = weight.nodes(t.i,xdomain)
    pt = c(pt,t.i)
    wt = c(wt,w.i)
  }
  quad = list(pt=pt,wt=wt)
  yqlist <- split(as.matrix(gene1), seq(nrow(gene1)))
  xqlist <- split(as.matrix(gene2), seq(nrow(gene2)))
  
  
  tryCatch({
    fit = flr2d(yqlist,xqlist,quad,method="v",alpha=0.3,id.basis=2*c(1:K),xdomain=xdomain)
    tgrid = seq(1,nTime, length = 101)
    est0 <- estimate.flr2d(fit,tgrid,se=TRUE)
    plot(tgrid,est0$fit,type="l",xlab = "t",ylab = "value", ylim = c(-0.5,0.5))
    #lines(tgrid,est_true, pch = 5, col = "red", type = "l", lty = 2)
    #lines(rho, pch = 5, col = "green", type = "l", lty = 2)
    lines(tgrid,est0$fit-1.96*est0$se.fit,lty = "dashed")
    lines(tgrid,est0$fit+1.96*est0$se.fit,lty = "dashed")
    title(main = paste0("Gene-",name1, " and Gene-", name2))
    fit_fit[[jj]] = t(est0$fit)
    names_name = c(names_name, paste0(name1,"_",name2))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  print(jj)
}

dev.off()


names(fit_fit) = names_name

# mucl <- data.frame(matrix(unlist(fit_fit), nrow=length(fit_fit), byrow=F))
# mucl <- as.data.frame(fit_fit)
mucl <- do.call(rbind.data.frame, fit_fit)
rownames(mucl) = names_name

mucl_list[[7]] = mucl
save(mucl_list, file = "mucl_all.RData")


######## reference 

gene_expression = get_gene_exp(expression = trajectory_expression, nGenes = nGenes)
names(gene_expression) = gene_name

#### correlation 
time = list()
for (ii in 1:nTime){
  temp = data_count[data_count$time == ii,]
  temp = t(temp[,c(1:(ncol(temp)-2))])
  corr = gcn(temp)
  time[[ii]] = corr
  print(ii)
}

#### Find cell means
gene_mean_cell = lapply(trajectory_expression, function(x){colMeans(x)})
gene_mean_cell = data.frame(matrix(unlist(gene_mean_cell), nrow=length(gene_mean_cell), byrow=TRUE))

###### register time
nTime_seq = seq(1:ncol(gene_mean_cell))
nfine   = 101
agefine = seq(min(nTime_seq), max(nTime_seq), length=nfine)

hgtf   = t(gene_mean_cell)
ncasef = dim(hgtf)[2]

growbasis = create.bspline.basis(nTime_seq, norder=3)
growfdPar = fdPar(growbasis,lambda = 0.1)
growthMon = smooth.basis(nTime_seq, hgtf, growfdPar)
hgtfhatfd = growthMon$fd

# un register curves
plot(hgtfhatfd, xlim=c(1,nTime), ylim=c(0,1.5), lty=1, lwd=2,
     cex=2, xlab="Time", ylab="Expression")
title(main = "Unregistered Cell Curves")

nwbasisCR = 6
norderCR  =  3
wbasisCR  = create.bspline.basis(c(1,nTime), nwbasisCR, norderCR)
Wfd0CR    = fd(matrix(0,nwbasisCR,ncasef),wbasisCR)
lambdaCR  = 0.5
WfdParCR  = fdPar(Wfd0CR, 1, lambdaCR)


registerlistCR = register.fd(yfd = hgtfhatfd,WfdParobj = WfdParCR)
regist = registerlistCR$regfd

plot(regist, xlim=c(1,nTime), ylim=c(0,1.5), lty=1, lwd=2,
     cex=2, xlab="Time", ylab="Expression")
title(main = "Registered Cell Curves")

biotime = eval.fd(c(1:nTime),registerlistCR$warpfd)

######################################################################

res = subsets(200,2)
mse_all = c()
name = colnames(count)
fit_fit = list()
names_name = c()
somePDFPath = "zebrafish_200.pdf"
pdf(file=somePDFPath)  
for (jj in 1:nrow(res)){
  ###### true ##########
  gene_gene = res[jj,]
  gene_name = c(name[gene_gene[1]], name[gene_gene[2]])
  rho = c()
  for (i in 1:nTime){
    temp = time[[i]]
    name1 = gene_name[1]
    name2 = gene_name[2]
    sub = subset(temp, rownames(temp) %in% name1)
    value = sub[, name2]
    rho = c(rho,value)
  }
  ss = data.frame(x = nTime_seq, y = rho)
  fitss = ssanova(y~x, data = ss,alpha = 0.4)
  hh = data.frame(x = seq(1,nTime,length = 101))
  est_true = predict(fitss, newdata = hh)
  ###### Estimate
  gene1 = gene_expression[[gene_gene[1]]]
  gene2 = gene_expression[[gene_gene[2]]]
  
  K = nTime #num of time points
  xdomain = c(0,20)
  
  #sigma = 0.5
  pt = c()
  wt = c()
  for (i in 1:ncol(biotime)){
    t.i = biotime[,i]
    w.i = weight.nodes(t.i,xdomain)
    pt = c(pt,t.i)
    wt = c(wt,w.i)
  }
  quad = list(pt=pt,wt=wt)
  yqlist <- split(as.matrix(gene1), seq(nrow(gene1)))
  xqlist <- split(as.matrix(gene2), seq(nrow(gene2)))
  
  
  tryCatch({
    fit = flr2d(yqlist,xqlist,quad,method="v",alpha=0.3,id.basis=2*c(1:K),xdomain=xdomain)
    tgrid = seq(1,nTime, length = 101)
    est0 <- estimate.flr2d(fit,tgrid,se=TRUE)
    plot(tgrid,est0$fit,type="l",xlab = "t",ylab = "value", ylim = c(-0.5,0.5))
    lines(tgrid,est_true, pch = 5, col = "red", type = "l", lty = 2)
    lines(rho, pch = 5, col = "green", type = "l", lty = 2)
    lines(tgrid,est0$fit-1.96*est0$se.fit,lty = "dashed")
    lines(tgrid,est0$fit+1.96*est0$se.fit,lty = "dashed")
    title(main = paste0("Gene-",name1, " and Gene-", name2))
    fit_fit[[jj]] = est0$fit
    names_name = c(names_name, paste0(name[gene_gene[1]],"_",name[gene_gene[2]]))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  print(gene_name)
}

dev.off()
names(fit_fit) = names_name

mucl <- data.frame(matrix(unlist(fit_fit), nrow=length(fit_fit), byrow=TRUE))
rownames(mucl) = names_name

save(mucl, file = "zebrafish_correlation.RData")

######### D ####







######## Smoothing spline clustering

#### Make networks 
library(igraph)

mat <- time[[11]]
# Keep only high correlations
mat[mat<0.15] <- 0
network <- graph_from_adjacency_matrix( mat, weighted=T, mode="undirected", diag=F)
plot(network, vertex.size=9, vertex.label.font=0.5)

net3 <- network(mat,matrix.type="adjacency", 
                loops=F, multiple=F, ignore.eval = F)

plot(net3)

es <- data.frame(onset=1:49, terminus=50, 
                 head=as.matrix(net3, matrix.type="adjacency")[,1],
                 tail=as.matrix(net3, matrix.type="adjacency")[,2])



############### Figure #############
load("~/Documents/big_data/single-cell/dy_network_06012021/cell_cycle_count.RData")
##### seurat ####

library(Seurat)
library(SeuratObject)

object <- CreateSeuratObject(counts = t(count), project = "simulation", 
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
time_label = c(rep(c(1:33),each = 80),rep(c(1:33),each = 80))
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

