# Load package
library(splatter)
library(scater)
library(reticulate)
library(ggplot2)

conda_list()
use_condaenv(condaenv = 'ot', required = TRUE)
np = import("numpy")
ot = import("ot")

## Making figure A and B
nGenes = 500
n = 200
params <- newSplatParams()
params <- setParams(params, update = list(batchCells = c((n/2), (n/2)), nGenes = nGenes, 
                                          group.prob = c(0.20, 0.20, 0.20, 0.20, 0.20), seed = 56893))
sim <- splatSimulate(params = params, method = "groups", verbose = FALSE)
sim <- logNormCounts(sim)
data = logcounts(sim)
time = sim@colData@listData[["Batch"]]
group = sim@colData@listData[["Group"]]

#optimal transport
batch = sim$Batch
type = as.character(sim$Group)
count = data.frame(t(data))
count$batch = batch
count$type = type

batch1 = count[count$batch == "Batch1",]
batch2 = count[count$batch == "Batch2",]

Xs = data.matrix(batch1[1:(length(batch1)-2)])
Xt = data.matrix(batch2[1:(length(batch2)-2)])

ot_emd = ot$da$EMDTransport()
ot_emd$fit(Xs=Xs, Xt=Xt)
coupling = ot_emd$coupling_

# connecting
index = apply(coupling, 1, which.max)
res = data.frame(from = seq(1:(n/2)), to = index)

result = res
result$b1 = batch1$type

re = batch2[c(res$to),]
result$b2 = re$type

result2 = result[-c(23,36,39,52,57,60,63,67,71),]

## plot a
trj = paste0("trajectory", seq(1:91))
index = sample(91, 40)

d1_count = batch1[result2$from,]
d1_mean = rowMeans(d1_count[,1:500])
#d1_mean = d1_count[,10]
d1_total = data.frame(cell = rownames(d1_count), expression = d1_mean, cell_type = d1_count$type, time = d1_count$batch)
d1_total$trajectory = trj
d1_total = d1_total[index,]

d2_count = batch2[result2$to,]
d2_mean = rowMeans(d2_count[,1:500])
#d2_mean = d2_count[,10]
d2_total = data.frame(cell = rownames(d2_count), expression = d2_mean, cell_type = d2_count$type, time = d2_count$batch)
d2_total$trajectory = trj
d2_total = d2_total[index,]

total = rbind(d1_total, d2_total)
total$cell_type = gsub("Group", "Type", total$cell_type)
total$time = gsub("Batch", "Time", total$time)

gra_a = ggplot(total, aes(x=time, y=expression, color=factor(cell_type))) + 
  geom_point()

gra_a

gra_a = ggplot(total, aes(x=time, y=expression)) + 
  geom_point(data = total, aes(x=time, y=expression), size = 2, color = "purple") + 
  geom_line(data = total, aes(x = time, y = expression, group = trajectory), size = 0.8, color = "red")

gra_a[["labels"]][["colour"]] = "Cell Type"
gra_a[["labels"]][["y"]] = "Expression"
gra_a[["labels"]][["x"]] = "Time Points"
gra_a + theme_bw()


# umap
library(Seurat)
library(SeuratObject)

object <- CreateSeuratObject(counts = data, project = "simulation", 
                             min.cells = 3, min.features = 10)
object

object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 400)
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


mucl = merge(total, countUMAP, by="row.names")

gra_b = ggplot(mucl, aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(data = mucl, aes(x=UMAP_1, y=UMAP_2, color = factor(cell_type)), size = 2) + 
  geom_line(data = mucl, aes(x=UMAP_1, y=UMAP_2, group = trajectory), size = 0.1)

gra_b[["labels"]][["colour"]] = "Cell Type"
gra_b[["labels"]][["y"]] = "UMAP2"
gra_b[["labels"]][["x"]] = "UMAP1"
gra_b + theme_bw()


# all_gra = ggplot() + 
#   geom_point(data = points,aes(x = id,y = value), size=2, color = "purple") + 
#   geom_line(data = plot1,aes(x = id,y = value,group = variable), size=0.8, color = "red") + 
#   geom_line(data = plot2,aes(x = x,y = y), size=1, color = "blue") + 
#   geom_line(data = plot3,aes(x = x,y = y), size=1, color = "black") + 
#   geom_line(data = plot4,aes(x = id,y = value,group = variable), size=0.8, color = "red") + 
#   geom_line(data = plot5,aes(x = x,y = y), size=1, color = "blue") + 
#   geom_line(data = plot6,aes(x = x,y = y), size=1, color = "black") + xlab("Time Points") + 
#   ylab("Expression") + ggtitle("Individual Cell Trajectories")
# 
# all_gra = all_gra + theme_bw()
# all_gra

##### Compare three method ########
rep = seq(1:4)
gene = seq(100, 500, by=20)

n = 700
accuracy_ot = matrix(0,nrow = 4, ncol = 21)
accuracy_pearson = matrix(0,nrow = 4, ncol = 21)
accuracy_eucl = matrix(0,nrow = 4, ncol = 21)
time_ot = matrix(0,nrow = 4, ncol = 21)
time_pearson = matrix(0,nrow = 4, ncol = 21)
time_eucl = matrix(0,nrow = 4, ncol = 21)

for (time in rep){
  set.seed(20*time)
  for (ngenes in 1:length(gene)){
    nGenes = gene[ngenes]
    params <- newSplatParams()
    params <- setParams(params, update = list(batchCells = c((n/2), (n/2)), nGenes = nGenes, 
                                              group.prob = c(0.20, 0.20, 0.20, 0.20, 0.20), seed = 20*time))
    sim <- splatSimulate(params = params, method = "groups", verbose = FALSE)
    sim <- logNormCounts(sim)
    data = logcounts(sim)
    batch = sim$Batch
    type = as.character(sim$Group)
    
    count = data.frame(t(data))
    count$batch = batch
    count$type = type
    
    batch1 = count[count$batch == "Batch1",]
    batch2 = count[count$batch == "Batch2",]
    
    # Pearson correlation
    start.time <- Sys.time()
    cellaim = c()
    for (nn in 1: nrow(batch1)){
      corr = c()
      temp = as.numeric(batch1[nn,c(1:nGenes)])
      bb2 = t(as.matrix(batch2[,c(1:nGenes)]))
      b = apply(bb2, 2, function(x) {
        cor.test(temp, x, method = "pearson")
      })
      
      res <- sapply(b, "[[", "estimate")
      corr = res
      corr[cellaim] = 0
      cellaim = c(cellaim, which.max(corr))
    }
    sub = batch2[cellaim,]
    data1 = data.frame(cell1 = seq(1:(n/2)), type1 = batch1$type, cell2 = cellaim, type2 = sub$type)
    coun = 0
    for (i in 1:nrow(data1)){
      if (data1[i,]$type1 == data1[i,]$type2){
        coun = coun +1
      }
      
    }
    acc = coun/(n/2)
    accuracy_pearson[time,ngenes] = acc
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time_pearson[time,ngenes] = time.taken
    
    ##### euclidean distance
    start.time <- Sys.time()
    cellaim = c()
    for (nn in 1: nrow(batch1)){
      corr = c()
      temp = as.numeric(batch1[nn,c(1:nGenes)])
      bb2 = as.matrix(batch2[,c(1:nGenes)])
      res = apply(bb2,1,function(x) dist(rbind(temp,x),method = "euclidean"))
      corr = res
      corr[cellaim] = 9999
      cellaim = c(cellaim, which.min(corr))
    }
    
    sub = batch2[cellaim,]
    
    data1 = data.frame(cell1 = seq(1:(n/2)), type1 = batch1$type, cell2 = cellaim, type2 = sub$type)
    
    coun = 0
    for (i in 1:nrow(data1)){
      if (data1[i,]$type1 == data1[i,]$type2){
        coun = coun +1
      }
      
    }
    acc_eucl = coun/(n/2)
    accuracy_eucl[time,ngenes] = acc_eucl
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time_eucl[time,ngenes] = time.taken
    
    # optimal transport using EMD distance
    start.time <- Sys.time()
    Xs = data.matrix(batch1[1:(length(batch1)-2)])
    Xt = data.matrix(batch2[1:(length(batch2)-2)])
    
    ot_emd = ot$da$EMDTransport()
    ot_emd$fit(Xs=Xs, Xt=Xt)
    coupling = ot_emd$coupling_
    
    # connecting
    index = apply(coupling, 1, which.max)
    res = data.frame(from = seq(1:(n/2)), to = index)
    
    result = res
    result$b1 = batch1$type
    
    re = batch2[c(res$to),]
    result$b2 = re$type
    
    coun = 0
    for (i in 1:nrow(res)){
      if (result[i,]$b1 == result[i,]$b2){
        coun = coun +1
      }
      
    }
  
    acc_ot = coun/(n/2)
    accuracy_ot[time,ngenes] = acc_ot
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time_ot[time,ngenes] = time.taken
    
    print(paste0("done time ",time," and gene number ", nGenes))
  }
}

mean_pearson = colMeans(accuracy_pearson)
mean_eucl = colMeans(accuracy_eucl)
mean_ot = colMeans(accuracy_ot)

cell_700 = data.frame(Pearson = mean_pearson, Euclidean = mean_eucl, OT = mean_ot)
cell_700_table = data.frame(rbind(mean_pearson, mean_eucl, mean_ot))
colnames(cell_700_table) = gene

cell_700_use = cell_700_table[c("100","200","300","400","500")]
#xtable(cell_200_table)

write.csv(cell_200_table, file = "cell_200.csv")
write.csv(cell_300_table, file = "cell_300.csv")
write.csv(cell_400_table, file = "cell_400.csv")
write.csv(cell_500_table, file = "cell_500.csv")
write.csv(cell_600_table, file = "cell_600.csv")
write.csv(cell_700_table, file = "cell_700.csv")

# mean_pearson = c(0.7266667, 0.7166667, 0.7075000, 0.7166667, 0.7175000, 0.7225000,
#                  0.7258333, 0.7541667, 0.7250000, 0.7091667, 0.7408333, 0.7333333, 
#                  0.7241667, 0.7450000, 0.7450000, 0.7233333, 0.7475000, 0.7316667, 
#                  0.7133333, 0.7258333, 0.7541667)
# mean_eucl = c(0.7333333, 0.7216667, 0.7225000, 0.7125000, 0.7200000, 0.7133333, 0.7175000, 
#               0.7483333, 0.6958333, 0.7233333, 0.7333333, 0.7175000, 0.7066667, 0.7600000, 
#               0.7350000, 0.7375000, 0.7533333, 0.7450000, 0.7200000, 0.7208333, 0.7550000)
# mean_ot = c(0.7700000, 0.7875000, 0.7900000, 0.7966667, 0.7975000, 0.8108333, 0.7900000, 0.8325000, 
#             0.8458333, 0.8341667, 0.8400000, 0.8366667, 0.8466667, 0.8516667, 0.8325000, 0.8333333, 
#             0.8333333, 0.8375000, 0.8641667, 0.8583333, 0.8691667)
# 
# cell_700_table = data.frame(rbind(mean_pearson, mean_eucl, mean_ot))
# colnames(cell_700_table) = gene
# write.csv(cell_700_table, file = "cell_700.csv")

### Plot ##

rplot1 = data.frame(x = seq(100, 500, by=20), y = mean_pearson)
rplot2 = data.frame(x = seq(100, 500, by=20), y = mean_eucl)
rplot3 = data.frame(x = seq(100, 500, by=20), y = mean_ot)

p3 = ggplot() + 
  geom_line(data = rplot1, aes(x = x, y = y, colour = "green")) +
  geom_line(data = rplot2, aes(x = x, y = y, colour = "red")) +
  geom_line(data = rplot3, aes(x = x, y = y, colour = "black")) +
  xlab('Number of Genes per Cell') + scale_color_discrete(name = "Legend", labels = c("CellST", "Euclidean", "Pearson")) +
  ylim(0.25,1) +
  ylab('Accuracy') + ggtitle("700 Cells")

p3 = p3 + theme_bw()
p3


#### time

rplot1 = data.frame(x = seq(100, 500, by=20), y = time_pear)
rplot2 = data.frame(x = seq(100, 500, by=20), y = time_ot)
rplot3 = data.frame(x = seq(100, 500, by=20), y = time_eucl)
#rplot4 = data.frame(x = seq(80, 500, by=20), y = accuracy_ken)

tim3 = ggplot() + 
  geom_line(data = rplot1, aes(x = x, y = y, colour = "green")) +
  geom_line(data = rplot2, aes(x = x, y = y, colour = "red")) +
  geom_line(data = rplot3, aes(x = x, y = y, colour = "black")) +
  #geom_line(data = rplot4, aes(x = x, y = y), color = "green") +
  xlab('Number of Genes per Cell') + scale_color_discrete(name = "Legend", labels = c("Euclidean", "Pearson", "CellOT")) +
  ylab('Time') + ggtitle("400 Cells Time")

tim3 = tim3 + theme_bw()
tim3

library(ggpubr)
ggarrange(p, p2, p3, p4,p5, p6, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
ggarrange(tim, tim2, tim3, tim4,tim5, tim6, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")


# sim <- splatSimulate(batchCells = c((n/2), (n/2)), nGenes = nGenes, group.prob = c(0.20, 0.20, 0.20, 0.20, 0.20),
#                      method = "groups", verbose = FALSE, seed = 20*time)

# params <- newSplatParams()
# params <- setParams(params, update = list(nGenes = 8000, mean.rate = 0.5))
# params <- setParams(params, update = list(batchCells = c((n/2), (n/2)), nGenes = nGenes, group.prob = c(0.20, 0.20, 0.20, 0.20, 0.20)
#                                           , seed = 20*time))
# 
# params
# 
# sim <- splatSimulate(params = params, method = "groups", verbose = FALSE)



#### Plot change y value
library(ggplot2)
rownames(cell_600) = cell_600$X; cell_600 = cell_600[,-1]
data = data.frame(t(cell_600))

mean_pearson = data$mean_pearson
mean_eucl = data$mean_eucl
mean_ot = data$mean_ot
#mean_pearson = subset(cell_200, rownames(cell_200) %in% "mean_pearson")
#mean_eucl = subset(cell_200, rownames(cell_200) %in% "mean_pearson")

rplot1 = data.frame(x = seq(100, 500, by=20), y = mean_pearson)
rplot2 = data.frame(x = seq(100, 500, by=20), y = mean_eucl)
rplot3 = data.frame(x = seq(100, 500, by=20), y = mean_ot)

p3 = ggplot() + 
  geom_line(data = rplot1, aes(x = x, y = y, colour = "green")) +
  geom_line(data = rplot2, aes(x = x, y = y, colour = "red")) +
  geom_line(data = rplot3, aes(x = x, y = y, colour = "black")) +
  xlab('Number of Genes per Cell') + scale_color_discrete(name = "Legend", labels = c("CellST", "Euclidean", "Pearson")) +
  ylim(0.5,0.9) +
  ylab('Accuracy') + ggtitle("600 Cells")

p3 = p3 + theme_bw()
p3


# ##### Compare OT with Pearson correlation ########
# n = 400
# accuracy_ot = matrix(0,nrow = 5, ncol = 21)
# accuracy_spear = c()
# accuracy_eucl = c()
# time_ot = c()
# time_pear = c()
# time_eucl = c()
# 
# for (ngenes in seq(100, 500, by=20)){
#   nGenes = ngenes
#   set.seed(Sys.time())
#   sim <- splatSimulate(batchCells = c((n/2), (n/2)), nGenes = nGenes, group.prob = c(0.20, 0.20, 0.20, 0.20, 0.20),
#                        method = "groups", verbose = FALSE)
#   sim <- logNormCounts(sim)
#   data = logcounts(sim)
#   batch = sim$Batch
#   type = as.character(sim$Group)
#   
#   count = data.frame(t(data))
#   count$batch = batch
#   count$type = type
#   
#   batch1 = count[count$batch == "Batch1",]
#   batch2 = count[count$batch == "Batch2",]
#   
#   # Pearson correlation
#   start.time <- Sys.time()
#   cellaim = c()
#   for (nn in 1: nrow(batch1)){
#     corr = c()
#     temp = as.numeric(batch1[nn,c(1:nGenes)])
#     bb2 = t(as.matrix(batch2[,c(1:nGenes)]))
#     b = apply(bb2, 2, function(x) {
#       cor.test(temp, x, method = "pearson")
#     })
#     
#     res <- sapply(b, "[[", "estimate")
#     corr = res
#     corr[cellaim] = 0
#     cellaim = c(cellaim, which.max(corr))
#   }
#   
#   sub = batch2[cellaim,]
#   
#   data1 = data.frame(cell1 = seq(1:(n/2)), type1 = batch1$type, cell2 = cellaim, type2 = sub$type)
#   
#   coun = 0
#   for (i in 1:nrow(data1)){
#     if (data1[i,]$type1 == data1[i,]$type2){
#       coun = coun +1
#     }
#     
#   }
#   print(ngenes)
#   acc = coun/(n/2)
#   accuracy = c(accuracy, acc)
#   end.time <- Sys.time()
#   time.taken <- end.time - start.time
#   time_pear = c(time_pear, time.taken)
#   
#   ##### euclidean distance
#   start.time <- Sys.time()
#   cellaim = c()
#   for (nn in 1: nrow(batch1)){
#     corr = c()
#     temp = as.numeric(batch1[nn,c(1:nGenes)])
#     bb2 = as.matrix(batch2[,c(1:nGenes)])
#     res = apply(bb2,1,function(x) dist(rbind(temp,x),method = "euclidean"))
#     corr = res
#     corr[cellaim] = 9999
#     cellaim = c(cellaim, which.min(corr))
#   }
#   
#   sub = batch2[cellaim,]
#   
#   data1 = data.frame(cell1 = seq(1:(n/2)), type1 = batch1$type, cell2 = cellaim, type2 = sub$type)
#   
#   coun = 0
#   for (i in 1:nrow(data1)){
#     if (data1[i,]$type1 == data1[i,]$type2){
#       coun = coun +1
#     }
#     
#   }
#   print(ngenes)
#   acc_eucl = coun/(n/2)
#   accuracy_eucl = c(accuracy_eucl, acc_eucl)
#   end.time <- Sys.time()
#   time.taken <- end.time - start.time
#   time_eucl = c(time_eucl, time.taken)
#   
#   # optimal transport using EMD distance
#   start.time <- Sys.time()
#   Xs = data.matrix(batch1[1:(length(batch1)-2)])
#   Xt = data.matrix(batch2[1:(length(batch2)-2)])
# 
#   ot_emd = ot$da$EMDTransport()
#   ot_emd$fit(Xs=Xs, Xt=Xt)
#   coupling = ot_emd$coupling_
#   
#   # connecting
#   index = apply(coupling, 1, which.max)
#   res = data.frame(from = seq(1:(n/2)), to = index)
#   
#   result = res
#   result$b1 = batch1$type
#   
#   re = batch2[c(res$to),]
#   result$b2 = re$type
#   
#   coun = 0
#   for (i in 1:nrow(res)){
#     if (result[i,]$b1 == result[i,]$b2){
#       coun = coun +1
#     }
#     
#   }
#   
#   print(ngenes)
#   acc_ot = coun/(n/2)
#   accuracy_ot = c(accuracy_ot, acc_ot)
#   end.time <- Sys.time()
#   time.taken <- end.time - start.time
#   time_ot = c(time_ot,time.taken)
# }


