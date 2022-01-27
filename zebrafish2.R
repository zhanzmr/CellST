library(ndtv)
library(gss)
library(intergraph)
library(networkDynamic)
library(igraph) 

setwd("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development")
load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/mucl_all.RData")
load("~/Documents/big_data/single-cell/new_dy_network_07082021/zebrafish_development/zebrafish1_total.RData")
gene_name = colnames(count)

mucl = mucl_list[[7]]

network = list() 

for (ii in 1:101){
  
  temp = mucl[,ii]
  class(temp) <- 'dist'
  attr(temp,'Size') <- length(gene_name)
  temp_corr = as.matrix(temp)
  colnames(temp_corr) = gene_name
  rownames(temp_corr) = gene_name
  
  
  temp_corr[temp_corr < 0.7] <- 0
  nn = network(temp_corr, directed = FALSE)
  network[[ii]] = nn
  print(ii)
}
names(network) = seq(1,101)


# dyn_network<-networkDynamic(network.list=network, vertex.pid="vertex.names")
# 
# compute.animation(dyn_network, animation.mode = "kamadakawai",
#                   slice.par=list(start=0, end=100, interval=1,
#                                  aggregate.dur=1, rule='any'))
# render.d3movie(dyn_network,displaylabels=T)

#file:///private/var/folders/lm/5f3_tpnj56z6wtvtsx85mznr0000gn/T/RtmpbQQzL3/filea3ec62951037.html

########### Gene Community Detection #######################
## loop for community Detection
group_louvain_label = list()

for (ii in 1:101){
  net = network[[ii]]
  g2 <- asIgraph(net)
  lc <- cluster_louvain(g2)
  m = as.numeric(membership(lc))
  data = data.frame(gene = gene_name, group = m)
  group_louvain_label[[ii]] = data
  print(ii)
}
names(group_louvain_label) = seq(1,101)

somePDFPath = "cluster7.pdf"
pdf(file=somePDFPath) 
for (ii in 1:101){
  net = network[[ii]]
  g2 <- asIgraph(net)
  
  lc <- cluster_louvain(g2)
  cc = lc$membership
  mm = membership(lc)
  mm[10] <- "red"
  plot(lc, g2, vertex.label = NA, col = factor(mm))
  print(ii)
}
dev.off()

#V(g2)$vertex.names
net = network[[1]]
g2 <- asIgraph(net)

lc <- cluster_louvain(g2)
mm = membership(lc)
plot(lc, g2, vertex.label = NA)
#V(g2)$vertex.names

test = V(g2)$vertex.names


######### Connecting community 
########### Loop to connect communities

community_connect = list()

for (kk in 1:100){
  temp_net = group_louvain_label[[kk]]
  temp2_net = group_louvain_label[[kk+1]]
  n = nlevels(factor(temp_net$group))
  data_connect = data.frame(source = NULL, target = NULL)
  
  for (ii in 1:n){
    num = c()
    subset1 = temp_net[temp_net$group == ii, ]
    v1 = subset1$gene
    n2 = nlevels(factor(temp2_net$group))
    for (jj in 1:n2){
      subset2 = temp2_net[temp2_net$group == jj, ]
      v2 = subset2$gene
      inter = intersect(v1,v2)
      if (length(inter) > 1){
        num = c(num,length(inter))
      } else {
        num = c(num, -1)
      }
    }
    max_gene = which.max(num)
    if (num[max_gene] > 1){
      pp = data.frame(source = ii, target = max_gene)
      data_connect = rbind(data_connect, pp)
    }
  }
  community_connect[[kk]] = data_connect
  print(kk)
}



##### gene expression for average connected communities

# temp1 = community_connect[[1]]
# colnames(temp1) = c("use", "connect")
# temp2 = community_connect[[2]]
# colnames(temp2) = c("connect", "use")
# 
# community_total = data.frame(time= temp1$use )
# index = temp1$connect
# temp2 = temp2[temp2$connect == index, ]
# col = temp2$index
# community_total = cbind(cell_name_traj,col)
# 
# temp1 = temp2
# 
# 
# for (n in 2:(nTime-1)){
#   temp2 = count[count$time == n+1,]
#   temp_coup = couple[[n]]
#   temp_coup = temp_coup[c(index),]
#   
#   index = apply(temp_coup, 1, which.max)
#   temp2 = temp2[c(index),]
#   col = temp2$index
#   cell_name_traj = cbind(cell_name_traj,col)
#   
#   temp1 = temp2
# }


community_total = community_connect[[1]]

for (ii in 2:100 ) {
  temp_use = community_connect[[ii]]
  index = community_total[,ii]
  
  for (jj in 1:length(index)) {
    num = index[jj]
    temp = temp_use[temp_use$source == num, ]
    
    if (nrow(temp) > 0){
      community_total[jj, ii+1] = temp$target
    } else {
      community_total[jj, ii+1] = NA
    }
  }
}

names(gene_expression) = gene_name


com_label_list = as.integer(community_total[1,])

group = com_label_list[[ii]]
lou = group_louvain_label[[ii]]
intersect_gene = lou[lou$group == group,]$gene

for (ii in 2:101){
  group = com_label_list[[ii]]
  lou = group_louvain_label[[ii]]
  gene = lou[lou$group == group,]$gene
  intersect_gene = intersect(intersect_gene, gene)
  
}


smoothed_gene_curves = data.frame()
for (ii in 1:96){
  temp = gene_expression[[ii]]
  ss = colMeans(temp)
  cc = seq(1:17)
  fit.random <- ssanova(ss~cc, alpha=0.5)
  x = data.frame(cc = seq(1, 17, length = 101))
  pred <- predict(fit.random, newdata = x, se.fit = T )
  smoothed_gene_curves = rbind(smoothed_gene_curves, pred$fit)
}

colnames(smoothed_gene_curves) = seq(1,101)
rownames(smoothed_gene_curves) = gene_name

matplot(t(smoothed_gene_curves), type = "l", ylab = "Expression", xlab = "Time", main = "Average expression for all genes")



com_label_list = as.integer(community_total[3,])
community_expression = c()
for (ii in 1:101){
  group = com_label_list[[ii]]
  lou = group_louvain_label[[ii]]
  gene = lou[lou$group == group,]$gene
  exp = mean(smoothed_gene_curves[gene,1])
  community_expression = c(community_expression, exp)
  
}

plot(seq(1,101), community_expression, type = "l", xlab = "Time", main = "Gene Community Average Expression")

write.table(community_expression, "myFile.txt", row.names = F, col.names = F)
myFile <- read.table("~/Documents/big_data/single-cell/new_dy_network_07082021/mouse_embryogenesis/pathway/cell_cycle/myFile.txt", quote="\"", comment.char="")

test = myFile$V1
plot(seq(1,101), test, type = "l", xlab = "Time", main = "Gene Community Average Expression")

df = data.frame(Time = seq(0,17,length.out = 101), Expression = test)

ggplot(data=df, aes(x=Time, y=Expression)) +
  geom_line()+
  geom_point() + theme_bw() + ylim(0,0.5) + ggtitle("Gene Community Average Expression")

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Sugar maple', 'White ash', 'Black walnut',
                            'Red oak', 'Eastern hemlock'), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c('orange', 'red', 'green', 'blue', 'purple'))
mtext("Species", at=0.2, cex=2)

