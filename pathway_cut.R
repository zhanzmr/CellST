library(ggplot2)
library(dplyr)
library(gtools)
library(tidyverse)
library(pheatmap)
library(gss)
library(reshape2)
library(reticulate)

conda_list()
use_condaenv(condaenv = 'ot', required = TRUE)
# import python module
np = import("numpy")
ot = import("ot")

set.seed(234987)
setwd("cell_cycle/Cell_cycle_smoothed/")

n = 1956
file <- mixedsort(list.files(path = "Cell_cycle_data/", full.names = T))

total = data.frame()
for (i in 1:17){
  temp = read.csv(file[i])
  sample_cell = sample_n(temp,n)
  total = rbind(total,sample_cell)
  print(i)
}

total$time = total$time + 1


#### Optimal Transport ####

couple = list()
for ( jj in 1:16){
  temp1 = total[total$time == jj,]
  temp2 = total[total$time == jj+1,]
  Xs = data.matrix(temp1[2:(length(temp1)-1)])
  Xt = data.matrix(temp2[2:(length(temp2)-1)])
  # EMD Transport
  ot_emd = ot$da$EMDTransport()
  ot_emd$fit(Xs=Xs, Xt=Xt)
  coupling = ot_emd$coupling_
  couple[[jj]] = coupling
  print(jj)
} 

### connecting cells
temp1 = total[total$time == 1,]
temp2 = total[total$time == 2,]
temp_coup = couple[[1]]

cell_name_traj = data.frame(time= temp1$index )
index = apply(temp_coup, 1, which.max)
temp2 = temp2[c(index),]
col = temp2$index
cell_name_traj = cbind(cell_name_traj,col)

temp1 = temp2


for (n in 2:16){
  temp2 = total[total$time == n+1,]
  temp_coup = couple[[n]]
  temp_coup = temp_coup[c(index),]
  
  index = apply(temp_coup, 1, which.max)
  temp2 = temp2[c(index),]
  col = temp2$index
  cell_name_traj = cbind(cell_name_traj,col)
  
  temp1 = temp2
  print(n)
}

x = seq(1:17)
col = paste0("time",x)
colnames(cell_name_traj) = col

####################################################
n_col = ncol(total)
### output for WOT
day = data.frame(day = total$time)
rownames(day) = total$index
pathway_index = total$index

#write(pathway_index,file = paste0(pathway_name,".txt"))
#write.csv(total[,-c(n_col)], file = paste0(pathway_name,"_count.csv"), row.names = F)
#write.csv(day, file = paste0(pathway_name,"_day.csv"))




### Get individual cell trajectories expression
expression = list()
#somePDFPath = paste0(pathway_name,"_heatmap.pdf")
#pdf(file=somePDFPath)

for (row in 1:1956){
  
  row_num = row
  name = as.character(cell_name_traj[row_num,1])
  value = as.numeric(t(total[total$index == name,]))
  value = value[c(2:(length(value)-1))]
  gene_value = as.data.frame(value)
  
  for (nn in 2:17){
    name = as.character(cell_name_traj[row_num,nn])
    value = as.numeric(t(total[total$index == name,]))
    value = value[c(2:(length(value)-1))]
    value = as.data.frame(value)
    gene_value = cbind(gene_value,value)
  }
  
  gga = as.matrix(gene_value)
  colnames(gga) = paste0("time_", seq(1:17))
  expression[[row]] = gga
  #pheatmap(gga, kmeans_k = 35, cluster_cols = F, show_rownames = F)
  print(row)
}
#dev.off()

#### ssanova 
### read wot result ###

wot_compare_trends <- read.delim("Cell_cycle_trends.txt")

geneset_mclu_tran = t(as.matrix(wot_compare_trends[,-1]))
data <- as.data.frame(geneset_mclu_tran)
col_mean = colMeans(data, na.rm = FALSE, dims = 1)

### plot CellOT and wot

pred_list = list()
plot_list = list()

somePDFPath = "cell_cycle_CellOT_compare.pdf"
pdf(file=somePDFPath)
for (nn in 1:length(expression)){
  data = expression[[nn]]
  geneset_mclu_tran = t(data)
  data <- as.data.frame(geneset_mclu_tran)
  data$id <- 1:nrow(data) 
  plot_data <- melt(data,id.var="id")
  yy <- plot_data$value
  cc <- plot_data$id
  id <- as.factor(plot_data$variable)
  
  fit.random <- ssanova(yy~cc)
  x = data.frame(cc = seq(1, 17, length = 300))
  x2 = data.frame(cc = seq(1:17))
  pred2 <- predict(fit.random, newdata = x, se.fit = T )
  pred_point <- predict(fit.random, newdata = x2, se.fit = T )
  pred_list[[nn]] = pred2$fit
  plot_list[[nn]] = pred_point$fit
  plot(x$cc,pred2$fit, type = "l", xlab = "Time_points", ylab = "Gene Set Expression", main = "CellOT compare")
  lines(x$cc,pred2$fit, col= "red", lwd = 2)
  lines(seq(1:17),col_mean, col= "Blue", lwd = 2)
  print(nn)
}
dev.off()


#### plot all curves
#index = c(4,9,23,28,32,44,46,47,58,65,78)
index = c(4,9,23,28)
subset_pred = pred_list[index]
selected_cell = expression[index]

real_cell_mean = data.frame()
for (ii in 1:length(selected_cell)){
  temp = selected_cell[[ii]]
  res = colMeans(temp)
  real_cell_mean = rbind(real_cell_mean, res)
}
colnames(real_cell_mean) = paste0("time_", seq(1:17))
dd2 <- data.frame(t(real_cell_mean))
dd2$id <- 1:nrow(dd2) 
dot1_real <- melt(dd2,id.var="id")


dd  <-  data.frame(matrix(unlist(subset_pred), nrow = length(subset_pred), byrow = T))
dd <- data.frame(t(dd))
dd$id <- 1:nrow(dd) 
plot1 <- melt(dd,id.var="id")
plot1$id = rep(x$cc,4)


plot2 <- data.frame(x = seq(1:17), y = col_mean)

subset_point = plot_list[index]
dd2  <-  data.frame(matrix(unlist(subset_point), nrow = length(subset_point), byrow = T))
dd2 <- data.frame(t(dd2))
dd2$id <- 1:nrow(dd2) 
dot1 <- melt(dd2,id.var="id")

#dot1 <- plot1[sample(nrow(plot1), 500), ]

# compare = ggplot() + 
#   geom_line(data = plot1,aes(x = id,y = value,group = variable, colour = "red"), size=0.3) + 
#   geom_line(data = plot2,aes(x = x,y = y, colour = "blue"), size=1) + 
#   geom_point(data = dot1,aes(x = id,y = value,group = variable, colour = "black"), size=0.6) +
#   ylim(-0.15, 1.3) + 
#   labs(title = paste0("Comparsion"), x = "Time", y = "Expression", color = "Legend") +
#   scale_color_manual(labels = c("Cell Expression","WOT","CellST"), values = c("blue", "red","black"))
# 
# compare = compare + theme_bw()
# compare

compare = ggplot() + 
  geom_point(data = dot1,aes(x = id,y = value, group = variable), size=2, color = "purple") + 
  geom_point(data = dot1_real,aes(x = id,y = value), size=2, color = "grey") +
  geom_line(data = plot1,aes(x = id,y = value,group = variable), size=0.8, color = "red") + 
  geom_line(data = plot2,aes(x = x,y = y), size=1, color = "blue") + 
  xlab("Time Points") + ylab("Expression") + ylim(-0.15, 1.3)

compare = compare + theme_bw() + theme(legend.position="none")
compare

gra1 = gra1 + theme(legend.position="none")

library(ggpubr)
ggarrange(gra1, compare, ncol = 2, nrow = 1)


