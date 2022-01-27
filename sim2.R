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

# set.seed(34895)
# load("~/Documents/big_data/single_cell/new_dy_network_07082021/sim2/sim2_for_tscan.RData")
# import python module
# np = import("numpy")
# ot = import("ot")

## Two Pathways
simple_mod <- suppressMessages(inSilicoCellModel(initialNum=80, runTime=72,
                                                 density=0.1, outputIncrement=6, randSeed=277))

# first pathway
mitosisGeneNames <- paste("m_", seq(1:100), sep="")
mitosisExpression <- function(model, cell, time)
{
  ifelse(getCellPhase(model, time, cell) == "M", 1, 0)
}

pwyMitosis <- new("Pathway", genes=mitosisGeneNames,
                  expressionScale=mitosisExpression)


# second pathway
contactInhibitionGeneNames <- paste("ci_", seq(1:100), sep="")
contactInhibitionExpression <- function(model, cell, time)
{
  getLocalDensity(model, time, cell, 3.3)
}
pwyContactInhibition <- new("Pathway", genes=contactInhibitionGeneNames,
                            expressionScale=contactInhibitionExpression)

# create simulated data set
allGenes <- c(mitosisGeneNames, contactInhibitionGeneNames)
geneMeans <- 2 + rexp(length(allGenes), 1/20)
data <- t(pmax(sapply(geneMeans, rnorm, n=25, sd=2), 0))
rownames(data) <- allGenes

# calibrate pathways
pwyMitosis <- calibratePathway(pwyMitosis, data)
pwyContactInhibition <- calibratePathway(pwyContactInhibition, data)

# plot pathway
params <- new("GeneExpressionParams")
params@randSeed <- 288 # control this for reporducibility
params@nCells <- 30 # sample 30 cells at each time point to measure activity
params@sampleFreq <- 6 # measure activity every 6 hours

pwyMitosis@expressionScale = function(model, cell, time)
{
  window <- c(max(time - 2, 0), min(time + 2, model@runTime))
  a1 <- getAxisLength(model, window[1], cell)
  a2 <- getAxisLength(model, window[2], cell)
  if (is.na(a1)) a1 <- 0 # in case cell was just born
  return(ifelse(a2 < a1, 1, 0))
}
pwys <- c(pwyMitosis, pwyContactInhibition)
pwyActivity <- inSilicoGeneExpression(simple_mod, pwys, params)$pathways

# mitosis
plot(seq(0,72,6), pwyActivity[[1]], type="l", col="orange", ylim=c(0,1))
# contact inhibition
lines(seq(0,72,6), pwyActivity[[2]], col="blue")

# pathway activities (microarray)

params@randSeed <- 288 
params@nCells <- 80 # sample 30 cells at each time point to measure activity
params@sampleFreq <- 6 # measure activity every 6 hours
params@RNAseq <- FALSE # generate microarray data
params@singleCell <- FALSE # generate bulk data
params@perError <- 0.1 # parameter for simulated noise

pwys <- c(pwyMitosis, pwyContactInhibition)
ge <- inSilicoGeneExpression(simple_mod, pwys, params)$expression

colnames(ge) = paste0("time_", seq(1:13))
pheatmap(ge, kmeans_k = 35, cluster_cols = F)

### RNA-seq data
params@randSeed <- 288
params@nCells <- 80 # sample 30 cells at each time point to measure activity
params@sampleFreq <- 6 # measure activity every 6 hours
params@RNAseq <- TRUE # 
params@singleCell <- FALSE # generate bulk data
params@perError <- 0 # parameter for simulated noise

pwys <- c(pwyMitosis, pwyContactInhibition)
ge_rna <- inSilicoGeneExpression(simple_mod, pwys, params)$expression

### Show RNA-seq
#matplot(t(ge_rna), type="l", lty=1, lwd=2, xlab = "Time_points", 
#        ylab = "Gene Expression", main = "Gene Trajectories")

ce <- SingleCellExperiment(list(counts=ge_rna))
ce <- logNormCounts(ce)

ge_log = logcounts(ce)

colnames(ge_log) = paste0("time_", seq(1:13))

gra_true_rna = pheatmap(ge_log, kmeans_k = 35, cluster_cols = F,show_rownames = F)
gra_true_rna

### Comupute Mean curves ####
## first pathway
ge_mitosisGene = ge_log[mitosisGeneNames,]
geneset_mclu_tran = t(ge_mitosisGene)
data <- as.data.frame(geneset_mclu_tran)
data$id <- 1:nrow(data) 
plot_data <- melt(data,id.var="id")

yy <- plot_data$value
cc <- plot_data$id
id <- as.factor(plot_data$variable)

fit.random <- ssanova(yy~cc)
x = data.frame(cc = seq(1, 13, length = 300))
pred_mitosis <- predict(fit.random, newdata = x, se.fit = T )

plot(x$cc,pred_mitosis$fit, type = "l", xlab = "Time_points", 
     ylab = "Gene Set Expression", main = "Mitosis Pathway Mean Curve")

#### second pathway
ge_contact = ge_log[contactInhibitionGeneNames,]
geneset_mclu_tran = t(ge_contact)
data <- as.data.frame(geneset_mclu_tran)
data$id <- 1:nrow(data)
plot_data <- melt(data,id.var="id")

yy <- plot_data$value
cc <- plot_data$id
id <- as.factor(plot_data$variable)

fit.random <- ssanova(yy~cc)
x = data.frame(cc = seq(1, 13, length = 300))
pred_contact <- predict(fit.random, newdata = x, se.fit = T )

plot(x$cc,pred_contact$fit, type = "l", xlab = "Time_points", 
     ylab = "Gene Set Expression", main = "Contact Pathway Mean Curve")

plot(x$cc,pred_contact$fit, type = "l", xlab = "Time_points", 
     ylab = "Gene Set Expression", main = "Simulated Two Pathways Mean Trajectory (Truth)", ylim=c(0.5,2))
lines(x$cc,pred_mitosis$fit, col= "red", lwd = 2)
lines(x$cc,pred_contact$fit, col= "blue", lwd = 2)


#### Optimal Transport
params@randSeed <- 288 
params@RNAseq <- TRUE
params@singleCell <- TRUE
params@dropoutPresent <- F
pwys <- c(pwyMitosis, pwyContactInhibition)
ge_single <- inSilicoGeneExpression(simple_mod, pwys, params)$expression


# log normal
sce <- SingleCellExperiment(list(counts=ge_single))
sce <- logNormCounts(sce)

Count = logcounts(sce)
count = data.frame(t(Count))
time = rep(1:13,each = 80)
count$time = time
count$index = rownames(count)

count_mitosisGene = count[,c(mitosisGeneNames,"time","index")]
count_contactInhibition = count[,c(contactInhibitionGeneNames,"time","index")]

### optimal transport
# first pathway
couple_mitosisGene = list()
for ( ii in 1:12){
  temp1 = count_mitosisGene[count_mitosisGene$time == ii,]
  temp2 = count_mitosisGene[count_mitosisGene$time == ii+1,]
  Xs = data.matrix(temp1[1:(length(temp1)-2)])
  Xt = data.matrix(temp2[1:(length(temp2)-2)])
  # EMD Transport
  ot_emd = ot$da$EMDTransport()
  ot_emd$fit(Xs=Xs, Xt=Xt)
  coupling = ot_emd$coupling_
  couple_mitosisGene[[ii]] = coupling
} 

### connecting cells
temp1 = count_mitosisGene[count_mitosisGene$time == 1,]
temp2 = count_mitosisGene[count_mitosisGene$time == 2,]
temp_coup = couple_mitosisGene[[1]]

traj_mitosisGene = data.frame(time= temp1$index )
index = apply(temp_coup, 1, which.max)
temp2 = temp2[c(index),]
col = temp2$index
traj_mitosisGene = cbind(traj_mitosisGene,col)

temp1 = temp2


for (n in 2:12){
  temp2 = count_mitosisGene[count_mitosisGene$time == n+1,]
  temp_coup = couple_mitosisGene[[n]]
  temp_coup = temp_coup[c(index),]
  
  index = apply(temp_coup, 1, which.max)
  temp2 = temp2[c(index),]
  col = temp2$index
  traj_mitosisGene = cbind(traj_mitosisGene,col)
  
  temp1 = temp2
  print(n)
}

## Second pathway

couple_contactInhibition = list()
for ( ii in 1:12){
  temp1 = count_contactInhibition[count_contactInhibition$time == ii,]
  temp2 = count_contactInhibition[count_contactInhibition$time == ii+1,]
  Xs = data.matrix(temp1[1:(length(temp1)-2)])
  Xt = data.matrix(temp2[1:(length(temp2)-2)])
  # EMD Transport
  ot_emd = ot$da$EMDTransport()
  ot_emd$fit(Xs=Xs, Xt=Xt)
  coupling = ot_emd$coupling_
  couple_contactInhibition[[ii]] = coupling
} 

### connecting cells
temp1 = count_contactInhibition[count_contactInhibition$time == 1,]
temp2 = count_contactInhibition[count_contactInhibition$time == 2,]
temp_coup = couple_contactInhibition[[1]]

traj_contactInhibition = data.frame(time= temp1$index )
index = apply(temp_coup, 1, which.max)
temp2 = temp2[c(index),]
col = temp2$index
traj_contactInhibition = cbind(traj_contactInhibition,col)

temp1 = temp2


for (n in 2:12){
  temp2 = count_contactInhibition[count_contactInhibition$time == n+1,]
  temp_coup = couple_contactInhibition[[n]]
  temp_coup = temp_coup[c(index),]
  
  index = apply(temp_coup, 1, which.max)
  temp2 = temp2[c(index),]
  col = temp2$index
  traj_contactInhibition = cbind(traj_contactInhibition,col)
  
  temp1 = temp2
  print(n)
}


## build individual trajectories
expression_mitosisGene = list()
expression_contactInhibition = list()

# collect gene expression
for (row in 1:80){
  
  # first pathway
  row_num = row
  name = as.character(traj_mitosisGene[row_num,1])
  value = as.numeric(t(count_mitosisGene[count_mitosisGene$index == name,]))
  value = value[c(1:(length(value)-2))]
  gene_value = as.data.frame(value)
  
  for (nn in 2:13){
    name = as.character(traj_mitosisGene[row_num,nn])
    value = as.numeric(t(count_mitosisGene[count_mitosisGene$index == name,]))
    value = value[c(1:(length(value)-2))]
    value = as.data.frame(value)
    gene_value = cbind(gene_value,value)
  }
  
  gga = as.matrix(gene_value)
  colnames(gga) = paste0("time_", seq(1:13))
  expression_mitosisGene[[row]] = gga
  
  # second pathway
  row_num = row
  name = as.character(traj_contactInhibition[row_num,1])
  value = as.numeric(t(count_contactInhibition[count_contactInhibition$index == name,]))
  value = value[c(1:(length(value)-2))]
  gene_value = as.data.frame(value)
  
  for (nn in 2:13){
    name = as.character(traj_contactInhibition[row_num,nn])
    value = as.numeric(t(count_contactInhibition[count_contactInhibition$index == name,]))
    value = value[c(1:(length(value)-2))]
    value = as.data.frame(value)
    gene_value = cbind(gene_value,value)
  }
  
  gga = as.matrix(gene_value)
  colnames(gga) = paste0("time_", seq(1:13))
  expression_contactInhibition[[row]] = gga
  
  print(row)
}

# ssanova mean curves

pred2_mitosis = list()
pred2_contactInhibition = list()

plot_mitosis = list()
plot_contactInhibition = list()

# First pathway
somePDFPath = "CellOT_compare_sim2_mitosis.pdf"
pdf(file=somePDFPath)
for (nn in 1:length(expression_mitosisGene)){
  
  data = expression_mitosisGene[[nn]]
  geneset_mclu_tran = t(data)
  data <- as.data.frame(geneset_mclu_tran)
  data$id <- 1:nrow(data) 
  plot_data <- melt(data,id.var="id")
  yy <- plot_data$value
  cc <- plot_data$id
  id <- as.factor(plot_data$variable)
  
  fit.random <- ssanova(yy~cc)
  x = data.frame(cc = seq(1, 13, length = 300))
  x2 = data.frame(cc = seq(1:13))
  pred2 <- predict(fit.random, newdata = x, se.fit = T )
  pred_point <- predict(fit.random, newdata = x2, se.fit = T )
  pred2_mitosis[[nn]] = pred2$fit
  plot_mitosis[[nn]] = pred_point$fit
  
  plot(x$cc,pred2$fit, type = "l", xlab = "Time_points", ylab = "Gene Set Expression", main = "Mitosis Mean Curve", ylim=c(0.5,2))
  lines(x$cc,pred2$fit, col= "red", lwd = 2)
  lines(x$cc,pred_mitosis$fit, col= "black", lwd = 2)
  lines(x$cc,pred_mitosis$fit+1.96*pred_mitosis$se,col=5)
  lines(x$cc,pred_mitosis$fit-1.96*pred_mitosis$se,col=5)
  print(nn)
}
dev.off()


# Second pathway
somePDFPath = "CellOT_compare_sim2_contactInhibition.pdf"
pdf(file=somePDFPath)
for (nn in 1:length(expression_contactInhibition)){
  
  data = expression_contactInhibition[[nn]]
  geneset_mclu_tran = t(data)
  data <- as.data.frame(geneset_mclu_tran)
  data$id <- 1:nrow(data) 
  plot_data <- melt(data,id.var="id")
  yy <- plot_data$value
  cc <- plot_data$id
  id <- as.factor(plot_data$variable)
  
  fit.random <- ssanova(yy~cc)
  x = data.frame(cc = seq(1, 13, length = 300))
  x2 = data.frame(cc = seq(1:13))
  pred2 <- predict(fit.random, newdata = x, se.fit = T )
  pred_point <- predict(fit.random, newdata = x2, se.fit = T )
  pred2_contactInhibition[[nn]] = pred2$fit
  plot_contactInhibition[[nn]] = pred_point$fit
  
  plot(x$cc,pred2$fit, type = "l", xlab = "Time_points", ylab = "Gene Set Expression", main = "contactInhibition Mean Curve", ylim=c(0.5,2.3))
  lines(x$cc,pred2$fit, col= "red", lwd = 2)
  lines(x$cc,pred_contact$fit, col= "black", lwd = 2)
  lines(x$cc,pred_contact$fit+1.96*pred_contact$se,col=5)
  lines(x$cc,pred_contact$fit-1.96*pred_contact$se,col=5)
  print(nn)
}
dev.off()

#### plot real cell and estimated cell
#index = c(4,8,16,17,19,27,30,32)
index = c(4,8,19)
subset_point = plot_mitosis[index]
subset_pred = pred2_mitosis[index]
selected_cell_mitosis = expression_mitosisGene[index]

dd2  <-  data.frame(matrix(unlist(subset_point), nrow = length(subset_point), byrow = T))
dd2 <- data.frame(t(dd2))
dd2$id <- 1:nrow(dd2) 
dot1 <- melt(dd2,id.var="id")

real_cell_mean_mitosis = data.frame()
for (ii in 1:length(selected_cell_mitosis)){
  temp = selected_cell_mitosis[[ii]]
  res = colMeans(temp)
  real_cell_mean_mitosis = rbind(real_cell_mean_mitosis, res)
}
colnames(real_cell_mean_mitosis) = paste0("time_", seq(1:13))
dd2 <- data.frame(t(real_cell_mean_mitosis))
dd2$id <- 1:nrow(dd2) 
dot1_real <- melt(dd2,id.var="id")

dd  <-  data.frame(matrix(unlist(subset_pred), nrow = length(subset_pred), byrow = T))
dd <- data.frame(t(dd))
dd$id <- 1:nrow(dd)
plot1 <- melt(dd,id.var="id")
plot1$id = rep(x$cc,3)

plot3 <- data.frame(x = x$cc, y = pred_mitosis$fit)


#### second pathway
#index2 = c(5,7,12,16,17,30,44,51)
index2 = c(5,7,16)
subset_pred = pred2_contactInhibition[index2]
subset_point = plot_contactInhibition[index2]
selected_cell_contactInhibition = expression_contactInhibition[index2]

dd  <-  data.frame(matrix(unlist(subset_pred), nrow = length(subset_pred), byrow = T))
dd <- data.frame(t(dd))
dd$id <- 1:nrow(dd) 
plot4 <- melt(dd,id.var="id")
plot4$id = rep(x$cc,3)

plot6 <- data.frame(x = x$cc, y = pred_contact$fit)

subset_point = plot_contactInhibition[index2]
dd2  <-  data.frame(matrix(unlist(subset_point), nrow = length(subset_point), byrow = T))
dd2 <- data.frame(t(dd2))
dd2$id <- 1:nrow(dd2) 
dot2 <- melt(dd2,id.var="id")

real_cell_mean_contactInhibition = data.frame()
for (ii in 1:length(selected_cell_contactInhibition)){
  temp = selected_cell_contactInhibition[[ii]]
  res = colMeans(temp)
  real_cell_mean_contactInhibition = rbind(real_cell_mean_contactInhibition, res)
}
colnames(real_cell_mean_contactInhibition) = paste0("time_", seq(1:13))
dd2 <- data.frame(t(real_cell_mean_contactInhibition))
dd2$id <- 1:nrow(dd2) 
dot2_real <- melt(dd2,id.var="id")

points = rbind(dot1,dot2)

points2 = rbind(dot1_real,dot2_real)


gra_dd = ggplot() + 
  geom_point(data = points,aes(x = id,y = value), size=2, color = "purple") + 
  geom_line(data = plot1,aes(x = id,y = value,group = variable), size=0.8, color = "red") + 
  #geom_line(data = plot2,aes(x = x,y = y), size=1, color = "blue") + 
  geom_point(data = points2,aes(x = id,y = value), size=2, color = "grey") + 
  geom_line(data = plot3,aes(x = x,y = y), size=1, color = "black") + 
  geom_line(data = plot4,aes(x = id,y = value,group = variable), size=0.8, color = "red") + 
  #geom_line(data = plot5,aes(x = x,y = y), size=1, color = "blue") + 
  geom_line(data = plot6,aes(x = x,y = y), size=1, color = "black") + xlab("Time Points") + 
  ylab("Expression")

gra_dd = gra_dd + theme_bw()
gra_dd


### read wot result ###
## First pathway
# wot_compare_trends <- read.delim("mitosisGene_trends.txt")
# 
# geneset_mclu_tran = t(as.matrix(wot_compare_trends[,-1]))
# data <- as.data.frame(geneset_mclu_tran)
# col_mean_mitosisGene= colMeans(data, na.rm = FALSE, dims = 1)
# 
# ## Second pathway
# wot_compare_trends <- read.delim("contactInhibition_trends.txt")
# 
# geneset_mclu_tran = t(as.matrix(wot_compare_trends[,-1]))
# data <- as.data.frame(geneset_mclu_tran)
# col_mean_contactInhibition = colMeans(data, na.rm = FALSE, dims = 1)

#matplot(log(wot_compare_trends), type="l", lty=1, lwd=2, xlab = "Time_points", 
#        ylab = "Gene Expression", main = "Gene Trajectories")



# #### plot first pathway
# index = c(46,42,39,27,77,71,7,6)
# subset_pred = pred2_mitosis[index]
# 
# dd  <-  data.frame(matrix(unlist(subset_pred), nrow = length(subset_pred), byrow = T))
# dd <- data.frame(t(dd))
# dd$id <- 1:nrow(dd) 
# plot1 <- melt(dd,id.var="id")
# plot1$id = rep(x$cc,8)
# 
# 
# plot2 <- data.frame(x = seq(1:13), y = col_mean_mitosisGene)
# plot3 <- data.frame(x = x$cc, y = pred_mitosis$fit)
# 
# subset_point = plot_mitosis[index]
# dd2  <-  data.frame(matrix(unlist(subset_point), nrow = length(subset_point), byrow = T))
# dd2 <- data.frame(t(dd2))
# dd2$id <- 1:nrow(dd2) 
# dot1 <- melt(dd2,id.var="id")
# 
# 
# compare = ggplot() + 
#   geom_line(data = plot1,aes(x = id,y = value,group = variable, colour = "red"), size=0.3) + 
#   geom_line(data = plot2,aes(x = x,y = y, colour = "blue"), size=1) + 
#   geom_point(data = dot1,aes(x = id,y = value,group = variable, colour = "black"), size=0.6) +
#   geom_line(data = plot3,aes(x = x,y = y,colour = "purple"), size=1) +
#   labs(title = "Simulation Comparsion(Two Pathways)", x = "Time", y = "Expression", color = "Legend") +
#   scale_y_continuous(limits = c(0.5, 2.3)) +
#   scale_color_manual(labels = c("Cell Mean Expression","WOT","True Pathway expression", "CellOT"), values = c("blue", "red","black","purple"))
# 
# gra_compare = compare + theme_bw()
# gra_compare
# 
# 
# 
# #### plot second pathway
# index = c(1,17,39,33,50,51,55,75)
# subset_pred = pred2_contactInhibition[index]
# 
# dd  <-  data.frame(matrix(unlist(subset_pred), nrow = length(subset_pred), byrow = T))
# dd <- data.frame(t(dd))
# dd$id <- 1:nrow(dd) 
# plot1 <- melt(dd,id.var="id")
# plot1$id = rep(x$cc,8)
# 
# 
# plot2 <- data.frame(x = seq(1:13), y = col_mean_contactInhibition)
# plot3 <- data.frame(x = x$cc, y = pred_contact$fit)
# 
# subset_point = plot_contactInhibition[index]
# dd2  <-  data.frame(matrix(unlist(subset_point), nrow = length(subset_point), byrow = T))
# dd2 <- data.frame(t(dd2))
# dd2$id <- 1:nrow(dd2) 
# dot1 <- melt(dd2,id.var="id")
# 
# 
# compare = ggplot() + 
#   geom_line(data = plot1,aes(x = id,y = value,group = variable, colour = "red"), size=0.3) + 
#   geom_line(data = plot2,aes(x = x,y = y, colour = "blue"), size=1) + 
#   geom_point(data = dot1,aes(x = id,y = value,group = variable, colour = "black"), size=0.6) +
#   geom_line(data = plot3,aes(x = x,y = y,colour = "purple"), size=1) +
#   labs(title = "Simulation Comparsion(Two Pathways)", x = "Time", y = "Expression", color = "Legend") +
#   scale_y_continuous(limits = c(0.5, 2.3)) +
#   scale_color_manual(labels = c("Cell Mean Expression","WOT","True Pathway expression", "CellOT"), values = c("blue", "red","black","purple"))
# 
# gra_compare = compare + theme_bw()
# gra_compare




### Plot Two pathway together
index = c(4,8,16,17,19,27,30,32)
subset_pred = pred2_mitosis[index]

dd  <-  data.frame(matrix(unlist(subset_pred), nrow = length(subset_pred), byrow = T))
dd <- data.frame(t(dd))
dd$id <- 1:nrow(dd) 
plot1 <- melt(dd,id.var="id")
plot1$id = rep(x$cc,8)


plot2 <- data.frame(x = seq(1:13), y = col_mean_mitosisGene)
plot3 <- data.frame(x = x$cc, y = pred_mitosis$fit)

subset_point = plot_mitosis[index]
dd2  <-  data.frame(matrix(unlist(subset_point), nrow = length(subset_point), byrow = T))
dd2 <- data.frame(t(dd2))
dd2$id <- 1:nrow(dd2) 
dot1 <- melt(dd2,id.var="id")

index2 = c(5,7,12,16,17,30,44,51)
subset_pred = pred2_contactInhibition[index2]

dd  <-  data.frame(matrix(unlist(subset_pred), nrow = length(subset_pred), byrow = T))
dd <- data.frame(t(dd))
dd$id <- 1:nrow(dd) 
plot4 <- melt(dd,id.var="id")
plot4$id = rep(x$cc,8)


plot5 <- data.frame(x = seq(1:13), y = col_mean_contactInhibition)
plot6 <- data.frame(x = x$cc, y = pred_contact$fit)

subset_point = plot_contactInhibition[index2]
dd2  <-  data.frame(matrix(unlist(subset_point), nrow = length(subset_point), byrow = T))
dd2 <- data.frame(t(dd2))
dd2$id <- 1:nrow(dd2) 
dot2 <- melt(dd2,id.var="id")

points = rbind(dot1,dot2)

gra1 = ggplot() + 
  geom_point(data = points,aes(x = id,y = value), size=2, color = "purple") + 
  geom_line(data = plot1,aes(x = id,y = value,group = variable), size=0.8, color = "red") + 
  geom_line(data = plot2,aes(x = x,y = y), size=1, color = "blue") + 
  geom_line(data = plot3,aes(x = x,y = y), size=1, color = "black") + 
  geom_line(data = plot4,aes(x = id,y = value,group = variable), size=0.8, color = "red") + 
  geom_line(data = plot5,aes(x = x,y = y), size=1, color = "blue") + 
  geom_line(data = plot6,aes(x = x,y = y), size=1, color = "black") + xlab("Time Points") + 
  ylab("Expression")

gra1 = gra1 + theme_bw()
gra1



# all_gra = ggplot() + 
#   geom_line(data = plot1,aes(x = id,y = value,group = variable, color = "red"), size=0.3) + 
#   geom_line(data = plot2,aes(x = x,y = y, colour = "blue"), size=1) + 
#   geom_point(data = dot1,aes(x = id,y = value,group = variable, colour = "black"), size=3) +
#   geom_line(data = plot3,aes(x = x,y = y,colour = "purple"), size=1) +
#   geom_line(data = plot4,aes(x = id,y = value,group = variable, colour = "red"), size=0.3) + 
#   geom_line(data = plot5,aes(x = x,y = y, colour = "blue"), size=1) + 
#   geom_point(data = dot2,aes(x = id,y = value,group = variable, colour = "black"), size=3) +
#   geom_line(data = plot6,aes(x = x,y = y,colour = "purple"), size=1) + 
#  labs(title = "Simulation Comparsion (Two Pathways)", x = "Time", y = "Expression", color = "Legend") +
#  scale_y_continuous(limits = c(0.5, 2.3)) +
#  scale_color_manual(labels = c("Gene Expression in Cell","Waddington-OT Method","Benchmark Pathway expressions", "CellST Individual Cell Trajectories"), values = c("orange", "blue","black","red"))
# 
# all_gra = all_gra + theme_bw()
# all_gra


### Two together for heatmap


### optimal transport
couple = list()
for ( ii in 1:12){
  temp1 = count[count$time == ii,]
  temp2 = count[count$time == ii+1,]
  Xs = data.matrix(temp1[1:(length(temp1)-2)])
  Xt = data.matrix(temp2[1:(length(temp2)-2)])
  # EMD Transport
  ot_emd = ot$da$EMDTransport()
  ot_emd$fit(Xs=Xs, Xt=Xt)
  coupling = ot_emd$coupling_
  couple[[ii]] = coupling
} 

### connecting cells
temp1 = count[count$time == 1,]
temp2 = count[count$time == 2,]
temp_coup = couple[[1]]

cell_name_traj = data.frame(time= temp1$index )
index = apply(temp_coup, 1, which.max)
temp2 = temp2[c(index),]
col = temp2$index
cell_name_traj = cbind(cell_name_traj,col)

temp1 = temp2


for (n in 2:12){
  temp2 = count[count$time == n+1,]
  temp_coup = couple[[n]]
  temp_coup = temp_coup[c(index),]
  
  index = apply(temp_coup, 1, which.max)
  temp2 = temp2[c(index),]
  col = temp2$index
  cell_name_traj = cbind(cell_name_traj,col)
  
  temp1 = temp2
  print(n)
}

# Plot Heatmap
gra_plot_list = list()
somePDFPath = "CellOT_heatmap.pdf"
pdf(file=somePDFPath)

for (row in 1:2){
  
  row_num = row
  name = as.character(cell_name_traj[row_num,1])
  value = data.frame(as.numeric(t(count[count$index == name,])))
  value = value[c(1:(nrow(value)-2)),]
  gene_value = as.data.frame(value)
  
  for (nn in 2:13){
    name = as.character(cell_name_traj[row_num,nn])
    value = data.frame(as.numeric(t(count[count$index == name,])))
    value = value[c(1:(nrow(value)-2)),]
    value = as.data.frame(value)
    gene_value = cbind(gene_value,value)
    print(name)
  }
  
  gga = as.matrix(gene_value)
  colnames(gga) = paste0("time_", seq(1:13))
  gra = pheatmap(gga, kmeans_k = 35, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F)
  gra_plot_list[[row]] = gra
}
dev.off()

heat = gga
gra = pheatmap(heat, kmeans_k = 35, cluster_cols = F, cluster_rows = T, show_rownames = F, show_colnames = F)



