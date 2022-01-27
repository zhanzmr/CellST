library(ggplot2)
library(dplyr)
library(gtools)
library(tidyverse)
library(pheatmap)
library(gss)
library(reshape2)

############# Plot UMAP Graph ######################
setwd("~/Documents/big_data/single_cell/new_dy_network_07082021/mouse_embryogenesis/pathway")
### Plot UMAP visualization
umap = read.csv("all_umap.csv")
colnames(umap)[1] <- "id"
cell_time <- read.csv("cell_time.csv")
umap_expression <- merge(umap,cell_time,by="id")
subset <- umap_expression[umap_expression$time <= 16,]

# total = data.frame()
# for (i in 1:17){
#   temp = read.csv(file[i])
#   sample_cell = sample_n(temp,n)
#   total = rbind(total,sample_cell)
#   print(i)
# }
# 
# name <- total$index
# rownames(subset) = subset$id
# cell_subset <- subset[name,]

plot_cell <- sample_n(subset, 700)

ggplot(plot_cell, aes(Umap1, Umap2, colour = factor(time))) + geom_point() + 
  theme_bw() + theme(legend.position="none")

#### Plot trajectries
load("~/Documents/big_data/single_cell/new_dy_network_07082021/mouse_embryogenesis/pathway/cell_cycle/r_envir/cell_cycle.RData")
cell_name_traj = trajectory_name
tra_num = 10
sub_traj = sample_n(cell_name_traj, tra_num)
nn = as_vector(sub_traj)
cell_subset <- subset[nn,]

ggplot(plot_cell, aes(Umap1, Umap2, colour = factor(time))) + geom_point() + 
  theme_bw() + theme(legend.position="none")

subset_vector = stack(sub_traj)$values
subset_umap <- subset[subset$id %in% subset_vector, ]

traject = seq(1:tra_num)
traject_umap <- paste("traject ", traject, sep="")
umap_group = rep(traject_umap,17)
subset_umap$group = umap_group

##add line
line_connect1 <- subset_umap[subset_umap$time == 0 | subset_umap$time == 1, ]
line_connect2 <- subset_umap[subset_umap$time == 1 | subset_umap$time == 2, ]
line_connect3 <- subset_umap[subset_umap$time == 2 | subset_umap$time == 3, ]
line_connect4 <- subset_umap[subset_umap$time == 3 | subset_umap$time == 4, ]
line_connect5 <- subset_umap[subset_umap$time == 4 | subset_umap$time == 5, ]
line_connect6 <- subset_umap[subset_umap$time == 5 | subset_umap$time == 6, ]
line_connect7 <- subset_umap[subset_umap$time == 6 | subset_umap$time == 7, ]
line_connect8 <- subset_umap[subset_umap$time == 7 | subset_umap$time == 8, ]
line_connect9 <- subset_umap[subset_umap$time == 8 | subset_umap$time == 9, ]
line_connect10 <- subset_umap[subset_umap$time == 9 | subset_umap$time == 10, ]
line_connect11 <- subset_umap[subset_umap$time == 10 | subset_umap$time == 11, ]
line_connect12 <- subset_umap[subset_umap$time == 11 | subset_umap$time == 12, ]
line_connect13 <- subset_umap[subset_umap$time == 12 | subset_umap$time == 13, ]
line_connect14 <- subset_umap[subset_umap$time == 13 | subset_umap$time == 14, ]
line_connect15 <- subset_umap[subset_umap$time == 14 | subset_umap$time == 15, ]
line_connect16 <- subset_umap[subset_umap$time == 15 | subset_umap$time == 16, ]
line_connect17 <- subset_umap[subset_umap$time == 16 | subset_umap$time == 17, ]

gra = ggplot(subset_umap, aes(Umap1, Umap2))+
  geom_point(data = subset_umap, aes(x = Umap1, y = Umap2, color = factor(time))) + 
  geom_line(data = line_connect1, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect2, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect3, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect4, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect5, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect6, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect7, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect8, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect9, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect10, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect11, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect12, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect13, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect14, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect15, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect16, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  geom_line(data = line_connect17, aes(x = Umap1, y = Umap2, group = group), size = 0.1) +
  xlab("Umap1") + ylab("Umap2")

gra[["labels"]][["colour"]] = "Experimental Time"
gra1 = gra + theme_bw() + theme(legend.position="bottom")
gra1