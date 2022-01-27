suppressMessages(library(reticulate))
### use Conda environment for optimal transport
use_condaenv(condaenv = "myenv", conda = "/Users/mengruizhang/miniconda3/condabin/conda")

CellST<- function(count, time_label,nTime) {
  np = import("numpy")
  ot = import("ot")
  count$time = time_label
  count$index = rownames(count)
  couple = list()
  for ( ii in 1:(nTime-1)){
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
  
  
  for (n in 2:(nTime-1)){
    temp2 = count[count$time == n+1,]
    temp_coup = couple[[n]]
    temp_coup = temp_coup[c(index),]
    
    index = apply(temp_coup, 1, which.max)
    temp2 = temp2[c(index),]
    col = temp2$index
    cell_name_traj = cbind(cell_name_traj,col)
    
    temp1 = temp2
  }
  return(cell_name_traj)
}

get_traj_exp <- function(cell_name_traj, all, nTime, nGenes, gene_name){
  
  expression = list()
  for (row in 1:nrow(cell_name_traj)){
    row_num = row
    name = as.character(cell_name_traj[row_num,1])
    value = data.frame(as.numeric(t(all[all$index == name,])))
    value = value[c(1:(nrow(value)-2)),]
    gene_value = as.data.frame(value)
    
    for (nn in 2:nTime){
      name = as.character(cell_name_traj[row_num,nn])
      value = data.frame(as.numeric(t(all[all$index == name,])))
      value = value[c(1:(nrow(value)-2)),]
      value = as.data.frame(value)
      gene_value = cbind(gene_value,value)
    }
    
    gga = as.matrix(gene_value)
    colnames(gga) = paste0("time_", seq(1:nTime))
    rownames(gga) = gene_name
    expression[[row]] = gga
  }
  return(expression)
}


get_gene_exp <- function(expression, nGenes, gene_name) {
  gene_name = gene_name
  gene_expression_list = list()
  
  for (row in 1:nGenes){
    gene = gene_name[row]
    temp = lapply(expression, function(x){x[rownames(x) %in% gene]})
    df = data.frame(matrix(unlist(temp), nrow=length(temp), byrow=TRUE))
    gene_expression_list[[row]] = df
  }
  return(gene_expression_list)
}


get_true_correlation <- function(sim, datalist, nTime){
  time = list()
  for (t in 1:nTime){
    group2_gcnlist = lapply(datalist, 
                            function(data){gcn(data[,which(colData(sim)$Path == t)], CPM2 = TRUE)})
    time[[t]] = group2_gcnlist$`simulated truth`
  }
  return(time)
}










