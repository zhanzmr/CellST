#' Optimal transport across different time points
#'
#' @description \code{cell_name_traj()} accomplishes a cell lineage tracing method, which aligns
#' two individual cells between any two adjacent time points via the optimal
#' transport technique(Earth mover's distance is used in this function), and returns
#' unique trajectories for individual cells.
#'
#' Note that users need to set up a conda environment with package "pot" installed
#' before using this function.
#'
#' @param data a dataframe with \eqn{N*t} rows and \eqn{d+2} columns, where \eqn{N} is
#' the number of cells in one time point(the number of cells need to be equal
#' among different time points),
#' \eqn{t} is the number of time points, \eqn{d} is the number of genes. The first column
#' is the cell names and the last column is the time points which contains continuous
#' integers beginning from 1. The column names of the first and the last columns need to be set as
#' 'index' and 'time', respectively.
#' This dataframe could be regraded as a result of \code{rbind()} of cell-gene
#' dataset ordered by \emph{increasing} time points.
#'
#' @return a dataframe with \eqn{N} rows and \eqn{t} columns, giving \eqn{N} unique
#' trajectories for individual cells, each trajectory has length \eqn{t}.
#'
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import gtools
#' @import tidyverse
#' @import pheatmap
#' @import gss
#' @import reshape2
#' @import reticulate
#' @import fda
#' @import CombMSC
#' @import igraph
#' @import network
#' @import ndtv
#' @import intergraph
#' @import networkDynamic
#' @import ESCO
#' @import progress
#' @import purrr
#'
#'
cell_name_traj <- function(data){
  ot = import("ot")
  couple = list()
  t = max(data$time)
  pb <- progress_bar$new(
    format = 'Processing [:bar] :percent eta: :eta',
    total = (t-1), clear = FALSE, width = 80)
  for ( jj in 1:(t-1)){
    pb$tick()
    temp1 = data[data$time == jj,]
    temp2 = data[data$time == jj+1,]
    Xs = data.matrix(temp1[2:(length(temp1)-1)])
    Xt = data.matrix(temp2[2:(length(temp2)-1)])
    # EMD Transport
    ot_emd = ot$da$EMDTransport()
    ot_emd$fit(Xs=Xs, Xt=Xt)
    coupling = ot_emd$coupling_
    couple[[jj]] = coupling
    #print(jj)
  }

  ### connecting cells
  temp1 = data[data$time == 1,]
  temp2 = data[data$time == 2,]
  temp_coup = couple[[1]]

  cell_name_traj1 = data.frame(time= temp1$index )
  index = apply(temp_coup, 1, which.max)
  temp2 = temp2[c(index),]
  col = temp2$index
  cell_name_traj1 = cbind(cell_name_traj1,col)

  temp1 = temp2
  if(t > 2){
    for (n in 2:(t-1)){
      temp2 = data[data$time == n+1,]
      temp_coup = couple[[n]]
      temp_coup = temp_coup[c(index),]

      index = apply(temp_coup, 1, which.max)
      temp2 = temp2[c(index),]
      col = temp2$index
      cell_name_traj1 = cbind(cell_name_traj1,col)

      temp1 = temp2
      #print(n)
    }
  }


  x = seq(1:t)
  col = paste0("time",x)
  colnames(cell_name_traj1) = col
  return(cell_name_traj1)
}
