% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_network.R
\name{build_network}
\alias{build_network}
\title{Build network}
\usage{
build_network(mucl, gene_name, thres = 0.4)
}
\arguments{
\item{mucl}{a dataframe containing strength of pairwise-gene linear relationship in
different predicted time points. This dataframe could be generated by function
 \code{\link{dynamic_linear_relationship}}.}

\item{gene_name}{a vector of character containing gene names. It is recommended using
function \code{\link{dynamic_linear_relationship}} to generate \code{mucl} so that the order of
gene names of this parameter  is just the order in your \code{data}. See
\code{\link{dynamic_linear_relationship}} for details.}

\item{thres}{a threshold determine whether the edge exist in the network. When the
correlation between two genes is smaller than the threshold, the edge between these
genes will not exist.}
}
\value{
a list containing objects of class 'network'. Each element in
the list represents networks in different predicted time points.
}
\description{
\code{build_network()} establishes network structures between genes
in different predicted time points.
}
