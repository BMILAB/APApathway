#' @title Obtain the score of each pathway through the gene expression matrix
#' @description This function first randomly walks the gene expression matrix in
#'            the protein-protein interaction network, so that each gene in the
#'            PPI network gets a score. Afterwards, effective genes can be
#'            screened according to the set parameters. Finally, according to
#'            the pathway data in the pathway database,the score is matched with
#'            each gene after the score is obtained, and the score of each
#'            pathway is finally obtained.
#' @param x A expression count matrix. The rows correspond to genes and
#'            the columns correspond to cells.
#' @param network A adjacency matrix contation gene-gene interaction network.
#' @param gamma A number between 0 and 1 (default: 0.5).
#' @param pathway pathway data
#' @param select A threshold value, the absolute value of the score is less than
#'             the value is regarded as an invalid value and will not participate
#'             in the calculation
#' @export
#' @return The score of each path in the pathway data set

data("string")
data("LUCS")
data("KEGG")
GetScore<-function(x=LUCS.data,network=string, gamma=0.5,select=0.005,pathway=KEGG){
  ep.data <- matrix(rep(0,nrow(network)*dim(x)[2]),
                    nrow=nrow(network))
  rownames(ep.data) <- rownames(network)
  colnames(ep.data) <- colnames(x)

  gene.overlap <- intersect(rownames(ep.data),
                            rownames(x))
  ep.data[gene.overlap,] <- x[gene.overlap,]
  Anorm <- network/Matrix::rowSums(network)
  eye <- diag(dim(network)[1])
  AA <- Matrix::t(eye - gamma*Anorm)
  BB <- (1-gamma) * ep.data
  smooth.x <- solve(AA, BB)
  smooth.x<-as.matrix(smooth.x)
  smooth.select = data.frame(row.names = rownames(smooth.x)[abs(smooth.x)>select],
                                   level = smooth.x[abs(smooth.x) >select])
  TEST<-as.matrix(apply(pathway, 1,function(x,data){
    x<-as.matrix(x)
    num<-sum(!is.na(x))
    data1<-x[2:num,1]
    data1<-as.matrix(data1)
    index <- match(as.character(data1[,1]),rownames(data))
    index<-as.matrix(index)
    index<-index[!is.na(index[,1]),]
    total<-smooth.select[index,]
    total<-as.matrix(total)
    score<-sum(total[,1])
    score<-score/nrow(total)
  } ,data=smooth.select))
  pathway.score = cbind(as.matrix(pathway[,1]),TEST)
  return(pathway.score)
}

