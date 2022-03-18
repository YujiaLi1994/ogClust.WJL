#' @title Simulation of survival outcome in the paper
#' @param n number of total samples (the first n/3 samples belong to cluster 1, next n/3 belongs to cluster 2 and the last n/3 belong to cluster 3)
#' @param beta1 coefficent of first covariate
#' @param beta2 coefficent of second covariate
#' @param q number of true informative genes
#' @param q1 number of genes for random clusters
#' @param q2 number of genes correlated with covariates
#' @param q3 number of noise genes
#' @param c1 effect size of outcome separation
#' @param var_g control the correlation of the q3 genes with covariates
#' @param mu effect size of gene separation
#' @param mu1 effect size of random clusters
#' @param sigma_y standard deviation for outcome
#' @param censor time point of right censor, the outcome larger than it will be censored.
#' @return
#' result_list: A list with four components.
#' \itemize{
#' \item{x: }{simulated covariate matrix}
#' \item{y: }{simulated outcome}
#' \item{y.ind: }{censor indication}
#' \item{G: }{simulated gene expression data}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' Data<-Sim_surv(n=99,beta1 = 0.5,beta2 = 0.5,q = 50,q1 = 300,q2 = 50,q3 = 1600,c1 = 2,var_g = 2,mu = 1.5,mu1 = 1.7,sigma_y = 0.5,censor = 100)
#' }


Sim_surv<-function(n,beta1,beta2,q,q1,q2,q3,c1,var_g,mu,mu1,sigma_y,censor){
  #simulatex
  x1<-rnorm(n,1,0.5)
  x2<-rnorm(n,1,0.5)
  x<-cbind(x1,x2)
  #simulate y
  y<-rep(NA,n)
  y[1:(n/3)]<-0+beta1*x1[1:(n/3)]+beta2*x2[1:(n/3)]+rlogis(n/3, location = 0, scale = 1)*sigma_y
  y[(n/3+1):(2*n/3)]<-c1+beta1*x1[(n/3+1):(2*n/3)]+beta2*x2[(n/3+1):(2*n/3)]+rlogis(n/3, location = 0, scale = 1)*sigma_y
  y[(2*n/3+1):n]<-2*c1+beta1*x1[(2*n/3+1):(n)]+beta2*x2[(2*n/3+1):(n)]+rlogis(n/3, location = 0, scale = 1)*sigma_y
  y<-exp(y)
  y.ind<-ifelse(y>censor,0,1)
  y<-ifelse(y>censor,censor,y)
  X1 = rbind(matrix(rnorm((n/3) * q), ncol = q),
             matrix(rnorm((n/3) * q), ncol = q),
             matrix(rnorm((n/3) * q), ncol = q))
  X2 = rbind(matrix(rnorm((n/3) * q1), ncol = q1),
             matrix(rnorm((n/3) * q1), ncol = q1),
             matrix(rnorm((n/3) * q1), ncol = q1))
  X3 = rbind(matrix(rnorm((n/3) * q2), ncol = q2),
             matrix(rnorm((n/3) * q2), ncol = q2),
             matrix(rnorm((n/3) * q2), ncol = q2))
  X.noise = rbind(matrix(rnorm((n/3) * q3), ncol = q3),
                  matrix(rnorm((n/3) * q3), ncol = q3),
                  matrix(rnorm((n/3) * q3), ncol = q3))
  X1[1:(n/3), 1:q] <- X1[1:(n/3), 1:q] - mu
  X1[(n/3+1):(2*n/3), 1:q] <- X1[(n/3+1):(2*n/3), 1:q]
  X1[(2*n/3+1):n, 1:q] <- X1[(2*n/3+1):n, 1:q] + mu
  index<-sample(1:n)
  X2[index[1:(n/3)], ] <- X2[index[1:(n/3)],]-mu1
  X2[index[(n/3+1):(2*n/3)], ] <- X2[index[(n/3+1):(2*n/3)], ]
  X2[index[(2*n/3+1):n], ] <- X2[index[(2*n/3+1):n], ]+mu1
  for(i in 1:ncol(X3)){
    X3[,i]<-X3[,i]+rnorm(1,0,var_g)*x1
  }
  X<-cbind(X1,X2,X3,X.noise)
  data<-list(x=x,y=y,y.ind=y.ind,G=X)
  return(data)
}
