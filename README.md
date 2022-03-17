# ogClust.WJL

This package includes the following parts. 
* ogClust() and ogClust_Surv are the outcome-guided clustering method with weighted joint likelihood approach for continuous and survival outcome respectively. 
* Rsquare() is the function to calculate the Rsquare of genes and continuous outcome. Predict_test() are the function to predict the cluster label for the testing/validation data based on the ogClust model from training/discovery data and the gene expression from testing/validation data. 
* Sim_continous() and Sim_surv() are the functions to simulate datasets with continous or survival outcome. 
* region_lambda() and region_lambda_surv() are the functions to generate efficient grids of lambda for a given K and outcome weight w. 
* Four real datsets (two lung disease datasets and two triple negative breast cancer datasets) are also included in this package.



## Installing
ogClust.WJL package files are in `ogClust.WJL/` folder, You can install by copying and paste the following code into R

```
devtools::install_github("YujiaLi1994/ogClust.WJL/ogClust.WJL")
```

## Sample code for continous outcomes

* filter genes by marginal screening
```r
t.stat<-rep(NA,nrow(GSE47460_GPL6480$expression))
for(i in 1:nrow(GSE47460_GPL6480$expression)){
  
  data<-data.frame(y=GSE47460_GPL6480$outcome,x=GSE47460_GPL6480$expression[i,])
  mod<-lm(y ~ x, data = data)
  mod<-summary(mod)
  t.stat[i]<-mod$coefficients[2,3]
}
Index<-order(abs(t.stat),decreasing = T)
x<-GSE47460_GPL6480$covariates
G<-GSE47460_GPL6480$expression[Index[1:2000],]
y<-GSE47460_GPL6480$outcome
```
* Get the initial cluster center by K-means
```r
mod.kmeans<-kmeans(t(G),centers = 3,nstart = 20)
center<-t(mod.kmeans$centers)
x<-as.matrix(x)
```

* Set the input parameters
```r
s_G<-200
w<-0.727
w<-(s_G*w)/(s_G*w+1-w)
lambda<-8
K<-3
mod = ogClust(x=x,G=G,y=y,c_center=center,
                        lambda=lambda,v_int=NULL,pi_int=NULL,K=K,max_iter=200,w_outcome=w,w_G=1-w,z_int=NULL)
apply(mod$result_list$z,1,which.max)#cluster label
which(apply(mod$result_list$mu,1,function(x){length(unique(x))})>1)#get the position of informative genes

``` 

* Get the Rsquare of genes and outcome
```r
cluster<-apply(mod$result_list$z,1,which.max)
Rsquare(cluster,Y = y,G = G,X = x)
``` 

* Predict the cluster label for the validation data
```r
G.test<-GSE47460_GPL14550$expression
X.test<-GSE47460_GPL14550$covariates
index<-match(rownames(G),rownames(G.test))
G.test<-G.test[index,]
mod.predict<-predict_test(mod,K = 3,D.test =t(G.test),X1 =X.test, p = nrow(G))

```

## Example code for survival outcomes
* Simulate a dataset of survival outcome using Sim_surv()
```r
Data<-Sim_surv(n=99,beta1 = 0.5,beta2 = 0.5,q = 50,q1 = 300,q2 = 50,q3 = 1600,c1 = 2,var_g = 2,mu = 1.5,mu1 = 1.7,sigma_y = 0.5,censor = 100)
G<-Data$G
X<-Data$x
Y<-Data$y
Y.ind<-Data$y.ind
```
* Implement ogClust_Surv()
```r
mod.kmeans<-kmeans(G,centers = 3,nstart = 20)
center<-t(mod.kmeans$centers)
region_lambda_surv(lambda1=15,lambda2 = 0,iteration = 10,Y = Y,G = G,X = X,center =center ,w = 0.5,K = 3,delta = Y.ind)
```
## Sample code for getting efficient grids of lambdas for fixed w and K

* continuous outcome
```r
Data<-Sim_continuous(n=99,beta1 = 1,beta2 = 1,q = 50,q1 = 50,q2 = 50,q3 = 1850,c1 = 3,var_g = 2,mu = 1.5,mu1 = 1.7,sigma_y = 1)
G<-Data$G
X<-Data$x
Y<-Data$y
mod.kmeans<-kmeans(G,centers = 3,nstart = 20)
center<-t(mod.kmeans$centers)
region_lambda(lambda1=15,lambda2 = 0,iteration = 10,Y = Y,G = G,X = X,center =center ,w = 0.5,K = 3)
```

* survival outcome
```r
Data<-Sim_surv(n=99,beta1 = 0.5,beta2 = 0.5,q = 50,q1 = 300,q2 = 50,q3 = 1600,c1 = 2,var_g = 2,mu = 1.5,mu1 = 1.7,sigma_y = 0.5,censor = 100)
G<-Data$G
X<-Data$x
Y<-Data$y
Y.ind<-Data$y.ind
mod.kmeans<-kmeans(G,centers = 3,nstart = 20)
center<-t(mod.kmeans$centers)
region_lambda_surv(lambda1=15,lambda2 = 0,iteration = 10,Y = Y,G = G,X = X,center =center ,w = 0.5,K = 3,delta = Y.ind)

```

## Access four real datasets used in ogClust.WJL paper
All the datasets are imported and accessible after library(ogClust.WJL).
```r
#Data sets in package ‘ogClust.WJL’:
# GSE47460_GPL6480              Microarray data of COPD or ILD lung disease after preprocessing(GSE47460)
# GSE47460_GPL14550             Microarray data of COPD or ILD lung disease after preprocessing(GSE47460)
# Data.Metabric                 Triple negative breast cancer patients in METABRIC breast cancer dataset
# Data.ScanB                   Triple negative breast cancer patients in SCAN-B breast cancer dataset(GSE60789)
```
