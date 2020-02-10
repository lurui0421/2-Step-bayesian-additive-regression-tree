### Fix rJava issue of installing 
install.packages("rJava")
library(rJava)
install.packages('parallel')
library("parallel")

### simulate data set use to test functions
test_parallel<-function(N,p)
      {mat<-data.frame(matrix(rep(rnorm(N),10),ncol=10))
      colnames(mat)<-paste("X",c(1:10),sep="")
      attach(mat)
      Y1<-X1+X2+X3+X4+X5+X6+rnorm(N)
      Y2<-X1+X2*X6+X3+X4+X5+X6+X7+X8+X9+X10+rnorm(N)
      treat<-rbinom(N,1,p)
      Y<-treat*Y1+(1-treat)*Y2
      treat<-ifelse(treat==0,1,2)
      data<-data.frame(mat,treat,Y1,Y2,Y)
      detach(mat)
 return(data)
}

### BART with parallel computing 
parallel_bart<-function(data_list,mc.cores = 4){
ATEs <- mclapply(X= data_test, FUN = BART,mc.cores =mc.cores)
mean_ATE<-mean(unlist(ATEs))
var_ATE<-var(unlist(ATEs))
result<-c(mean_ATE,var_ATE)
names(result)<-c("mean_ATE","var_ATE")
return(result)}
