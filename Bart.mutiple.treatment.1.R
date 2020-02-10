##################################################################################################################
##                                              Bart for Mutiple Treatment 
##################################################################################################################
### Author: Rui Lu 
#### Mutiple Treatment Bart Functions 
## Vanilla BART 
install.packages("dbarts")
library("dbarts")
install.packages("BART")
library("BART")
## Bart Functions 
Vanilla.bart <-function(data,cov,treat)
    { 
      ### xcov is the covariate set 
      ### treat is the treatment indicator
      ### Step 1: Fit Bart Model
      xt<-cbind(cov,treat)
      bart_mod<-bart(x.train = xt,y.train = data$Y_obs,keeptrees = TRUE)
      n.class<-length(levels(as.factor(treat)))
      y.potent<-matrix(,nrow=nrow(data),ncol=n.class)
      ### Step 2:  Impute potencial outcomes
      for( i in 1:n.class){
           xp<-xt
           xp[,grep("treat",names(xp))]<-i
           bart_pred<-predict(bart_mod,xp)
           y.potent[,i]<-rowMeans(bart_pred)
           }
      return(y.potent)}

#### Test for Vanilla bart
fit<-Vanilla.bart(data,cov=data[,17:24],treat=treatment)  
### The Propensity Score Bart 
install.packages("twang")
library("twang")
installed.packages("BART")
library("BART")
### The P.S Bart Function
ps.bart<-function(data,cov,treat,ps.method)
    { ### Method is the methods to fit propensity score: linear,GBM or BART
      xt<-cbind(cov,treat)
      n.class<-length(levels(as.factor(data$treat)))
      f<-as.formula(paste("treat~",paste0("X",1:ncol(xcov),collapse = "+"),sep=""))
      ### Fit propensity score model
      if(ps.method=="linear"){
      ### Mutinominal logistic regression 
              muti.logit<- nnet::multinom(formula=f,data=data)
              p.score<-fitted(muti.logit)
              p.score<-as.data.frame(p.score)
              names(p.score)<-paste("ps",0:(n.class-1),sep="")}
      else if (ps.method=="twang"){
      ### Mutinominal GBM  ####
              data$treat<-as.factor(data$treat)
              fit<-mnps(f,data=data,estimand="ATE",stop.method = "es.mean")
              p.score<-data.frame(fit$psList[[1]]$ps$es.mean.ATE,fit$psList[[2]]$ps$es.mean.ATE, fit$psList[[3]]$ps$es.mean.ATE)
              names(p.score)<-paste("ps",0:(n.class-1),sep="")}
      else if (ps.method=="BART"){
      ###Mutinominal Bart ###
                p.bart<-mbart(xcov,treat, nskip=1, ndpost=1,ntype=3)
                ps<-predict(p.bart,xcov)
                idx = seq(1, 3*length(treat), by=3)
                ps1=ps$prob.test.mean[idx]
                ps2=ps$prob.test.mean[idx+1]
                ps3=ps$prob.test.mean[idx+2]
                p.score<-data.frame(ps1,ps2,ps3)
                ## Check rowSums(p.score)
                names(p.score)<-paste("ps",0:(n.class-1),sep="")}
                xt<-cbind(p.score,xt)
               ### Fit Bart Model
      bart_mod<-bart(x.train = xt,y.train = data$Y_obs,keeptrees = TRUE)
      y.potent<-matrix(,nrow=nrow(data),ncol=n.class)
      ### Impute potencial outcome
      for( i in 1:n.class){
                          xp<-xt
                          xp[,grep("treat",names(xp))]<-i
                          bart_pred<-predict(bart_mod,xp)
                          y.potent[,i]<-rowMeans(bart_pred)
      }
      return(y.potent) }

## The BART.CV function 
bart.cv<-function(data,cov,treat,k)
        {## Step 1: corss validate K
            xt<-cbind(xcov,treat)
            xval <- dbarts::xbart(xt,data$Y_obs,k = k, n.trees=200)
            bestk <- k[which.min(apply(xval, 2L, mean))]
          ## Step 2: Fit Bart Model use best K
            bart_mod<-bart(x.train = xt,y.train = data$Y_obs,keeptrees = TRUE,k=bestk)
            n.class<-length(levels(as.factor(treat)))
            y.potent<-matrix(,nrow=nrow(data),ncol=n.class)
          ## Step 3: Impute potencial outcome
        for( i in 1:n.class){
            xp<-xt
            xp[,grep("treat",names(xp))]<-i
            bart_pred<-predict(bart_mod,xp)
            y.potent[,i]<-rowMeans(bart_pred)
          }
          return(y.potent) }

### Use Bartcause in stead
install.packages("devtools")
devtools::install_github("vdorie/bartCause")
library("bartCause")
head(test.data.1[[1]])
test.data.1[[1]]$Y_obs

treatment<-lca.assginment(data.list.linear.large[[i]],method="Modal")
treatment.1<-ifelse(treatment==1,1,)
bartc(test.data.1[[1]]$Y_obs,


### Attaching Bart for latent class analysis 
Bart.latent<-function(data,cov,treat,bart.method="Vanilla.bart",ps.method="twang")
            {  if (bart.method=="Vanilla.bart"){
                 outcome<-Vanilla.bart(data,cov,treat)
                 ATE.1.0<-mean(outcome[,2]-outcome[,1])
                 ATE.2.1<-mean(outcome[,3]-outcome[,2])
                 ATE.2.0<-mean(outcome[,3]-outcome[,1])
                 result<-c(ATE.1.0,ATE.2.0,ATE.2.1)
              }
              else if (bart.method=="Bart.CV"){
                 k <- c(0.5, 1, 2, 4, 8,10,20)
                 outcome<-bart.cv(data,xcov,treat)
                 ATE.1.0<-mean(outcome[,2]-outcome[,1])
                 ATE.2.1<-mean(outcome[,3]-outcome[,2])
                 ATE.2.0<-mean(outcome[,3]-outcome[,1])
                 result<-c(ATE.1.0,ATE.2.0,ATE.2.1)
              }
              else {
                 outcome<-ps.bart(data,cov,data$treat,ps.method=ps.method)
                 ATE.1.0<-mean(outcome[,2]-outcome[,1])
                 ATE.2.0<-mean(outcome[,3]-outcome[,1])
                 ATE.2.1<-mean(outcome[,3]-outcome[,2])
                 result<-c(ATE.1.0,ATE.2.0,ATE.2.1) }
                 names(result)<-c("ATE.1.0","ATE.2.0","ATE.2.1")
                 return(result)
            }

### Function to get result from bart model 
Bart.result<-function(data,
                      cov,
                      treat,
                      bart.method="Vanilla.bart",
                      ps.method="twang",
                      assign.method="pseduo",
                      K=20)
            {if(assign.method=="modal")
              {   result<-Bart.latent(data,cov,treat=treatment,bart.method=bart.method,ps.method=ps.method) }
             else{result<-matrix(rep(NA,3*K),ncol=3,nrow=K,byrow=TRUE)
              for (i in 1:20){
                        result[i,]<-Bart.latent(data,cov,treat=treatment[,i],bart.method=bart.method,ps.method=ps.method)}
                        result<-apply(result,2,mean) }
                        return(result)}



