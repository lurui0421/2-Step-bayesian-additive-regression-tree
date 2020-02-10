library("twang")
data<-Latent.GDP.2(1000,
                   15,6,
                   0.5,0.4,
                   0.5,0.2,
                   0.4,0.3,
                   3,
                   -0.41,1,
                   -0.51,1,
                   1,1,
                   1,1,
                   1,1,
                   c(0,5,10),"linear")

data<-data[[1]]

latent.modal.sim<-function(data,num.tx=3)
{
  ## modal class assignment
  f.8cov<-cbind(U1, U2, U3, U4, U5)~X1+X2+X3+X4+X5+X6+X7+X8 ##linear combination of the x on latent class
  nes.8cov<-poLCA(f.8cov, data, nclass=num.tx, maxiter=50000, verbose = FALSE)
  prob.start<-poLCA.reorder(nes.8cov$probs.start, order(nes.8cov$P, decreasing = TRUE))
  nes.8cov<-poLCA(f.8cov, data, nclass=3, maxiter=50000, probs.start = prob.start,verbose=FALSE, graph=TRUE)
}

modal<-latent.modal.sim(data)
modal$predclass

twostep<-function(data,ps.method="linear")
{
  modal<-latent.modal.sim(data)
  data$modal_c<-modal$predclass
  
  ###############################################################
  ## ESTIMATE PROPENSITY SCORES: Modal assignment
  ###############################################################
## modal class Indicators for LATENT CLASS 1, LATENT CLASS 2, LATENT CLASS 3
 
  data$m1 <-as.numeric(data$modal_c==1)
  data$m2 <-as.numeric(data$modal_c==2)
  data$m3 <-as.numeric(data$modal_c==3)
  
if (ps.method=="twang")
  {ps.1m.model<-ps(m1~X1+X2+X3+X4+X5+X6+X7+X8,data=data,n.tree=5000,interaction.depth=2,shrinkage=0.01,perm.test.iters=0,stop.method=c('es.mean'),estimand="ATE",verbose=FALSE)
    ps.2m.model<-ps(m2~X1+X2+X3+X4+X5+X6+X7+X8,data=data,n.tree=5000,interaction.depth=2,shrinkage=0.01,perm.test.iters=0,stop.method=c('es.mean'),estimand="ATE",verbose=FALSE)
    ps.3m.model<-ps(m3~X1+X2+X3+X4+X5+X6+X7+X8,data=data,n.tree=5000,interaction.depth=2,shrinkage=0.01,perm.test.iters=0,stop.method=c('es.mean'),estimand="ATE",verbose=FALSE)
    ps.1m<-ps.1m.model$ps
    ps.2m<-ps.2m.model$ps
    ps.3m<-ps.3m.model$ps}
 { ## Run 3 (binary) propensity score models...
    ps.1m <- glm(m1~X1+X2+X3+X4+X5+X6+X7+X8, data=data, family=binomial)
    ps.2m <- glm(m2~X1+X2+X3+X4+X5+X6+X7+X8, data=data, family=binomial)
    ps.3m <- glm(m3~X1+X2+X3+X4+X5+X6+X7+X8, data=data, family=binomial)
    ## Extract estimated PS from binary PS models
    ps1m <- predict(ps.1m, type = "response")
    ps2m <- predict(ps.2m, type = "response")
    ps3m <- predict(ps.3m, type = "response")}
  
  ## Extract estimated PS from binary PS models

  ## Calculate ATE weights from binary PS models
  w.1m <- ifelse(data$modal_c==1, 1/ps1m, 1/(1-ps1m))
  w.2m <- ifelse(data$modal_c==2, 1/ps2m, 1/(1-ps2m))
  w.3m <- ifelse(data$modal_c==3, 1/ps3m, 1/(1-ps3m))
  ## Set final estimated ATE weight to ATE weight calc from observed LC
  data$w.fm[data$modal_c==1] <- w.1m[data$modal_c==1]
  data$w.fm[data$modal_c==2] <- w.2m[data$modal_c==2]
  data$w.fm[data$modal_c==3] <- w.3m[data$modal_c==3]
  ## Set weights above 98th percentile to 98th percentile value
  data$w.fm[data$w.fm>quantile(data$w.fm,probs=.98)]<-quantile(data$w.fm,probs=.98)
  
  ###############################################################
  ## ESTIMATE TX EFFECTS: Modal + Propensity Scores
  ###############################################################
  
  # outcome regression for Y ~ Tx (Tx is factor variable)
  data$modal_c <- factor(data$modal_c, levels=c(1,2,3))
  m.design <- svydesign(id=~1, weights=~w.fm, data=data)
  out1m<-svyglm(Y_obs ~ C(modal_c,contr.treatment(3,base=1)),
                design=m.design,data=data)
  out2m <- svyglm(Y_obs ~C(modal_c,contr.treatment(3,base=2)),
                  design=m.design,data=data)
  ## Class 2 v Class 1
  m.2v1 <- summary(out1m)$coefficient[2,1]
  
  ## Class 3 v Class 1
  m.3v1 <- summary(out1m)$coefficient[3,1]
  
  ## Class 3 v Class 2
  m.3v2 <- summary(out2m)$coefficient[3,1]
  
  output <- NULL
  output<-as.numeric(c(m.2v1,
                       m.3v1,
                       m.3v2))
   return(output) }





### 2-step method for original data ###

step2.result.origin<-matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100)
{
  step2.result.origin[i,]<-twostep(data=data.list.origin[[i]])
}

apply(step1.result,2,mean)
apply(step2.result,2,mean)
apply(step3.result,2,mean)

### 2-step method for new data ###

step2.result<-matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100)
{
  step2.result[i,]<-twostep(data=data.list[[i]][[1]])
}

### 2-step method for linear small confounding effect ###

step2.result.linear.small<-matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100)
{
  step2.result.linear.small[i,]<-twostep(data=data.list.linear.small[[i]])
}

### 2-step method for linear medium confounding effect ###

step2.result.linear.medium<-matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100)
{
  step2.result.linear.medium[i,]<-twostep(data=data.list.linear.medium[[i]])
}

### 2-step method for linear large confounding effect ###

step2.result.linear.large<-matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100)
{
  step2.result.linear.large[i,]<-twostep(data=data.list.linear.large[[i]])
}

### 2-step method for nonlinear small confounding effect ###

step2.result.nonlinear.small<-matrix(NA, nrow = 100, ncol = 3)

for(i in 93:100)
{
  step2.result.nonlinear.small[i,]<-twostep(data=data.list.nonlinear.small[[i]],ps.method = 'twang')
}


### 2-step method for nonlinear medium confounding effect ###

step2.result.nonlinear.medium<-matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100)
{
  step2.result.nonlinear.medium[i,]<-twostep(data=data.list.nonlinear.medium[[i]],ps.method = 'twang')
}

### 2-step method for nonlinear large confounding effect ###

step2.result.nonlinear.large<-matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100)
{
  step2.result.nonlinear.large[i,]<-twostep(data=data.list.nonlinear.large[[i]],ps.method = 'twang')
}

apply(abs(step2.result.linear.small),2,mean)-c(2,4,2)
apply(abs(step2.result.linear.medium),2,mean)-c(2,4,2)
apply(abs(step2.result.linear.large),2,mean)-c(2,4,2)
apply(abs(step2.result.nonlinear.small),2,mean)-c(2,4,2)
apply(abs(step2.result.nonlinear.medium),2,mean)-c(2,4,2)
apply(abs(step2.result.nonlinear.large),2,mean)-c(2,4,2)

step2.result.linear.small.bias<-apply(step2.result.linear.small,2,mean)-c(2,4,2)
step2.result.linear.medium.bias<-apply(step2.result.linear.medium,2,mean)-c(2,4,2)
step2.result.linear.large.bias<-apply(step2.result.linear.large,2,mean)-c(2,4,2)
step2.result.nonlinear.small.bias<-apply(step2.result.nonlinear.small,2,mean)-c(2,4,2)
step2.result.nonlinear.medium.bias<-apply(step2.result.nonlinear.medium,2,mean)-c(2,4,2)
step2.result.nonlinear.large.bias<-apply(step2.result.nonlinear.large,2,mean)-c(2,4,2)

step2.result.linear.small.sd<-apply(step2.result.linear.small,2,sd)
step2.result.linear.medium.sd<-apply(step2.result.linear.medium,2,sd)
step2.result.linear.large.sd<-apply(step2.result.linear.large,2,sd)
step2.result.nonlinear.small.sd<-apply(step2.result.nonlinear.small,2,sd)
step2.result.nonlinear.medium.sd<-apply(step2.result.nonlinear.medium,2,sd)
step2.result.nonlinear.large.sd<-apply(step2.result.nonlinear.large,2,sd)


#compare the predicted class and the real class
diff.comp.modal<-cbind(data$trueclass, modal$predclass,rep(NA,1000))
diff.comp.modal[,3]<-ifelse(abs(diff.comp.modal[,1]-diff.comp.modal[,2])>0, 1, 0)
colnames(diff.comp.modal)<-c("true class", "predict class", "difference")
sum(diff.comp.modal[,3])
