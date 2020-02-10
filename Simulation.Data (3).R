### data simulation 
### Author Rui Lu 
### Simulation function
### Simulation according to stuart's function 

DGP.latent.1<-function(N, tx.prb=c(0.33,0.33,0.34), # Probability to generate latent class 
              alpha1=4,alpha2=3,alpha3=4, # Correlation with treatment assignment 
              beta1=1,beta2=1,beta3=1,    # Correlation with responsse variables
              num.cov=8,num.item=15,# Number of observed covariates and indicator variables 
              response="linear", # Response surface is linear or nonlinear
              p=c(0.6,0.5,0.5)) # Indicator probability for class 1, class 2 and class 3 
            {
              ##Start with potencial covariates
                X1 <- rnorm(N)## Instrument 
                X2 <- rnorm(N)## Instrument
                X3 <- rnorm(N)## Confounder
                X4 <- rnorm(N)## Confounder
                X5 <- rnorm(N)## Confounder
                X6 <- rnorm(N)## Confounder
                X7 <- rnorm(N)## Risk Factor
                X8 <- rnorm(N) ## Risk Factor
              ## Create treatment index 
                tx1.ind <- tx2.ind <- tx3.ind <- tx <- rep(NA,N)
                Treat.gen <- rmultinom(N, size=1, prob=tx.prb)## True treatment class
                Treat.assign<-matrix(Treat.gen,ncol=3,nrow=N,byrow=TRUE)
                tx1.ind <- Treat.assign[,1]==1
                tx2.ind <- Treat.assign[,2]==1
                tx3.ind <- Treat.assign[,3]==1
                tx[tx1.ind] <- 1
                tx[tx2.ind] <- 2
                tx[tx3.ind] <- 3
                X.matrix <- cbind(X1,X2,X3,X4,X5,X6,X7,X8)
            ## Strengthen of the latent class with x
               C1.matrix <- X.matrix * alpha1
               C2.matrix <- X.matrix * alpha2
               C3.matrix <- X.matrix * alpha3
            ## Observed covariates 
              X.obs <- matrix(nrow=N,ncol=num.cov)
              X.obs[tx1.ind,] <- C1.matrix[tx1.ind,]
              X.obs[tx2.ind,] <- C2.matrix[tx2.ind,]
              X.obs[tx3.ind,] <- C3.matrix[tx3.ind,]
              colnames(X.obs) <- paste("X",1:num.cov,sep="")
            ### Generate response variables
            ### Linear and Nonlinear types of variable 
            if(response=="linear") {
            Y0<-beta1*X3+beta1*X4+beta1*X5+beta1*X6+beta1*X7+beta1*X8+rnorm(N)
            Y1<-2+beta1*X3+beta2*X4+beta2*X5+beta2*X6+beta2*X7+beta2*X8+rnorm(N)
            Y2<-4+beta3*X3+beta3*X4+beta3*X5+beta3*X6+beta3*X7+beta3*X8+rnorm(N)} 
            else{
            Y0<-beta1*X3+beta1*X4+beta1*X5+beta1*X6^2+beta1*X7^2+beta1*X8+X2*X5+rnorm(N)
            Y1<-2+beta2*X3*X4+beta2*X5+beta2*6*X7+beta2*X8+rnorm(N) 
            Y2<-4+beta3*X3+beta3*X4^2+beta3*X5+beta3*X6+beta3*X7*X8+rnorm(N)}
            Y_obs<-ifelse(tx==1,Y0,ifelse(tx==2,Y1,Y2))
            ### Generate latent class indicators
            U.matrix <- matrix(nrow = N, ncol= num.item)
            ### Create "Potential" latent class indicators for C1, C2, C3
            item1.matrix <- array(rbinom(n=N*num.item,1,prob=p[1]), dim=c(N,num.item)) 
            item2.matrix <- array(rbinom(n=N*num.item,1,prob=p[2]), dim=c(N,num.item)) 
            item3.matrix <- array(rbinom(n=N*num.item,1,prob=p[3]), dim=c(N,num.item)) 
            ### Convert 0 --> 2 for LCCA package
            item1.matrix[item1.matrix==0] <- 2
            item2.matrix[item2.matrix==0] <- 2
            item3.matrix[item3.matrix==0] <- 2
            ### Observe indicators that correspond to true tx class
            U.matrix[tx1.ind,] <- item1.matrix[tx1.ind,]
            U.matrix[tx2.ind,] <- item2.matrix[tx2.ind,]
            U.matrix[tx3.ind,] <- item3.matrix[tx3.ind,]
            colnames(U.matrix) <-paste("U",1:num.item,sep="")
            data <- data.frame(cbind(tx,X.obs,U.matrix,Y0,Y1,Y2,Y_obs))}

###Simulation according to poLCA package 
install.packages('poLCA')
library("poLCA")

### Data Simulation 
Latent.GDP.2<-function(N,      # Number of observations 
                       ndv=15, # Number of indicator variable 
                       niv=6,  # Number of variable related to treatment: 6+2 
                       p11=0.5,# Conditional probability P(U|C)
                       p12=0.3,
                       p21=0.4,
                       p22=0.1,
                       p31=0.2,
                       p32=0.2,
                       nclass=3,# Number of class fix to be 3
                       alpha1.ins,alpha1.conf,# Association between covariate and treatment 
                       alpha2.ins,alpha2.conf,# Association between covariate and treatment 
                       beta1,beta2,beta3,beta4,beta5,beta6,# Association between covariate and outcome 
                       mu=c(1,1.5,2),# The Treatment effect size 
                       response="linear") # Response surface
     { 
     ### Coeffcient of the x variables 
     b<-matrix(c(c(-1,alpha1.ins,alpha1.ins,alpha1.ins,alpha1.conf,alpha1.conf,alpha1.conf),
                c(-1,alpha2.ins,alpha2.ins,alpha2.ins,alpha2.conf,alpha2.conf,alpha2.conf)),
                 ncol=nclass-1,nrow=niv+1,byrow=FALSE)
     ### The conditional probability P(U|C) where C is the latent class and U is the indcator variable 
     ### replicate 15
     probs <- list(
      matrix( c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
                ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
                 ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
               ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
              ncol = 3, byrow = TRUE),
       matrix(c(p11,p12,1-p11-p12,
                p21,p22,1-p21-p22,
                p31,p32,1-p31-p32),
                ncol = 3, byrow = TRUE))
    ### Simulate the data 
        sim.data <- poLCA.simdata(N=N,probs=probs,ndv=ndv,niv=niv,nclass=nclass,b=b)
        names(sim.data$dat)[1:ndv]<-paste("U",1:ndv,sep="")
        latent.data<-cbind(trueclass=sim.data$trueclass,sim.data$dat)
        X7<-rnorm(1000)
        X8<-rnorm(1000)
        attach(latent.data)
    ### Response surfaces 
        if(response=="linear") {
        Y0<-mu[1]+beta1*X3+beta2*X4+beta3*X5+beta4*X6+beta5*X7+beta6*X8+rnorm(N)
        Y1<-mu[2]+beta1*X3+beta2*X4+beta3*X5+beta4*X6+beta5*X7+beta6*X8+rnorm(N)
        Y2<-mu[3]+beta1*X3+beta2*X4+beta3*X5+beta4*X6+beta5*X7+beta6*X8+rnorm(N)} 
        else{
        Y0<-mu[1]+beta1*X3+beta2*X4+beta3*X5+beta4*X6^2+beta5*X7^2+beta6*X8^2+rnorm(N)
        Y1<-mu[2]+beta1*X3*X4+beta2*X5+beta3*X6*X7+beta4*X8+beta5*X7+beta6*X8+rnorm(N) 
        Y2<-mu[3]+beta1*X3^2+beta2*X4^2+beta3*X5+beta3*X4*X5+beta5*X7^2+beta6*X8*X3+rnorm(N)}
        Y_obs<-ifelse(trueclass==1,Y0,ifelse(trueclass==2,Y1,Y2))
        detach(latent.data)
        data<-data.frame(latent.data,X7,X8,Y0,Y1,Y2,Y_obs)
        data.all<-list(data,alpha=round(sim.data$P,2))
        return(data.all)}

### Testing Simulation Function 
test.data.1<-Latent.GDP.2(1000,
                        15,6,
                        0.5,0.4,
                        0.5,0.2,
                        0.4,0.3,
                        3,
                        0.5,1,
                        0.5,1,
                        1,1,
                        1,1,
                        1,1,
                        c(0,2,4),"linear")


data.list<-vector('list',100)

for(i in 1:100)
{
  data.list[[i]]<-Latent.GDP.2(1000,
                                 15,6,
                                 0.5,0.4,
                                 0.5,0.2,
                                 0.4,0.3,
                                 3,
                                 2,2,
                                 2,2,
                                 1,1,
                                 1,1,
                                 1,1,
                                 c(0,2,4),"linear")
}


test.data.1<-Latent.GDP.2(1000,
                          15,6,
                          0.5,0.4,
                          0.5,0.2,
                          0.4,0.3,
                          3,
                          0.5,1,
                          0.5,1,
                          1,1,
                          1,1,
                          1,1,
                          c(0,2,4),"linear")

number<-c(sum(test.data.1[[1]]$trueclass==1),sum(test.data.1[[1]]$trueclass==2),sum(test.data.1[[1]]$trueclass==3))
number



class<-matrix(NA,nrow = 100,ncol = 3)
for(i in 1:100)
{
  class[i,]<-c(sum(data.list[[i]][[1]]$trueclass==1),sum(data.list[[i]][[1]]$trueclass==2),sum(data.list[[i]][[1]]$trueclass==3))
}



test.data.2<-Latent.GDP.2(1000,
                         15,6,
                         0.5,0.4,
                         0.5,0.2,
                         0.4,0.3,
                         3,
                         0.5,1,
                         0.5,1,
                         1,1,
                         1,1,
                         1,1,
                         c(0,2,4),"nonlinear")


### Test the linearity of the dat aset 
summary(lm(Y_obs~as.factor(trueclass)+X1+X2+X3+X4+X5+X6+X7+X8,test.data.1[[1]])) 
summary(lm(Y_obs~as.factor(trueclass)+X1+X2+X3+X4+X5+X6+X7+X8,test.data.2[[1]])) 

### get the latent class 

f3a <- cbind(U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,U12,U13,U14,U15)~X1+X2+X3+X4+X5+X6+X7+X8
lc3a <- poLCA(f3a,data,nclass=3,maxiter = 100000,verbose=FALSE)
prob.start<-poLCA.reorder(lc3a$probs.start,order(lc3a$P,decreasing = TRUE))
lc3a<- poLCA(f3a,data,nclass=3,maxiter = 100000,probs.start =prob.start )
length(which(data$trueclass-lc3a$predclass!=0))

pred.treat<-lc3a$predclass
data$treat<-pred.treat
cov<-data[,17:24]
fit<-Vanilla.bart(data=data,xcov=cov,data$treat)

ATE.1.0<-mean(fit[,2]-fit[,1])
ATE.1.0
ATE.2.1<-mean(fit[,3]-fit[,2])
ATE.2.1
ATE.2.0<-mean(fit[,3]-fit[,1])
ATE.2.0


mean(test.data.2$Y1-test.data.2$Y0)
mean(test.data.2$Y2-test.data.2$Y1)
mean(test.data.2$Y2-test.data.2$Y0)

