### 1-step LCCA function

onestep.lcca <- function(data, num.tx=3, num.item=15, iter=100000){
                ### Calculate true treatment effects
                diff.21<-mean(data$Y1-data$Y0)
                diff.31<-mean(data$Y2-data$Y0)
                diff.32<-mean(data$Y2-data$Y1)
                
  ### Adjusted LCCA model
  fit.lcca <- lcca(formula.treatment=cbind(U1,U2,U3,U4,U5)~X1+X2+X3+X4+X5+X6+X7+X8,
                   formula.outcome=Y_obs~X1+X2+X3+X4+X5+X6+X7+X8,data=data,
                   nclass=num.tx,iter.max=iter, flatten.rhos=1,flatten.gamma=1,
                   stabilize.alphas=1, se.method="STANDARD")
  ### Resolve class switching
  cond.mean <- NULL
  # Sum conditional probs (rhos) within each class to determine class names
  cond.mean[1]<-sum(fit.lcca$theta[1:num.item])
  cond.mean[2]<-sum(fit.lcca$theta[(num.item+1):(2*num.item)])
  cond.mean[3]<-sum(fit.lcca$theta[(2*num.item+1):(3*num.item)])

  # Check if classes are well differentiated
  test <- pairwise.diffs(cond.mean)
  if(sum(abs(test)<0.25)>0) print("Warning: Indistinguishable Classes")
  #Reorder classes and summarize lcca
  ord <- order(cond.mean)
  fit.lcca.ord <- permute.class(fit.lcca, ord)
  ### Save Average Treatment Effects, Bias, SE, and MSE
  os.2v1<-summary(fit.lcca.ord)[94][[1]][1]
  
  os.3v1<-summary(fit.lcca.ord)[94][[1]][2]
  
  os.3v2<-summary(fit.lcca.ord)[94][[1]][3]
  
  # Save Class Prevalences
  p1<-fit.lcca.ord$marg.gamma[1]
  p2<-fit.lcca.ord$marg.gamma[2]
  p3<-fit.lcca.ord$marg.gamma[3]
  output <- NULL
  output<-as.numeric(c(os.2v1,
                       os.3v1,
                       os.3v2))
  return(output)
}

### Pairwise difference function
pairwise.diffs <- function(x)
{ # create matrix of combination pairs
  prs <- cbind(rep(1:length(x), each = length(x)), 1:length(x))
  # drop ones that compare same classes
  drops <- NULL
  for(i in 1:(length(x))^2 )
  { if (prs[i,1]==prs[i,2]) drops <- c(drops, i)
  }
  new <- prs[-drops,]
  # do pairwise differences
  result <- x[new[,1]] - x[new[,2], drop = FALSE]
}


onestep.lcca(data=test)
