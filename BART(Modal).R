### Import packages that are needed ###
Y<-factor(data$H4EC1)
Y<-as.numeric(Y)
names(data_modal)
cov<-data[,c(21:34)]
names(cov)<-paste("X",1:14,sep="")
modal_treat<-as.numeric(data_modal$Modal)
pseudo_treat<-as.numeric(data_pseudo$Modal)
test<-bartc(response=Y,treatment=modal_treat, confounders=cov, commonSup.rule="sd")
summary(test)



### This function get the upper bound of the reponse surfaces
get_cutoff<-function(x){
  cutoff<-max(x)+sd(x)
  return(cutoff)}


install.packages("devtools")
devtools::install_github("vdorie/bartCause")

### This function will help to delete the non-overlap for predicted observations
check_olp<-function(data,FUN=get_cutoff)
           {mark<-apply(data,1,FUN)
            ### This loop will go over each observations to compare it with upper bound
            for (i in 1:nrow(data)){
              for (j in 1:ncol(data)){
                if(data[i,j]>as.vector(mark[i]))
                {data[i,j]<-NA
                print("warning:non-overlap case detected")}}
              }
             return(na.omit(data))}

### pairwise difference function
pairwise.diffs <- function(x)
{ stopifnot(is.matrix(x))
  # Create column combination pairs
  prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
  col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  # Do pairwise differences 
  result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
   # Set colnames
  if(is.null(colnames(x)))
  colnames(x) <- 1:ncol(x)
  colnames(result) <- paste(colnames(x)[col.diffs[, 1]], ".vs.",  colnames(x)[col.diffs[, 2]], sep = "")
  return(result)
  }


### Vanilla BART function 
BART<-function(data,FUN=check_olp,size)
    { y_obs<-data$Y
      treat<-data$treat
      cov<-data[,c(grep("X",names(data)))]
      xt<-cbind(cov,treat)
      ### Step 1: Fit a whole BART model 
      bart_mod<-bart(x.train = xt,y.train = y_obs,keeptrees = TRUE)
      plot(bart_mod)
      n.class<-length(levels(as.factor(treat)))
      y.pred<-list()
      ### Step 2: Impute counter-factuals based on observed outcome
      for(i in 1:n.class){
        xp<-xt
        xp[,grep("treat",names(xp))]<-i
        bart_pred<-predict(bart_mod,xp)
        y.pred[[i]]<-bart_pred}
      ### Step 3: Assess the overlap by using Hill and Su, 2013 1.st SD rule
      trim_data<-lapply(y.pred,FUN)## non-overlap is trimmed
      ### Step 4: get sereis of ATES based on assignments 
      result<-matrix(rep(NA,1000*size),ncol=size)
      for (i in 1:1000){
        if(n.class==2){
                    mat<-matrix(c(trim_data[[1]][i,],trim_data[[2]][i,]),ncol=n.class)}
        else if(n.class==3)
                    {mat<-matrix(c(trim_data[[1]][i,],trim_data[[2]][i,],trim_data[[3]][i,]),ncol=n.class)}
                    diff<-pairwise.diffs(mat)
                    result[i,]<-apply(diff,2,mean)
                   colnames(result)=colnames(diff)}
      result<-as.data.frame(result)
      result<-result[!duplicated(as.list(result))]
      return(result)}
#### model_BART
Modal_BART<-function(data,size){
  ATEs <- BART(data,FUN=check_olp,size)
  ## Estimate mean ATES and variance ATE based on assignments 
  mean_ATE<-apply(ATEs,2,mean)
  names(mean_ATE)<-paste("mean:",names(mean_ATE))
  sd_ATE<-apply(ATEs,2,sd)
  names(sd_ATE)<-paste("sd:",names(sd_ATE))
  return(c(mean_ATE,sd_ATE))
         }
### Real data analysis 
## Clean and prepare the data set 
data<-read.csv(file=file.choose(),header = TRUE)
data_modal <- read.csv(file=file.choose(),header = TRUE) 
data_pseudo<-read.csv(file=file.choose(),header=TRUE) 
data_pseudo_1<-read.csv(file=file.choose(),header=TRUE) 
Y<-factor(data$H4EC1)
Y<-as.numeric(Y)
names(data_modal)
cov<-data[,c(21:34)]
names(cov)<-paste("X",1:14,sep="")
modal_treat<-as.numeric(data_modal$Modal)
pseudo_treat<-as.numeric(data_pseudo$Modal)
pseudo_treat_1<-as.numeric(data_pseudo_1$Modal)
merge_data_modal<-data.frame(cov,treat=modal_treat,Y=Y)
merge_data_pseudo<-data.frame(cov=cov,treat=pseudo_treat,Y=Y)
merge_data_pseudo_1<-data.frame(cov=cov,treat=pseudo_treat_1,Y=Y)
### Run the model 
Modal_BART(merge_data_modal,size=1)
Modal_BART(merge_data_pseudo_1,size=1)


## Get confidence interval 
Test<-BART(merge_data_modal,size=1)
quantile(Test$`1.vs.2`, c(0.25, .95)) 

Modal_BART(merge_data_pseudo,size=3)

Test_1<-BART(merge_data_pseudo,size=3)
quantile(Test_1$`1.vs.2`, c(0.25, .95)) 
quantile(Test_1$`1.vs.3`, c(0.25, .95)) 
quantile(Test_1$`2.vs.3`, c(0.25, .95)) 




