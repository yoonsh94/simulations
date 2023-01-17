
library(pacman)
pacman::p_load(tidyverse,magrittr,flexsurv,dplyr,tidyr,tidyverse,data.table,plyr,TMB,survey,lme4)

sim_data <- function(sigmaA.true, sigmaB.true){
  set.seed(3105)
  # Number of subjects
  N <- 10000 

  # Unconditional true parameters (Generalised)
  alpha.true<- c(4,1,1,1) 
  beta.true  <- c(2,1,1,1)
  gamma.true <- c(1.5, 1, 1)
  lambda1.true <- 0.7 
  lambda2.true <- 0.5
  lambda3.true <- 0.5
  delta.true <- 0.5

  # Replicates 
  mean.repl <- 10 # or 20?
  replicates <- rpois(N,mean.repl)+3
  idK <- sort(rep(1:N, replicates))
  
  Ntotal <-   sum(replicates)
  
  X1 <- rbinom(Ntotal, 1, 0.5)
  X2 <- rnorm(Ntotal, X1, 0.4) 
  X3 <- rnorm(Ntotal, X1, 0.3) #X1 is now correlated with X2 and X3
  X <- cbind(1,X1,X2,X3)
  nparam_a <- ncol(X) + 1 # 4 alphas + logsigma_a
  nparam_b <- ncol(X) + 3 # 4 betas + logsigma_a + lambda1 + delta
  
  nalpha <- ncol(X)
  nbeta <- ncol(X)
  
  #for both binary and normal case?
  Z1 <- rbinom(N, 1, 0.4) #no repeated measure for Z.  Z1 is directly usable in cpp 
  Z2 <- rnorm(N, Z1, 0.4) 
  Z3 <- rnorm(N, Z1, 0.3)
  Z.1 <- Z1[idK]
  Z.2 <- Z2[idK]
  Z.3 <- Z3[idK]
  Z  <- cbind(Z.1, Z.2, Z.3)
  Z_one <- cbind(Z1, Z2, Z3)
  nparam_c <- ncol(Z) + 12 # 3 gammas + lambda2 + lambda3 + 10 h's
  
  ngamma <- ncol(Z)
  #partially missing variable available only for phase 2 sample
  #nth<-1
  #nth_missing <- cumsum(c(nth,replicates))[1:N]
  
  #aux1<-rbinom(Ntotal,1,ifelse(X1== 1 | X2 ==1, 0.9,0.1)) # binomial case
  #aux2<-rbinom(Ntotal,1,ifelse(X1+X2>=1, 0.95,0.05)) # binomial case
  #aux3<- rbinom(Ntotal, 5, ifelse(X1==1, 0.9,0.1)) #discrete case
  aux1<-rbinom(Ntotal,1,ifelse(X1== 1 | X2 ==1, 0.9,0.1)) # binomial case
  aux2<-rbinom(Ntotal,1,ifelse(X1+X2>=1, 0.95,0.05)) # binomial case
  aux3<- rnorm(Ntotal, X1, 0.2) #discrete case
  aux <- cbind(aux1,aux2,aux3)
  
  #aux vars for survival model
  Z_aux1 <- rbinom(N, 1, ifelse(Z1==1, 1,0)) #binomial
  Z_aux2 <- rbinom(N, 5, ifelse(Z1==1, 0.9,0.1)) #categ.
  #Z_aux3 <- rbinom(N, 1, ifelse(Z1==1, 0.95,0.05))
  Z_aux3 <- rnorm(N, Z1, 0.3) #continuous
  Z_aux <- cbind(Z_aux1, Z_aux2, Z_aux3) 
  
  numbaux <- 3 #number of aux vars
  
  
  # Random effects
  a<-rnorm(N,0,sigmaA.true)
  b<-rnorm(N,0,sigmaB.true)
  a_All <- a[idK]  #repeating baseline random effects
  b_All <- b[idK]
  
  
  # Conditional true parameters   
  eta.true <- X %*% alpha.true +a_All
  mu.true <- exp(X %*% beta.true + lambda1.true*a_All + b_All)
  sigmaiij2.true <- exp(delta.true)
  
  # Prepare a data frame with identifiers and covariates
  df <- data.frame(Identifier=idK,X=X[,-c(1)]) 
  df$YNoZero  <-  rbinom(Ntotal, 1, 1/(1+exp(-eta.true)))
  
  df$cost <- rgamma(Ntotal,shape=sigmaiij2.true, scale=mu.true/sigmaiij2.true)
  df <- df %>% mutate(cost2 =  cost*YNoZero )
  #df <- df %>% mutate(cost2 =  cost) 
  
  #survival part
  
  l0 <- 2  # hazard is exponential 
  U2     <-   runif(N, 0,1)[idK]
  idKu  <- unique(idK)
  survt0 <-survt <- (-log(U2)/(l0*exp(Z %*% gamma.true + lambda2.true*a_All + lambda3.true*b_All)))
  
  # censoring times:
  survt.censor <- rexp(N, rate= 2)[idK]
  # censoring:
  di    <- as.numeric(survt0 <= survt.censor);
  survt <- pmin(survt0, survt.censor)  #survt is the minimum of survival time and censoring time
  df$survt <- survt
  df$di <- di
  
  #df$cost  <- df$cost* df$di
  
  df$idK <- idK
  df$num <- unlist(sapply(1:N, function(v) 1:replicates[v]))
  
  #auxiliary
  df <- cbind(df,aux)
  
  dfS <- df[!duplicated(df$idK, fromLast = TRUE),]
  
  survt1<-dfS$survt #one survt per indiv
  di1 <- dfS$di
  replicates <- dfS$num #replicates after removing obs after event
  
  
  maxrep<-max(replicates)
  numbX<-ncol(X)-1 #excluding the first column of 1's
  
  
  setDT(df)
  
  #dcast(df, Identifier~num, value.var = c("cost2","X1","X2","X3"))
  df1 <- cbind(dcast(df, Identifier~num, value.var = c("cost2",paste0("X.X",1:numbX),paste0("aux",1:numbaux))),dfS[,c("survt","di","num","YNoZero")]) #complete wide-form dataset with nrow=N
  
  
  start<-cumsum(c(1,replicates))
  end<-cumsum(replicates)
  #modular now
  df1$X <- rowMeans(sapply(1:numbX, function(i) sapply(1:N, function(v)
    mean(get(paste0("X",i))[start[v]:end[v]])
  )))
  
  
  ncoldf <-ncol(df1)
  nX1col <- (1:ncoldf)[colnames(df1) %in% paste('cost2', 1:maxrep, sep='_')]
  dfcost<-as.matrix(df1[,nX1col , with = FALSE])
  costindsum<-list()
  costsum<-list()
  for(i in 1:N){
    costindsum[[i]] <- sum(dfcost[i,][1:replicates[i]] > 0)
    costsum[[i]] <- sum(dfcost[i,][1:replicates[i]] )
  }
  df1$costindexsum <- unlist(costindsum) #df1 is the dataset to be sampled. This has column-wise longitudinal measures
  df1$costsum <- unlist(costsum)
  df1 <- cbind(df1,Z_one,Z_aux) #Zs and aux_Zs (subject-level) integrated to the data. 
  
  returnList <- list("df1"=df1, "maxrep"=maxrep,"nparam_a"=nparam_a, "nparam_b"=nparam_b, 
                     "nparam_c"=nparam_c,"nalpha"=nalpha,"nbeta"=nbeta,"ngamma"=ngamma,"numbaux"=numbaux,"numbX"=numbX)
  return(returnList)
}

# data to be sampled and some objects required in the later stage
sim_data_list<-sim_data(1,1)
df1<-sim_data_list$df1
maxrep<-sim_data_list$maxrep
nparam_a<-sim_data_list$nparam_a
nparam_b<-sim_data_list$nparam_b
nparam_c<-sim_data_list$nparam_c
nalpha <- sim_data_list$nalpha
nbeta <- sim_data_list$nbeta
ngamma <- sim_data_list$ngamma
numbaux<-sim_data_list$numbaux
numbX<-sim_data_list$numbX


#Sampling from population
n.strata <- 20 #each stratum has 10 clusters and 500 obs
n.cluster <- 10 #for each cluster within a stratum there are 50 obs.  A total of 200 clusters.

n.per.cluster <- 10 # n.per.cluster * 2 * n.strata = sample size   #this adjusts the weight
sample_size <- n.per.cluster * 2 * n.strata


total.cluster <- n.strata * n.cluster
df1$strata<-rep(seq(1:n.strata), each=N/n.strata) 
df1$cluster<-rep(seq(1:n.cluster),N/n.cluster)

samplingfunc <- function(iter,data,n.per.cluster){

  
  phase2.df1 = list()
  phase2.df2 = list()
  datalist1=list()
  datalist2=list() 
  X.rowlist<-list()
  costlist<-list()
  datares<-list()
  for(i in 1:n.strata){
    strata.df <- subset(data, data$strata == i) #each stratum has 500 
    cluster.label <- sample(c(1:n.cluster), 2, prob=rep(1/n.cluster,n.cluster))
    
    cluster.df1 <- subset(strata.df, strata.df$cluster == cluster.label[1])
    cluster.df2 <- subset(strata.df, strata.df$cluster == cluster.label[2])
    
    cluster.df1$size <- (0.25+0.5*cluster.df1$X)*0.5  #size measure
    cluster.df1$pps <- cluster.df1$size / sum(cluster.df1$size) #selection prob
    #cluster.df1$wght <- 1/cluster.df1$prob    #sampling wght
    
    
    cluster.df2$size<- (0.25+0.5*cluster.df2$X)*0.5 #size measure
    cluster.df2$pps <- cluster.df2$size / sum(cluster.df2$size) #pps
    #cluster.df2$wght <- 1/cluster.df2$prob   #sampling wght
    
    indiv.selected <- sample(c(1:(N/total.cluster)),size=n.per.cluster,prob=c(cluster.df1$prob),
                             replace=FALSE)
    #select n.per.cluster individuals based on selection probs
    cluster.df1.indiv <- cluster.df1[indiv.selected,]
    cluster.df1.indiv$prob <- 0.2 * (n.per.cluster*cluster.df1.indiv$size)/ sum(cluster.df1$size)
    cluster.df1.indiv$weight <- 1/cluster.df1.indiv$prob
    
    indiv.selected2 <- sample(c(1:(N/total.cluster)),size=n.per.cluster,prob=c(cluster.df2$prob),replace=FALSE)
    cluster.df2.indiv <- cluster.df2[indiv.selected2,]
    cluster.df2.indiv$prob <- 0.2 * (n.per.cluster*cluster.df2.indiv$size)/ sum(cluster.df2$size)
    cluster.df2.indiv$weight <- 1/cluster.df2.indiv$prob
    
    # Xrow subset using identifier cbind(data[ident,],Xlist[[1]][ident,])
    
    datalist1[[i]] <- cluster.df1.indiv
    datalist2[[i]] <- cluster.df2.indiv
    #sample further from cluster.df1 and df.2 (subsample or phase 2 data)
    
    
    phase2.selected1 <- sample(c(1:n.per.cluster), size=5, prob=rep(1/n.per.cluster, n.per.cluster),replace=FALSE)
    phase2.selected2 <- sample(c(1:n.per.cluster), size=5, prob=rep(1/n.per.cluster, n.per.cluster),replace=FALSE)
    
    phase2.data1 <- cluster.df1.indiv[phase2.selected1,]
    phase2.data2 <- cluster.df2.indiv[phase2.selected2,]
    phase2.data1$weight2 <- phase2.data1$weight * 2
    phase2.data2$weight2 <- phase2.data2$weight * 2
    phase2.df1[[i]] <- phase2.data1
    phase2.df2[[i]] <- phase2.data2
    
    
  }
  
  data1 <- do.call("rbind", datalist1)
  data2 <- do.call("rbind", datalist2)
  
  # get phase 2 sample here  by simply  selecting individuals from the complete phase 1 data
  
  phase2_1 <- do.call("rbind", phase2.df1)
  phase2_2 <- do.call("rbind", phase2.df2) 
  phase2 <- rbind(phase2_1, phase2_2)
  data.final <-  rbind(data1,data2)
  
  datares<-list(data.final)
  phase2<-list(phase2)
  
  list(datares=datares,phase2=phase2)
  #returnList <- list("phase1"=data.final,"phase2"=phase2)
  #return(returnList)
}

# sampling  n times to generate dataset
par_est <- list()
n.iter<-1000
sampledata<-lapply(1:n.iter, function(v) samplingfunc(v,df1,n.per.cluster=10))




# individual-laplace cpp templates

compile("la_general.cpp",framework="TMBad")
dyn.load(dynlib("la_general"))

compile("lb_general.cpp",framework="TMBad")
dyn.load(dynlib("lb_general"))

#no randeff
compile("lc_general.cpp")
dyn.load(dynlib("lc_general"))


#Variance estimation function for weighted 

var_est <- function(data,gr_all, nparam_a,nparam_b,nparam_c,weights,sample_size,n.indiv,num.cluster){
  param_list<-list()
  nparam <- nparam_a + nparam_b + nparam_c
  for(i in 1:nparam){
    
    param_list[[i]] <- gr_all[(sample_size*i+1):(samples_size*(i+1))]
  }
  mat_gr <- data.frame(do.call(cbind,param_list),data[[1]]$datares[[1]]$cluster,data[[1]]$datares[[1]]$strata, data[[1]]$datares[[1]]$weight)
  mat_gr1 <- weights*mat_gr[1:(sample_size/2),]
  mat_gr2 <- weights*mat_gr[(sample_size/2+1):sample_size,]
  #cpp template returns unweighted  gradient
  
  #variance estimation
  var_ind <- seq(1,sample_size+1,n.indiv)
  G <- matrix(0,nrow=nparam, ncol=nparam)
  for(i in 1:n.strata){ #n.strata = 100
    e_hl_1 <- colSums(mat_gr1[,1:nparam][var_ind[i]:var_ind[i+1]-1,]) #1st cluster sum in str h
    e_hl_2 <- colSums(mat_gr2[,1:nparam][var_ind[i]:var_ind[i+1]-1,]) #2nd cluster sum in str h
    
    e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
    G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                 + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
    
  }
  invhess <- solve(opt$hessian) #this is weighted Hessian
  se <- sqrt(diag(invhess %*% G %*% invhess)) 
  return(se)
}




# the entire estimation using a single function


adfun <- function(dat,n.iter,sample_size, nparam_a,nparam_b,nparam_c,nalpha,nbeta,ngamma, n.indiv, num.cluster, n.strata, type) { 
  
  # data passed to this function is a list of data frames (1000)
  # can have sampling function here  so that  sample size can be also modular (by + or - the number of subjects selected in each cluster) if needed later
  # make the function extract gradient functions only when the type is weighted or calibrated because gradient functions aren't required for unweighted analysis
  parest_a <- matrix(NA,nrow=n.iter,ncol=nparam_a)
  parest_b <- matrix(NA,nrow=n.iter,ncol=nparam_b)
  parest_c <- matrix(NA,nrow=n.iter,ncol=nparam_c)
  na_ind_a <- c()
  na_ind_b <- c()
  na_ind_c <- c()
  SE <- matrix(NA,nrow=n.iter,ncol=nparam_a+nparam_b+nparam_c) #separately or jointly?
  SE_a <- matrix(NA,nrow=n.iter,ncol=nparam_a)
  SE_b <- matrix(NA,nrow=n.iter,ncol=nparam_b)
  SE_c <- matrix(NA,nrow=n.iter,ncol=nparam_c)
  
  for(i in 1:n.iter){
    
    if (type == "UW"){
      xmat<-list()
      for(j in 1:numbX){
        assign(paste0("colname_X", j), (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste0('X.X', j, "_",1:maxrep)])
        assign(paste0("X",j,"row"),as.matrix(dat[[i]]$datares[[1]][,get(paste0("colname_X",j)) , with = FALSE]) )
        xmat[[j]] <- as.vector(t(get(paste0("X",j,"row"))))
      }
      
      xmat<-do.call("cbind",xmat)
      X_sample<-cbind(1,xmat[!rowSums(!is.finite(xmat)),])
      
      ncoldfZ <-ncol(dat[[i]]$datares[[1]])
      colname_costZ <- (1:ncoldfZ)[colnames(dat[[i]]$datares[[1]]) %in% paste0('Z', 1:numbX)] 
      #to get Z matrix. Ignore the name 'cost'
      zmat <- as.matrix(dat[[i]]$datares[[1]][,colname_costZ, with =FALSE])
      
      # cost
      
      ncoldf <-ncol(dat[[i]]$datares[[1]])
      colname_cost <- (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste('cost2', 1:maxrep, sep='_')]
      
      dfcost<-as.matrix(dat[[i]]$datares[[1]][,colname_cost , with = FALSE])
      survt1 <-   dat[[i]]$datares[[1]]$survt
      di1 <-  dat[[i]]$datares[[1]]$di
      costindexsum<- dat[[i]]$datares[[1]]$costindexsum
      costsum<- dat[[i]]$datares[[1]]$costsum
      
      weight<-rep(1,sample_size) #unweighted so weight = 1
      
      qq <-10  #10 intervals
      Q.partition <- c(0,quantile(survt1,1:qq/qq)) # length 11, need to be calculated for entire sample each iteration
      cumhaz_int <- c(Q.partition[2],(Q.partition[3:(qq+1)]-Q.partition[(2):(qq)]))
      index_vec<-sapply(1:sample_size, function(v) min((1:length(Q.partition)-2)[survt1[v]<=Q.partition]))
      #(1:length(Q.partition))[survt1[1]<=Q.partition]) 
      
      replicates <- dat[[i]]$datares[[1]]$num  # num. repeated measurements per subject
      maxrep_tmb <- max(replicates)
      start<-c(0,cumsum(replicates))
      
      
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=weight,a1=0, start=start)
      parameters <- list(sigma=c(0,rep(1,nalpha)))  # sigma=c(0,rep(1,nparam_a-1))
      
      obj <- MakeADFun(data,parameters,DLL="la_general", method="L-BFGS-B") 
      #L-BFGS-B always works for lA when it doesnt work without it, LB works fine for weights
      obj$hessian <- TRUE
      opt<-suppressWarnings(do.call("optim",obj))
      parest_a[i,] <- opt$par #parameter estimates
      na_ind_a[i] <- is.na(opt$value)
      
      if(is.na(opt$value)){
        next
      }
      
      rep<-sdreport(obj)
      repsum <- summary(rep, "report") #estimated a's
      a_est<-repsum[1:sample_size]
      la_gr <- repsum[(sample_size+1):(sample_size*nparam_a),] # excluding a's, only the gr fn
      SE_a[i,] <- sqrt(diag(solve(opt$hessian)))  # check whether logsigma need to be exponentiated 
      
      #model B individual laplace
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=weight,a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      
      parameters <- list(sigma=c(0,1,1,rep(1,nbeta))) #make this modular  use eg nparam_b and rep
      
      obj <- MakeADFun(data,parameters,DLL="lb_general")
      obj$hessian <- TRUE
      opt<-suppressWarnings(do.call("optim",obj))
      parest_b[i,] <- opt$par
      na_ind_b[i] <- is.na(opt$value)
      
      if(is.na(opt$value)){
        next
      }
      
      rep<-sdreport(obj)
      repsum <- summary(rep, "report") #estimated b's
      b_est<-repsum[1:sample_size]
      lb_gr <- repsum[(sample_size+1):(sample_size*nparam_b),]
      
      
      SE_b[i,] <- sqrt(diag(solve(opt$hessian)))  # check whether logsigma need to be exponentiated 
      #lC
      
      data<-list(a1=a_est,b1=b_est, Z0=zmat,N=sample_size,Qpartition=Q.partition,survt1=survt1,di1=di1,cumhaz_int=cumhaz_int,index_vec=index_vec,weights=weight)
      
      parameters <- list(theta=c(rep(1,nparam_c))) #list(gamma=rep(1,ncol(Z)),lambda2=1,lambda3=1,h1=1,h2=1,h3=1,h4=1,h5=1,h6=1,h7=1,h8=1,h9=1,h10=1)
      
      obj <- MakeADFun(data,parameters,DLL="lc_general")
      obj$hessian <- TRUE
      opt <- suppressWarnings(do.call("optim",obj))
      parest_c[i,] <- opt$par
      na_ind_c[i] <- is.na(opt$value)
      rep<-sdreport(obj)
      repsum <- summary(rep, "report")
      lc_gr <- repsum[(sample_size+1):(sample_size*nparam_c),]
      
      gr_all <- rbind(la_gr,lb_gr,lc_gr)
      
      
      SE_c[i,] <- sqrt(diag(solve(opt$hessian)))  # check whether logsigma need to be exponentiated 
      SE[i,] <- c(SE_a[i],SE_b[i],SE_c[i]) #assigning a joint SE vector to first row of SE matrix
      
      
      
      
      
    }
    
    if (type == "W"){
      xmat<-list()
      for(j in 1:numbX){
        assign(paste0("colname_X", j), (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste0('X.X', j, "_",1:maxrep)])
        assign(paste0("X",j,"row"),as.matrix(dat[[i]]$datares[[1]][,get(paste0("colname_X",j)) , with = FALSE]) )
        xmat[[j]] <- as.vector(t(get(paste0("X",j,"row"))))
      }
      
      xmat<-do.call("cbind",xmat)
      X_sample<-cbind(1,xmat[!rowSums(!is.finite(xmat)),])
      
      ncoldfZ <-ncol(dat[[i]]$datares[[1]])
      colname_costZ <- (1:ncoldfZ)[colnames(dat[[i]]$datares[[1]]) %in% paste0('Z', 1:numbX)] 
      #to get Z matrix. Ignore the name 'cost'
      zmat <- as.matrix(dat[[i]]$datares[[1]][,colname_costZ, with =FALSE])
      
      # cost
      
      ncoldf <-ncol(dat[[i]]$datares[[1]])
      colname_cost <- (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste('cost2', 1:maxrep, sep='_')]
      
      dfcost<-as.matrix(dat[[i]]$datares[[1]][,colname_cost , with = FALSE])
      survt1 <-   dat[[i]]$datares[[1]]$survt
      di1 <-  dat[[i]]$datares[[1]]$di
      costindexsum<- dat[[i]]$datares[[1]]$costindexsum
      costsum<- dat[[i]]$datares[[1]]$costsum
      
      weight<- dat[[i]]$datares[[1]]$weight
      
      
      qq <-10  #10 intervals
      Q.partition <- c(0,quantile(survt1,1:qq/qq)) # length 11, need to be calculated for entire sample each iteration
      cumhaz_int <- c(Q.partition[2],(Q.partition[3:(qq+1)]-Q.partition[(2):(qq)]))
      index_vec<-sapply(1:sample_size, function(v) min((1:length(Q.partition)-2)[survt1[v]<=Q.partition]))
      #(1:length(Q.partition))[survt1[1]<=Q.partition]) 
      
      replicates <- dat[[i]]$datares[[1]]$num  # num. repeated measurements per subject
      maxrep_tmb <- max(replicates)
      start<-c(0,cumsum(replicates))
      
      
      
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=weight,a1=0, start=start)
      parameters <- list(sigma=c(0,rep(1,nalpha)))  # sigma=c(0,rep(1,nparam_a-1))
      
      obj <- MakeADFun(data,parameters,DLL="la_general", method="L-BFGS-B") 
      #L-BFGS-B always works for lA when it doesnt work without it, LB works fine for weights
      obj$hessian <- TRUE
      opt<-suppressWarnings(do.call("optim",obj))
      parest_a[i,] <- opt$par #parameter estimates
      na_ind_a[i] <- is.na(opt$value)
      
      if(is.na(opt$value)){
        next
      }
      
      rep<-sdreport(obj)
      repsum <- summary(rep, "report") #estimated a's
      a_est<-repsum[1:sample_size]
      la_gr <- repsum[(sample_size+1):(sample_size*nparam_a),] # excluding a's, only the gr fn
      
      
      #model B individual laplace
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=weight,a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      
      parameters <- list(sigma=c(0,1,1,rep(1,nbeta))) #make this modular  use eg nparam_b and rep
      
      obj <- MakeADFun(data,parameters,DLL="lb_general")
      obj$hessian <- TRUE
      opt<-suppressWarnings(do.call("optim",obj))
      parest_b[i,] <- opt$par
      na_ind_b[i] <- is.na(opt$value)
      
      if(is.na(opt$value)){
        next
      }
      
      rep<-sdreport(obj)
      repsum <- summary(rep, "report") #estimated b's
      b_est<-repsum[1:sample_size]
      lb_gr <- repsum[(sample_size+1):(sample_size*nparam_b),]
      #need some kind of indicator to  know whether the estimation is successful  (Eg value=?) in case grad is NaN
      
      #lC
      
      data<-list(a1=a_est,b1=b_est, Z0=zmat,N=sample_size,Qpartition=Q.partition,survt1=survt1,di1=di1,cumhaz_int=cumhaz_int,index_vec=index_vec,weights=weight)
      
      parameters <- list(theta=c(rep(1,nparam_c))) #list(gamma=rep(1,ncol(Z)),lambda2=1,lambda3=1,h1=1,h2=1,h3=1,h4=1,h5=1,h6=1,h7=1,h8=1,h9=1,h10=1)
      
      obj <- MakeADFun(data,parameters,DLL="lc_general")
      obj$hessian <- TRUE
      opt <- suppressWarnings(do.call("optim",obj))
      parest_c[i,] <- opt$par
      na_ind_c[i] <- is.na(opt$value)
      rep<-sdreport(obj)
      repsum <- summary(rep, "report")
      lc_gr <- repsum[(sample_size+1):(sample_size*nparam_c),]
      
      gr_all <- rbind(la_gr,lb_gr,lc_gr)
      
      
      #variance estimation using the G matrix and inverted Hessian
      SE[i,] <- var_est(data,gr_all, nparam_a,nparam_b,nparam_c,weights,sample_size,n.indiv,num.cluster)
      
      #end of weighted likelihood
    }
    
    if (type == "CAL"){
      
      xmat<-list()
      for(j in 1:numbX){
        assign(paste0("colname_X", j), (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste0('X.X', j, "_",1:maxrep)])
        assign(paste0("X",j,"row"),as.matrix(dat[[i]]$datares[[1]][,get(paste0("colname_X",j)) , with = FALSE]) )
        xmat[[j]] <- as.vector(t(get(paste0("X",j,"row"))))
      }
      
      xmat<-do.call("cbind",xmat)
      X_sample<-cbind(1,xmat[!rowSums(!is.finite(xmat)),])
      
      ncoldfZ <-ncol(dat[[i]]$datares[[1]])
      colname_costZ <- (1:ncoldfZ)[colnames(dat[[i]]$datares[[1]]) %in% paste0('Z', 1:numbX)]
      zmat <- as.matrix(dat[[i]]$datares[[1]][,colname_costZ, with =FALSE])
      
      # cost
      
      ncoldf <-ncol(dat[[i]]$datares[[1]])
      colname_cost <- (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste('cost2', 1:maxrep, sep='_')]
      
      dfcost<-as.matrix(dat[[i]]$datares[[1]][,colname_cost , with = FALSE])
      survt1 <-   dat[[i]]$datares[[1]]$survt
      di1 <-  dat[[i]]$datares[[1]]$di
      costindexsum<- dat[[i]]$datares[[1]]$costindexsum
      costsum<- dat[[i]]$datares[[1]]$costsum
      weight<- dat[[i]]$datares[[1]]$weight
      
      qq <-10  #10 intervals
      Q.partition <- c(0,quantile(survt1,1:qq/qq)) # length 11, need to be calculated for entire sample each iteration
      cumhaz_int <- c(Q.partition[2],(Q.partition[3:(qq+1)]-Q.partition[(2):(qq)]))
      index_vec<-sapply(1:sample_size, function(v) min((1:length(Q.partition)-2)[survt1[v]<=Q.partition]))
      #(1:length(Q.partition))[survt1[1]<=Q.partition]) 
      
      replicates <- dat[[i]]$datares[[1]]$num  # num. repeated measurements per subject
      maxrep_tmb <- max(replicates)
      start<-c(0,cumsum(replicates))
      
      
      
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=rep(1,sample_size),a1=0, start=start)
      parameters <- list(sigma=c(0,rep(1,nalpha)))  # sigma=c(0,rep(1,nparam_a-1))
      
      obj <- MakeADFun(data,parameters,DLL="la_general")
      obj$hessian <- TRUE
      opt<-suppressWarnings(do.call("optim",obj))
      #opt$par #parameter estimates
      rep<-sdreport(obj)
      repsum <- summary(rep, "report") #estimated a's
      a_est<-repsum[1:sample_size]
      #la_gr <- repsum[(sample_size+1):(sample_size*nparam_a)] # excluding a's, only the gr fn
      
      
      #model B individual laplace
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=rep(1,sample_size),a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      
      parameters <- list(sigma=c(0,1,1,rep(1,nbeta))) #make this modular  use eg nparam_b and rep
      
      obj <- MakeADFun(data,parameters,DLL="lb_general")
      obj$hessian <- TRUE
      opt<-suppressWarnings(do.call("optim",obj))
      #opt$par
      
      rep<-sdreport(obj)
      repsum <- summary(rep, "report") #estimated b's
      b_est<-repsum[1:sample_size]
      #lb_gr <- repsum[(sample_size+1):(sample_size*nparam_b)]
      #need some kind of indicator to  know whether the estimation is successful  (Eg value=?) in case grad is NaN
      
      #lC
      #using the randeffs estimated above, fit a imputed cox model
      lc_des <- svydesign(id=~1, weights=~weight2, strata=NULL, data=dat[[i]]$phase2[[1]])
      z_nam <- paste0("Z",missing_z) # Z1
      z_auxnam <- paste0("Z_aux",1:numbaux)
      fmla <- as.formula(paste(get("z_nam"),"~",paste(z_auxnam,collapse="+")))
      lc_glm <- svyglm(fmla, family=binomial, design=lc_des) #predicton model using phase2 data
      
      lc_fit <- summary(lc_glm) #use this to impute pmv
      ncoldf <-ncol(dat[[i]]$datares[[1]])
      colname_Z_aux <- (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste0('Z_aux', 1:numbaux)]
      Z_auxmat <- as.matrix(dat[[i]]$datares[[1]][,colname_Z_aux , with = FALSE])
      Z_auxmat <- cbind(1,Z_auxmat)
      xb <- Z_auxmat  %*% lc_fit$coefficients[,1]
      dat[[i]]$datares[[1]]$imp_Z1 <- ifelse(exp(xb)/(1+exp(xb))>0.5,1,0) 
      Z_nam <- c("imp_Z1",paste0("Z",2:numbX),"offset(a_est)", "offset(b_est)")
      fmla <- as.formula(paste("Surv(survt,di)~",paste(Z_nam,collapse="+")))  #use a's and b's as offsets? 
      lc_imp <- coxph(fmla, data=dat[[i]]$datares[[1]])
      
      inf_fun <- resid(lc_imp,  type='dfbeta')  #this is the IF for cox model
      
      colnames(inf_fun) <- paste0("if",1:ncol(inf_fun))
      des0 <- svydesign(id=~0, weights=~weight, data=dat[[i]]$datares[[1]])
      calformula <- make.formula(colnames(inf_fun))
      des0if <- update(des0, if1 = inf_fun[,1],
                       if2=inf_fun[,2], if3=inf_fun[,3])
      
      cal <- calibrate(des0if, formula=~if1+if2+if3, pop=c(N,0,0,0)) #put other ifs?
      cal_w_c <- weights(cal)
      
      #fit the 'actual' survival model using calibrated weights
      data<-list(a1=a_est,b1=b_est, Z0=zmat,N=sample_size,Qpartition=Q.partition,survt1=survt1,di1=di1,cumhaz_int=cumhaz_int,index_vec=index_vec,
                 weights=cal_w_c) #put c-weights
      
      parameters <- list(theta=c(rep(1,nparam_c))) #list(gamma=rep(1,ncol(Z)),lambda2=1,lambda3=1,h1=1,h2=1,h3=1,h4=1,h5=1,h6=1,h7=1,h8=1,h9=1,h10=1)
      
      obj <- MakeADFun(data,parameters,DLL="lc_general")
      obj$hessian <- TRUE
      opt <- suppressWarnings(do.call("optim",obj))
      parest_c[i,] <- opt$par
      rep<-sdreport(obj)
      repsum <- summary(rep, "report")
      lc_gr <- repsum[(sample_size*12+1):(sample_size*nparam_c)] #use this to calculate calibrated variance
      #exclude 10  hazs, cannot calibrate lambdas as randeffs are offsets
      
      #calibrated var est- fix!
      gr_gam1 <- lc_gr[1:sample_size];  #these are not weighed yet due to c++ code 
      gr_gam2 <- lc_gr[(sample_size+1):(2*sample_size)]
      gr_gam3 <- lc_gr[(2*sample_size+1):(3*sample_size)]
      mat_gr <- cbind(gr_gam1,gr_gam2,gr_gam3)
      mat_gr1 <- mat_gr[1:(sample_size/2),];
      mat_gr2 <- mat_gr[(sample_size/2+1):sample_size,]
      
      fit_gam1 <- lm(gr_gam1~inf_fun[,1]) #IF also need to be weighted??
      fit_gam2 <- lm(gr_gam2~inf_fun[,2])
      fit_gam3 <- lm(gr_gam3~inf_fun[,3])
      
      score_gam1 <- resid(fit_gam1, type="response")  #or fit_gam1$fitted.values ?
      score_gam2 <- resid(fit_gam2, type="response")
      score_gam3 <- resid(fit_gam3, type="response")
      
      score_est <- cbind(score_gam1,score_gam2,score_gam3)
      score_est1 <- score_est[1:(sample_size/2),] #for cluster 1 sample
      score_est2 <- score_est[(sample_size/2+1):sample_size,] #for cluster 2 sample
      
      #inf_fun_tst <- cbind(score_alp0, score_alp1, score_alp2, score_alp3)
      #inf_fun1 <- cbind(inf_fun_tst[ind,1],inf_fun_tst[ind,2],inf_fun_tst[ind,3],inf_fun_tst[ind,4])
      ngamma=3
      
      #variance estimation
      #n.indiv <- 10
      var_ind <- seq(1,sample_size/2+1,n.indiv) #n.indiv is  a parameter given when calling the function
      num.cluster <- 2 #num. cluster is fixed
      G <- matrix(0,nrow=ngamma, ncol=ngamma)
      for(i in 1:n.strata){ #n.strata = 20
        e_hl_1 <- colSums(score_est1[var_ind[i]:var_ind[i+1]-1,])  #1st cluster sum in str h
        e_hl_2 <- colSums(score_est2[var_ind[i]:var_ind[i+1]-1,]) #2nd cluster sum in str h
        
        e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
        G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                     + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
        
      }
      #invhess <- solve(opt$hessian) #make sure this is the weighted hessian
      #invhess_gam <- invhess[-c(1:12),-c(1:12)]
      #sqrt(diag(invhess_gam %*% G %*% invhess_gam))
      
      invhess <- solve(opt$hessian[-c(1:12),-c(1:12)])
      SE[i,] <- sqrt(diag(invhess %*% G %*% invhess))
    }
    
  } #end of for loop
  #return output.
  par <- data.frame(cbind(parest_a,parest_b,parest_c,na_ind_a,na_ind_b,na_ind_c))
  #for weighted  cbind parest_a.. etc inside the loop
  #so that this code works unilaterally for all types of analysis
  par_SE <- list(par=par, SE=SE)
  return(par_SE)
}