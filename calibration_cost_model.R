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
  
  #X1 <- rbinom(Ntotal, 1, 0.3)
  X1 <- rbinom(Ntotal, 1, 0.5)
  X2 <- rnorm(Ntotal, X1, 0.4) 
  X3 <- rnorm(Ntotal, X1, 0.3)
  #X1 <- rnorm(Ntotal, 2+0.5*X2+0.4*X3, 1) #X1 is now correlated with X2 and X3
  X <- cbind(1,X1,X2,X3)
  nparam_a <- ncol(X) + 1 # 4 alphas + logsigma_a
  nparam_b <- ncol(X) + 3 # 4 betas + logsigma_a + lambda1 + delta
  
  nalpha <- ncol(X)
  nbeta <- ncol(X)
  
  #for both binary and normal case?
  Z1 <- rbinom(N, 1, 0.4) #no repeated measure for Z.  Z1 is directly usable in cpp 
  Z2 <- rnorm(N, Z1, 0.4) 
  Z3 <- rnorm(N, Z1, 0.3)
  #Z1 <- rnorm(N, 1 + 0.5*Z2 + 0.6*Z3 ,1) #Z1 is now correlated with Z2 and Z3
  Z.1 <- Z1[idK] #longitudinal Z (needed for generating survival times)
  Z.2 <- Z2[idK]
  Z.3 <- Z3[idK]
  Z  <- cbind(Z.1, Z.2, Z.3)
  Z_one <- cbind(Z1, Z2, Z3)
  nparam_c <- ncol(Z) + 12 # 3 gammas + lambda2 + lambda3 + 10 h's
  
  ngamma <- ncol(Z)
  #partially missing variable available only for phase 2 sample
  #nth<-1
  #nth_missing <- cumsum(c(nth,replicates))[1:N]
  
  aux1<-rbinom(Ntotal,1,ifelse(X1== 1 | X2 ==1, 0.9,0.1)) # binomial case
  aux2<-rbinom(Ntotal,1,ifelse(X1+X2>=1, 0.95,0.05)) # binomial case
  aux3<- rnorm(Ntotal, X1, 0.2) #discrete case
  aux <- cbind(aux1,aux2,aux3)
  
  #aux vars for survival model
  Z_aux1 <- rbinom(N, 1, ifelse(Z1==1, 1,0)) #binomial
  Z_aux2 <- rbinom(N, 5, ifelse(Z1==1, 0.9,0.1)) #categ.
  #Z_aux3 <- rbinom(N, 1, ifelse(Z1==1, 0.95,0.05))
  Z_aux3 <- rnorm(N, Z1, 0.3) #continuous
  Z_aux <- cbind(Z_aux1, Z_aux2, Z_aux3) #matrix of aux vars for Z
  
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
sim_data_list<-sim_data(400,1,1)
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

#n.per.cluster <- 10 # n.per.cluster * 2 * n.strata = sample size   #this adjusts the weight
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
    #select 4 individuals based on selection probs
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

}



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
      
      # weighted design object using phase 2 data
      la_des <- svydesign(id=~1, weights=~weight2, strata=NULL,data=dat[[i]]$phase2[[1]])
      
      xnam <- paste0("X.X",missing_x, "_", nth)
      auxnam <- paste0("aux",1:numbaux)
      fmla <- as.formula(paste(get("xnam"),"~",paste(auxnam,collapse="+")))
      
      #weighted logistic model using phase 2 data
      la_glm <- svyglm(fmla, family=binomial, design=la_des)
      la_fit <- summary(la_glm)
      
      ncoldf <-ncol(data)
      colname_aux <- (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste0('aux', 1:numbaux, "_", nth)]   #extracting auxiliary variables used to predict nth measurement of some X variable
      auxmat<-as.matrix(dat[[i]]$datares[[1]][,colname_aux , with = FALSE])
      
      #matrix of auxiliary variables
      auxmat1 <- cbind(1,auxmat)
      xb <- auxmat1  %*% la_fit$coefficients[,1]
      data$imp_pmv <- ifelse(exp(xb)/(1+exp(xb))>0.5,1,0) #imputed vals for X. what if = 0.5?
      
      xmat<-list()
      maxrep<-max(dat[[i]]$datares[[1]]$num)
      for(i in 1:numbX){ 
        assign(paste0("colname_X", i), (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste0('X.X', i, "_",1:maxrep)])
        assign(paste0("X",i,"row"),as.matrix(dat[[i]]$datares[[1]][,get(paste0("colname_X",i)) , with = FALSE]) )
        xmat[[i]] <- as.vector(t(get(paste0("X",i,"row"))))
        
      }
      
      xmat<-do.call("cbind",xmat)
      X_sample<-cbind(1,xmat[!rowSums(!is.finite(xmat)),])
      
      ind<-cumsum(c(nth,dat[[i]]$datares[[1]]$num))[1:sample_size]
      X_sample[,(missing_x+1)][ind] <- imp_pmv  #replacing original pmv with its imputed vals
      
      ncoldf <-ncol(dat[[i]]$datares[[1]])
      colname_cost <- (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste('cost2', 1:maxrep, sep='_')]
      costmat<-as.matrix(dat[[i]]$datares[[1]][,colname_cost , with = FALSE])
      costmat1<-as.vector(t(costmat))
      cost2<-costmat1[is.finite(costmat1)] #not a matrix
      YNoZero_sample<-ifelse(cost2>0,1,0) #this is used for R glm
      
      idK_sample<-sort(rep(1:sample_size, dat[[i]]$datares[[1]]$num))
      
      xn<-paste0("X_sample[,",2:(numbX+1),"]") #excluding the intercept
      xvar_a<-c(xn,"(1|idK_sample)")
      la_fml<-as.formula(paste("YNoZero_sample~",paste(xvar_a,collapse="+")))
      
      #fit the interest model using sample
      
      la_imp<-glmer(la_fml, family=binomial)
      
      esta<-unlist(ranef(la_imp)) #est a's from la_imp  and then re-fitting the same model using these as offsets
      
      la_fml_offset <- as.formula(paste("YNoZero_sample~",paste(xn,collapse="+")))
      
      la_imp_offset<-glm(la_fml_offset + offset=esta[idK_sample], family=binomial) #fitting standard glm using estimated a's as offsets,  now this model  is used to get IFs
      modmat <- model.matrix(la_imp_offset)
      Ihat <- (t(modmat) %*% (modmat * fitted.values(la_imp_offset) * (1 - fitted.values(la_imp_offset))))/ nrow(modmat)
      score <- modmat * resid(la_imp_offset, type="response")
      inf_fun_a <- score %*% solve(Ihat)
      
      
      randeff_a <- unlist(getME(la_imp, "b"))
      
      #influence function for lA
      
      colnames(inf_fun_a) <- paste0("if",1:ncol(inf_fun_a))
      
      #u' for variance estimation
      #score.fit <- lm(score~inf_fun)
      #sc <- resid(score.fit, type="response")
      
      des0 <- svydesign(id=~0, weights=~weight[idK_sample], data=dat[[i]]$datares[[1]]) #define a design  with longitudinal weights
      calformula <- make.formula(colnames(inf_fun_a))
      
      des0if <- update(des0, if1 = inf_fun_a[,1],
                       if2=inf_fun_a[,2], if3=inf_fun_a[,3],
                       if4=inf_fun_a[,4])   #make this modular
      cal <- calibrate(des0if, formula=~if1+if2+if3+if4, pop=c(N,0,0,0,0))
      cal_w_a <- weights(cal) #longitudinal c-weights.  alternatively  get weights jointly with lB IFs
      
      idKdf <- data.frame(cal_w_a,idK_sample)
      cal_w_a <- aggregate(cal_w_a~idK_sample,idKdf,mean)$cal_w_a #subject-level weights
      
      #lB   no need to impute again as lB uses same  X predictors as lA
      
      # fit glmer with random effects and then fit glm  with offset random effects. delete where cost2 is 0 (the likelihood is 0 for zero cost, so deleting doesn't cause any issue)
      
      
      nozero_mat<-cbind(cost2,randeff_a[idK_sample],idK,X_sample)
      nozero_mat2<- testmat[-c(which(cost2==0)),]
      
      #redefining objects  w/o  zeroes
      cost2_nozero <- nozero_mat2[,1]
      randeff_a_nozero <- nozero_mat2[,2] #longitudinal without zeros
      idK_sample_nozero <- nozero_mat2[,3]
      X_sample_nozero <- nozero_mat2[,4:(ncol(nozero_mat2))]
      
      xn1<-paste0("X_sample_nozero[,",2:(numbX+1),"]")
      xvar_b<-c(xn1,"randeff_a_nozero", "(1|idK_sample_nozero)")
      lb_fml <- as.formula(paste("cost2_nozero~",paste(xvar_b,collapse="+")))
      
      #fit.b<-glmer(testmat2[,1]    ~ testmat2[,4]+  testmat2[,2] + (1|testmat2[,3]), data=df , Gamma) 
      
      lb_imp <- glmer(lb_fml, family=gamma)
      randeff_b <- unlist(getME(lb_imp, "b"))[idK] 
      
      lb_fml_offset <- as.formula(paste("cost2_nozero~",paste(xn1,collapse="+")))
      
      lb_imp_offset <- glm(lb_fml_offset + randeff_a_nozero + offset(randeff_b), family=gamma)
      
      #inf_fun <- model.matrix(lb_imp) * resid(lb_imp, type="response") #dfbeta for lb
      modmat <- model.matrix(lb_imp_offset)
      Ihat <- solve(vcov(lb_imp_offset)) #information matrix
      
      inf_fun_b <- modmat * resid(lb_imp_offset, type="response") %*% solve(Ihat)
      
      
      colnames(inf_fun_b) <- paste0("if",1:ncol(inf_fun_b))
      weight <- dat[[i]]$datares[[1]]$weight[-c(which(cost2==0))]
      #des0 <- svydesign(id=~0, weights=~weight, data=dat[[i]]$datares[[1]])
      des0 <- svydesign(id=~0, weights=~weight)
      calformula <- make.formula(colnames(inf_fun))
      
      
      des0if <- update(des0, if1 = inf_fun_b[,1],
                       if2=inf_fun_b[,2], if3=inf_fun_b[,3],
                       if4=inf_fun_b[,4])   #make this modular
      cal <- calibrate(des0if, formula=~if1+if2+if3+if4, pop=c(N,0,0,0,0))
      
      cal_w_b <- weights(cal)
      idKdf <- data.frame(cal_w_b,idK_sample)
      cal_w_b <- aggregate(cal_w_b~idK_sample,idKdf,mean)$cal_w_b #subject-level weights
      
      #now  calibrated weights are obtained for lA and lB
      # 1) code  MakeADFun  with original data and supply these weights
      xmat<-list()
      for(j in 1:numbX){
        assign(paste0("colname_X", j), (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste0('X.X', j, "_",1:maxrep)])
        assign(paste0("X",j,"row"),as.matrix(dat[[i]]$datares[[1]][,get(paste0("colname_X",j)) , with = FALSE]) )
        xmat[[j]] <- as.vector(t(get(paste0("X",j,"row"))))
      }
      
      xmat<-do.call("cbind",xmat)
      X_sample<-cbind(1,xmat[!rowSums(!is.finite(xmat)),])
      
      #ncoldfZ <-ncol(dat[[i]]$datares[[1]])
      #colname_costZ <- (1:ncoldfZ)[colnames(dat[[i]]$datares[[1]]) %in% paste0('Z', 1:numbX)] 
      #zmat <- as.matrix(dat[[i]]$datares[[1]][,colname_costZ, with =FALSE])
      
      # cost
      
      ncoldf <-ncol(dat[[i]]$datares[[1]])
      colname_cost <- (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste('cost2', 1:maxrep, sep='_')]
      
      dfcost<-as.matrix(dat[[i]]$datares[[1]][,colname_cost , with = FALSE])
      survt1 <-   dat[[i]]$datares[[1]]$survt
      di1 <-  dat[[i]]$datares[[1]]$di
      costindexsum<- dat[[i]]$datares[[1]]$costindexsum
      costsum<- dat[[i]]$datares[[1]]$costsum
      
      #weight<- dat[[i]]$datares[[1]]$weight
      
      
      qq <-10  #10 intervals
      Q.partition <- c(0,quantile(survt1,1:qq/qq)) # length 11, need to be calculated for entire sample each iteration
      cumhaz_int <- c(Q.partition[2],(Q.partition[3:(qq+1)]-Q.partition[(2):(qq)]))
      index_vec<-sapply(1:sample_size, function(v) min((1:length(Q.partition)-2)[survt1[v]<=Q.partition]))
      
      replicates <- dat[[i]]$datares[[1]]$num  # num. repeated measurements per subject
      maxrep_tmb <- max(replicates)
      start<-c(0,cumsum(replicates))
      
      
      
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=cal_w_a,a1=0, start=start)
      parameters <- list(sigma=c(0,rep(1,nalpha)))  # sigma=c(0,rep(1,nparam_a-1))
      
      obj <- MakeADFun(data,parameters,DLL="la_general", method="L-BFGS-B") 
      #L-BFGS-B always works for lA when it doesnt work without it, LB works fine for weights
      obj$hessian <- TRUE
      opt<-suppressWarnings(do.call("optim",obj))
      parest_a[i,] <- opt$par 
      na_ind_a[i] <- is.na(opt$value)
      
      if(is.na(opt$value)){
        next
      }
      
      rep<-sdreport(obj)
      repsum <- summary(rep, "report") #estimated a's
      a_est<-repsum[1:sample_size] # estimated random effects based on calibrated weights
      la_gr <- repsum[(sample_size+1):(sample_size*nparam_a),] # excluding a's, only the gr fn
      la_gr <- la_gr * cal_w_a
      
      #alpha gradients
      
      #gr_alp1 <- la_gr[1:sample_size]
      #gr_alp2 <- la_gr[(sample_size+1):(2*sample_size)]
      #gr_alp3 <- la_gr[(2*sample_size+1):(3*sample_size)]
      
      for(i in 1:nalpha){
        assign(paste0("gr_alp",i),la_gr[(sample_size*(i-1)+1):(i*sample_size)])
      }
      
      mat_gr <- cbind(gr_alp1,gr_alp2,gr_alp3,gr_alp4)
      mat_gr1 <- mat_gr[1:(sample_size/2),];
      mat_gr2 <- mat_gr[(sample_size/2+1):sample_size,]
      
      #residuals of  u' - u
      score_alp1 <- resid(lm(gr_alp1~inf_fun_a[,1]),type="response") #IF also need to be weighted??
      score_alp2 <- resid(lm(gr_alp2~inf_fun_a[,2]),type="response")
      score_alp3 <- resid(lm(gr_alp3~inf_fun_a[,3]),type="response")
      score_alp4 <- resid(lm(gr_alp4~inf_fun_a[,4]),type="response")
      
      score_est <- cbind(score_alp1,score_alp2,score_alp3,score_alp4)
      score_est1 <- score_est[1:(sample_size/2),] #for cluster 1 sample
      score_est2 <- score_est[(sample_size/2+1):sample_size,] #for cluster 2 sample
      
      var_ind <- seq(1,sample_size/2+1,n.indiv) #n.indiv is the number of subj selected in each cluster and is a parameter given when calling the function
      num.cluster <- 2 # num.cluster is fixed
      G <- matrix(0,nrow=nalpha, ncol=nalpha)
      for(i in 1:n.strata){ #n.strata = 20
        e_hl_1 <- colSums(score_est1[var_ind[i]:var_ind[i+1]-1,])  #1st cluster sum in str h
        # score_est is a matrix of parameter scores where each row has scores for each subject. the scores for subjects in the same cluster are summed using colSums(), for each parameter. 
        e_hl_2 <- colSums(score_est2[var_ind[i]:var_ind[i+1]-1,]) #2nd cluster sum in str h
        
        e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
        G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                     + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
        
      }
      
      invhess <- solve(opt$hessian[-c(1),-c(1)]) #excluding logsigma_A
      SE_a[i,] <- sqrt(diag(invhess %*% G %*% invhess))
      
      #model B individual laplace
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=cal_w_b,a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      
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
      lb_gr <- lb_gr * cal_w_b
      
      for(i in 1:nbeta){
        assign(paste0("gr_beta",i),lb_gr[(sample_size*(i-1)+1):(i*sample_size)])
      }
      
      mat_gr <- cbind(gr_beta1,gr_beta2,gr_beta3,gr_beta4)
      mat_gr1 <- mat_gr[1:(sample_size/2),];
      mat_gr2 <- mat_gr[(sample_size/2+1):sample_size,]
      
      #residuals of  u' - u
      score_beta1 <- resid(lm(gr_beta1~inf_fun_a[,1]),type="response") #IF also need to be weighted??
      score_beta2 <- resid(lm(gr_beta2~inf_fun_a[,2]),type="response")
      score_beta3 <- resid(lm(gr_beta3~inf_fun_a[,3]),type="response")
      score_beta4 <- resid(lm(gr_beta4~inf_fun_a[,4]),type="response")
      
      score_est <- cbind(score_beta1,score_beta2,score_beta3,score_beta4)
      score_est1 <- score_est[1:(sample_size/2),] #for cluster 1 sample
      score_est2 <- score_est[(sample_size/2+1):sample_size,] #for cluster 2 sample
      
      var_ind <- seq(1,sample_size/2+1,n.indiv) #n.indiv is the number of subj selected in each cluster and is a parameter given when calling the function
      num.cluster <- 2 # num.cluster is fixed
      G <- matrix(0,nrow=nbeta, ncol=nbeta)
      for(i in 1:n.strata){ #n.strata = 20
        e_hl_1 <- colSums(score_est1[var_ind[i]:var_ind[i+1]-1,])  #1st cluster sum in str h
        # score_est is a matrix of parameter scores where each row has scores for each subject. the scores for subjects in the same cluster are summed using colSums(), for each parameter. 
        e_hl_2 <- colSums(score_est2[var_ind[i]:var_ind[i+1]-1,]) #2nd cluster sum in str h
        
        e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
        G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                     + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
        
      }
      
      invhess <- solve(opt$hessian[-c(1:3),-c(1:3)]) #excluding logsigma_b,delta and lambda1
      SE_b[i,] <- sqrt(diag(invhess %*% G %*% invhess))
      
    }
    
  } #end of for loop
  #return output.
  
  #because estimates are saved every iteration,  can calculate emp SE from the output
  par <- data.frame(cbind(parest_a,parest_b,parest_c,na_ind_a,na_ind_b,na_ind_c))
  #for weighted  cbind parest_a.. etc inside the loop
  #so that this code works unilaterally for all types of analysis
  par_SE <- list(par=par, SE=SE)
  return(par_SE)
}

# sampling  n times to generate dataset

#par_est <- list()
n.iter<-1000
sampledata<-lapply(1:n.iter, function(v) samplingfunc(v,df1,n.per.cluster=20))
