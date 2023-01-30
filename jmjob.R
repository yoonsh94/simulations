library(pacman)
pacman::p_load(tidyverse,magrittr,flexsurv,dplyr,tidyr,data.table,plyr,TMB,survey,lme4)

compile("la_general.cpp",framework="TMBad")
dyn.load(dynlib("la_general"))

compile("lb_general.cpp",framework="TMBad")
dyn.load(dynlib("lb_general"))

compile("lc_general.cpp")
dyn.load(dynlib("lc_general"))

sim_data <- function(N){
  set.seed(3105)
  # Number of subjects
  N <- 10000 
  
  # Unconditional true parameters (Generalised)
  alpha.true<- c(2.5,1,1,1,1) 
  beta.true  <- c(2,1,1,1,1)
  gamma.true <- c(1, 1, 1)
  lambda1.true <- 0.7 
  lambda2.true <- 0.5
  lambda3.true <- 0.5
  delta.true <- 0.5
  
  # Replicates 
  mean.repl <- 20 # or 20?
  replicates <- rpois(N,mean.repl)+3
  idK <- sort(rep(1:N, replicates))
  
  Ntotal <-   sum(replicates)
  
  X1 <- rbinom(Ntotal, 1, 0.5)
  X2 <- rnorm(Ntotal, X1, 0.4) 
  X3 <- rnorm(Ntotal, X1, 0.3) #X1 is now correlated with X2 and X3
  #X <- cbind(1,X1,X2,X3)
  
  #for both binary and normal case?
  Z1 <- rbinom(N, 1, 0.5) #no repeated measure for Z.  Z1 is directly usable in cpp 
  Z2 <- rnorm(N, Z1, 0.9) #change this to 3
  Z3 <- rnorm(N, Z1, 0.8) #change this to 2
  Z.1 <- Z1[idK]
  X4 <- Z.2 <- Z2[idK] #this is shared with cost models
  Z.3 <- Z3[idK]
  Z  <- cbind(Z.1, Z.2, Z.3)
  Z_one <- cbind(Z1, Z2, Z3)
  
  
  X <- cbind(1,X1,X2,X3,X4) #subject-lv var added 
  nalpha <- ncol(X)
  nbeta <- ncol(X)
  
  nparam_a <- ncol(X) + 1 # 5 alphas + logsigma_a
  nparam_b <- ncol(X) + 3 # 5 betas + logsigma_a + lambda1 + delta
  
  nparam_c <- ncol(Z) + 12 # 3 gammas + lambda2 + lambda3 + 10 h's
  
  ngamma <- ncol(Z)
  #partially missing variable available only for phase 2 sample
  #nth<-1, nth_missing <- cumsum(c(nth,replicates))[1:N]
  
  aux1<-rbinom(Ntotal,1,ifelse(X1== 1 | X2 ==1, 0.9,0.1)) # binomial case
  aux2<-rbinom(Ntotal,1,ifelse(X1+X2>=1, 0.95,0.05)) # binomial case
  aux3<- rnorm(Ntotal, X1, 0.2) #discrete case
  aux <- cbind(aux1,aux2,aux3)
  
  #aux vars for survival model
  Z_aux1 <- rbinom(N, 1, ifelse(Z1==1, 1,0)) #binomial
  Z_aux2 <- rbinom(N, 5, ifelse(Z1==1, 0.9,0.1)) #categ.
  Z_aux3 <- rnorm(N, Z1, 0.3) #continuous
  Z_aux <- cbind(Z_aux1, Z_aux2, Z_aux3) 
  
  numbaux <- 3 #number of aux vars
  
  
  df <- data.frame(Identifier=idK,X=X[,-c(1)]) 
  df$idK <- idK
  df$num <- unlist(sapply(1:N, function(v) 1:replicates[v]))
  
  #auxiliary
  df <- cbind(df,aux)
  
  #dfS <- df[!duplicated(df$idK, fromLast = TRUE),]
  
  #replicates <- dfS$num #replicates after removing obs after event
  
  
  maxrep<-max(replicates)
  numbX<-ncol(X)-1 #excluding the first column of 1's
  
  
  setDT(df)
  
  #dcast(df, Identifier~num, value.var = c("cost2","X1","X2","X3"))
  df1 <- cbind(dcast(df, Identifier~num, value.var = c(paste0("X.X",1:numbX),paste0("aux",1:numbaux)))) #complete wide-form dataset with nrow=N
  
  df1$num <- replicates
  start<-cumsum(c(1,replicates))
  end<-cumsum(replicates)
  #modular now
  df1$X <- rowMeans(sapply(1:numbX, function(i) sapply(1:N, function(v)
    mean(get(paste0("X",i))[start[v]:end[v]])
  )))
  
  
  ncoldf <-ncol(df1)
  
  df1 <- cbind(df1,Z_one,Z_aux) #Zs and aux_Zs (subject-level) integrated to the data. 
  
  returnList <- list("df1"=df1, "maxrep"=maxrep,"nparam_a"=nparam_a, "nparam_b"=nparam_b, 
                     "nparam_c"=nparam_c,"nalpha"=nalpha,"nbeta"=nbeta,"ngamma"=ngamma,
                     "numbaux"=numbaux,"numbX"=numbX,"ncoldf"=ncoldf)
  
  
  return(returnList)
}

samplingfunc <- function(iter,data,n.per.cluster, n.per.cluster.p2){
  
  set.seed(iter) #this ensures data is same  for all simulation scenarios
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
    cluster.df1.indiv$prob <- (2/n.cluster) * (n.per.cluster*cluster.df1.indiv$size)/ sum(cluster.df1$size)
    cluster.df1.indiv$weight <- 1/cluster.df1.indiv$prob
    
    indiv.selected2 <- sample(c(1:(N/total.cluster)),size=n.per.cluster,prob=c(cluster.df2$prob),replace=FALSE)
    cluster.df2.indiv <- cluster.df2[indiv.selected2,]
    cluster.df2.indiv$prob <- (2/n.cluster) * (n.per.cluster*cluster.df2.indiv$size)/ sum(cluster.df2$size)
    cluster.df2.indiv$weight <- 1/cluster.df2.indiv$prob
    
    # Xrow subset using identifier cbind(data[ident,],Xlist[[1]][ident,])
    
    datalist1[[i]] <- cluster.df1.indiv
    datalist2[[i]] <- cluster.df2.indiv
    #sample further from cluster.df1 and df.2 (subsample or phase 2 data)
    
    
    phase2.selected1 <- sample(c(1:n.per.cluster), size=n.per.cluster.p2, prob=rep(1/n.per.cluster, n.per.cluster),replace=FALSE)
    phase2.selected2 <- sample(c(1:n.per.cluster), size=n.per.cluster.p2, prob=rep(1/n.per.cluster, n.per.cluster),replace=FALSE)
    
    cluster.df1.indiv$p2ind <- 0
    cluster.df1.indiv$p2ind[phase2.selected1]  <- rep(1,n.per.cluster.p2)  
    
    cluster.df2.indiv$p2ind <- 0
    cluster.df2.indiv$p2ind[phase2.selected2]  <- rep(1,n.per.cluster.p2)  
    
    datalist1[[i]] <- cluster.df1.indiv
    datalist2[[i]] <- cluster.df2.indiv
    
    phase2.data1 <- cluster.df1.indiv[phase2.selected1,]
    phase2.data2 <- cluster.df2.indiv[phase2.selected2,]
    n1 <- n.per.cluster*2*n.strata
    n2 <- n.per.cluster.p2 *2*n.strata
    
    phase2.data1$weight2 <- phase2.data1$weight * (n1/n2) #this need to be fixed
    phase2.data2$weight2 <- phase2.data2$weight * (n1/n2) # something like * (n_1 / n_2)
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
  xmat<-list()
  for(j in 1:numbX){
    assign(paste0("colname_X", j), (1:ncoldf)[colnames(data.final) %in% paste0('X.X', j, "_",1:maxrep)])
    assign(paste0("X",j,"row"),as.matrix(data.final[,get(paste0("colname_X",j)) , with = FALSE]) )
    xmat[[j]] <- as.vector(t(get(paste0("X",j,"row"))))
  }
  
  xmat<-do.call("cbind",xmat)
  X_sample<-cbind(1,xmat[!rowSums(!is.finite(xmat)),])
  phase2<-list(phase2)
  
  list(datares=datares,phase2=phase2,X_sample=X_sample)
  #returnList <- list("phase1"=data.final,"phase2"=phase2)
  #return(returnList)
}


#variance estimation function for weighted analysis

var_est <- function(dat,opta,optb,optc,gr_a,gr_b,gr_c, nparam_a,nparam_b,nparam_c,weights,sample_size,n.indiv,num.cluster,i){
  param_list_a<-list()
  param_list_b<-list()
  param_list_c<-list()
  for(a in 1:nparam_a){
    
    param_list_a[[a]] <- gr_a[(sample_size*(a-1)+1):(sample_size*a)]
  }
  mat_gr <- data.frame(do.call(cbind,param_list_a),dat[[i]]$datares[[1]]$cluster,dat[[i]]$datares[[1]]$strata, dat[[i]]$datares[[1]]$weight)
  mat_gr1 <- mat_gr[1:(sample_size/2),]
  mat_gr2 <- mat_gr[(sample_size/2+1):sample_size,]
  
  #variance estimation
  var_ind <- seq(1,sample_size+1,n.indiv)
  G <- matrix(0,nrow=nparam_a, ncol=nparam_a)
  for(astr in 1:n.strata){ #n.strata = 100
    e_hl_1 <- colSums(mat_gr1[,1:nparam_a][var_ind[astr]:var_ind[astr+1]-1,]) #1st cluster sum in str h
    e_hl_2 <- colSums(mat_gr2[,1:nparam_a][var_ind[astr]:var_ind[astr+1]-1,]) #2nd cluster sum in str h
    
    e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
    G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                 + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
    
  }
  
  invhessa <- solve(opta$hessian) #this is weighted Hessian
  se_a <- sqrt(diag(invhessa %*% G %*% invhessa)) 
  
  for(b in 1:nparam_b){   
    #logsigma_b, delta, 4 betas
    param_list_b[[b]] <- gr_b[(sample_size*(b-1)+1):(sample_size*b)]
  }
  mat_gr <- data.frame(do.call(cbind,param_list_b),dat[[i]]$datares[[1]]$cluster,dat[[i]]$datares[[1]]$strata, dat[[i]]$datares[[1]]$weight)
  mat_gr1 <- mat_gr[1:(sample_size/2),]
  mat_gr2 <- mat_gr[(sample_size/2+1):sample_size,]
  
  #variance estimation
  var_ind <- seq(1,sample_size+1,n.indiv)
  G <- matrix(0,nrow = nparam_b, ncol = nparam_b)
  for(bstr in 1:n.strata){ #n.strata = 100
    e_hl_1 <- colSums(mat_gr1[,1:(nparam_b)][var_ind[bstr]:var_ind[bstr+1]-1,]) #1st cluster sum in str h
    e_hl_2 <- colSums(mat_gr2[,1:(nparam_b)][var_ind[bstr]:var_ind[bstr+1]-1,]) #2nd cluster sum in str h
    
    e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
    G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                 + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
    
  }
  invhessb <- solve(optb$hessian) #this is weighted Hessian
  se_b <- sqrt(diag(invhessb %*% G %*% invhessb)) 
  
  for(c in 1:(nparam_c-10)){   
    #15 param in total minus 10 hazs
    param_list_c[[c]] <- gr_c[(sample_size*(c-1)+1):(sample_size*c)]
  }
  mat_gr <- data.frame(do.call(cbind,param_list_c),dat[[i]]$datares[[1]]$cluster,dat[[i]]$datares[[1]]$strata, dat[[i]]$datares[[1]]$weight)
  mat_gr1 <- mat_gr[1:(sample_size/2),]
  mat_gr2 <- mat_gr[(sample_size/2+1):sample_size,]
  
  #variance estimation
  var_ind <- seq(1,sample_size+1,n.indiv)
  G <- matrix(0,nrow = nparam_c-10, ncol = nparam_c-10)
  for(cstr in 1:n.strata){ #n.strata = 100
    e_hl_1 <- colSums(mat_gr1[,1:(nparam_c-10)][var_ind[cstr]:var_ind[cstr+1]-1,]) #1st cluster sum in str h
    e_hl_2 <- colSums(mat_gr2[,1:(nparam_c-10)][var_ind[cstr]:var_ind[cstr+1]-1,]) #2nd cluster sum in str h
    
    e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
    G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                 + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
    
  }
  invhessc <- solve(optc$hessian[-c(1:10),-c(1:10)]) #this is weighted Hessian
  se_c <- sqrt(diag(invhessc %*% G %*% invhessc)) 
  
  se <- c(se_a,se_b,se_c) # nalpha + (nparam_b-1) + (nparam-10)
  return(se)
}


adfun <- function(dat,n.iter,sample_size, nparam_a,nparam_b,nparam_c,nalpha,nbeta,ngamma, n.indiv, num.cluster, numbX, n.strata, sigmaA.true, sigmaB.true, type) { 
  # data passed to this function is a list of data frames (1000)
  parest_a <- matrix(NA,nrow=n.iter,ncol=nparam_a)
  parest_b <- matrix(NA,nrow=n.iter,ncol=nparam_b)
  parest_c <- matrix(NA,nrow=n.iter,ncol=nparam_c)
  na_ind_a <- c()
  na_ind_a_calS <- c()
  na_ind_b <- c()
  na_ind_b_calS <- c()
  na_ind_c <- c()
  
  SE_a <- matrix(NA,nrow=n.iter,ncol=nparam_a) # exc. a's
  SE_b <- matrix(NA,nrow=n.iter,ncol=nparam_b) # exc. b's
  SE_c <- matrix(NA,nrow=n.iter,ncol=nparam_c-10) # exc. hazs
  SE <- matrix(NA,nrow=n.iter,ncol=nparam_a+nparam_b+nparam_c-10)
  SE_cals <-matrix(NA,nrow=n.iter,ncol = nparam_a +nparam_b + nparam_c-10)
  SE_calc <- matrix(NA,nrow=n.iter,ncol = nparam_a +nparam_b + nparam_c-10)
  print(sigmaA.true)
  for(i in 1:n.iter){
    #set.seed(n.iter)
    ncoldf <-ncol(dat[[i]]$datares[[1]]) 
    maxrep <- max(dat[[i]]$datares[[1]]$num)
    # X matrix
    
    X_sample <- dat[[i]]$X_sample
    
    # Z matrix : Note that Z is a fixed covariate and this is transforming it to a matrix
    replicates <- dat[[i]]$datares[[1]]$num
    idK <- sort(rep(1:sample_size, dat[[i]]$datares[[1]]$num))
    
    ncoldfZ <-ncol(dat[[i]]$datares[[1]])
    colname_Z <- (1:ncoldfZ)[colnames(dat[[i]]$datares[[1]]) %in% paste0('Z', 1:numbX)] 
    #to get Z matrix. '
    Z_sample <- as.matrix(dat[[i]]$datares[[1]][,colname_Z, with = FALSE]) #subject level
    Z.1 <- Z_sample[,1][idK] # In the thesis, this is the target var for cox model
    Z.2 <- Z_sample[,2][idK] # aux var  (with extra aux var)
    Z.3 <- Z_sample[,3][idK] # aux var2 (with extra aux var)
    Z  <- cbind(Z.1, Z.2, Z.3)
    
    
    #random effect
    
    a<-rnorm(sample_size,0,sigmaA.true) #make this a param
    b<-rnorm(sample_size,0,sigmaB.true)
    idK <- sort(rep(1:sample_size, dat[[i]]$datares[[1]]$num))
    a_All <- a[idK]  #repeating baseline random effects
    b_All <- b[idK]
    ntotal <-   sum(dat[[i]]$datares[[1]]$num)
    
    # Conditional true parameters   
    alpha.true<- c(2.5, 1, 1, 1, 1) 
    beta.true  <- c(2,1,1,1,1)
    gamma.true <- c(1,1,1)
    
    #nalpha <- length(alpha.true)
    #nbeta <- length(beta.true)
    #ngamma <- length(gamma.true)
    lambda1.true <- 0.7 
    lambda2.true <- 0.5
    lambda3.true <- 0.5
    delta.true <- 0.5
    #need X and Z  matrix here
    eta.true <- X_sample %*% alpha.true + a_All
    mu.true <- exp(X_sample %*% beta.true + lambda1.true*a_All + b_All)
    sigmaiij2.true <- exp(delta.true)
    
    # Prepare a longitudinal data frame with identifiers and covariates that will eventually be transformed to wide-form with new survt, cost etc 
    df <- data.frame(Identifier=idK,X_sample=X_sample[,-c(1)]) #longitudinal df
    df$YNoZero  <-  rbinom(ntotal, 1, 1/(1+exp(-eta.true)))
    
    df$cost <- rgamma(ntotal,shape=sigmaiij2.true, scale=mu.true/sigmaiij2.true)
    df <- df %>% mutate(cost2 =  cost*YNoZero ) # df$cost2
    
    #survival part
    
    l0 <- 2  # hazard is exponential 
    U2     <-   runif(sample_size, 0,1)[idK]
    idKu  <- unique(idK)
    survt0 <-survt <- (-log(U2)/(l0*exp(Z %*% gamma.true + lambda2.true*a_All + lambda3.true*b_All)))
    
    # censoring times:
    survt.censor <- rexp(sample_size, rate= 2)[idK]
    # censoring:
    di    <- as.numeric(survt0 <= survt.censor);
    survt <- pmin(survt0, survt.censor)  #survt is the minimum of survival time and censoring time
    df$survt <- survt
    df$di <- di
    df$idK <- idK
    df$num <- unlist(sapply(1:sample_size, function(v) 1:replicates[v]))
    
    #this removes rows after event occurrence
    dfS <- df[!duplicated(df$idK, fromLast = TRUE),] 
    
    survt1<-dfS$survt 
    di1 <- dfS$di 
    replicates <- dfS$num #replicates after removing obs after event
    
    qq <-10  #10 intervals
    Q.partition <- c(0,quantile(survt1,1:qq/qq)) # length 11, need to be calculated for entire sample each iteration
    cumhaz_int <- c(Q.partition[2],(Q.partition[3:(qq+1)]-Q.partition[(2):(qq)]))
    index_vec<-sapply(1:sample_size, function(v) min((1:length(Q.partition)-2)[survt1[v]<=Q.partition]))
    
    
    maxrep<-max(replicates) #maximum number of replicates needed for obtaining wide-form data
    numbX<-ncol(X)-1 #excluding the first column of 1's (intercept)
    
    setDT(df)
    
    # wide-form data (YNoZero removed. this can be obtained using dfcost and ifelse)
    #df1 <- cbind(dcast(df, Identifier~num, value.var = c("cost2")),dfS[,c("survt","di","num")]) 
    df1 <- cbind(dcast(df, Identifier~num, value.var = c("cost2"))) 
    
    ncoldf <-ncol(df1)
    colnames(df1) <- c("Identifier", paste0("cost2_",1:maxrep))
    nX1col <- (1:ncoldf)[colnames(df1) %in% paste('cost2', 1:maxrep, sep='_')]
    dfcost<-as.matrix(df1[,nX1col , with = FALSE]) # cost matrix generated each iteration
    costindsum<-list()
    costsum<-list()
    
    for(n in 1:sample_size){ #do not use 'i' here
      costindsum[[n]] <- sum(dfcost[n,][1:replicates[n]] > 0)
      costsum[[n]] <- sum(dfcost[n,][1:replicates[n]] )
    }
    
    costindexsum <- unlist(costindsum) 
    costsum <- unlist(costsum)
    
    if (type == "UW"){
      
      weight<-rep(1,sample_size) # unweighted 
      
      
      
      replicates <- dat[[i]]$datares[[1]]$num  # num. repeated measurements each subj
      #maxrep_tmb <- max(replicates)
      start<-c(0,cumsum(replicates))
      
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=weight,a1=0, start=start)
      parameters <- list(sigma=c(0,rep(1,nalpha)))  # sigma=c(0,rep(1,nparam_a-1))
      
      obj <- MakeADFun(data,parameters,DLL="la_general", method="L-BFGS-B",silent=T) 
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
      SE_a[i,] <- sqrt(diag(solve(opt$hessian)))  # check whether logsigma need to be exponentiated 
      
      #model B individual laplace
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=weight,a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      
      parameters <- list(sigma=c(0,1,1,rep(1,nbeta))) 
      
      obj <- MakeADFun(data,parameters,DLL="lb_general",method="L-BFGS-B",silent=T)
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
      
      
      SE_b[i,] <- sqrt(diag(solve(opt$hessian)))  # check whether logsigma need to be exponentiated 
      
      #lC
      
      data<-list(a1=a_est,b1=b_est, Z0=Z_sample,N=sample_size,Qpartition=Q.partition,survt1=survt1,di1=di1,cumhaz_int=cumhaz_int,index_vec=index_vec,weights=weight)
      
      parameters <- list(theta=c(rep(1,nparam_c)))     
      obj <- MakeADFun(data,parameters,DLL="lc_general",silent=T)
      obj$hessian <- TRUE
      opt <- suppressWarnings(do.call("optim",obj))
      parest_c[i,] <- opt$par
      na_ind_c[i] <- is.na(opt$value)
      rep<-sdreport(obj)
      repsum <- summary(rep, "report")
      
      SE_c[i,] <- sqrt(diag(solve(opt$hessian[-c(1:10),-c(1:10)])))  # check whether logsigma need to be exponentiated 
      SE[i,] <- c(SE_a[i,],SE_b[i,],SE_c[i,]) 
      
    }
    
    if (type == "W"){
      weight<- dat[[i]]$datares[[1]]$weight
      replicates <- dat[[i]]$datares[[1]]$num  # num. repeated measurements per subject
      start<-c(0,cumsum(replicates))
      
      
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=weight,a1=0, start=start)
      parameters <- list(sigma=c(0,rep(1,nalpha)))  # sigma=c(0,rep(1,nparam_a-1))
      
      obja <- MakeADFun(data,parameters,DLL="la_general", method="L-BFGS-B",silent=T) 
      #L-BFGS-B always works for lA when it doesnt work without it, LB works fine for weights
      obja$hessian <- TRUE
      opta<-suppressWarnings(do.call("optim",obja))
      parest_a[i,] <- opta$par 
      na_ind_a[i] <- is.na(opta$value)
      
      if(is.na(opta$value)){
        next
      }
      
      rep<-sdreport(obja)
      repsum <- summary(rep, "report") #estimated a's
      a_est<-repsum[1:sample_size]
      la_gr <- repsum[(sample_size+1):(sample_size*(nparam_a+1))] * weight # excluding a's, only the gr fn
      
      
      #model B individual laplace
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=weight,a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      
      parameters <- list(sigma=c(0,1,1,rep(1,nbeta))) 
      
      objb <- MakeADFun(data,parameters,DLL="lb_general", method="L-BFGS-B", silent=T)
      objb$hessian <- TRUE
      optb<-suppressWarnings(do.call("optim",objb))
      parest_b[i,] <- optb$par
      na_ind_b[i] <- is.na(optb$value)
      
      if(is.na(optb$value)){
        next
      }
      
      rep<-sdreport(objb)
      repsum <- summary(rep, "report") #estimated b's
      b_est<-repsum[1:sample_size]
      lb_gr <- repsum[(sample_size+1):(sample_size*(nparam_b+1))] * weight
      
      
      #lC
      
      data<-list(a1=a_est,b1=b_est, Z0=Z_sample,N=sample_size,Qpartition=Q.partition,survt1=survt1,di1=di1,cumhaz_int=cumhaz_int,index_vec=index_vec,weights=weight)
      
      parameters <- list(theta=c(rep(1,nparam_c))) 
      
      objc <- MakeADFun(data,parameters,DLL="lc_general",silent=T)
      objc$hessian <- TRUE
      optc <- suppressWarnings(do.call("optim",objc))
      parest_c[i,] <- optc$par
      na_ind_c[i] <- is.na(optc$value)
      rep<-sdreport(objc)
      repsum <- summary(rep, "report")
      lc_gr <- repsum[(sample_size*10+1):(sample_size*nparam_c)] * weight
      
      #variance estimation using the G matrix and inverted Hessian
      SE[i,] <- var_est(dat,opta,optb,optc,la_gr,lb_gr,lc_gr,
                        nparam_a,nparam_b,nparam_c,weights,sample_size,n.indiv,num.cluster,i)
      
      #end of weighted likelihood
    }
    
    if (type == "CAL"){
      
      #weight<- dat[[i]]$datares[[1]]$weight
      
      #replicates <- dat[[i]]$datares[[1]]$num  # num. repeated measurements per subject
      #start<-c(0,cumsum(replicates))
      
      ##### estimating a's #####
      #data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=rep(1,sample_size),a1=0, start=start) #not weighted  because these are for randeffs estimation for phase-1 imputed model, which should be unweighted
      #parameters <- list(sigma=c(0,rep(1,nalpha)))  
      #obj <- MakeADFun(data,parameters,DLL="la_general",method="L-BFGS-B",silent=T)
      #obj$hessian <- TRUE
      #opt<-suppressWarnings(do.call("optim",obj))
      #na_ind_a[i] <- is.na(opt$value)
      #rep<-sdreport(obj)
      #repsum <- summary(rep, "report") #estimated a's
      #a_est<-repsum[1:sample_size]
      
      
      ##### estimating b's #####
      #data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=rep(1,sample_size),a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      #parameters <- list(sigma=c(0,1,1,rep(1,nbeta))) 
      #obj <- MakeADFun(data,parameters,DLL="lb_general",method="L-BFGS-B",silent=T)
      #obj$hessian <- TRUE
      #opt<-suppressWarnings(do.call("optim",obj))
      #na_ind_b[i] <- is.na(opt$value)
      #rep<-sdreport(obj)
      #repsum <- summary(rep, "report") #estimated b's
      #b_est<-repsum[1:sample_size]
      
      
      weight<- dat[[i]]$datares[[1]]$weight
      
      replicates <- dat[[i]]$datares[[1]]$num  # num. repeated measurements per subject
      start<-c(0,cumsum(replicates))
      
      
      
      nX1col <- (1:ncoldf)[colnames(df1) %in% paste('cost2', 1:maxrep, sep='_')]
      
      ####
      
      colname_cost <- nX1col
      costmat<-as.matrix(df1[,colname_cost , with = FALSE])
      costmat1<-as.vector(t(costmat))
      cost2<-costmat1[is.finite(costmat1)] #not a matrix
      YNoZero_sample<-ifelse(cost2>0,1,0) #this is used for R glm
      
      idK_sample<-sort(rep(1:sample_size, dat[[i]]$datares[[1]]$num))
      
      xn<-paste0("X_sample[,",2:(numbX+1),"]") #excluding the intercept
      xvar_a<-c(xn,"(1|idK_sample)")
      la_fml<-as.formula(paste("YNoZero_sample~",paste(xvar_a,collapse="+")))
      
      #fit the interest model using sample
      
      la_fit<-glmer(la_fml, family=binomial)
      
      #esta<-unlist(ranef(la_fit)) #est a's from la_imp  and then re-fitting the same model using these as offsets
      
      a_est <- unlist(getME(la_fit, "b")) #subject-lev a
      
      #lB   no need to impute again as lB uses same  X predictors as lA
      
      nozero_mat<-cbind(cost2,a_est[idK_sample],idK,X_sample)
      nozero_mat2<- nozero_mat[-c(which(cost2==0)),]
      
      a_est <- a_est[,1]
      
      #redefining objects  w/o  zeroes
      cost2_nozero <- nozero_mat2[,1]
      randeff_a_nozero <- nozero_mat2[,2] #longitudinal without zeros
      idK_sample_nozero <- nozero_mat2[,3]
      X_sample_nozero <- nozero_mat2[,4:(ncol(nozero_mat2))]
      
      xn1<-paste0("X_sample_nozero[,",2:(numbX+1),"]")
      xvar_b<-c(xn1,"randeff_a_nozero", "(1|idK_sample_nozero)")
      lb_fml <- as.formula(paste("cost2_nozero~",paste(xvar_b,collapse="+")))
      
      
      lb_imp <- glmer(lb_fml, family=Gamma(link = "log"))
      b_est <- unlist(getME(lb_imp, "b"))[,1]
      
      
      
      ### using estimated randeffs, fit a imputed cox model ###
      
      # Phase-two prediction model
      lc_des <- svydesign(id=~1, weights=~weight2, strata=NULL, data=dat[[i]]$phase2[[1]])
      z_nam <- paste0("Z",1)
      #z_auxnam <- paste0("Z_aux",1:numbaux)
      z_auxnam <- paste0("Z_aux",1:3)
      z_nam2 <- paste0("Z",2:ncol(Z))
      fmla <- as.formula(paste(get("z_nam"),"~",paste(c(z_nam2, z_auxnam),collapse="+")))
      lc_glm <- svyglm(fmla, family=binomial, design=lc_des) 
      
      lc_fit <- summary(lc_glm) # use this to impute partially missing var
      ncoldf <-ncol(dat[[i]]$datares[[1]])
      colname_Z_aux <- (1:ncoldf)[colnames(dat[[i]]$datares[[1]]) %in% paste0('Z_aux', 1:numbaux)]
      Z_auxmat <- as.matrix(dat[[i]]$datares[[1]][,colname_Z_aux , with = FALSE])
      Z_auxmat <- cbind(1,Z_sample[,2],Z_sample[,3],Z_auxmat)
      xb <- Z_auxmat  %*% lc_fit$coefficients[,1]
      
      #incorporating imputed Z1 to matrix of Z so this can be used in the final model fitting
      dat[[i]]$datares[[1]]$imp_Z1 <- Z_sample[,1] <- ifelse(exp(xb)/(1+exp(xb))>0.5,1,0)
      # assigned to Z_sample so that  imputed value is used for final model fitting
      
      #using original values for phase 2 sample
      
      
      for(p in 1:sample_size){
        if (dat[[i]]$datares[[1]]$p2ind[p] ==1){
          dat[[i]]$datares[[1]]$imp_Z1[p] <- subset(dat[[i]]$phase2[[1]], dat[[i]]$phase2[[1]]$Identifier == dat[[i]]$datares[[1]]$Identifier[p])$Z1 
        }
      }
      
      
      dat[[i]]$datares[[1]]$survt <- survt1
      dat[[i]]$datares[[1]]$di <- di1
      
      Z_nam <- c("imp_Z1",paste0("Z",2:numbX),"offset(a_est)", "offset(b_est)")
      fmla <- as.formula(paste("Surv(survt,di)~",paste(Z_nam,collapse="+")))  #randeffs as offset
      lc_imp <- coxph(fmla, data=dat[[i]]$datares[[1]])
      
      inf_fun <- resid(lc_imp,  type='dfbeta')  #this is the IF for cox model
      
      colnames(inf_fun) <- paste0("if",1:ncol(inf_fun))
      des0 <- svydesign(id=~0, weights=~weight, data=dat[[i]]$datares[[1]])
      calformula <- make.formula(colnames(inf_fun))
      des0if <- update(des0, if1 = inf_fun[,1],
                       if2=inf_fun[,2], if3=inf_fun[,3])
      
      cal <- calibrate(des0if, formula=~if1+if2+if3, pop=c(N,0,0,0)) #put other ifs?
      cal_w_c <- weights(cal) #calibrated weights 
      
      
      #fitting the cost models using the weights calibrated to the survival influence functions
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=cal_w_c,a1=0, start=start)
      parameters <- list(sigma=c(0,rep(1,nalpha)))  # sigma=c(0,rep(1,nparam_a-1))
      
      obja <- MakeADFun(data,parameters,DLL="la_general",method="L-BFGS-B",silent=T)
      obja$hessian <- TRUE
      opta<-suppressWarnings(do.call("optim",obja))
      na_ind_a_calS[i] <- is.na(opta$value)
      parest_a[i,]<- opta$par #parameter estimates
      rep<-sdreport(obja)
      repsum <- summary(rep, "report") #estimated a's
      a_est<-repsum[1:sample_size]
      la_gr <- repsum[(sample_size+1):(sample_size*(nparam_a+1))]   
      
      #model B individual laplace
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=cal_w_c, a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      
      parameters <- list(sigma=c(0,1,1,rep(1,nbeta))) #make this modular  use eg nparam_b and rep
      
      objb <- MakeADFun(data,parameters,DLL="lb_general",method="L-BFGS-B",silent=T)
      objb$hessian <- TRUE
      optb<-suppressWarnings(do.call("optim",objb))
      na_ind_b_calS[i] <- is.na(optb$value)
      parest_b[i,] <-optb$par
      
      rep<-sdreport(objb)
      repsum <- summary(rep, "report") #estimated b's
      b_est<-repsum[1:sample_size]
      lb_gr <- repsum[(sample_size+1):(sample_size*(nparam_b+1))] 
      
      #fitting the survival model using calibrated weights
      data<-list(a1=a_est,b1=b_est, Z0=Z_sample,N=sample_size,Qpartition=Q.partition,survt1=survt1,di1=di1,cumhaz_int=cumhaz_int,index_vec=index_vec,
                 weights=cal_w_c) #put c-weights
      
      parameters <- list(theta=c(rep(1,nparam_c))) 
      obj <- MakeADFun(data,parameters,DLL="lc_general",silent=T)
      obj$hessian <- TRUE
      opt <- suppressWarnings(do.call("optim",obj))
      parest_c[i,] <- opt$par
      rep<-sdreport(obj)
      repsum <- summary(rep, "report")
      lc_gr <- repsum[(sample_size*10+1):(sample_size*nparam_c)] 
      
      #this is only if lambdas are needed whose SEs can be calc. as in Xu et al
      #lc_gr_lambda <- repsum[(sample_size*10+1):(sample_size*nparam_c)] * cal_w_c
      # SE_lambdas <- var_est_c(...)[1:2]
      
      #calibrated var est
      gr_logsigma_a <- la_gr[1 : sample_size]
      gr_alp0 <- la_gr[(sample_size*1+1) : (sample_size*2)] 
      gr_alp1 <- la_gr[(sample_size*2+1) : (sample_size*3)]
      gr_alp2 <- la_gr[(sample_size*3+1) : (sample_size*4)]
      gr_alp3 <- la_gr[(sample_size*4+1) : (sample_size*5)]
      gr_alp4 <- la_gr[(sample_size*5+1) : (sample_size*6)] #subject level 
      
      gr_logsigma_b <- lb_gr[1  : sample_size]
      gr_delta <- lb_gr[(sample_size*1+1)  : (sample_size*2)]
      gr_lambda1 <- lb_gr[(sample_size*2+1)  : (sample_size*3)]
      gr_beta0 <- lb_gr[(sample_size*3+1)  : (sample_size*4)]
      gr_beta1 <- lb_gr[(sample_size*4+1)  : (sample_size*5)]
      gr_beta2 <- lb_gr[(sample_size*5+1)  : (sample_size*6)]
      gr_beta3 <- lb_gr[(sample_size*6+1)  : (sample_size*7)]
      gr_beta4 <- lb_gr[(sample_size*7+1)  : (sample_size*8)]
      
      gr_lambda2 <- lc_gr[1 :sample_size]
      gr_lambda3 <- lc_gr[(sample_size*1+1) : (sample_size*2)]
      gr_gam1 <- lc_gr[(sample_size*2+1) :(sample_size*3)]   
      gr_gam2 <- lc_gr[(sample_size*3+1) :(sample_size*4)]
      gr_gam3 <- lc_gr[(sample_size*4+1) :(sample_size*5)]
      
      #mat_gr <- cbind(gr_alp1, gr_alp2, gr_alp3, gr_beta1, gr_beta2, gr_beta3, gr_gam1,gr_gam2,gr_gam3)
      #mat_gr1 <- mat_gr[1:(sample_size/2),];
      #mat_gr2 <- mat_gr[(sample_size/2+1):sample_size,]
      
      
      ##### lA #####
      score_logsigma_a <- cal_w_c * resid(lm(gr_logsigma_a ~ inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_alp0 <- cal_w_c * resid(lm(gr_alp0~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      score_alp1 <- cal_w_c * resid(lm(gr_alp1~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      score_alp2 <- cal_w_c * resid(lm(gr_alp2~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      score_alp3 <- cal_w_c * resid(lm(gr_alp3~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_alp4 <- cal_w_c * resid(lm(gr_alp4~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_est <- cbind(score_logsigma_a, score_alp0, score_alp1,score_alp2,score_alp3,score_alp4)
      score_est1 <- score_est[1:(sample_size/2),] #for cluster 1 sample
      score_est2 <- score_est[(sample_size/2+1):sample_size,] #for cluster 2 sample
      
      var_ind <- seq(1,sample_size/2+1,n.indiv) #n.indiv is  a parameter given when calling the function
      num.cluster <- 2 #num. cluster is fixed
      G <- matrix(0,nrow=nparam_a, ncol=nparam_a)
      for(s in 1:n.strata){ #n.strata = 20
        e_hl_1 <- colSums(score_est1[var_ind[s]:var_ind[s+1]-1,])  #1st cluster sum in str h
        e_hl_2 <- colSums(score_est2[var_ind[s]:var_ind[s+1]-1,]) #2nd cluster sum in str h
        
        e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
        G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                     + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
        
      }
      
      invhess_a <- solve(opta$hessian) #no sigma_a and alpha_0
      SE_a <- sqrt(diag(invhess_a %*% G %*% invhess_a))
      
      
      
      ##### lB #####
      score_logsigma_b <- cal_w_c * resid(lm(gr_logsigma_b~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_delta <- cal_w_c * resid(lm(gr_delta~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_lambda1 <- cal_w_c * resid(lm(gr_lambda1~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_beta0 <- cal_w_c * resid(lm(gr_beta0~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_beta1 <- cal_w_c * resid(lm(gr_beta1~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_beta2 <- cal_w_c * resid(lm(gr_beta2~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_beta3 <- cal_w_c * resid(lm(gr_beta3~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_beta4 <- cal_w_c * resid(lm(gr_beta4~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_est <- cbind(score_logsigma_b, score_delta, score_lambda1, score_beta0, score_beta1,score_beta2,score_beta3, score_beta4)
      
      score_est1 <- score_est[1:(sample_size/2),] #for cluster 1 sample
      score_est2 <- score_est[(sample_size/2+1):sample_size,] #for cluster 2 sample
      
      #variance estimation for lC
      var_ind <- seq(1,sample_size/2+1,n.indiv) #n.indiv is  a parameter given when calling the function
      num.cluster <- 2 #num. cluster is fixed
      G <- matrix(0,nrow=nparam_b, ncol=nparam_b)
      for(s in 1:n.strata){ #n.strata = 20
        e_hl_1 <- colSums(score_est1[var_ind[s]:var_ind[s+1]-1,])  #1st cluster sum in str h
        e_hl_2 <- colSums(score_est2[var_ind[s]:var_ind[s+1]-1,]) #2nd cluster sum in str h
        
        e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
        G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                     + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
        
      }
      
      invhess_b <- solve(optb$hessian) #no sigma_b, lambda1, delta and beta_0
      SE_b <- sqrt(diag(invhess_b %*% G %*% invhess_b))
      
      ##### lC #####
      score_lambda2 <- cal_w_c * resid(lm(gr_lambda2~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      score_lambda3 <- cal_w_c * resid(lm(gr_lambda3~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      score_gam1 <- cal_w_c * resid(lm(gr_gam1~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      score_gam2 <- cal_w_c * resid(lm(gr_gam2~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      score_gam3 <- cal_w_c * resid(lm(gr_gam3~inf_fun[,1]+inf_fun[,2]+inf_fun[,3]), type="response")
      
      score_est <- cbind(score_lambda2, score_lambda3, score_gam1, score_gam2,score_gam3)
      score_est1 <- score_est[1:(sample_size/2),] #for cluster 1 sample
      score_est2 <- score_est[(sample_size/2+1):sample_size,] #for cluster 2 sample
      
      
      #variance estimation for lC
      var_ind <- seq(1,sample_size/2+1,n.indiv) #n.indiv is  a parameter given when calling the function
      num.cluster <- 2 #num. cluster is fixed
      G <- matrix(0,nrow=ngamma+2, ncol=ngamma+2) #gamma has no intercept
      for(s in 1:n.strata){ #n.strata = 20
        e_hl_1 <- colSums(score_est1[var_ind[s]:var_ind[s+1]-1,])  #1st cluster sum in str h
        e_hl_2 <- colSums(score_est2[var_ind[s]:var_ind[s+1]-1,]) #2nd cluster sum in str h
        
        e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
        G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                     + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
        
      }
      
      invhess <- solve(opt$hessian[-c(1:10),-c(1:10)]) #no lambdas as randeffs are offsets
      SE_c <- sqrt(diag(invhess %*% G %*% invhess))
      
      SE_cals[i,] <- c(SE_a, SE_b, SE_c)
    }
    
    if (type == "CAL-C"){
      
      # weighted design object using phase 2 data
      weight <- dat[[i]]$datares[[1]]$weight  #to calibrate the weights
      la_des <- svydesign(id=~1, weights=~weight2, strata=NULL,data=dat[[i]]$phase2[[1]])
      
      #xnam <- paste0("X.X",missing_x, "_", nth)
      xnam <- paste0("X.X","_",1:maxrep)
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
      
      la_imp_offset<-glm(la_fml_offset + offset(esta[idK_sample]), family=binomial) #fitting standard glm using estimated a's as offsets,  now this model  is used to get IFs
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
      
      obj <- MakeADFun(data,parameters,DLL="la_general", method="L-BFGS-B",silent = T) 
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
      for(s in 1:n.strata){ #n.strata = 20
        e_hl_1 <- colSums(score_est1[var_ind[s]:var_ind[s+1]-1,])  #1st cluster sum in str h
        # score_est is a matrix of parameter scores where each row has scores for each subject. the scores for subjects in the same cluster are summed using colSums(), for each parameter. 
        e_hl_2 <- colSums(score_est2[var_ind[s]:var_ind[s+1]-1,]) #2nd cluster sum in str h
        
        e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
        G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                     + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
        
      }
      
      invhess <- solve(opt$hessian) #excluding logsigma_A
      SE_a[i,] <- sqrt(diag(invhess %*% G %*% invhess))
      
      #model B individual laplace
      data <- list(Xrow=X_sample, N=sample_size, replicates=replicates, dfcost=dfcost, weights=cal_w_b,a1=a_est, b1=0, costindexsum=costindexsum, costsum=costsum, start=start)
      
      parameters <- list(sigma=c(0,1,1,rep(1,nbeta))) #make this modular  use eg nparam_b and rep
      
      obj <- MakeADFun(data,parameters,DLL="lb_general",method="L-BFGS-B",silent=T)
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
      score_beta1 <- resid(lm(gr_beta1~inf_fun_a[,1]),type="response") 
      score_beta2 <- resid(lm(gr_beta2~inf_fun_a[,2]),type="response")
      score_beta3 <- resid(lm(gr_beta3~inf_fun_a[,3]),type="response")
      score_beta4 <- resid(lm(gr_beta4~inf_fun_a[,4]),type="response")
      
      score_est <- cbind(score_beta1,score_beta2,score_beta3,score_beta4)
      score_est1 <- score_est[1:(sample_size/2),] #for cluster 1 sample
      score_est2 <- score_est[(sample_size/2+1):sample_size,] #for cluster 2 sample
      
      var_ind <- seq(1,sample_size/2+1,n.indiv) 
      num.cluster <- 2 # num.cluster is fixed
      G <- matrix(0,nrow=nbeta, ncol=nbeta)
      for(s in 1:n.strata){ #n.strata = 20
        e_hl_1 <- colSums(score_est1[var_ind[s]:var_ind[s+1]-1,])  #1st cluster sum in str h
        # score_est is a matrix of parameter scores where each row has scores for each subject. the scores for subjects in the same cluster are summed using colSums(), for each parameter. 
        e_hl_2 <- colSums(score_est2[var_ind[s]:var_ind[s+1]-1,]) #2nd cluster sum in str h
        
        e_h <- (1/num.cluster) * (e_hl_1 + e_hl_2) #colsum returns a vector
        G <- G + (num.cluster / (num.cluster -1)) * ((as.matrix(e_hl_1 - e_h) %*% t(as.matrix(e_hl_1 - e_h)))
                                                     + (as.matrix(e_hl_2 - e_h) %*% t(as.matrix(e_hl_2 - e_h))))
        
      }
      
      invhess <- solve(opt$hessian) #excluding logsigma_b,delta and lambda1
      SE_b[i,] <- sqrt(diag(invhess %*% G %*% invhess))
      
    }
    if (i %% 10 == 0) cat("Iteration:" , i, "is done\n")
    
  } # end of for loop
  
  # return output.
  par <- data.frame(cbind(parest_a,parest_b,parest_c,na_ind_a,na_ind_b,na_ind_c,na_ind_a_calS,na_ind_b_calS))
  
  par_SE <- list(par=par, SE=SE, SE_cals = SE_cals)
  return(par_SE)
}



sigmaA<-rep(c(0,1,2), 9)
sigmaB<-rep(c(0,1,2), 9)
p2<-rep(c(5,5,5,10,10,10,15,15,15),3)
type<-rep(c("UW","W","CAL"),each=9)
jobid_mat<-cbind(sigmaA,sigmaB,p2,type)

sim_data_list<-sim_data(10000) #population data 

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
ncoldf<-sim_data_list$ncoldf

jobid = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

N<-10000
n.strata <- 20 #each stratum has 10 clusters and 500 obs
n.cluster <- 5 #10 #for each cluster within a stratum there are 50 obs.  A total of 200 clusters.

n.per.cluster <- 25 #10 # n.per.cluster * 2 * n.strata = sample size   #this adjusts the weight
sample_size <- n.per.cluster * 2 * n.strata

total.cluster <- n.strata * n.cluster
df1$strata<-rep(seq(1:n.strata), each=N/n.strata) 
df1$cluster<-rep(seq(1:n.cluster),N/n.cluster)

n.iter<-1000
sampledata<-lapply(1:n.iter, function(v) samplingfunc(v,df1,n.per.cluster=n.per.cluster, n.per.cluster.p2 = as.numeric(jobid_mat[jobid,3]))) 

sim_cal<-adfun(dat=sampledata,
               n.iter=n.iter,
               sample_size=sample_size, 
               nparam_a=nparam_a,
               nparam_b=nparam_b,
               nparam_c=nparam_c,
               nalpha=nalpha,nbeta=nbeta,ngamma=ngamma, 
               n.indiv=n.per.cluster, #number of subjects selected each cluster 
               num.cluster=2,  # fixed to 2. this is the number of clusters each stratum  from which samples are selected, not the number of clusters each stratum has
               numbX=numbX,
               n.strata=n.strata, sigmaA.true=as.numeric(jobid_mat[jobid,1]), sigmaB.true=as.numeric(jobid_mat[jobid,2]),
               type=jobid_mat[jobid,4])


saveRDS(sim_cal, paste("a",as.numeric(jobid_mat[jobid,1]),"_","b",
                       as.numeric(jobid_mat[jobid,2]),"_","p2","_",
                       as.numeric(jobid_mat[jobid,3]),"_",jobid_mat[jobid,4], ".rda", sep=""))