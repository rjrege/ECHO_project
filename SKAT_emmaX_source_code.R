
function (formula, data = NULL, K = NULL, Kin.File = NULL, ngrids = 100, 
          llim = -10, ulim = 10, esp = 1e-10, Is.GetEigenResult = FALSE) 
{
  obj1 <- model.frame(formula, na.action = na.omit, data)
  obj2 <- model.frame(formula, na.action = na.pass, data)
  n1 <- dim(obj2)[1]
  n <- dim(obj1)[1]
  id_include <- SKAT_Null_Model_Get_Includes(obj1, obj2)
  X <- model.matrix(formula, data = data)
  X <- Get_SKAT_Residuals.Get_X1(X)
  y <- model.response(obj1)
  q <- ncol(X)
  if (n1 - n > 0) {
    MSG <- sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!", 
                   n1 - n)
    warning(MSG, call. = FALSE)
  }
  stopifnot(nrow(X) == n)
  if (is.null(K) && is.null(Kin.File)) {
    stop("both K and Kin.File are NULL!")
  }
  else if (is.null(K)) {
    if (!file.exists(Kin.File)) {
      msg <- sprintf("File %s is not exist!", Kin.File)
      stop(msg)
    }
    cat("Read ", Kin.File, "\n")
    K = as.matrix(read.delim(Kin.File, header = FALSE, colClasses = "numeric", 
                             sep = "\t"))
    t <- nrow(K)
    cat("Read complete. ", Kin.File, " has ", t, " rows! \n")
  }
  if (!is.matrix(K)) {
    stop("K is not a matrix!")
  }
  if (n1 - n > 0) {
    K <- K[id_include, id_include]
  }
  t <- nrow(K)
  stopifnot(ncol(K) == t)
  eig.R <- emma.eigen.R.wo.Z(K, X)
  etas <- crossprod(eig.R$vectors, y)
  logdelta <- (0:ngrids)/ngrids * (ulim - llim) + llim
  m <- length(logdelta)
  delta <- exp(logdelta)
  Lambdas <- matrix(eig.R$values, n - q, m) + matrix(delta, 
                                                     n - q, m, byrow = TRUE)
  Etasq <- matrix(etas * etas, n - q, m)
  LL <- 0.5 * ((n - q) * (log((n - q)/(2 * pi)) - 1 - log(colSums(Etasq/Lambdas))) - 
                 colSums(log(Lambdas)))
  dLL <- 0.5 * delta * ((n - q) * colSums(Etasq/(Lambdas * 
                                                   Lambdas))/colSums(Etasq/Lambdas) - colSums(1/Lambdas))
  optlogdelta <- vector(length = 0)
  optLL <- vector(length = 0)
  if (dLL[1] < esp) {
    optlogdelta <- append(optlogdelta, llim)
    optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim, 
                                                   eig.R$values, etas))
  }
  if (dLL[m - 1] > 0 - esp) {
    optlogdelta <- append(optlogdelta, ulim)
    optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim, 
                                                   eig.R$values, etas))
  }
  for (i in 1:(m - 1)) {
    if ((dLL[i] * dLL[i + 1] < 0 - esp * esp) && (dLL[i] > 
                                                  0) && (dLL[i + 1] < 0)) {
      r <- uniroot(emma.delta.REML.dLL.wo.Z, lower = logdelta[i], 
                   upper = logdelta[i + 1], lambda = eig.R$values, 
                   etas = etas)
      optlogdelta <- append(optlogdelta, r$root)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root, 
                                                     eig.R$values, etas))
    }
  }
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  maxva <- sum(etas * etas/(eig.R$values + maxdelta))/(n - 
                                                         q)
  maxve <- maxva * maxdelta
  va = maxva
  ve = maxve
  V = va * K
  diag(V) = diag(V) + ve
  if (ve/va < 0.01) {
    eig.V <- eigen(V, symmetric = TRUE)
    idx.lambda <- which(abs(eig.V$values) > mean(abs(eig.V$values))/10000)
    V_inv <- eig.V$vectors[, idx.lambda] %*% (t(eig.V$vectors[, 
                                                              idx.lambda]) * (1/eig.V$values[idx.lambda]))
  }
  else {
    V_inv <- solve(V)
  }
  XVX = t(X) %*% (V_inv %*% X)
  XVX_inv <- solve(XVX)
  XV_inv = t(X) %*% V_inv
  P = V_inv - t(XV_inv) %*% (XVX_inv %*% XV_inv)
  res = P %*% y
  re <- list(LL = maxLL, va = va, ve = ve, P = P, res = res, 
             id_include = id_include)
  if (Is.GetEigenResult) {
    re$eig.R = eig.R
  }
  class(re) <- "SKAT_NULL_Model_EMMAX"
  return(re)
}

Get_SKAT_Residuals.Get_X1 = function(X1){
  
  qr1<-qr(X1)
  q1<-ncol(X1)
  if(qr1$rank < q1){
    
    X1.svd<-svd(X1)
    X1 = X1.svd$u	
  } 
  
  return(X1)
  
}


Get_SKAT_Residuals.linear = function(formula, data, n.Resampling, type.Resampling, id_include ){
  
  
  mod = lm(formula, data=data)
  X1<-model.matrix(formula,data=data)
  X1<-Get_SKAT_Residuals.Get_X1(X1)
  
  s2 = summary(mod)$sigma**2
  res = mod$resid
  n1<-length(res)
  res.out<-NULL
  
  if(n.Resampling > 0){
    
    if(type.Resampling=="permutation"){
      res.out<-res %x% t(rep(1,n.Resampling))
      res.out<-apply(res.out,2,sample)
    } else if(type.Resampling=="bootstrap"){
      res.out<-matrix(rnorm(n1*n.Resampling,mean=0,sd=sqrt(s2)),ncol=n.Resampling)
      X1_inv<-solve(t(X1) %*% X1)
      res.out<- res.out - (X1 %*% X1_inv) %*% (t(X1) %*% res.out)
    } else if(type.Resampling=="perturbation"){
      res.out<-matrix(rnorm(n1*n.Resampling,mean=0,sd=1),ncol=n.Resampling)
      res.out<-res.out * res
      stop("Error: Perturbation is no more provided!")
    } else {
      stop("Error: Wrong resampling method!")
    }
  }
  
  return(list(res=res, X1=X1,res.out=res.out,out_type="C", 
              n.Resampling=n.Resampling, type.Resampling=type.Resampling,
              id_include=id_include, s2=s2))
}


Get_SKAT_Residuals.logistic = function(formula, data, n.Resampling, type.Resampling,id_include){
  
  
  mod = lm(formula, data)
  X1<-model.matrix(formula,data=data)
  X1<-Get_SKAT_Residuals.Get_X1(X1)
  
  glmfit= glm(formula, data=data, family = "binomial")
  betas = glmfit$coef
  mu    = glmfit$fitted.values
  eta   = glmfit$linear.predictors
  n.case = sum(glmfit$y)
  
  pi_1 = mu*(1-mu)
  res = glmfit$y- mu
  n1<-length(res)
  res.out<-NULL
  
  if(n.Resampling > 0){
    if(type.Resampling=="bootstrap.fast"){
      
      res.out<-Get_Resampling_Bin(n.case, mu, n.Resampling)
      if(is.null(res.out)){
        type.Resampling="bootstrap"
      }
      res.out<-res.out - mu
    } 
    
    if(type.Resampling=="permutation"){
      res.out1<-res %x% t(rep(1,n.Resampling))
      res.out<-apply(res.out1,2,sample)
    } else if(type.Resampling=="bootstrap"){
      mu1<-mu/sum(mu)	# all prob
      res.out<-matrix(rep(0,n.Resampling*n1),ncol=n.Resampling)
      for(i in 1:n.Resampling){
        #id_case<-sample(1:n1,n.case,prob=mu1)
        #res.out[id_case,i]<-1
        #res.out[,i]<-rbinom(n1,1,mu)
        
        res.out1<-rbinom(n1,1,mu)
        res.out2<-rbinom(n1,1,mu)
        
        id_case1<-which(res.out1 ==1)
        id_case2<-which(res.out2 ==1)
        
        id_c1<-intersect(id_case1,id_case2)
        id_c2<-union(setdiff(id_case1,id_case2),setdiff(id_case2,id_case1))
        if(n.case <= length(id_c1)){
          id_case<-sample(id_c1,n.case)
        }else if (n.case > length(id_c1) && n.case <= length(id_c1)+length(id_c2)){
          id_c3<-sample(id_c2,n.case - length(id_c1))
          id_case<-c(id_c1,id_c3)
        }else {					
          id_case3<-union(id_c1,id_c2)
          id_c4<-setdiff(1:n1,id_case3)
          n.needed<-n.case - length(id_case3)
          
          id_c5<-sample(id_c4,n.needed,prob=mu[id_c4])
          id_case<-union(id_case3,id_c5)
        }
        
        res.out[id_case,i]<-1
      }	
      res.out<-res.out - mu
    } else if(type.Resampling=="perturbation"){
      res.out<-matrix(rnorm(n1*n.Resampling,mean=0,sd=1),ncol=n.Resampling)
      res.out<-res.out * res
      stop("Error: Perturbation is no more provided!")
    } else {
      if(is.null(res.out)){
        stop("Error: Wrong resampling method!")
      }
    }
    
  }
  
  return(list(res=res, X1=X1,res.out=res.out,out_type="D", 
              n.Resampling=n.Resampling, type.Resampling=type.Resampling,
              id_include=id_include, mu=mu,pi_1=pi_1))
  
}

#
#	type	:  	
#		: permu - permutation
#		: bootstrap - bootstrap
#		: 
#
SKAT_Null_Model = function(formula, data=NULL, out_type="C", n.Resampling=0, type.Resampling="bootstrap", Adjustment=TRUE){
  
  SKAT_MAIN_Check_OutType(out_type)
  
  
  # check missing 
  obj1<-model.frame(formula,na.action = na.omit,data)
  obj2<-model.frame(formula,na.action = na.pass,data)
  
  n<-dim(obj2)[1]
  n1<-dim(obj1)[1]
  id_include<-SKAT_Null_Model_Get_Includes(obj1,obj2)
  n1<-length(id_include)
  
  # Check whether n < 2000 and out_type="D", apply the adjustment 
  # if No_Adjustment = FALSE
  if(n1< 2000 && out_type=="D" && Adjustment){
    MSG<-sprintf("Sample size (non-missing y and X) = %d, which is < 2000. The small sample adjustment is applied!\n",n )
    cat(MSG)
    n.Resampling.kurtosis=10000
    #if(n > 1000){
    #	n.Resampling.kurtosis = floor(10000 - (n-1000) * 5)	
    #} 
    #if(n.Resampling.kurtosis < 5000){
    #	n.Resampling.kurtosis = 5000
    #}
    
    
    re<-SKAT_Null_Model_MomentAdjust(formula, data, n.Resampling, type.Resampling=type.Resampling, is_kurtosis_adj=TRUE, n.Resampling.kurtosis=n.Resampling.kurtosis)
    return(re)
  }
  
  
  if(n - n1 > 0){
    MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n - n1)
    warning(MSG,call.=FALSE)
  }
  
  if(out_type=="C"){
    re<-Get_SKAT_Residuals.linear(formula, data, n.Resampling, type.Resampling, id_include )
  } else {
    re<-Get_SKAT_Residuals.logistic (formula, data, n.Resampling, type.Resampling, id_include )
  }
  
  
  class(re)<-"SKAT_NULL_Model"
  re$n.all<-n
  return(re)
  
}


SKAT_Null_Model_MomentAdjust = function(formula, data=NULL, n.Resampling=0, type.Resampling="bootstrap", is_kurtosis_adj=TRUE, n.Resampling.kurtosis=10000){
  
  
  # check missing 
  obj1<-model.frame(formula,na.action = na.omit,data)
  obj2<-model.frame(formula,na.action = na.pass,data)
  
  n<-dim(obj2)[1]
  n1<-dim(obj1)[1]
  id_include<-SKAT_Null_Model_Get_Includes(obj1,obj2)
  
  
  if(n - n1 > 0){
    MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n - n1)
    warning(MSG,call.=FALSE)
  }
  
  re1<-Get_SKAT_Residuals.logistic (formula, data, n.Resampling, type.Resampling, id_include )
  re1$n.all<-n
  re2<-NULL
  
  if(is_kurtosis_adj == TRUE){
    re2<-Get_SKAT_Residuals.logistic (formula, data, n.Resampling.kurtosis, type.Resampling, id_include )
  }
  
  class(re1)<-"SKAT_NULL_Model"
  re<-list(re1=re1, re2=re2, is_kurtosis_adj= is_kurtosis_adj, type = "binary")
  
  
  class(re)<-"SKAT_NULL_Model_ADJ"
  return(re)
  
}


SKAT_Null_Model_Get_Includes<-function( obj_omit, obj_pass){
  
  ID1<-rownames(obj_omit)
  ID2<-rownames(obj_pass)
  
  d1<-data.frame(ID=ID1)
  d2<-data.frame(ID=ID2, idx=1:length(ID2))
  
  d3<-merge(d1, d2,by.x="ID", by.y="ID")
  id_include = sort(d3$idx)
  
  return(id_include)
}


emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}


emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

#######################
#
# Changed by SLEE

SKAT.emma.eigen.R.wo.Z <- function(K, X) {
  
  n <- nrow(X)
  q <- ncol(X)
  XX1<-X %*% (solve(crossprod(X,X)))
  K1<-t(X) %*% K
  K2<-K %*% X
  
  #Mat<-K - XX1 %*% K1 - K2 %*% t(XX1) + XX1 %*% ((K1 %*% X) %*% t(XX1)) - XX1 %*% t(X)
  Mat<-K - K2 %*% t(XX1) - XX1 %*% (K1  +  ((K1 %*% X) %*% t(XX1)) -  t(X))
  
  diag(Mat) = diag(Mat) +1
  eig <- eigen(Mat,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

#
#	X : covariates
#	Z = NULL
#	Updated by SLEE 07/21/2017
SKAT_NULL_emmaX<- function(formula, data=NULL, K=NULL, Kin.File=NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, Is.GetEigenResult=FALSE) {
  
  # check missing 
  obj1<-model.frame(formula,na.action = na.omit,data)
  obj2<-model.frame(formula,na.action = na.pass,data)
  
  n1<-dim(obj2)[1]
  n<-dim(obj1)[1]
  id_include<-SKAT_Null_Model_Get_Includes(obj1,obj2)
  X<-model.matrix(formula,data=data)
  X<-Get_SKAT_Residuals.Get_X1(X)
  
  y<-model.response(obj1)
  q <- ncol(X)
  # n and n1 are opposite (compare to SKAT)
  if(n1 - n > 0){
    MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n1 - n)
    warning(MSG,call.=FALSE)
  }
  
  stopifnot(nrow(X) == n)
  
  ########################################################
  # Read kinship 
  
  if(is.null(K) && is.null(Kin.File)){
    stop("both K and Kin.File are NULL!")
  } else if(is.null(K)){
    if(!file.exists(Kin.File)){
      msg<-sprintf("File %s is not exist!", Kin.File)
      stop(msg)
    }
    cat("Read ", Kin.File, "\n")
    K = as.matrix(read.delim(Kin.File, header=FALSE, colClasses="numeric", sep = "\t"))
    t <- nrow(K)
    cat("Read complete. ", Kin.File, " has ", t, " rows! \n")
  }
  
  if(!is.matrix(K)){
    stop("K is not a matrix!")
  }
  
  if(n1 - n > 0){
    K<-K[id_include, id_include]
  }
  t <- nrow(K)
  stopifnot(ncol(K) == t)
  
  #######################################################
  # Estimate parameters
  
  eig.R <- emma.eigen.R.wo.Z(K,X)
  #eig.R<-SKAT.emma.eigen.R.wo.Z(K,X) # This code sometimes makes an error!
  
  etas <- crossprod(eig.R$vectors,y)
  #etas<-crossprod(eig.R$vectors,res)
  
  
  logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
  m <- length(logdelta)
  delta <- exp(logdelta)
  
  Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
  Etasq <- matrix(etas*etas,n-q,m)
  LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
  dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
  
  optlogdelta <- vector(length=0)
  optLL <- vector(length=0)
  if ( dLL[1] < esp ) {
    optlogdelta <- append(optlogdelta, llim)
    optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
  }
  if ( dLL[m-1] > 0-esp ) {
    optlogdelta <- append(optlogdelta, ulim)
    optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
  }
  
  for( i in 1:(m-1) ){
    if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ){
      r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
      optlogdelta <- append(optlogdelta, r$root)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
    }
  }
  
  #############################################################################
  # variance term : maxva K + maxve I
  # eig.R$vectors $*$ diag(maxva * eig.R$values  + maxdelta) %*% t(eig.R$vectors)
  
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)    # additive effect, genetic
  
  maxve <- maxva*maxdelta	# noise
  
  #############################################################################
  # Get V^-1 (y - X \beta)
  # = V-1 y - V-1 X (X'V-1X)-1 X'V-1 y
  # = P y
  # P = V-1 -  V-1 X (X'V-1X)-1 X'V-1
  
  va=maxva
  ve=maxve 
  
  #if(Is.EPACTS){
  #	idx.lambda<-which(abs(eig.R$values) > mean(abs(eig.R$values)) / 10000)
  #  
  #	lambda_inv<-1/(va * eig.R$values[idx.lambda] + ve)
  #
  #	P = eig.R$vectors[,idx.lambda] %*% (t(eig.R$vectors[,idx.lambda]) * (lambda_inv))
  #	res = P %*% y
  #	
  
  
  V = va * K
  diag(V) = diag(V) + ve
  
  # numerical reason
  # added
  if(ve/va < 0.01){
    eig.V<-eigen(V, symmetric=TRUE)
    idx.lambda<-which(abs(eig.V$values) > mean(abs(eig.V$values)) / 10000)
    V_inv<- eig.V$vectors[,idx.lambda] %*% (t(eig.V$vectors[,idx.lambda]) * (1/eig.V$values[idx.lambda]) )
    
  } else {
    V_inv<- solve(V)
  }
  
  XVX = t(X) %*% (V_inv %*% X)
  XVX_inv<-solve(XVX)
  XV_inv = t(X) %*% V_inv
  P = V_inv -  t(XV_inv) %*% (XVX_inv %*% XV_inv)
  res = P %*% y
  
  re<-list( LL=maxLL, va=va, ve=ve, P=P, res=res, id_include=id_include)
  if(Is.GetEigenResult){
    re$eig.R=eig.R
  }
  
  class(re)<-"SKAT_NULL_Model_EMMAX"
  return (re)
}


SKAT_emmaX = function( Z, obj, kernel= "linear.weighted", method="davies", weights.beta=c(1,25), weights=NULL, impute.method="fixed"
                       , r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, max_maf=max_maf, estimate_MAF=1, SetID = NULL){
  
  
  if(!Check_Class(obj, "SKAT_NULL_Model_EMMAX")){
    stop("ERROR: obj is not an returned object from SKAT.emmaX.null")
  } 
  #if(method=="optimal.adj"){
  #	stop("SKAT-O is not implemented for SKAT_EmmaX yet!")
  #}
  #if(length(r.corr) > 1){
  #	stop("SKAT-O is not implemented for SKAT_EmmaX yet! r.corr should be a scalar.")
  #}
  
  
  m = ncol(Z)
  n = nrow(Z)
  
  # Added by SLEE 4/24/2017
  out.method<-SKAT_Check_Method(method,r.corr, n=n, m=m)
  method=out.method$method
  r.corr=out.method$r.corr
  IsMeta=out.method$IsMeta
  
  if(method=="optimal.adj"){
    IsMeta=TRUE
  }
  
  SKAT_Check_RCorr(kernel, r.corr)
  #
  
  #####################################
  # Check genotypes and parameters
  out.z<-SKAT_MAIN_Check_Z(Z, n, obj$id_include, SetID, weights, weights.beta, impute.method, is_check_genotype
                           , is_dosage, missing_cutoff, max_maf=max_maf, estimate_MAF=estimate_MAF)
  if(out.z$return ==1){
    out.z$param$n.marker<-m
    return(out.z)
  }
  Z = out.z$Z.test
  weights = out.z$weights
  res = obj$res
  
  if(!IsMeta){
    re = SKAT_emmaX_work(Z=Z, obj=obj, kernel=kernel, method=method, weights=weights, r.corr=r.corr)
  } else {
    re = SKAT_RunFrom_MetaSKAT(res=obj$res,Z=Z, X1=NULL, kernel=kernel, weights=weights
                               , P0=obj$P, out_type="V", method=method, res.out=NULL, n.Resampling=0, r.corr=r.corr)
    
  }
  
  re$param$n.marker<-m
  re$param$n.marker.test<-ncol(Z)
  re$test.snp.mac<-SingleSNP_INFO(out.z$Z.test)
  
  return(re)
}

SKAT_emmaX_work = function( res, Z, obj, kernel, method, weights=NULL, r.corr=0){
  
  m = ncol(Z)
  n = nrow(Z)
  
  # Weighted Linear Kernel 
  if (any(kernel == "linear.weighted")) {
    Z = t(t(Z) * (weights))
  }
  
  if(r.corr == 1){
    Z<-cbind(rowSums(Z))
  } else if(r.corr > 0){
    
    p.m<-dim(Z)[2]	
    R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
    L<-chol(R.M,pivot=TRUE)
    Z<- Z %*% t(L) 
  }
  
  
  # get Q
  Q.Temp = t(obj$res)%*%Z
  Q = Q.Temp %*% t(Q.Temp)/2
  Q.res=NULL
  
  # Get Z' P0 Z
  W.1= t(Z) %*% (obj$P %*% Z) # t(Z) P0 Z
  
  if( method == "liu.mod" ){
    
    out<-Get_Liu_PVal.MOD(Q, W.1, Q.res)    
    
  } else if( method == "davies" ){
    
    out<-Get_Davies_PVal(Q, W.1, Q.res)    
    
  } else {
    stop("Invalid Method!")
  }
  
  
  re<-list(p.value = out$p.value, Test.Type = method, Q = Q, param=out$param ) 
  
  
  return(re)
}



#
# x is either y or SKAT_NULL_Model 
#

SKAT_emmaX.SSD.OneSet = function(SSD.INFO, SetID, obj, ..., obj.SNPWeight=NULL){
  
  id1<-which(SSD.INFO$SetInfo$SetID == SetID)
  if(length(id1) == 0){
    MSG<-sprintf("Error: cannot find set id [%s] from SSD!", SetID)
    stop(MSG)
  }	
  SetIndex<-SSD.INFO$SetInfo$SetIndex[id1]
  
  re = SKAT_emmaX.SSD.OneSet_SetIndex(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=obj.SNPWeight)
  
  return(re)
}

SKAT_emmaX.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, obj, ..., obj.SNPWeight=NULL){
  
  re1 = SKAT.SSD.GetSNP_Weight(SSD.INFO, SetIndex, obj.SNPWeight=obj.SNPWeight)
  SetID = SSD.INFO$SetInfo$SetID[SetIndex]
  if(!re1$Is.weights){
    re<-SKAT_emmaX(re1$Z, obj, ...)
  } else {
    
    re<-SKAT_emmaX(re1$Z, obj, weights=re1$weights, ...)
  }
  
  return(re)
  
}

#
# Only SKAT_Null_Model obj can be used
#
SKAT_emmaX.SSD.All = function(SSD.INFO, obj, ..., obj.SNPWeight=NULL){
  
  N.Set<-SSD.INFO$nSets
  OUT.Pvalue<-rep(NA,N.Set)
  OUT.Marker<-rep(NA,N.Set)
  OUT.Marker.Test<-rep(NA,N.Set)
  OUT.Error<-rep(-1,N.Set)
  OUT.Pvalue.Resampling<-NULL
  
  Is.Resampling = FALSE
  n.Resampling = 0
  
  pb <- txtProgressBar(min=0, max=N.Set, style=3)
  for(i in 1:N.Set){
    
    try1 = try(SKAT_emmaX.SSD.OneSet_SetIndex(SSD.INFO=SSD.INFO, SetIndex=i, obj=obj, ..., obj.SNPWeight=obj.SNPWeight))
    if(!Is_TryError(try1)){
      
      err.msg<-geterrmessage()
      msg<-sprintf("Error to run SKAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
      warning(msg,call.=FALSE)
      
    } else {
      
      OUT.Pvalue[i]<-re$p.value
      OUT.Marker[i]<-re$param$n.marker
      OUT.Marker.Test[i]<-re$param$n.marker.test
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test)
  re<-list(results=out.tbl)
  class(re)<-"SKAT_SSD_ALL"
  
  return(re)	
}