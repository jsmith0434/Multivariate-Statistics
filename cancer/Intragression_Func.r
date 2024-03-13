##========================================================================================##
# Codes for nonparametric estimation of the row space of the dependence kernel in the model
# E[Y] = PhiM, where E[Y] is no longer modeled through a link function
# Assumption: Independence, smoothness, bounded parameter domain, and higher order moments

## Copyright:  Xiongzhi Chen and John D Storey
## Distributed under GNU license without explicit or implicit responsibility for the consequences of
## using these codes. These codes should not be used for commercial purposes.


#######################################################################################
#===== section 1: Modified nScree from nFactors package to estimate # of factors  =====#
#######################################################################################

nScreeA = function (eig = NULL, x = eig, aparallel = NULL, cor = TRUE,
    model = "components", criteria = NULL, ...)
{
    eig <- eigenComputes(x, cor = cor, model = model, ...)

    if (any(is.infinite(eig)) | any(is.na(eig)) | any(is.nan(eig)))
          cat("-- nScree operates ill-posed matrix","\n")

    if (is.null(aparallel))
        aparallel <- rep(1, length(eig))
    nk <- length(eig)
    k <- 1:nk
    proportion <- eig/sum(eig)
    cumulative <- proportion
    if (is.null(criteria))
        criteria <- mean(eig)
    for (i in 2:nk) cumulative[i] = cumulative[i - 1] + proportion[i]
    proportion[proportion < 0] <- 0
    cond1 <- TRUE
    cond2 <- TRUE
    i <- 0
    pred.eig <- af <- rep(NA, nk)                    #i<nk changed into nk-1
    while ((cond1 == TRUE) && (cond2 == TRUE) && (i < nk-1)) {
        i <- i + 1
       # cat("--nScree accessing",i+1, "and",c(i+1, nk), "\n")
        ind <- k[c(i+1, nk)]; # cat("--nscree: ind", ind,"\n")
        eigA = eig[c(i+1, nk)]

        # edited lines
        rto = eigA/ind
        if (any(is.infinite(rto)) | any(is.na(rto)) | any(is.nan(rto)))
              cat("--nScree:singular index or selected eigvenvalues", "\n")

        if (max(rto) == min(rto)) {
            cf = min(rto)
        } else {
        vp.p <- lm(eigA ~ ind)
        cf = coef(vp.p) }

        vp.prec <- pred.eig[i] <- sum(c(1, i) * cf)
        # end of edited lines

        ## original lines
       # vp.p <- lm(eig[c(i + 1, nk)] ~ ind)
#        vp.prec <- pred.eig[i] <- sum(c(1, i) * coef(vp.p))
        ## end of original lines

        cond1 <- (eig[i] >= vp.prec)
        cond2 <- (eig[i] >= aparallel[i])
        nc <- i - 1
    }
    tag <- 1
    for (j in 2:(nk - 1)) {
        if (eig[j - 1] >= aparallel[j - 1]) {
            af[j] <- (eig[j + 1] - 2 * eig[j]) + eig[j - 1]
        }
    }
    if (model == "components")
        p.vec <- which(eig >= aparallel, TRUE)
    else p.vec <- which((eig - aparallel) >= 0 & eig >= criteria)
    npar <- sum(p.vec == (1:length(p.vec)))
    nkaiser <- sum(eig >= rep(criteria, nk))
    naf <- which(af == max(af, na.rm = TRUE), TRUE) - 1
    for (i in (nc + 1):(nk - 2)) {
        ind <- k[c(i + 1, nk)]
        vp.p <- lm(eig[c(i + 1, nk)] ~ ind)
        vp.prec <- pred.eig[i] <- sum(c(1, i) * coef(vp.p))
    }
    for (j in 2:(nk - 1)) af[j] <- (eig[j + 1] - 2 * eig[j]) +
        eig[j - 1]
    coc <- rep("", nk)
    coc[nc] = "(< OC)"
    caf <- rep("", nk)
    caf[naf] = "(< AF)"
    result <- (list(Components = data.frame(noc = nc, naf = naf,
        nparallel = npar, nkaiser = nkaiser), Analysis = data.frame(Eigenvalues = eig,
        Prop = proportion, Cumu = cumulative, Par.Analysis = aparallel,
        Pred.eig = pred.eig, OC = coc, Acc.factor = af, AF = caf),
        Model = model))
    class(result) <- "nScree"
    return(result)
} # end of scree function

########################################################################################
##=== Section 2: Estimate # of factors by thresholding scaled eigenvalues        ========#
########################################################################################

ModDimfunc <- function(dat, SecondPar, dist, ita = -1/3, aL=10000, nsub=100, plot=FALSE){
  
  if (nsub < 100) {
    stop("--Number of variables should be at least 100")
  }
  
 dims <- dim(dat)
 a <- seq(0,30,length=aL)   # range adjustable ; a is a vector
 b <- sample(1:dims[1])
 n <- floor(dims[1]/nsub)
 if ( n ==0) 
  stop("--Number of rows in data too small")
  
 rhat <- matrix(0,nrow=aL,ncol=nsub)
 
 # obtain Dhat first
 Vhat = EstVarsFunc(dat, SecondPar, dist)
 Dhat = diag(colMeans(Vhat))
                     	
 for(j in 1:nsub){
   dats <- dat[b[1:(j*n)],]
   
   # R^hat with one D^hat
       Rhat = (j*n)^(-1)*t(dats)%*%dats - Dhat
       RsvdD = svd(Rhat)$d 
   # end of R^hat
   
   v <- c(rep(T, aL), rep(F, dims[2]))
   v <- v[order(c(a*(j*n)^ita*dims[2],RsvdD), decreasing = T)]  
   u <- 1:length(v)
   w <- 1:aL
   rhat[,j] <- rev((u[v==TRUE]-w))
 }
 ss <- rowVars(rhat)

 peak <- which.max(ss)
 start <- which.max(c(rep(1e5,peak),ss[(peak+1):aL]) < 0.5*(ss[peak] + ss[aL]))
 finish <- which.max(ss*c(rep(0,start),rep(1,aL-start)) > 0.75*(ss[peak] + ss[aL]))
 if(finish==1){finish <- aL}

 est <- modefunc(rhat[start:finish,nsub])
 if(plot){
   par(mar=c(5,5,5,5))
   plot(a*dims[2]*(j*nsub)^ita,rhat[,nsub],type="l",lwd=3,col="blue",xlab="Threshold",ylab="",ylim=c(0,dims[2]),cex.axis=1.5,yaxt="n")
   lines(a*dims[2]*(j*nsub)^ita,ss/max(ss)*dims[2],type="l",col="red",lwd=3)
   axis(4,at=seq(0,dims[2],length=5),labels=round(seq(0,max(ss),length=5),2),cex.axis=1.5,col.axis="red")
   axis(2,at=seq(0,dims[2],length=5),labels=round(seq(0,dims[2],length=5),2),cex.axis=1.5,col.axis="blue")
   lines(rep(a[start]*dims[2]*(j*nsub)^ita,aL),seq(0,dims[2],length=aL),col="green",lwd=3,lty=2)
   lines(rep(a[finish]*dims[2]*(j*nsub)^ita,aL),seq(0,dims[2],length=aL),col="green",lwd=3,lty=2)
   lines(a*dims[2]*(j*nsub)^ita,rep(est,aL),lty=2,lwd=2,col="black")
   mtext("Estimated # of Factors",side=2,line=3,col="blue",cex=1.5)
   mtext("Variance of Estimate",side=4,line=3,col="red",cex=1.5)
   legend(a[nsub]*dims[2]*(j*nsub)^ita,dims[2],legend=c("2nd Stability Interval",paste("Estimated Value =",est)),col=c("green","black"),lwd=c(3,3),lty=c(2,2),bg="white")	
 }
 return(est)
}


## Borrowed from genefilter
rowVars <- function(x, ...) 
{
    sqr = function(x) x * x
    n = rowSums(!is.na(x))
    n[n <= 1] = NA
    return(rowSums(sqr(x - rowMeans(x, ...)), ...)/(n - 1))
}

modefunc <- function(x){
	return(as.numeric(names(sort(-table(x)))[1]))
}

## Estimate variaces for use with abov "ModDimFunc"

EstVarsFunc = function(Y = NULL, v = NULL, Dist= NULL)
{
  
  if (Dist != "Normal" & Dist != "Poisson") {
    if(is.null(v)) {
      stop("---Please specify the second (fixed) parameter of the distribution.")
    }
    if (v == 0) {
      stop("---The second (fixed) parameter of the distribution cannot be zero.")  
    } 
  }
  
  if (Dist == "Normal") {   Y0 = matrix(1,nrow(Y),ncol(Y))    }
  
  if (Dist == "Poisson") {    Y0 = Y     }
  
  if (Dist == "Binomial") {    Y0 = (v*Y - Y^2)/(v-1) }
  
  if (Dist == "NegativeBinomial") {    Y0 = (v*Y + Y^2)/(v+1)     }
  
  if (Dist == "Gamma") {     Y0 = Y^2/(v+1)     }
  
  return(Y0)
}  # end 

###############################################################
#===== section 4*:  compute R = k^{-1} Y^T Y - D^{\hat}   =====#
###############################################################

GetMainMatrixFunc = function(Y = NULL, v = NULL, Dist= NULL) {
  
  k = dim(Y)[1]
  Y0 = EstVarsFunc(Y, v, Dist) # get variances
  
  # compute Dhat and Rhat matrix
  Dhat = diag(colMeans(Y0))
  Rhat = k^(-1)*t(Y)%*%Y - Dhat 
  
  # return R matrix
  return(list("Rhat"=Rhat, "Dhat"= Dhat))
} 


###############################################################
#===== section 5: Estimate # of factors   =====#
###############################################################

EstRankFunc = function(Y = NULL, v= NULL, Dist=NULL, rEstMethod = "Ratio") {
  
  k = dim(Y)[1]; n = dim(Y)[2]
  # if estimate r by modified leek
  if (rEstMethod == "Scaling")  {
    rEst = ModDimfunc(Y,v,Dist)    
  } else {
    
    Rhat=  GetMainMatrixFunc(Y, v, Dist)$Rhat
    Rsvd = svd(Rhat);   RsvdD = Rsvd$d ;  RsvdV = Rsvd$v
  
  if (any(is.infinite(RsvdD)) | any(is.na(RsvdD)) | any(is.nan(RsvdD)))
    stop("--Ill-posed matrix for estimation of latent space")
  
  # estimate r by nScreeA
  if (rEstMethod == "Scree") {
    screeRs = nScreeA(RsvdD)
    # use majority vote to get r1
    rcad = screeRs$Components
    rcad3 = as.vector(as.matrix(rcad))
    tr = as.matrix(table(rcad3)) # cases-by-1
    ridx = which.max(tr)
    rScree = as.numeric(rownames(tr))[ridx]  # Desai uses $naf
    rEst = max(0, min(rScree,n))
  } 
  
  # est r by ratio of eigenvalues
  if (rEstMethod == "Ratio") {
    rs = RsvdD[1:(floor(0.5*n)-1)]/RsvdD[2:floor(0.5*n)]
    rEst = which.min(rs)+1
  } 
  
} # end of else

  if (rEst ==0) {
    warning("---Estimated r is zero")  }
  
  # return
  return("rEst" = rEst)
}


############################################################
##=== Section 6: Estimate row space of M         ========#
##########################################################

EstSpanMFunc = function(Y= NULL, v = NULL, Dist = NULL, r = NULL, rEstMeth = "Ratio")
  {
       k = dim(Y)[1]; n = dim(Y)[2]
       
       RhatDhat = GetMainMatrixFunc(Y, v, Dist)
       Rhat = RhatDhat$Rhat 
         
       if (is.null(r)){
         rEst = EstRankFunc(Y, v,Dist,rEstMeth)         
       } else {
         rEst = r
       }
      
       Rsvd = svd(Rhat);  RsvdV = Rsvd$v
       Mh = t(RsvdV[,1:rEst])
     # return 
     return(list("Mhat" = Mh, "rEst" = rEst))

}   # end of function 


################################################
##=== Section 8:  compute distance        ========#
################################################

# compute distance between M and Mh with rEst
NormProj_Func = function(Mhat = NULL, M = NULL) {
      n = ncol(Mhat); rEst = nrow(Mhat)
      Mh = Mhat
      HMh = t(Mh)%*%Mh
      HM =  t(M)%*%ginv(M%*%t(M))%*%M
      resMh = (diag(n) - HMh)%*%t(M)
      resV = (diag(n) - HM) %*%t(Mh)
      M2Norm = sqrt(sum(resMh^2)+sum(resV^2))/sqrt(n*rEst)

     return(M2Norm)
}

# compute directional measure
NormAProj_Func = function(Mhat = NULL, M = NULL) {
      n = ncol(Mhat); rEst = nrow(Mhat)
      Mh = Mhat
      HMh = t(Mh)%*%Mh
      resMh = (diag(n) - HMh)%*%t(M)
      M2Norm = sqrt(sum(resMh^2))/sqrt(n*rEst)

     return(M2Norm)
}

# another directional measure with range in [0,1]
NormBProj_Func = function(Mhat = NULL, M = NULL) {
  n = ncol(Mhat); rEst = nrow(Mhat)
  Mh = Mhat
  HMh = t(Mh)%*%Mh
  resMh = (diag(n) - HMh)%*%t(M)
  M2Norm = sqrt(sum(resMh^2))/sqrt(sum(M^2))
  
  return(M2Norm)
}