#'fit a conditional logistic regression with lasso 
#'
#' Fit a sequence of conditional logistic regression with lasso penalty,for small to large sized samples
#' 
#' @param X Input matrix, of dimension nobs x nvars; each row is an observation vector
#' @param y Binary response variable, with 1 for cases and 0 for controls
#' @param strata Vector with stratum membership of each observation
#' @param fraction Sequence of lambda values
#' @param nbfraction The number of lambda values - default is 100
#' @param nopenalize List of coefficients not to penalize starting at 0
#' @param BACK If TRUE, use Backtracking-line search -default is TRUE
#' @param BIC If TRUE compute BIC-default is FALSE
#' @param standardize Logical flag for x variable standardization, prior to fitting the model sequence.
#' @param maxit Maximum number of iterations of outer loop - default is 100
#' @param maxitB Maximum number of iterations  in Backtracking-line search - default is 100
#' @param thr Threshold for convergence in lassoshooting. Default value is 1e-10. Iterations stop when max absolute parameter change is less than thr
#' @param tol Threshold for convergence-default value is 1e-10
#' @param epsilon ratio of smallest to largest value of regularisation parameter  at which we find parameter estimates
#' @param trace If TRUE the algorithm will print out information as iterations  proceed -default is TRUE
#' @param log If TRUE, fraction are spaced uniformily on the log scale
#' @param adaptive If TRUE adaptive lasso is fitted-default is FALSE
#' @param separate If TRUE, the weights in adaptive lasso are build serarately using univariate models. Default is FALSE,
#' weights are build using multivariate model
#' @param ols If TRUE, weights less than 1 in adaptive lasso are set to 1. Default is FALSE
#' @param p.fact Weights for adaptive lasso
#' @param remove If TRUE, invariable covariates are removed-default is FALSE
#' @return An object of type \code{clogitLasso} which is a list with the following
#' components:
#' 
#' \item{beta}{nbfraction-by-ncol matrix of estimated coefficients. First row has all 0s}
#' 
#' \item{fraction}{A sequence of regularisation parameters at which we obtained the fits}
#' 
#' \item{nz}{A vector of length nbfraction containing the number of nonzero parameter estimates for
#'  the fit at the corresponding regularisation parameter}
#'  
#'  \item{w}{}
#'   
#'  \item{bic}{A vector of length nbfraction containing BIC  for
#'  the fit at the corresponding regularisation parameter}
#'  @author Marius Kwemou and Marta Avalos
#'  @details The sequence of models implied by fraction  is fit by  IRLS (iteratively 
#'  reweighted least squares) algorithm. 
#'  by coordinate descent with warm starts and sequential strong rules
#'  @references Avalos, M., Pouyes, H., Grandvalet, Y., Orriols, L., & Lagarde, E. (2015). \emph{Sparse conditional logistic 
#'  regression for analyzing large-scale matched data from epidemiological studies: a simple algorithm.} BMC bioinformatics, 16(6), 1.
#'  @importFrom lassoshooting lassoshooting
#'  @importFrom stats glm var aggregate
#'  @examples
#'  \dontrun{
#'  #generate data}
#'  y = rep(c(1,0), 100)
#'  X = matrix (rnorm(20000, 0, 1), ncol = 100) # pure noise
#'  strata = sort(rep(1:100, 2))
#'  
#'  fitLasso = clogitLasso(X,y,strata,log=TRUE)
#'  @export


clogitLasso=function(X,
                     y,
                     strata,
                     fraction=NULL,
                     nbfraction=100,
                     nopenalize=NULL,
                     BACK=TRUE,
                     BIC=FALSE,
                     standardize=FALSE,
                     maxit=100,
                     maxitB=500,
                     thr=1e-10,
                     tol=1e-10,
                     epsilon=0.0001,
                     trace=TRUE,
                     log=FALSE,
                     adaptive=FALSE,
                     separate=FALSE,
                     ols=FALSE,
                     p.fact,
                     remove=FALSE){
    
    if (length(y) != nrow(X)) 
        stop("Please ensure that each observation has predictors and response")
    
    # Matrix of differences 
    
    x = X[y==1,]-X[y==0,]
    
    
    if(missing(strata)){
        stop("'strata' is missing")
    }else{
        
        y = aggregate(as.data.frame(y),by = list(strata),diff)$y
        
        if (any(y != 1 )) 
            stop("Response vector should be 1 case and 1 control in each strata, starting by the case")
    }
    
    
    
    #Algorithme IRLS-LASSOSHOOTING
    #require(lassoshooting)
    x=as.matrix(x)
    n=dim(x)[1]
    
    #invariable covariates are removed
    sds=stand(x,n,unit=F)    
    if(remove){remove<-(sds==0)}
    else{remove<-(sds<=-1000)}    #don't remove invariable covariates
    x<-x[,!remove]
    m=dim(x)[2]
    
    # stabilize/standardize. 
    if (standardize) {
        x <- x / matrix(sds[!remove], n, m, byrow=T)
    }  
    
    #Calculate regularization parameter for the 'lassoshooting'
    if (is.null(fraction)){
        if(missing(nbfraction)){nbfraction=100}
        
        
        fraction=frac(epsilon=epsilon,log=log,x=x,n=n,m=m,nbfraction=nbfraction)
    }else{nbfraction=length(fraction)}
    
    
    if(adaptive){
        if(missing(p.fact)){p.fact=ponderation(m=m,separate=separate,n=n,x=x,ols=ols)}
        else{p.fact=p.fact[!remove]}
    }
    
    #Estimation of coefficient for each value of fraction 
    nb_coef_non_nuls=c()
    beta=matrix(0,nbfraction,length(remove))
    nz=0
    W=rep(0,n)
    bic=rep(0,length(fraction))
    
    if(!all(remove)){
        betanew=matrix(0,nbfraction,m)
        nbf=nbfraction
        for (i in (1:nbf)){
            if(trace){if (mod(i,40)==0){cat("fraction ",i,"\n")}}
            if (i==1) {betaold=rep(0,m)} else {betaold=betanew[(i-1),]}
            fold=likelihood.diff(x,betaold)
            a=0
            
            while(a<maxit){
                a=a+1
                
                z=rep(0,n)
                
                z=x%*%betaold+1/sigmoid(x,betaold)
                lambda=sigmoid(x,betaold)*(1-sigmoid(x,betaold))
                X=sqrt(as.vector(lambda))*x
                Y=sqrt(as.vector(lambda))*z
                rm(z)
                
                if(adaptive==TRUE){
                    gamma=(lassoshooting(X=X,y=Y,lambda=fraction[i],penaltyweight=p.fact,thr=thr,nopenalize=nopenalize))$coefficient
                    rm(X);rm(Y)
                } else {
                    XtX=t(X)%*%X
                    Xty=t(X)%*%Y
                    rm(X);rm(Y)
                    gamma=(lassoshooting(XtX=XtX,Xty=Xty,lambda=fraction[i],thr=thr,nopenalize=nopenalize))$coefficient
                    rm(XtX);rm(Xty)
                }
                
                #Backtracking-line search
                if(BACK){
                    step=gamma-betaold
                    t=1 
                    delta=grad.diff(x,betaold)%*%step
                    for (l in (1:maxitB)){
                        gamma=betaold+t*step
                        if (likelihood.diff(x,gamma)<=(likelihood.diff(x,betaold)+0.3*t*delta)) break
                        t=0.9*t
                    }
                    betanew[i,]=(1-t)*betaold+t*gamma
                }else{betanew[i,]=betaold}
                
                fnew=likelihood.diff(x,betanew[i,])
                if (abs(fnew-fold)/abs(fnew)<tol) break
                betaold=betanew[i,]
                fold=fnew
                
                
                
            }
            
            nb_coef_non_nuls[i]=sum(betanew[i,]!=0)
        }
        dimnames(betanew)[2]=dimnames(x)[2]
        
        bic.fraction=c()
        if(BIC){
            for (i in seq(length(fraction)))
                bic.fraction[i]=(-2*likelihood.diff(x,betanew[i,]))+nb_coef_non_nuls[i]*log(n)
        }
        
        
        
        if (standardize) {
            for (i in seq(length(fraction)))
                betanew[i,betanew[i,] != 0] <- betanew[i,betanew[i,] != 0] / sds[betanew[i,] != 0]
        }
        
        beta[,!remove]=betanew
        nz=apply(betanew,1,function(x) sum(x!=0))
        W=lambda
        bic=bic.fraction
    }
    
    list(beta=beta,fraction=fraction,nz=nz,W=W,bic=bic)
}




#tools
"%+%"<- function(x,y) paste(x,y,sep="")

# sigmoid function
sigmoid=function(x,beta){1/(1+exp(-x%*%beta))}

#likelihood function
likelihood.diff=function(x,beta){sum(-log(1/(1+exp(-x%*%beta))))}

# gradient function  
grad.diff=function(x,beta){colSums(sweep(-x,MARGIN=1,STATS=(1-1/(1+exp(-x%*%beta))),FUN="*"))}

#Calcul of fraction
frac=function(epsilon,log,x,n,m,nbfraction){
    fracmax=max(abs((t(x)%*%rep(1/2,n))))
    fracmin=fracmax*epsilon
    if (log==TRUE){fraction=exp(seq(from=log(fracmax),to=log(fracmin),length.out=nbfraction))
    } else {fraction=seq(from=fracmax,to=fracmin,length.out=nbfraction)}
    
    return(fraction)
}

#weigths
ponderation=function(m,separate,n,x,ols){
    beta.ols<-as.vector(rep(0,m))
    if(separate){
        for(j in 1:m){beta.ols[j] <- glm(rep(1,n)~ x[,j]+0,family="binomial")$coefficients}
    }else{beta.ols <- glm(rep(1,n)~ x+0,family="binomial")$coefficients}				
    
    p.fact <- abs(beta.ols)^(-1)
    if(ols==FALSE){p.fact[p.fact<1]=1}
    
    return(p.fact)
}

# stabilize/standardize
stand=function(x,n,unit=T) {
    vars <- apply(x,2,var) * (n-1)/n
    if(unit){vars[vars == 0] <- 1}
    sds <- sqrt(vars)
    
    return(sds)
}  

mod<-function(x,m){
    t1<-floor(x/m)
    return(x-t1*m)
}

# function for aggregate

diff=function(x) {
    res = x[1]-x[2]
    return(res)
}

