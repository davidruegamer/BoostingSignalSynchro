############################################################################
### Function to generate random historical effects 
### mostly based on simulation code by Fabian Scheipl and Sarah Brockhaus 
### See e.g. online appendix of 
### Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
### The functional linear array model. Statistical Modelling, 15(3), 279-300.
############################################################################

genRanHist <- function(s, t, nrLevels = 12, # number of random effect levels
                       coef=NULL, seed=NULL, df=10, # degrees of freedom for splines
                       pen=c(2,2), # difference penalties in s and t direction
                       lambda=c(1,1) # smoothing parameter for s and t direction
)
{
  
  if(!is.null(seed)) set.seed(seed)
  require(splines)
  Bs <- bs(s, df=df, intercept = TRUE)
  Bt <- bs(t, df=df, intercept = TRUE)
  
  # Recursion for difference operator matrix
  makeDiffOp <- function(degree, dim){
    if(degree==0){
      return(diag(dim))  
    } else {
      return(diff(makeDiffOp(degree-1, dim)))
    }    
  }
  Pt <- lambda[1] * kronecker(crossprod(makeDiffOp(pen[1], df)), diag(df))
  Ps <- lambda[2] * kronecker(diag(df), crossprod(makeDiffOp(pen[2], df)))            
  P <- .1*diag(df^2) + Pt + Ps
  
  if(is.null(coef)){
    
    coef <- vector("list", length = nrLevels)
    
    for(i in 1:length(coef)){
      
      coef[[i]] <- matrix(solve(P, rnorm(df^2)), df, df) 
      
    }
    
    meanCoef <- Reduce("+",coef)/length(coef)
    coef <- lapply(coef,function(singleCoef)return(singleCoef - meanCoef))
    
  }
  
  if(!is.list(coef)) coef <- list(coef)
  ret <- vector("list", length = length(coef))
  
  for(i in 1:length(coef)){
    
    ret[[i]] <- Bs%*%coef[[i]]%*%t(Bt)
    
    rownames(ret[[i]]) <- round(s, 2)
    colnames(ret[[i]]) <- round(t, 2)
    
    for(k in 1:length(s)){
      for(j in 1:length(t)){
        if(s[k] > t[j]){
          ret[[i]][k, j] <- 0
        } 
      }
    }
    
    attr(ret[[i]], "coef") <- coef[[i]]
    
  }
  
  return(ret)
  
}