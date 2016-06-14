### function to create folds in the case in which
### correlated observation units are given
createRandomRespFolds <- function(ranVar, # the factor variable which corresponds to the observation units
                                  folds = length(levels(ranVar)), # number of folds
                                  sLength = NULL)
{
  
  if(!is.factor(ranVar)) stop("ranVar must be a factor.")
  
  if(length(levels(ranVar)) %% folds != 0) # if given folds are not a multiple of factor levels
  {
    
    # calculate folds
    folds <- ceiling(length(levels(ranVar)) / round(length(levels(ranVar)) / folds))
    warning(paste0("Number of levels of ranVar is not a multiple of folds. Settings folds to ",folds))
    
  }
  
  if(is.null(sLength)) 
    stop("For cvrisk, please supply #observations per trajectory 'sLength'.")
  
  leaveThisNrPerPerFoldOut <- ceiling(length(levels(ranVar))/folds)
  
  # if(forValidateFDboost){
  #   
  #   ppMat <- matrix(rep(ranVar,folds),ncol=folds)
  #   ppMat <- sapply(0:(folds-1),function(i)
  #     as.numeric(sapply(ppMat[,i+1],function(x)!x%in%levels(ranVar)[nrPPperFold*i+1:nrPPperFold])))
  #   
  # }else{
  
  ppVec <- rep(ranVar, sLength)
  # create fold matrix
  ppMat <- sapply(1:folds, function(i){ 
    
    kickThis <- ((i-1)*leaveThisNrPerPerFoldOut)+1:leaveThisNrPerPerFoldOut
    as.numeric(!ppVec%in%levels(ranVar)[kickThis])
    
  })
    
  # }
  
  return(ppMat)
    
  
}