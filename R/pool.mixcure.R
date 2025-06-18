
pool.mixcure<-function (object)
{
  require(mice)
  call <- match.call()
  if (!is.list(object))
    stop("Argument 'object' not a list", call. = FALSE)
  #object <- as.mira(object)
  m <- length(object$analyses)
  fa <- getfit(object, 1)
  if (m == 1) {
    warning("Number of multiple imputations m = 1. No pooling done.")
    return(fa)
  } else {
  var.length<-length(object$analyses[[1]]$coefficients)
  N<-var.length*2*m
  data.array <- array(rep(0,N), dim=c(var.length,2,m))
  for (i in 1:m) {
  data.array[,1,i] <- c(object$analyses[[i]]$coefficients)    #coefficients matrix is coef;
 # data.array[,2,i] <- sqrt(diag(solve(object$analyses[[i]]$cov)))    # se;
  data.array[,2,i] <- diag(object$analyses[[i]]$cov)  # var;

  }
  out.table          <- matrix(rep(0,var.length*2),ncol=2)
  out.table[,c(1,2)] <- apply(data.array,c(1,2),mean)

  results            <- matrix(rep(0,var.length*7),ncol=7)
  results[,1]        <- out.table[,1]   #coef
  results[,2]        <- exp(out.table[,1])   #exp(coef)
  results[,3]        <- sqrt(out.table[,2]+ (apply(data.array,c(1,2),sd)[,1])^2*(1+1/m)) #se
  results[,4]        <- out.table[,1]/results[,3]   #z
  results[,5]        <- pnorm(-abs(results[,4]))   #p-value
  results[,6]        <- results[,1]-1.96*results[,3]   #LCI95
  results[,7]        <- results[,1]+1.96*results[,3]   #UCI95

  rownames(results)  <- names(object$analyses[[1]]$coefficients)
  colnames(results)  <- c("coef","exp(coef)","se","z","p(z)","LCI95%","UCI95%")

  class(results) <- c("mipo", "data.frame")

  return(results)
  }
  }
