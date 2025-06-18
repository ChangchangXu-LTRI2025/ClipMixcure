################################
### Function for mixcure model##
#### penalized loglikelihoods ##
##########################################################
#### Last modified Dec31 2018 for 3 or more x variables ##
##########################################################

mixcure.penal.mi <- function(formula, data, init, pl, iterlim = 200) { 
require(splines)
require(survival)
require(abind)
require(R.utils)

  cal <- match.call()
  #########################################################################################
  mat.inv <- function(matx) {
    
    detm = det(matx)
    #2x2 matrix inverse;
    if (ncol(matx) == 2) {
      inv.matx = (1/detm) * matrix(c(matx[2,2],-matx[2,1],-matx[1,2],matx[1,1]), nrow = 2)
    }
    
    else {
      #For any n>2 dimension square matrix;  
      adjug.matx <- matrix(rep(0, ncol(matx)^2), nrow = nrow(matx))
      for (i in 1:nrow(matx)) {
        for (j in 1:ncol(matx)) {
          adjug.matx[i,j] <- (-1)^(i+j)*det(matx[-i,][,-j])
        }
      }
      inv.matx <- t(adjug.matx/detm)
    } 
    
    return(inv.matx)
  }
  
  #########################################################################################
  
  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];
  
  design.matrix <- model.matrix(formula, data = design.matrix);
  
  # index ranges of coefficients of glm and cox models
  index.cure.v <- 1 : ncol(design.matrix); 
  index.surv.v <- (ncol(design.matrix) + 1) : (2*length(index.cure.v))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.v)+1;
  
  #samp.s <- nrow(design.matrix)
  
  
  ####################################################
  ## nonlinear minimization algoritm to solve       ##
  ## penalized mixture cure loglikelihood functions ##
  ####################################################      
  
  loglik.mixture <- function(p, survt, design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {
    
   ####  parameter and variable dep parameters;
   #####
    theta = 1/(1+exp(-design.matrix%*%p[index.cure.var]))
    eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.var])
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap = theta*(1-theta)*(1-eta)+(1-theta)^2*eta*(1-eta)  # for est and PLCI
   # kap= (1-eta)*(1-theta)*(theta + eta)    # exp for est and PLCI 
     pi = exp(eps)*eps*eta^2
    # lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    # phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi

    #calculate loglikelihood for the unpenalized;
    cure.par <- p[1 : ncol(design.matrix) ];
    surv.par <- p[ (ncol(design.matrix) + 1) : (2*length(cure.par)) ];
    p.gamma <- p[ 2*length(cure.par) + 1 ];  #use original shape parameter instead of exp();
    
    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer; 
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) - 
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
    
 if (pl==F)    {
   loglik = loglikelihood
 } else {
  ####calculate inverse of info matrix by block matrix;
 
    n.elema = length(index.cure.var)^2
    a.sub1 <- matrix(rep(0,n.elema), nrow = length(index.cure.var))
    a.sub2 <- matrix(rep(0,n.elema), nrow = length(index.cure.var))
    
     for (i in c(index.cure.var)) {
      for (j in c(index.cure.var)) {
        a.sub1[i,j] <- sum((design.matrix[,i]*design.matrix[,j]*theta*(1-theta))[survt[, 2] == 1])
        a.sub2[i,j] <- sum((design.matrix[,i]*design.matrix[,j]*kap)[survt[, 2] == 0])
        }
    }
    info.a = a.sub1 + a.sub2
    
    design.xt <- cbind(design.matrix, log(survt[,1]))
    n.elemb <- length(index.cure.var)*(length(index.cure.var)+1)
    b.sub <- matrix(rep(0,n.elemb), nrow = length(index.surv.var))
    
    for (i in c(index.cure.var)) {
      for (j in c(index.cure.var,length(index.surv.var)+1)) {
        #b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*theta*(1-theta)*pi)[survt[, 2] == 0]) #for est
        #b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])  #for LRT
        b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*eps*(1-delta)*delta)[survt[, 2] == 0]) #alternative expression for est
        
              }
    }
    info.b = b.sub  #Upper right block of fisher.info;
    

    n.elemd <- (length(index.surv.var)+1)^2
    d.sub1 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
    d.sub2 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
    
    for (i in c(index.cure.var,length(index.surv.var)+1)) {
      for (j in c(index.cure.var,length(index.surv.var)+1)) {
        d.sub1[i,j] <- sum((design.xt[,i]*design.xt[,j]*eps)[survt[, 2] == 1])
        #d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta^2))[survt[, 2] == 0]) #for est, PLCI
        #d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*(delta*(1-delta))))[survt[, 2] == 0]) #for LRT, same as below
        d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0]) #for est, PLCI
        
        }
    }
    info.d = d.sub1 + d.sub2 + 
           matrix(c(rep(0, (n.elemd-1)),sum(survt[, 2] == 1)/(p[index.gamma]^2)),nrow = (length(index.surv.var)+1))
    
    
    info.d.inv = mat.inv(info.d)

#    fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
    #hessian.mat = -fisher.info
    
    # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
   info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)

        #determinant of hessian matrix;
    det.info = det(info.set0)*det(info.d)
 #   det.info = matrix.det(fisher.info)
      
      
      loglik = loglikelihood - 0.5*log(det.info)
 
 }
    
    #loglik = loglikelihood
    return(loglik)
    
  }
  
  ######END of loglik.mixture####################################
  
  
  # Parameter estimation under Ha (non-restricted likelihood)
  # maximize penalized or unpenalized loglikelihood by nlm; 
  maximizer0 <- nlm(
    f = loglik.mixture, p = init, survt=survt, design.matrix=design.matrix, 
    pl = pl, 
    iterlim = iterlim, hessian=TRUE);

  hessmat <- maximizer0$hessian
  if (det(hessmat) < 1e-05) 
    diag(hessmat) <- diag(hessmat) + 1e-06
  
var.mat <- solve(hessmat)

alpha.hat <- maximizer0$estimate[index.gamma];

loglik <- -maximizer0$minimum  #in loglik function loglik was calculated as minus of actual loglik value

###### Calculate restricted likelihood #######
# loglik.mixture.part <- function(p, survt,part.cure=F, design.matrix1, design.matrix0, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {  #design.matrix1-surv, design.matrix0-cure
#   
#   design.mtx.comb = cbind(design.matrix0,design.matrix1)
#   
#   
#   #parameter and variable dep parameters;
#  
#   if (k > length(index.cure.v)) {
#     if (ncol(design.matrix) <3)   {
#       theta = 1/(1+exp(-design.matrix[,index.cure.var]%*%as.matrix(p[index.cure.var])))
#       eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]*as.matrix(p[-c(index.cure.var,index.gamma-1)]))
#       
#     } else {
#     theta = 1/(1+exp(-design.matrix[,index.cure.var]%*%as.matrix(p[index.cure.var])))
#     eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[-c(index.cure.var,index.gamma-1)]))
#            }
#                                 } else {
#     if (ncol(design.matrix) <3) {                              
#     theta = 1/(1+exp(-design.matrix[,index.cure.var]*as.matrix(p[-c(index.surv.var-1,index.gamma-1)])))
#     eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.var-1]))
#     }
#                                        }
#   # theta = 1/(1+exp(-design.matrix0%*%as.matrix(p[index.cure.var])))
#   # eps = survt[,1]^(p[index.gamma])*exp(design.matrix1%*%as.matrix(p[index.surv.var]))
#   
#   eta = 1/((exp(eps)-1)*theta+1)
#   delta = 1/(theta/(1-theta)*exp(eps)+1)
#   kap = theta*(1-theta)*(1-eta)+(1-theta)^2*eta*(1-eta) #for LRT
#   #kap= (1-eta)*(1-theta)*(theta + eta)  #for est, PLCI
#   pi = exp(eps)*eps*eta^2
#   #lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
#   #phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
#   
#   ####################################################################################################
#   # Note: below constructs fisher info matrix; steps are divide into 4 blocks, 2 square blocks (A&D) #
#   # on upper left and lower right, 2 identical transposed blocks (B) on upper right and lower left;  #
#   # the idential B blocks are not identical in reduced models unless it's a global LRT, needs to be C#
#   ####################################################################################################
#   
#   #calculate loglikelihood for the unpenalized;
#   cure.par <- p[1 : ncol(design.matrix) ];
#    p.gamma <- p[index.gamma-1];  #use original shape parameter instead of exp();
#   
#   # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
#   loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
#                            + log(eps)-eps )[survt[, 2] == 1] ) -
#     sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
#   
#   if (pl == F) {loglik.part = loglikelihood} else {
#     
#     
#     max.len = max(length(index.cure.var),length(index.surv.var))
#     n.elema = max.len^2
#     a.sub1 <- matrix(rep(0,n.elema), nrow = max.len)
#     a.sub2 <- matrix(rep(0,n.elema), nrow = max.len)
#     
#     for (i in c(index.cure.var)) {
#       for (j in c(index.cure.var)) {
#         a.sub1[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*theta*(1-theta))[survt[, 2] == 1])
#         a.sub2[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*kap)[survt[, 2] == 0])
#       }
#     }
#     info.a = (a.sub1 + a.sub2)[index.cure.var,index.cure.var]
#     
#     ##info matrix block B
#     design.xt0 <- cbind(design.matrix0, log(survt[,1]))
#     n.elemb <- max.len*(max.len+1)
#     b.sub <- matrix(rep(0,n.elemb), nrow = max.len)
#     
#     for (i in c(index.cure.var)) {
#       for (j in c((index.surv.var-max.len), max.len+1)) {
#         #b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])  #equivalent to expression below
#         b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*eps*(1-delta)*delta)[survt[, 2] == 0])  #equivalent to expression below
#         #b.sub[i,j] <- -sum((design.matrix1[,i]*design.xt0[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])
#       }
#     }
#     info.b = b.sub[index.cure.var,c(index.surv.var-max.len,index.gamma-max.len)]
#     
#     design.xt1 <- cbind(design.matrix1, log(survt[,1]))
#     
#     n.elemd <- (max.len+1)^2
#     d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
#     d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
#     
#     for (i in c(index.surv.var-max.len, max.len +1)) {
#       for (j in c(index.surv.var-max.len, max.len +1)) {
#         d.sub1[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*eps)[survt[, 2] == 1])
#         #d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
#         #d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta^2))[survt[, 2] == 0])
#         d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*(delta*(1-delta))))[survt[, 2] == 0])
#         
#       }
#     }
#     d.sub = d.sub1 + d.sub2 +
#       matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-1]^2)),
#              nrow = (max.len + 1))
#     
#     info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]
#     
#     info.d.inv = mat.inv(info.d)
#     
#     # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
#     if (part.cure == T & length(cure.par)==2) {info.set0 = info.a-t(info.b)%*%info.d.inv%*%info.b} else
#     {info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)}
#     
#     det.info = det(info.set0)*det(info.d)
#     
#     loglik.part = loglikelihood - 0.5*log(det.info)
#     
#   }
#   
#   return(loglik.part)
# }


loglik.mixture.part <- function(p, survt, design.matrix1, design.matrix0,
                                index.cure.var=index.cure.v,
                                index.surv.var=index.surv.v, pl) {  #design.matrix1-surv, design.matrix0-cure
  
  design.mtx.comb = cbind(design.matrix0,design.matrix1)
  
  #parameter and variable dep parameters;
  if (k <= length(index.cure.v)) {
    theta = 1/(1+exp(-design.mtx.comb[,index.cure.var]%*%as.matrix(p[index.cure.v[-length(index.cure.v)]])))
    eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.var-1]))
  } else {
    theta = 1/(1+exp(-design.mtx.comb[,index.cure.var]%*%as.matrix(p[index.cure.var])))
    eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.v[-length(index.surv.v)]]))
  }
  eta = 1/((exp(eps)-1)*theta+1)
  delta = 1/(theta/(1-theta)*exp(eps)+1)
  kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta)
 # kap= (1-eta)*(1-theta)*(theta + eta)    # exp for est and PLCI 
  pi = exp(eps)*eps*eta^2
  lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
  phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
  
  
  max.len = max(length(index.cure.var),length(index.surv.var))
  n.elema = max.len^2
  a.sub1 <- matrix(rep(0,n.elema), nrow = max.len)
  a.sub2 <- matrix(rep(0,n.elema), nrow = max.len)
  
  for (i in c(index.cure.v)) {
    for (j in c(index.cure.v)) {
      a.sub1[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*theta*(1-theta))[survt[, 2] == 1])
      a.sub2[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*kap)[survt[, 2] == 0])
    }
  }
  info.a = (a.sub1 + a.sub2)
  
  
  ##info matrix block B
  design.xt0 <- cbind(design.matrix0, log(survt[,1]))
  n.elemb <- max.len*(max.len+1)
  b.sub <- matrix(rep(0,n.elemb), nrow = max.len)
  
  for (i in c(index.cure.v)) {
    for (j in c(1:length(index.surv.v), max.len+1)) {
      b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])
    }
  }
  info.b = b.sub
  
  design.xt1 <- cbind(design.matrix1, log(survt[,1]))
  
  n.elemd <- (max.len+1)^2
  d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
  d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
  
  for (i in c(index.cure.v, max.len +1)) {
    for (j in c(index.cure.v, max.len +1)) {
      d.sub1[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*eps)[survt[, 2] == 1])
      d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
    }
  }
  d.sub = d.sub1 + d.sub2 +
    matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-1]^2)),
           nrow = (max.len + 1))
  
  info.d = d.sub
  
  
  info.d.inv = mat.inv(info.d)
  
  #fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
  #hessian.mat = -fisher.info
  
  # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
  info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)
  
  #determinant of hessian matrix;
  det.info = det(info.set0)*det(info.d)
  
  #calculate loglikelihood for the unpenalized;
  cure.par <- p[index.cure.var];
  surv.par <- p[index.surv.var];
  p.gamma <- p[index.gamma-1];  #use original shape parameter instead of exp();
  
  # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
  loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                           +log(eps)-eps )[survt[, 2] == 1] ) -
    sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
  
  
  if (pl == FALSE)
  {
    loglik.part = loglikelihood
  }
  else if (pl == TRUE)
  {
    loglik.part = loglikelihood - 0.5*log(det.info)
  }
  
  return(loglik.part)
}


# dim.v <- ncol(design.matrix)
# ll.cure <- rep(0, dim.v)
# est.cure <- matrix(0, nrow = dim.v, ncol = (2*dim.v + 1))
# 
# for (k in index.cure.v[-1]) {
#   maximizer <- nlm(
#     f = loglik.mixture.part,
#     p = init[-k], part.cure=T,
#     survt = survt, design.matrix0 = design.matrix,
#     design.matrix1=design.matrix,
#     index.cure.var=index.cure.v[-k],
#     pl=pl,
#     iterlim = iterlim, hessian=F
#   );
#   
#   ll.cure[k] = -maximizer$minimum;
#   est.cure[k,] = insert(maximizer$estimate, ats=k, values=0)
# }
# 
# 
# ll.surv <- rep(0,ncol(design.matrix))
# est.surv <- matrix(0, nrow = dim.v, ncol = (2*dim.v + 1))
# 
# for (k in index.surv.v[-1]) {
#   # mle under the reduced (null) model for surv parameter;
#   is=k-length(index.cure.v)
#   # if (k==11|k==12) {init=c(1,rep(0,5),-5,rep(0,5),0.1)}
#   maximizer <- nlm(
#     f = loglik.mixture.part, p =  init[-k],   #p=c(-3,0.1,-1,-0.1,1),c(-2,0.1,0,-0.1,1) for univar.null.lr.pl;
#     survt = survt, design.matrix1 = design.matrix,
#     design.matrix0=design.matrix,
#     index.surv.var=index.surv.v[-is],
#     pl=pl,
#     iterlim = iterlim, hessian=FALSE
#   );
#   
#   ll.surv[is] = -maximizer$minimum;
#   est.surv[is,] = insert(maximizer$estimate, ats=k, values=0)
# }

dim.v <- ncol(design.matrix)
ll.cure <- rep(0,dim.v)
llr.cure <- rep(0,dim.v)
pval.cure <- rep(0,dim.v)
est.cure <- matrix(0, nrow = dim.v, ncol = (2*dim.v + 1))
# index.cure.v[-1] for no intercept calculation
for (k in index.cure.v) {
  maximizer <- nlm(
    f = loglik.mixture.part, p = init[-k],
    survt = survt, design.matrix0 = design.matrix,
    design.matrix1=design.matrix,
    index.cure.var=index.cure.v[-k], pl=pl,
    iterlim = iterlim, hessian=T
  );
  loglik.part = -maximizer$minimum;
  dif.ll = -2*(loglik.part-loglik);  #loglik is ll under Ha;
  pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
  ll.cure[k]<- loglik.part
  llr.cure[k]<- dif.ll
  pval.cure[k]<- pval
  est.cure[k,] = insert(maximizer$estimate, ats=k, values=0)
  if (det(maximizer$hessian) < 1e-05)
    diag(maximizer$hessian) <- diag(maximizer$hessian) + 1e-06
  
}



### loglikelihood calculation for each surv part variable;

ll.surv <- rep(0,ncol(design.matrix))
llr.surv <- rep(0,ncol(design.matrix))
pval.surv <- rep(0,ncol(design.matrix))
est.surv <- matrix(0, nrow = dim.v, ncol = (2*dim.v + 1))

for (k in index.surv.v) {
  is=k-length(index.cure.v)
  maximizer <- nlm(
    f = loglik.mixture.part, p = init[-k],
    survt = survt, design.matrix1 = design.matrix,
    design.matrix0=design.matrix,
    index.surv.var=index.surv.v[-is], pl=pl,
    iterlim = iterlim, hessian=FALSE
  );
  
  loglik.part = -maximizer$minimum;
  dif.ll = -2*(loglik.part-loglik);
  pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
  ll.surv[is]<- loglik.part
  llr.surv[is]<-dif.ll
  pval.surv[is]<-pval
  est.surv[is,] = insert(maximizer$estimate, ats=k, values=0)
}



######### Wald statistic inference ##########
  # confidence intervals for estimated coefficients
  z.score <- maximizer0$estimate / sqrt(diag(var.mat));
  
coef.table.cure <- cbind(
  'coef'        = maximizer0$estimate[index.cure.v],
  'exp(coef)'   = exp(maximizer0$estimate[index.cure.v]),
  'se(coef)'    = sqrt(diag(var.mat)[index.cure.v]),
  'z'           = z.score[index.cure.v],
  'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.cure.v]))),
  'LCI.95%' = maximizer0$estimate[index.cure.v] - 1.96 * sqrt(diag(var.mat)[index.cure.v]),
  'UCI.95%' = maximizer0$estimate[index.cure.v] + 1.96 * sqrt(diag(var.mat)[index.cure.v])
);
rownames(coef.table.cure) <- colnames(design.matrix);

coef.table.surv <- cbind(
  'coef'        = maximizer0$estimate[index.surv.v],
  'exp(coef)'   = exp(maximizer0$estimate[index.surv.v]),
  'se(coef)'    = sqrt(diag(var.mat)[index.surv.v]),
  'z'           = z.score[index.surv.v],
  'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.surv.v]))),
  'LCI.95%' = maximizer0$estimate[index.surv.v] - 1.96 * sqrt(diag(var.mat)[index.surv.v]),
  'UCI.95%' = maximizer0$estimate[index.surv.v] + 1.96 * sqrt(diag(var.mat)[index.surv.v])
);
rownames(coef.table.surv) <- colnames(design.matrix);

coef.table.alpha <- cbind(
  'coef'     = alpha.hat,
  'se(coef)' = sqrt(diag(var.mat)[index.gamma]),
  'z'        = z.score[index.gamma],
  'Pr(>|z|)' = 2 * (1 - pnorm(abs(z.score[index.gamma]))),
  'LCI.95%'  = maximizer0$estimate[index.gamma] - 1.96 * sqrt(diag(var.mat)[index.gamma]),
  'UCI.95%'  = maximizer0$estimate[index.gamma] + 1.96 * sqrt(diag(var.mat)[index.gamma]),
  'loglik' = -maximizer0$minimum 
);
rownames(coef.table.alpha) <- 'alpha';

#######################################
## Output tables from either method; ##
#######################################

colnames(var.mat) <- c(
  paste('cure.', colnames(design.matrix)), 
  paste('surv.', colnames(design.matrix)),
  'alpha'
);
rownames(var.mat) <- colnames(var.mat);

out <- list(call=cal,
  coefficients = c(coef.table.cure[,1], coef.table.surv[,1], coef.table.alpha[,1]),
 ci.lower = c(coef.table.cure[,6], coef.table.surv[,6], coef.table.alpha[,5]),
 ci.upper = c(coef.table.cure[,7], coef.table.surv[,7], coef.table.alpha[,6]),
 loglik = list(
   cure = ll.cure[-1],
   surv = ll.surv[-1],
   cure.est = est.cure,
   surv.est = est.surv,
   full = -maximizer0$minimum
 ),
 terms = colnames(var.mat),
 alpha.level = 0.05,
 cov = var.mat
);
class(out) <- c('mixcure', 'list');

return(out);

}


#### print.mixcure #############################################################
# DESCRIPTION
#   To print a mixcure object.
# INPUT
#   object : a mixcure object, which is an outcome of function mixcure.
#   digits : number of digits for printing, passed to print.default.
#   ...    : other parameters passed to print.default.
# OUTPUT
#   NULL.   

# print.mixcure <- function(object, digits = 3, ...) {
#   sep.line.cure   <- paste(c(rep('-', 37),  ' CURE ' , rep('-', 37)), collapse = '');
#   sep.line.surv   <- paste(c(rep('-', 37),  ' SURVIVAL ' , rep('-', 37)), collapse = '');
#   sep.line.alpha <- paste(c(rep('-', 36), ' ALPHA ', rep('-', 36)), collapse = '');
#   
#   message(sep.line.cure);
#   print.default(object$coefficients$cure,   digits = digits, ...);
#   
#   message(sep.line.surv);
#   print.default(object$coefficients$surv,   digits = digits, ...);
#   
#   message(sep.line.alpha);
#   print.default(object$coefficients$alpha, digits = digits, ...);
#   
#   return(NULL);
# };


#### coef.mixcure ##############################################################
# coef.mixcure <- function(object) {
#   coefs <- c(
#     object$coefficients$cure[, 'coef'], 
#     object$coefficients$surv[, 'coef'], 
#     object$coefficients$alpha[, 'coef']
#   );
#   names(coefs) <- c( paste('cure.', rownames(object$coefficients$cure)), 
#                      paste('surv.', rownames(object$coefficients$surv)),
#                      rownames(object$coefficients$alpha) );
#   
#   return(coefs);
# }
# 
# 
# vcov.mixcure <- function(object) {
#   return(object$cov);
# }


