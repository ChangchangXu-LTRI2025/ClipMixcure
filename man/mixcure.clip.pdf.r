###################################
#### generate mixcure.clip.pdf() ##
###################################

# formula = Surv(TSURV2, CENS == 0) ~ Her2;  beta = c(-1,0.1,-5,-0.1,0.1); pl=F
# 
# data = big.data[imputation.indicator == zz,]; 
# beta = beta[zz, ]; loglik = loglik[zz]; pos = pos; pl = F; b = z;


mixcure.clip.pdf <- function (formula, data, pl, iterlim = 200, pos, b, beta = NULL, loglik = NULL) 
{
  require(splines)
  require(survival)
  require(abind)

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
  
  
  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];
  
  design.matrix <- model.matrix(formula, data = design.matrix);
  
  # index ranges of coefficients of glm and cox models
  index.cure.v <- 1 : ncol(design.matrix);
  index.surv.v <- (ncol(design.matrix) + 1) : (2*length(index.cure.v))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.v)+1;
  
  n <- nrow(design.matrix)
  
  
  res <- matrix(0, 1, 3)
 
  init <- beta
  #init[pos] <- b
  
  loglik.mixture.profile <- function(p, survt, k, design.matrix1=design.matrix, design.matrix0=design.matrix, param.est, index.cure.var=index.cure.v,index.surv.var=index.surv.v, pl) {
    
    
    design.mtx.comb = cbind(design.matrix0,design.matrix1)
    ik = k-length(index.cure.var);
    
    #parameter and variable dep parameters;
    if (k > length(index.cure.v)) {
      theta = 1/(1+exp(-design.matrix[,index.cure.var]%*%as.matrix(p[index.cure.var])))
      eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var[-ik]]%*%as.matrix(p[-c(index.cure.var,index.gamma-1)])+design.mtx.comb[,k]*param.est)
      
      } else {
      theta = 1/(1+exp(-design.matrix[,index.cure.var[-k]]%*%as.matrix(p[-c(index.surv.var-1,index.gamma-1)])-design.mtx.comb[,k]*param.est))
      eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.var-1]))
      }
    
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta) #for LRT
    #kap= (1-eta)*(1-theta)*(theta + eta)  # for est,PLCI
    pi = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
    
    ####################################################################################################
    # Note: below constructs fisher info matrix; steps are divide into 4 blocks, 2 square blocks (A&D) #
    # on upper left and lower right, 2 identical transposed blocks (B) on upper right and lower left;  #
    # the idential B blocks are not identical in reduced models unless it's a global LRT, needs to be C#
    ####################################################################################################
    
    #calculate loglikelihood for the unpenalized;
    p.gamma <- p[index.gamma-1];  #use original shape parameter instead of exp();
    
    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) -
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
    
    if (pl==F) {
      loglik.part = loglikelihood
    } else {
      max.len = max(length(index.cure.var),length(index.surv.var))
      n.elema = max.len^2
      a.sub1 <- matrix(rep(0,n.elema), nrow = max.len)
      a.sub2 <- matrix(rep(0,n.elema), nrow = max.len)
      
      for (i in c(index.cure.var)) {
        for (j in c(index.cure.var)) {
          a.sub1[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*theta*(1-theta))[survt[, 2] == 1])
          a.sub2[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*kap)[survt[, 2] == 0])
        }
      }
      # if (k <= length(index.cure.var)) {info.a = (a.sub1 + a.sub2)[index.cure.var[-k],index.cure.var[-k]]} else
      info.a = (a.sub1 + a.sub2)[index.cure.var,index.cure.var]
      
      ##info matrix block B
      design.xt0 <- cbind(design.matrix0, log(survt[,1]))
      n.elemb <- max.len*(max.len+1)
      b.sub <- matrix(rep(0,n.elemb), nrow = max.len)
      
      for (i in c(index.cure.var)) {
        for (j in c(1:length(index.surv.var), max.len+1)) {
          b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*eps*(1-delta)*delta)[survt[, 2] == 0])
          #b.sub[i,j] <- -sum((design.matrix1[,i]*design.xt0[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])
          
        }
      }
      # if (k <= length(index.cure.var)) {info.b = b.sub[index.cure.var[-k],c(index.surv.var-max.len,index.gamma-max.len)]} else
      # {info.b = b.sub[index.cure.var,c(index.surv.var[-ik]-max.len,index.gamma-max.len)]}
      info.b = b.sub[index.cure.var,c(index.surv.var-max.len,index.gamma-max.len)]
      
      ###info matrix block d
      design.xt1 <- cbind(design.matrix1, log(survt[,1]))
      
      n.elemd <- (max.len+1)^2
      d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
      d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
      
      for (i in c(index.surv.var-max.len, max.len +1)) {
        for (j in c(index.surv.var-max.len, max.len +1)) {
          d.sub1[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*eps)[survt[, 2] == 1])
          d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
          #d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta^2))[survt[, 2] == 0])
        }
      }
      d.sub = d.sub1 + d.sub2 +
        matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-1]^2)),
               nrow = (max.len + 1))
      
      # if (k <= length(index.cure.var))
      #   {info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]} else
      #   {info.d = d.sub[c(index.surv.var[-ik]-max.len,index.gamma-max.len),c(index.surv.var[-ik]-max.len,index.gamma-max.len)]}
      
      info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]
      
      info.d.inv = mat.inv(info.d)
      
      fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
      #hessian.mat = -fisher.info
      
      # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
      info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)
      
      #determinant of hessian matrix;
      det.info = det(info.set0)*det(info.d)
      #det.info = matrix.det(fisher.info)
      
      loglik.part = loglikelihood - 0.5*log(det.info)
      } 
  
    return(loglik.part)
  }
  
  maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = init[-pos], survt=survt,
                          param.est = b, k = pos,
                          pl = pl, iterlim = iterlim, hessian=TRUE)
  red.loglik = -maximizer.temp1$minimum
  
 
  res[1, 1] <- b  #submitted parameter estimate, ie, lowerbound.lo      
  res[1, 2] <- 2 * abs(loglik - red.loglik)  #chisq value by 2*LR
  ##########################################################################
  res[1, 3] <- 1 - (1 - pchisq(res[, 2], 1))/2   #b corresponding pdf value as signed sqrt of chisq distributed LR
  ##########################################################################
  res[1, 3][res[, 1] < beta[pos]] <- 1 - res[, 3][res[, 1] < beta[pos]]
  
  results <- list(beta = res[1, 1], chisq = res[1, 2], pdf = res[1, 3])
  results
}