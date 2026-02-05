###################################
#### generate mixcure.clip.pdf() ##
###################################

# formula = Surv(TSURV2, CENS == 0) ~ Her2;  beta = c(-1,0.1,-5,-0.1,0.1); pl=F
#
# data = big.data[imputation.indicator == zz,];
# beta = beta[zz, ]; loglik = loglik[zz]; pos = pos; pl = F; b = z;

#' Compute PDF/CDF pieces used by the clipped-likelihood algorithm
#'
#' @noRd
#' @export
mixcure.clip.pdf <- function (formula, data, pl, iterlim = 200, pos, b, beta = NULL, loglik = NULL)
{
  require(splines)
  require(survival)
  require(abind)

  #########################################################################################

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

  loglik.mixture.profile <- function(
    p, survt, k,
    design.matrix1 = design.matrix, design.matrix0 = design.matrix,
    param.est,
    index.cure.var = index.cure.v, index.surv.var = index.surv.v,
    pl
  ) {
    t      <- survt[, 1]
    status <- survt[, 2]
    event  <- (status == 1L)
    cens   <- !event
    logt   <- log(t)

    design.mtx.comb <- cbind(design.matrix0, design.matrix1)
    ik <- k - length(index.cure.var)

    # index.gamma is assumed defined in caller environment
    gamma <- p[index.gamma - 1]
    if (!is.finite(gamma) || gamma <= 0) return(Inf)

    # ---- theta & eps ----
    if (k > length(index.cure.v)) {
      lp_cure <- drop(design.matrix0[, index.cure.var, drop = FALSE] %*%
                        as.matrix(p[index.cure.var]))
      theta <- plogis(lp_cure)

      lp_surv <- drop(design.mtx.comb[, index.surv.var[-ik], drop = FALSE] %*%
                        as.matrix(p[-c(index.cure.var, index.gamma - 1)]) +
                        design.mtx.comb[, k] * param.est)
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
    } else {
      lp_cure <- drop(design.matrix0[, index.cure.var[-k], drop = FALSE] %*%
                        as.matrix(p[-c(index.surv.var - 1, index.gamma - 1)]) -
                        design.mtx.comb[, k] * param.est)
      theta <- plogis(lp_cure)

      lp_surv <- drop(design.mtx.comb[, index.surv.var, drop = FALSE] %*%
                        as.matrix(p[index.surv.var - 1]))
      logeps <- gamma * logt + lp_surv
      eps <- exp(logeps)
    }

    # Stable rewrite: D = theta + (1-theta)*exp(-eps)
    emeps <- exp(-eps)
    D <- theta + (1 - theta) * emeps

    eta   <- emeps / D
    delta <- (1 - theta) * eta

    # IMPORTANT: keep your LRT kap
    kap <- theta * (1 - theta) * (1 - eta) - (1 - theta)^2 * eta * (1 - eta)

    # pi = exp(eps)*eps*eta^2 (rewrite safely)
    pi <- eps * emeps / (D * D)

    # ---- negative log-likelihood ----
    event_term <- log1p(-theta) + log(gamma) - logt + logeps - eps
    loglikelihood <- -sum(event_term[event]) - sum(log(D[cens]))

    if (isFALSE(pl)) return(loglikelihood)

    # =================================
    # Fast Fisher-block construction
    # =================================
    max.len <- max(length(index.cure.var), length(index.surv.var))

    X0 <- as.matrix(design.matrix0)
    X1 <- as.matrix(design.matrix1)

    # Use first max.len columns to match your i/j indexing (1..max.len)
    X0A <- X0[, 1:max.len, drop = FALSE]
    X1B <- X1[, 1:max.len, drop = FALSE]

    e <- as.numeric(event)
    c <- 1 - e

    # --- Block A: event weight theta(1-theta), cens weight kap
    wA <- e * (theta * (1 - theta)) + c * kap
    A_full <- crossprod(X0A, X0A * wA)
    info.a <- A_full[index.cure.var, index.cure.var, drop = FALSE]

    # --- Block B: -sum( design.matrix1[,i] * design.xt0[,j] * eps*(1-delta)*delta ) over cens
    Xt0 <- cbind(X0A, logt)  # (n x (max.len+1))
    wB <- c * (eps * (1 - delta) * delta)
    B_full <- -crossprod(X1B, Xt0 * wB)  # (max.len x (max.len+1))

    cols_B <- c(index.surv.var - max.len, index.gamma - max.len)
    info.b <- B_full[index.cure.var, cols_B, drop = FALSE]

    # --- Block D: event weight eps, cens weight (eps*delta - eps^2*delta + eps^2*delta^2)
    Xt1 <- cbind(X1B, logt)  # (n x (max.len+1))
    wd2 <- eps * delta - (eps^2) * delta + (eps^2) * (delta^2)
    wD <- e * eps + c * wd2

    D_full <- crossprod(Xt1, Xt1 * wD)
    D_full[max.len + 1, max.len + 1] <- D_full[max.len + 1, max.len + 1] + sum(event) / (gamma * gamma)

    rowscols_D <- c(index.surv.var - max.len, index.gamma - max.len)
    info.d <- D_full[rowscols_D, rowscols_D, drop = FALSE]

    # Schur complement without explicit inverse
    tmp <- solve(info.d, t(info.b))
    info.set0 <- info.a - info.b %*% tmp

    # Stable log(det(info.set0)*det(info.d))
    detS <- determinant(info.set0, logarithm = TRUE)
    detD <- determinant(info.d, logarithm = TRUE)
    if (detS$sign * detD$sign <= 0) return(Inf)

    logdet <- as.numeric(detS$modulus) + as.numeric(detD$modulus)
    loglikelihood - 0.5 * logdet
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
