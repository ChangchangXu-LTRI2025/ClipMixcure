###############################################################
###############################################################
####                                                       ####
#### FOR GENERATING CLIP BOUNDS FOR MI DATA UNDER MC MODEL ####
####                                                       ####
###############################################################
###############################################################

##feed mira object to obj, and specify TRUE or FALSE for pl;

#' Generate CLIP bounds for multiply-imputed mixture cure model fits
#'
#' This is a legacy entry point used internally by the ClipMixcure engine.
#' The meta-layer uses \code{clipmixcure_fit()} as the stable front-end.
#'
#' @noRd
#' @export
clip.mixcure <- function (obj = NULL, variable = NULL, pl =F , ci.level = c(0.025, 0.975),
                           pvalue = TRUE, bound.lo = NULL, bound.up = NULL, iternum=20) {

  ## Prefer mids object attached by clipmixcure_fit()
  data.imp <- attr(obj, "mids_data", exact = TRUE)

  ## Legacy fallback: object created in user's global environment
  if (!inherits(data.imp, "mids")) {
    data.imp <- eval(obj$call$data, envir = parent.frame())
  }

  if (!inherits(data.imp, "mids")) {
    stop("clip.mixcure(): could not recover a 'mids' object. ",
         "Provide a mira object created by with(mids, ...) or attach it via attr(obj,'mids_data').",
         call. = FALSE)
  }

  ## Robust number of imputations
  nimp <- data.imp$m
  if (is.null(nimp) || length(nimp) == 0L) {
    nimp <- length(obj$analyses)
  }

  ## Use mice::complete explicitly
  data <- lapply(seq_len(nimp), function(x) mice::complete(data.imp, action = x))

  fits <- obj$analyses
  formula <- as.formula(obj$call$expr[[2]])

  if (is.null(variable))
    variable <- fits[[1]]$terms
  nvar <- length(variable)

  ### generate central function inner.CLIP within CLIP.MIXCURE
  ############################################################

  inner.CLIP <- function(myvar) {
    variable <- myvar
    imputations <- nimp

    # res<- matrix(0, length(grid), 3)

    nperimp <- unlist(lapply(1:imputations, function(x) nrow(data[[x]])))
    variable.names <- head(fits[[1]]$terms,-1)
    imputation.indicator <- unlist(sapply(1:imputations, function(x) rep(x, nperimp[x]))[TRUE])   # transform of listed 20 imputed datasets to a vector of imp indicators
    mat.data <- do.call(rbind, data)                                                #check do.call; #append 20 imputed datasets
    big.data <- data.frame(mat.data, stringsAsFactors = F)
    colnames(big.data) <- colnames(data[[1]])
    k <- length(variable.names)
    # xyw <- matrix(0, sum(nperimp), k + 1)
    #xyw[, 1:k] <- model.matrix(formula, data = big.data)
    #xyw[, k + 1] <- model.response(model.frame(as.formula(formula), data = big.data))

    cov.name <- variable.names
    cov.name2 <- variable
    pos <- match(cov.name2, cov.name)

    if (is.null(bound.lo)) {
      lower.collect <- (unlist(lapply(1:imputations, function(x) fits[[x]]$ci.lower[pos]))) # the lower confidence limits of the parameter for all imputations
      lowerbound.lo <- min(lower.collect) #min of all lower CIs from MI
      upperbound.lo <- max(lower.collect) #max of all lower CIs from MI
    } else {
      lowerbound.lo <- min(bound.lo)
      upperbound.lo <- max(bound.lo)
    }
    if (is.null(bound.up)) {
      upper.collect <- (unlist(lapply(1:imputations, function(x) fits[[x]]$ci.upper[pos]))) # the upper confidence limits of the parameter for all imputations
      lowerbound.up <- min(upper.collect) #min of all upper CIs from MI
      upperbound.up <- max(upper.collect) #max of all upper CIs from MI
    } else {
      lowerbound.up <- min(bound.up)
      upperbound.up <- max(bound.up)
    }
    if (lowerbound.lo == upperbound.lo) {
      lowerbound.lo <- lowerbound.lo - 1/2
      upperbound.lo <- upperbound.lo + 1/2
    }
    if (lowerbound.up == upperbound.up) {
      lowerbound.up <- lowerbound.up - 1/2
      upperbound.up <- upperbound.up + 1/2
    }
    estimate <- mean(unlist(lapply(1:imputations, function(x) fits[[x]]$coefficients[pos])))  #final est is the mean of all imputed estimates
    iter <- numeric(0)

    loglik <- unlist(lapply(1:imputations, function(x) fits[[x]]$loglik$full))   # form a list of full model loglikelihood value of each imputed parameter; loglik[2] is full model loglik;
    beta <- t(matrix(unlist(lapply(1:imputations, function(x) fits[[x]]$coefficients)), (k+1), imputations))  # formatting coefficients of each imputation in a summary matrix

    ### generate mixcure.clip.pdf() externally
    ##########################################


    lpdf <- function(zz, z) mixcure.clip.pdf(formula = formula,
                                             data = big.data[imputation.indicator == zz,],
                                             beta = beta[zz, ], loglik = loglik[zz], pos = pos, pl = pl, b = z)$pdf
    f = function(z) mean(unlist(lapply(1:imputations, function(zz) {
      lpdf(zz, z)
    })))

    ################################## finding lower end of CI ####################################
    f.lower <- f(lowerbound.lo) - ci.level[1] #ci.level=c(0.025,0.975); #compare CDF at beta=lowerbound.lo and 0.025
    f.upper <- f(upperbound.lo) - ci.level[1] #compare CDF at beta=upperbound.lo and 0.025
    iter[1] <- 2

    # start iterations for CDF of beta=lowerbound.lo to converge to 0.025
    itwhile <- 0
    while (f.lower > 0 & (upperbound.lo - lowerbound.lo) > 0 & itwhile < iternum) {
      itwhile <- itwhile + 1
      ######################################################################
      # adjustment condition for lowerbound.lo to converge from right side #
      ######################################################################
      lowerbound.lo <- lowerbound.lo - (upperbound.lo - lowerbound.lo)/2
      f.lower <- f(lowerbound.lo) - ci.level[1]
      iter[1] <- iter[1] + 1
    }
    ci <- numeric(0)
    if (itwhile >= iternum & f.lower > 0)
      # stop("pool.pl can not find a lower boundary for the lower confidence limit.\n Try to increase number of imputations or supply boundaries by bound.lo=c(x,xx).\n")
    {
      ci[1] <- NA  ########## final value for lower end of CI
    } else {
    # start iterations for CDF of beta=upperbound.lo to converge to 0.025
    itwhile <- 0
    while (f.upper < 0 & (upperbound.lo - lowerbound.lo) > 0 & itwhile < iternum) {
      itwhile <- itwhile + 1
      ######################################################################
      # adjustment condition for upperbound.lo to converge from right side #
      ######################################################################
      upperbound.lo <- upperbound.lo + (upperbound.lo - lowerbound.lo)/2
      f.upper <- f(upperbound.lo) - ci.level[1]
      iter[1] <- iter[1] + 1
    }
    if (itwhile >= iternum & f.upper < 0)
      stop("pool.pl can not find an upper boundary for the lower confidence limit.\n Try to increase number of imputations or supply boundaries by bound.lo=c(x,xx).\n")

   # ci <- numeric(0)
    res.ci <- uniroot(f = function(z) {
      f(z) - ci.level[1]
    }, lower = lowerbound.lo, upper = upperbound.lo,
  #  extendInt = "yes",
    f.lower = f.lower, f.upper = f.upper
    )  #root finding for function f(z)=0.025
    ci[1] <- res.ci$root  ########## final value for lower end of CI
    iter[1] <- res.ci$iter + iter[1]
    }
    ################################## finding upper end of CI ####################################

    f.lower <- f(lowerbound.up) - ci.level[2]    #compare CDF at beta=lowerbound.up and 0.025
    f.upper <- f(upperbound.up) - ci.level[2]    #compare CDF at beta=upperbound.up and 0.025
    iter[2] <- 2

    # start iterations for CDF of beta=lowerbound.up to converge to 0.975
    itwhile <- 0
    while (f.upper < 0 & (upperbound.up - lowerbound.up) > 0 & itwhile < iternum) {
      itwhile <- itwhile + 1
      ######################################################################
      # adjustment condition for upperbound.up to converge from right side #
      ######################################################################
      upperbound.up <- upperbound.up + (upperbound.up - lowerbound.up)/2
      f.upper <- f(upperbound.up) - ci.level[2]
      iter[2] <- iter[2] + 1
    }
    if (itwhile >= iternum & f.upper < 0)
     # stop("pool.pl can not find an upper boundary for the upper confidence limit.\n Try to increase number of imputations or supply boundaries by bound.up=c(x,xx).\n")
    {ci[2] <- NA} else {

    # start iterations for CDF of beta=lowerbound.up to converge to 0.975
    itwhile <- 0
    while (f.lower > 0 & (upperbound.up - lowerbound.up) > 0 & itwhile < iternum) {
      itwhile <- itwhile + 1
      ######################################################################
      # adjustment condition for lowerbound.up to converge from right side #
      ######################################################################
      lowerbound.up <- lowerbound.up - (upperbound.up - lowerbound.up)/2
      f.lower <- f(lowerbound.up) - ci.level[2]
      iter[2] <- iter[2] + 1
    }
    if (itwhile >= iternum & f.lower > 0)
      stop("pool.pl can not find a lower boundary for the upper confidence limit.\n Try to increase number of imputations or supply boundaries by bound.up=c(x,xx).\n")

    res.ci <- uniroot(f = function(z) {
      f(z) - ci.level[2]
    }, lower = lowerbound.up, upper = upperbound.up,
  #  extendInt = "yes",
    f.lower = f.lower, f.upper = f.upper) #root finding for function f(z)=0.975
    ci[2] <- res.ci$root ########## final value for upper end of CI
    iter[2] <- res.ci$iter + iter[2]
    }

    if (pvalue == "TRUE") {
      pvalue1 <- f(0)
      pvalue <- 2 * min(pvalue1, (1 - pvalue1))
    } else pvalue <- NA
    res <- list(estimate = estimate, ci = ci, pvalue = pvalue,
                imputations = imputations, ci.level = ci.level, myvar = myvar,
                call = match.call(), bound.lo = c(lowerbound.lo,
                                                  upperbound.lo), bound.up = c(lowerbound.up, upperbound.up),
                iter = iter)
    return(res)
  }

  ##check
  estimate <- numeric(nvar)
  ci <- matrix(0, nvar, 2)
  pvalue.out <- numeric(nvar)
  bound.lo.out <- matrix(0, nvar, 2)
  bound.up.out <- matrix(0, nvar, 2)
  iter <- matrix(0, nvar, 2)
  for (i in 1:(nvar-1)) {
    res.tmp <- inner.CLIP(myvar = variable[i])
    estimate[i] <- res.tmp$estimate
    ci[i, ] <- res.tmp$ci
    pvalue.out[i] <- res.tmp$pvalue
    bound.lo.out[i, ] <- res.tmp$bound.lo
    bound.up.out[i, ] <- res.tmp$bound.up
    iter[i, ] <- res.tmp$iter
  }
  res <- list(variable = variable, estimate = estimate, ci = ci,
              pvalue = pvalue.out, imputations = res.tmp$imputations,
              ci.level = ci.level, bound.lo = bound.lo.out, bound.up = bound.up.out,
              iter = iter, call = match.call())

coef.table <- cbind(
  'coef'        = res$estimate,
  'exp(coef)'   = exp(res$estimate),
  'pval(LRT)' = res$pvalue,
  'LCI.95%' = res$ci[,1],
  'UCI.95%' = res$ci[,2]
);
rownames(coef.table) <- res$variable


attr(coef.table,'class') <- c("clip.mixcure")

return(coef.table)

}
