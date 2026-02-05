#' Adapter for meta-layer: fit ClipMixcure from a mc_spec
#'
#' This function is a stable entry point intended for use by the mixcure.meta
#' package. It should map the unified `spec` object to the existing ClipMixcure
#' fitting routine.
#'
#' @param spec A mc_spec object (list-like) containing formulas and data.
#' @param control A list of engine-specific controls.
#' @return An engine-specific fitted model object.
#' @export
clipmixcure_fit <- function(spec, control = list()) {
  # Minimal adapter used by the meta-layer.
  #
  # ClipMixcure's legacy core routines are written around a *single* survival
  # formula Surv(time, status) ~ ... and assume the same covariates are used in
  # both the incidence and latency parts of the mixture cure model.
  #
  # The meta-layer provides incidence/latency formulas separately; here we form
  # a combined survival formula using the union of covariate terms.
  stopifnot(is.list(spec))
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("clipmixcure_fit(): package 'survival' is required.", call. = FALSE)
  }
  if (is.null(spec$data) || !is.data.frame(spec$data)) {
    stop("clipmixcure_fit(): 'spec$data' must be a data.frame.", call. = FALSE)
  }

  dat <- spec$data
  time_var <- spec$time %||% "time"
  status_var <- spec$status %||% "status"

  if (!time_var %in% names(dat)) {
    stop("clipmixcure_fit(): time column not found in data: ", time_var, call. = FALSE)
  }
  if (!status_var %in% names(dat)) {
    stop("clipmixcure_fit(): status column not found in data: ", status_var, call. = FALSE)
  }

  # Build combined RHS from incidence + latency formulas.
  inc_terms <- character(0)
  lat_terms <- character(0)
  if (!is.null(spec$incidence) && inherits(spec$incidence, "formula")) {
    inc_terms <- attr(stats::terms(spec$incidence), "term.labels")
  }
  if (!is.null(spec$latency) && inherits(spec$latency, "formula")) {
    lat_terms <- attr(stats::terms(spec$latency), "term.labels")
  }
  rhs_terms <- unique(c(inc_terms, lat_terms))
  rhs <- if (length(rhs_terms) == 0) {
    "1"
  } else {
    paste(rhs_terms, collapse = " + ")
  }

  # Use survival::Surv(...) in the formula so we don't depend on attaching survival.
  formula <- stats::as.formula(
    paste0("survival::Surv(", time_var, ", ", status_var, ") ~ ", rhs)
  )

  pl <- isTRUE(control$pl)
  iterlim <- control$iterlim %||% 200

  # If no init provided, create a sensible default: zeros for regression
  # coefficients and gamma=1 for the Weibull shape.
  init <- control$init
  if (is.null(init)) {
    mf <- stats::model.frame(formula, data = dat, na.action = stats::na.omit)
    mm <- stats::model.matrix(formula, data = mf)
    k <- ncol(mm)
    init <- c(rep(0, 2 * k), 1)
  }

  mixcure.penal.mi(
    formula = formula,
    data = dat,
    init = init,
    pl = pl,
    iterlim = iterlim
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
