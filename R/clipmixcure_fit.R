#' Adapter for mixcuref.meta: fit from a mc_spec
#'
#' @param spec A mc_spec object
#' @param control Engine-specific options
#' @return An engine-specific fitted object
#' @export
clipmixcure_fit <- function(spec, control = list()) {

  method <- control$method %||% "penal_mi"

  # For ClipMixcure MI workflow we require a Surv(...) ~ ... formula in control
  fml <- control$formula
  if (is.null(fml) || !inherits(fml, "formula")) {
    stop("clipmixcure_fit(): control$formula must be a formula like Surv(Time, CENS==1) ~ ...",
         call. = FALSE)
  }

  init <- control$init
  pl   <- isTRUE(control$pl)

  dat <- spec$data

  if (identical(method, "penal_mi")) {
    # 1) Multiple imputation case: mids -> mice::with(...)
    if (inherits(dat, "mids")) {
      if (!requireNamespace("mice", quietly = TRUE)) {
        stop("clipmixcure_fit(): package 'mice' is required for mids objects.", call. = FALSE)
      }
      return(with(
        data = dat,
        expr = ClipMixcure::mixcure.penal.mi(fml, init = init, pl = pl)
      ))
    }

    # 2) Standard data.frame case
    if (is.data.frame(dat)) {
      return(ClipMixcure::mixcure.penal.mi(fml, init = init, pl = pl))
    }

    stop("clipmixcure_fit(): for method='penal_mi', spec$data must be a data.frame or a mids object.",
         call. = FALSE)
  }

  if (identical(method, "clip")) {
    # This assumes clip.mixcure() expects an MI-fit object (like your example)
    # So method='clip' should generally be called on a previous MI fit, not directly here.
    stop("clipmixcure_fit(): method='clip' is a post-fit step. Fit with method='penal_mi' first, then call clip.mixcure().",
         call. = FALSE)
  }

  stop("clipmixcure_fit(): unknown control$method = ", method, call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
