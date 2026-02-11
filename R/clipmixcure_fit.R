#' @export
clipmixcure_fit <- function(spec, control = list()) {

  method <- control$method
  if (is.null(method)) method <- "penal_mi"

  fml <- control$formula
  if (is.null(fml) || !inherits(fml, "formula")) {
    stop("clipmixcure_fit(): control$formula must be a formula like Surv(Time, CENS==1) ~ ...",
         call. = FALSE)
  }

  init <- control$init
  pl   <- isTRUE(control$pl)

  dat <- spec$data

  if (identical(method, "penal_mi")) {

    if (inherits(dat, "mids")) {
      if (!requireNamespace("mice", quietly = TRUE)) {
        stop("clipmixcure_fit(): package 'mice' is required for mids objects.", call. = FALSE)
      }

      res <- with(
        data = dat,
        expr = local({
          fml2 <- fml
          environment(fml2) <- environment()
          ClipMixcure::mixcure.penal.mi(fml2, init = init, pl = pl)
        })
      )

      ## KEY: make object self-contained for downstream clip.mixcure()
      attr(res, "mids_data") <- dat
      attr(res, "model_formula") <- fml

      return(res)
    }

    if (is.data.frame(dat)) {
      return(ClipMixcure::mixcure.penal.mi(fml, init = init, pl = pl))
    }

    stop("clipmixcure_fit(): for method='penal_mi', spec$data must be a data.frame or a mids object.",
         call. = FALSE)
  }

  stop("clipmixcure_fit(): unknown control$method = ", method, call. = FALSE)
}
