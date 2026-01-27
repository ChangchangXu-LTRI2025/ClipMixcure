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
  stop("clipmixcure_fit(): adapter not yet wired to ClipMixcure internal fit function.")
}
