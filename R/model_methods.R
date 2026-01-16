# Display Results ----

## Table of coefficients ----

VAR_TO_ENGLISH_DICT <- c(
  a1 = "a₁ (early)",
  a2 = "a₂ (late)",

  beta0_1 = "β₀₁ (early)",
  beta0_2 = "β₀₂ (late)",

  t_half_early = "t½ (early)",
  t_half_late  = "t½ (late)",

  eta_early   = "η (early)",
  eta_late    = "η (late)",

  gamma_early = "γ (early)",
  gamma_late  = "γ (late)",

  sigma = "σ"
)

#' Print method for multi-mix model objects
#'
#' This S3 method prints a nicely formatted table of estimated parameters
#' from a `multimix_model` or `multimix_model_lite` object.
#'
#' @param x An object of class `multimix_model` or `multimix_model_lite`,
#'          typically returned by `multimix()` or its lite variant.
#' @param ... Additional arguments (currently ignored) to allow method dispatch compatibility.
#'
#' @return Invisibly returns the original object `x`.
#' @export
#' @method print multimix_model
print.multimix_model <- function(x, ...) {
  est <- x$est

  est_long <- data.frame(
    parameter = names(est),
    estimate  = unname(est),
    stringsAsFactors = FALSE
  )

  est_long$parameter <- ifelse(
    est_long$parameter %in% names(VAR_TO_ENGLISH_DICT),
    VAR_TO_ENGLISH_DICT[est_long$parameter],
    est_long$parameter
  )

  rownames(est_long) <- est_long$parameter
  est_long$parameter <- NULL

  print(
    round(est_long, 2),
    ...
  )

  invisible(x)
}


# Export the lite class method
#' @export
#' @method print multimix_model_lite
print.multimix_model_lite <- print.multimix_model


# Plotting functions ----
#' Plot drug probability trajectories for a multi-mix model
#'
#' This S3 method plots the estimated probabilities of drug administration
#' over time for each drug class, using the multimix model object.
#'
#' @param x An object of class `multimix_model`, typically returned by `fit_model_with_retries()`.
#' @param ... Additional arguments passed to underlying plotting functions (currently ignored).
#'
#' @return `ggplot` object
#' @export
#' @method plot multimix_model
plot.multimix_model <- function(x, ...) {
  plot_drug_probabilities(x, ...)
}

#' Plot simplified drug probability trajectories for a multi-mix model
#'
#' This S3 method plots the estimated probabilities of drug administration
#' over time for each drug class, using a lighter version of the multimix model
#'
#' @param x An object of class `multimix_model`, typically returned by `multimix()`.
#' @param ... Additional arguments passed to underlying plotting functions (currently ignored).
#'
#' @return `ggplot` object
#' @export
#' @method plot multimix_model_lite
plot.multimix_model_lite <- function(x, ...) {
  plot_drug_probabilities_lite(x, ...)
}

