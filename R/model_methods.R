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
#' from a `multi_mix_model` or `multi_mix_model_lite` object.
#'
#' @param x An object of class `multi_mix_model` or `multi_mix_model_lite`,
#'          typically returned by `fit_multi_mix_model()` or its lite variant.
#' @param ... Additional arguments (currently ignored) to allow method dispatch compatibility.
#'
#' @return Invisibly returns the original object `x`.
#' @export
#' @method print multi_mix_model
print.multi_mix_model <- function(x, ...) {
  est <- x$est

  est_long <- data.frame(
    parameter = names(est),
    estimate  = unname(est),
    stringsAsFactors = FALSE
  )

  est_wide <- est_long |>
    mutate(
      parameter = if_else(
        parameter %in% names(VAR_TO_ENGLISH_DICT),
        VAR_TO_ENGLISH_DICT[parameter],
        parameter
      )
    ) |>
    pivot_wider(
      names_from  = model,
      values_from = estimate
    )

  gt_tbl <- est_wide |>
    gt(rowname_col = "parameter") |>
    fmt_number(
      columns = everything(),
      decimals = 2
    ) |>
    tab_header(
      title = "Estimated Model Parameters"
    ) |>
    opt_table_outline()

  print(gt_tbl)
  invisible(x)
}

# Export the lite class method
#' @export
#' @method print multi_mix_model_lite
print.multi_mix_model_lite <- print.multi_mix_model


# Plotting functions ----
#' Plot drug probability trajectories for a multi-mix model
#'
#' This S3 method plots the estimated probabilities of drug administration
#' over time for each drug class, using the multimix model object.
#'
#' @param x An object of class `multi_mix_model`, typically returned by `fit_model_with_retries()`.
#' @param ... Additional arguments passed to underlying plotting functions (currently ignored).
#'
#' @return `ggplot` object
#' @export
#' @method plot multi_mix_model
plot.multi_mix_model <- function(x, ...) {
  plot_drug_probabilities(x, ...)
}

#' Plot simplified drug probability trajectories for a multi-mix model
#'
#' This S3 method plots the estimated probabilities of drug administration
#' over time for each drug class, using a lighter version of the multimix model
#'
#' @param x An object of class `multi_mix_model`, typically returned by `fit_multi_mix_model()`.
#' @param ... Additional arguments passed to underlying plotting functions (currently ignored).
#'
#' @return `ggplot` object
#' @export
#' @method plot multi_mix_model_lite
plot.multi_mix_model_lite <- function(x, ...) {
  plot_drug_probabilities_lite(x, ...)
}

