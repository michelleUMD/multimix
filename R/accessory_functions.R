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

# Plotting functions ----
plot.multi_mix_model <- function(x, ...) {
  plot_drug_probabilities(x)
}

plot.multi_mix_model_lite <- function(x, ...) {
  plot_drug_probabilities_lite(x)
}

