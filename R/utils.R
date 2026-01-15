# Example initial parameters ----
default_init_example <- c(
  beta0_1 = 10,
  beta0_2 = 0,
  a1 = 1,
  a2 = 2,
  log_sigma1 = 1,
  log_sigma2 = 1,
  log_t_half_early = log(3),
  log_t_half_late = log(10),
  eta_early = 8,
  eta_late = 8,
  gamma_early = 0,
  gamma_late = 0
)

lower_bounds_example <- c(
  beta0_1 = -10,
  beta0_2 = -10,
  a1 = -5,
  a2 = -5,
  log_sigma1 = log(1e-4),
  log_sigma2 = log(1e-4),
  log_t_half_early = log(0.01),
  log_t_half_late = log(0.01),
  eta_early = -5,
  eta_late = -5,
  gamma_early = -5,
  gamma_late = -5
)

upper_bounds_example <- c(
  beta0_1 = 5,
  beta0_2 = 5,
  a1 = 4,
  a2 = 4,
  log_sigma1 = log(10),
  log_sigma2 = log(10),
  log_t_half_early = log(10),
  log_t_half_late = log(20),
  eta_early = 5,
  eta_late = 5,
  gamma_early = 5,
  gamma_late = 5
)


default_init_example_lite <- c(
  beta0_1 = 10,
  beta0_2 = 0,
  a1 = 1,
  a2 = 2,
  log_sigma = 1,
  log_t_half_early = log(3),
  log_t_half_late = log(10),
  eta_early = 8,
  eta_late = 8,
  gamma_early = 0,
  gamma_late = 0
)

lower_bounds_example_lite <- c(
  beta0_1 = -10,
  beta0_2 = -10,
  a1 = -5,
  a2 = -5,
  log_sigma = log(1e-4),
  log_t_half_early = log(0.01),
  log_t_half_late = log(0.01),
  eta_early = -2,
  eta_late = -2,
  gamma_early = -2,
  gamma_late = -2
)

upper_bounds_example_lite <- c(
  beta0_1 = 5,
  beta0_2 = 5,
  a1 = 4,
  a2 = 4,
  log_sigma = log(10),
  log_t_half_early = log(10),
  log_t_half_late = log(20),
  eta_early = 2,
  eta_late = 2,
  gamma_early = 2,
  gamma_late = 2
)

generate_random_init_from_bounds <- function(
    lower_bounds,
    upper_bounds
) {

  repeat {
    init <- mapply(
      function(lo, hi) runif(1, lo, hi),
      lower_bounds,
      upper_bounds
    )

    # Make sure not both are negative (undefined)
    eta_early   <- init["eta_early"]
    eta_late    <- init["eta_late"]
    gamma_early <- init["gamma_early"]
    gamma_late  <- init["gamma_late"]

    # reject invalid draws
    if (
      (eta_early < 0 && gamma_early < 0) ||
      (eta_late  < 0 && gamma_late  < 0 ||
       abs(eta_early) < 0.1 || abs(eta_late) < 0.1)
    ) {
      next
    }

    return(init)
  }

  init
}

# This is just to show the df structure
# Should be in long format
# Columns are: Subject_ID, Time, Binary_outcome
get_df_long <- function(df) {
  df_long <- df %>%
    pivot_longer(
      cols = -time,
      names_to = "patient",
      values_to = "on_drug"
    ) %>%
    mutate(
      patient = factor(patient),
      on_drug = as.integer(on_drug)
    ) %>%
    filter(!is.na(on_drug))

  df_long <- data.frame(
    Subject_ID = df_long$patient,
    Time = df_long$time + 0.03,
    Binary_outcome = df_long$on_drug
  ) %>%
    arrange(Subject_ID, Time)

  rownames(df_long) <- NULL

  return(df_long)
}
