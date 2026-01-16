plot_drug_probabilities <- function(est_list, tgrid = seq(0.03, 12*5, length.out = 400),
                                    title = "Temporal Decomposition Logistic Mixed Effects model") {

  df_long <- est_list$df_long
  u_hat <- est_list$u_hat  # now a matrix with columns u1, u2
  est <- est_list$est

  # Extract parameters
  a1      <- as.numeric(est["a1"])
  a2      <- as.numeric(est["a2"])
  beta0_1 <- as.numeric(est["beta0_1"])
  beta0_2 <- as.numeric(est["beta0_2"])
  t_half_early <- as.numeric(est["t_half_early"])
  t_half_late  <- as.numeric(est["t_half_late"])
  eta_early    <- as.numeric(est["eta_early"])
  eta_late     <- as.numeric(est["eta_late"])
  gamma_early  <- as.numeric(est["gamma_early"])
  gamma_late   <- as.numeric(est["gamma_late"])

  # Conditional odds function with 2 random effects
  conditional_odds <- function(t, u1, u2) {
    Omega1 <- exp(beta0_1 + a1 * u1)
    Omega2 <- exp(beta0_2 + a2 * u2)
    Omega1 * get_early_phase(t, t_half_early, eta_early, gamma_early) +
      Omega2 * get_late_phase(t, t_half_late, eta_late, gamma_late)
  }

  subjects <- unique(df_long$Subject_ID)

  # Compute fitted probabilities for all subjects
  sub_list <- lapply(subjects, function(i) {
    u_i <- u_hat[i, ]   # vector c(u1, u2)
    pi_i <- conditional_odds(tgrid, u_i[1], u_i[2]) /
      (1 + conditional_odds(tgrid, u_i[1], u_i[2]))
    data.frame(Subject_ID = i, Time = tgrid, Fitted = pi_i)
  })

  sub_df <- bind_rows(sub_list)

  # Average and 95% CI
  avg_ci_df <- sub_df %>%
    group_by(Time) %>%
    summarise(
      Average_Fitted = mean(Fitted),
      Lower_CI = pmax(0, quantile(Fitted, 0.025)),
      Upper_CI = pmin(1, quantile(Fitted, 0.975)),
      .groups = "drop"
    )

  ggplot() +
    geom_line(data = sub_df, aes(x = Time, y = Fitted, group = Subject_ID),
              alpha = 0.5, color = "grey", linetype = "dashed") +
    geom_line(data = avg_ci_df, aes(x = Time, y = Average_Fitted),
              color = "blue", size = 1) +
    scale_y_continuous(name = "Probability", limits = c(0, 1)) +
    scale_x_continuous(name = "Time (months)", limits = c(0, 60)) +
    labs(title = title) +
    theme_minimal(base_size = 14)
}

# Plotting functions ----
plot_drug_probabilities_lite <- function(est_list, tgrid = seq(0.03, 12*5, length.out = 400),
                                    title = "Temporal Decomposition Logistic Mixed Effects model") {

  df_long <- est_list$df_long
  u_hat <- est_list$u_hat
  est <- est_list$est

  # Extract parameters as numeric
  a1      <- as.numeric(est["a1"])
  a2      <- as.numeric(est["a2"])
  beta0_1      <- as.numeric(est["beta0_1"])
  beta0_2      <- as.numeric(est["beta0_2"])
  t_half_early <- as.numeric(est["t_half_early"])
  t_half_late  <- as.numeric(est["t_half_late"])
  eta_early    <- as.numeric(est["eta_early"])
  eta_late     <- as.numeric(est["eta_late"])
  gamma_early  <- as.numeric(est["gamma_early"])
  gamma_late   <- as.numeric(est["gamma_late"])

  # Conditional odds function
  conditional_odds <- function(t, u_i) {
    Omega1 <- exp(beta0_1 + a1 * u_i)
    Omega2 <- exp(beta0_2 + a2 * u_i)
    Omega1 * get_early_phase(t, t_half_early, eta_early, gamma_early) +
      Omega2 * get_late_phase(t, t_half_late, eta_late, gamma_late)
  }

  subjects <- unique(df_long$Subject_ID)

  # Compute fitted probabilities for all subjects
  sub_list <- lapply(subjects, function(i) {
    u_i <- as.numeric(u_hat[i])
    pi_i <- conditional_odds(tgrid, u_i) / (1 + conditional_odds(tgrid, u_i))
    data.frame(Subject_ID = i, Time = tgrid, Fitted = pi_i)
  })

  sub_df <- bind_rows(sub_list)

  # Average and 95% CI
  avg_ci_df <- sub_df %>%
    group_by(Time) %>%
    summarise(
      Average_Fitted = mean(Fitted),
      Lower_CI = pmax(0, quantile(Fitted, 0.025)),
      Upper_CI = pmin(1, quantile(Fitted, 0.975)),
      .groups = "drop"
    )

  ggplot() +
    geom_line(data = sub_df, aes(x = Time, y = Fitted, group = Subject_ID),
              alpha = 0.5, color = "grey", linetype = "dashed") +
    geom_line(data = avg_ci_df, aes(x = Time, y = Average_Fitted),
              color = "blue", size = 1) +
    scale_y_continuous(name = "Probability", limits = c(0, 1)) +
    scale_x_continuous(name = "Time (months)", limits = c(0, 60)) +
    labs(title = title) +
    theme_minimal(base_size = 14)
}




