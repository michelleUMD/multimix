plot_drug_probabilities <- function(
    est_list,
    tgrid = NULL,
    n_time_quartiles = 10,
    title = "Temporal Decomposition Logistic Mixed Effects model"
) {

  df_long <- est_list$df_long
  u_hat   <- est_list$u_hat
  est     <- est_list$est

  if(is.null(tgrid)) {
    tgrid <- seq(0, max(df_long$Time), length.out = max(df_long$Time) * 10)
  }

  # small time shift to avoid t = 0 issues in the model
  t_shift <- 0.03
  t_plot  <- tgrid
  t_model <- tgrid + t_shift

  # Extract parameters
  a1 <- as.numeric(est["a1"])
  a2 <- as.numeric(est["a2"])

  beta0_1 <- as.numeric(est["beta0_1"])
  beta0_2 <- as.numeric(est["beta0_2"])

  t_half_early <- as.numeric(est["t_half_early"])
  t_half_late  <- as.numeric(est["t_half_late"])

  eta_early   <- as.numeric(est["eta_early"])
  eta_late    <- as.numeric(est["eta_late"])
  gamma_early <- as.numeric(est["gamma_early"])
  gamma_late  <- as.numeric(est["gamma_late"])

  # Conditional odds with two random effects
  conditional_odds <- function(t, u1, u2) {
    Omega1 <- exp(beta0_1 + a1 * u1)
    Omega2 <- exp(beta0_2 + a2 * u2)

    Omega1 * get_early_phase(t, t_half_early, eta_early, gamma_early) +
      Omega2 * get_late_phase(t, t_half_late, eta_late, gamma_late)
  }

  subjects <- unique(df_long$Subject_ID)

  # Subject-level fitted probabilities
  sub_list <- lapply(subjects, function(i) {
    u_i <- u_hat[i, ]

    odds <- conditional_odds(t_model, u_i[1], u_i[2])
    pi_i <- odds / (1 + odds)

    data.frame(
      Subject_ID = i,
      Time       = t_plot,
      Fitted     = pi_i
    )
  })

  sub_df <- bind_rows(sub_list)

  # Average and 95% CI across subjects
  avg_ci_df <- sub_df %>%
    group_by(Time) %>%
    summarise(
      Average_Fitted = mean(Fitted),
      Lower_CI = pmax(0, quantile(Fitted, 0.025)),
      Upper_CI = pmin(1, quantile(Fitted, 0.975)),
      .groups = "drop"
    )

  # Quartile points from raw data
  quartile_points <- df_long %>%
    mutate(
      Time = as.numeric(Time),
      time_quartile = cut(
        Time,
        breaks = quantile(
          Time,
          probs = seq(0, 1, length.out = n_time_quartiles + 1),
          na.rm = TRUE
        ),
        include.lowest = TRUE
      )
    ) %>%
    group_by(time_quartile) %>%
    summarise(
      Time = mean(Time, na.rm = TRUE),
      Prop_TRUE = mean(Binary_outcome, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot() +
    geom_line(
      data = sub_df,
      aes(x = Time, y = Fitted, group = Subject_ID),
      alpha = 0.4,
      color = "grey",
      linetype = "dashed"
    ) +
    geom_line(
      data = avg_ci_df,
      aes(x = Time, y = Average_Fitted),
      color = "blue",
      linewidth = 1
    ) +
    geom_point(
      data = quartile_points,
      aes(x = Time, y = Prop_TRUE),
    ) +
    scale_y_continuous(name = "Probability", limits = c(0, 1)) +
    scale_x_continuous(name = "Time (months)", limits = c(0, 60)) +
    labs(title = title) +
    theme_minimal(base_size = 14)
}

plot_drug_probabilities_lite <- function(
    est_list,
    tgrid = NULL,
    n_time_quartiles = 10,
    title = "Temporal Decomposition Logistic Mixed Effects model"
) {

  df_long <- est_list$df_long
  u_hat   <- est_list$u_hat
  est     <- est_list$est

  # Default time grid from data
  if (is.null(tgrid)) {
    tgrid <- seq(
      0,
      max(df_long$Time, na.rm = TRUE),
      length.out = max(df_long$Time, na.rm = TRUE) * 10
    )
  }

  # small time shift to avoid t = 0 issues in the model
  t_shift <- 0.03
  t_plot  <- tgrid
  t_model <- tgrid + t_shift

  # Extract parameters
  a1 <- as.numeric(est["a1"])
  a2 <- as.numeric(est["a2"])

  beta0_1 <- as.numeric(est["beta0_1"])
  beta0_2 <- as.numeric(est["beta0_2"])

  t_half_early <- as.numeric(est["t_half_early"])
  t_half_late  <- as.numeric(est["t_half_late"])

  eta_early   <- as.numeric(est["eta_early"])
  eta_late    <- as.numeric(est["eta_late"])
  gamma_early <- as.numeric(est["gamma_early"])
  gamma_late  <- as.numeric(est["gamma_late"])

  # Conditional odds (single random effect)
  conditional_odds <- function(t, u_i) {
    Omega1 <- exp(beta0_1 + a1 * u_i)
    Omega2 <- exp(beta0_2 + a2 * u_i)

    Omega1 * get_early_phase(t, t_half_early, eta_early, gamma_early) +
      Omega2 * get_late_phase(t, t_half_late, eta_late, gamma_late)
  }

  subjects <- unique(df_long$Subject_ID)

  # Subject-level fitted probabilities
  sub_list <- lapply(subjects, function(i) {
    u_i <- as.numeric(u_hat[i])

    odds <- conditional_odds(t_model, u_i)
    pi_i <- odds / (1 + odds)

    data.frame(
      Subject_ID = i,
      Time       = t_plot,
      Fitted     = pi_i
    )
  })

  sub_df <- bind_rows(sub_list)

  # Average and 95% CI across subjects
  avg_ci_df <- sub_df %>%
    group_by(Time) %>%
    summarise(
      Average_Fitted = mean(Fitted),
      Lower_CI = pmax(0, quantile(Fitted, 0.025)),
      Upper_CI = pmin(1, quantile(Fitted, 0.975)),
      .groups = "drop"
    )

  # Quartile point from raw data
  quartile_points <- df_long %>%
    mutate(
      Time = as.numeric(Time),
      time_quartile = cut(
        Time,
        breaks = quantile(
          Time,
          probs = seq(0, 1, length.out = n_time_quartiles + 1),
          na.rm = TRUE
        ),
        include.lowest = TRUE
      )
    ) %>%
    group_by(time_quartile) %>%
    summarise(
      Time = mean(Time, na.rm = TRUE),
      Prop_TRUE = mean(Binary_outcome, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot() +
    geom_line(
      data = sub_df,
      aes(x = Time, y = Fitted, group = Subject_ID),
      alpha = 0.4,
      color = "grey",
      linetype = "dashed"
    ) +
    geom_line(
      data = avg_ci_df,
      aes(x = Time, y = Average_Fitted),
      color = "blue",
      linewidth = 1
    ) +
    geom_point(
      data = quartile_points,
      aes(x = Time, y = Prop_TRUE),
      size = 2.5,
      alpha = 0.8
    ) +
    scale_y_continuous(name = "Probability", limits = c(0, 1)) +
    scale_x_continuous(name = "Time (months)", limits = c(0, 60)) +
    labs(title = title) +
    theme_minimal(base_size = 14)
}
