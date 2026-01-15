# Model ----

fit_model <- function(df_long,
                      nGH = 40,
                      fixed_pars = list(),
                      default_init = default_init_example) {

  # Gauss–Hermite quadrature ----
  gh <- gauss.quad.prob(nGH, dist = "normal")
  nodes <- gh$nodes
  weights <- gh$weights

  # Distinguish parameters to optimize vs fixed
  full_param_names <- names(default_init_example)
  train_pars <- setdiff(full_param_names, names(fixed_pars))
  init_pars  <- default_init[train_pars]

  # Compute negative log likelihood for optimization
  negLogLik <- function(pars, data, nodes, weights, fixed_pars = list(), eps = 1e-300) {

    # Combine fixed and optimized parameters
    full_pars <- pars
    for (p in names(fixed_pars)) full_pars[[p]] <- fixed_pars[[p]]
    last_pars <<- full_pars

    beta0_1 <- full_pars["beta0_1"]
    beta0_2 <- full_pars["beta0_2"]
    a1       <- full_pars["a1"]
    a2       <- full_pars["a2"]

    # Must be > 0, so using exp and log
    sigma1 <- exp(full_pars["log_sigma1"])
    sigma2 <- exp(full_pars["log_sigma2"])
    t_half_early <- exp(full_pars["log_t_half_early"])
    t_half_late  <- exp(full_pars["log_t_half_late"])

    eta_early    <- full_pars["eta_early"]
    eta_late     <- full_pars["eta_late"]
    gamma_early  <- full_pars["gamma_early"]
    gamma_late   <- full_pars["gamma_late"]

    if(eta_early < 0 & gamma_early < 0 | eta_late < 0 & gamma_late < 0) {
      stop("Generic function undefined for eta and gamma both < 0")
    }

    Lambda_1 <- function(t) get_early_phase(t, t_half_early, eta_early, gamma_early)
    Lambda_2 <- function(t) get_late_phase(t, t_half_late, eta_late, gamma_late)

    conditional_odds <- function(t, u1, u2) {
      exp(beta0_1 + a1 * u1) * Lambda_1(t) +
        exp(beta0_2 + a2 * u2) * Lambda_2(t)
    }

    subjects <- unique(data$Subject_ID)
    total_loglik <- 0

    for (i in subjects) {
      subdat <- data[data$Subject_ID == i, ]
      y_i <- subdat$Binary_outcome
      t_i <- subdat$Time

      # 2D Gauss-Hermite quadrature for two random effects
      logL_k <- sapply(nodes, function(k1) {
        sapply(nodes, function(k2) {
          u1_k <- sigma1 * k1
          u2_k <- sigma2 * k2
          odds <- conditional_odds(t_i, u1_k, u2_k)
          pi_k <- odds / (1 + odds)
          pi_k <- pmin(pmax(pi_k, 1e-12), 1 - 1e-12)
          sum(y_i * log(pi_k) + (1 - y_i) * log(1 - pi_k))
        })
      })

      # Integrate using weights
      Li <- sum(weights %o% weights * exp(logL_k))  # outer product of weights
      total_loglik <- total_loglik + log(Li + eps)
    }

    return(-total_loglik)
  }

  # Optimization ----
  last_pars <- NULL

  cat("Starting Nelder–Mead\n")
  opt1 <- tryCatch(
    optim(
      init_pars, negLogLik,
      data = df_long,
      nodes = nodes,
      weights = weights,
      fixed_pars = fixed_pars,
      method = "Nelder-Mead",
      control = list(maxit = 5000)
    ),
    error = function(e) {
      message("Nelder-Mead failed: ", conditionMessage(e))
      message("Last parameters evaluated:")
      print(last_pars)
      stop(e)
    }
  )

  cat("After Nelder–Mead:\n")
  print(opt1$par)
  cat("LogLik =", -opt1$value, "\n\n")

  cat("Starting L-BFGS-B refinement\n")
  opt2 <- tryCatch(
    optim(
      opt1$par, negLogLik,
      data = df_long,
      nodes = nodes,
      weights = weights,
      fixed_pars = fixed_pars,
      method = "L-BFGS-B",
      hessian = TRUE,
      control = list(maxit = 2000)
    ),
    error = function(e) {
      message("L-BFGS-B failed: ", conditionMessage(e))
      message("Last parameters evaluated:")
      print(last_pars)
      stop(e)
    }
  )

  cat("After L-BFGS-B:\n")
  print(opt2$par)
  cat("LogLik =", -opt2$value, "\n\n")

  # Combine optimized params with fixed params ----
  est <- c(opt2$par, unlist(fixed_pars, use.names = TRUE))
  names(est) <- c(train_pars, names(fixed_pars))
  est <- est[full_param_names]
  est <- as.numeric(est)
  names(est) <- full_param_names

  # Exponentiate all "log_" variables
  log_params <- grep("^log_", names(est), value = TRUE)
  for (lp in log_params) {
    plain_name <- sub("^log_", "", lp)
    est[plain_name] <- exp(est[lp])
    est <- est[!names(est) %in% lp]
  }

  # Empirical Bayes estimates for two random effects ----
  conditional_odds_fun <- function(t, u1, u2) {
    exp(est["beta0_1"] + est["a1"] * u1) * get_early_phase(t, est["t_half_early"], est["eta_early"], est["gamma_early"]) +
      exp(est["beta0_2"] + est["a2"] * u2) * get_late_phase(t, est["t_half_late"], est["eta_late"], est["gamma_late"])
  }

  subjects <- unique(df_long$Subject_ID)

  estimate_u_i <- function(y, t) {
    obj <- function(u) {
      u1 <- u[1]
      u2 <- u[2]
      pi <- conditional_odds_fun(t, u1, u2) / (1 + conditional_odds_fun(t, u1, u2))
      -sum(y * log(pi) + (1 - y) * log(1 - pi)) +
        u1^2 / (2 * est["sigma1"]^2) + u2^2 / (2 * est["sigma2"]^2)
    }
    optim(c(0, 0), obj, method = "BFGS")$par
  }

  u_hat <- t(sapply(subjects, function(i) {
    sub <- subset(df_long, Subject_ID == i)
    estimate_u_i(sub$Binary_outcome, sub$Time)
  }))
  rownames(u_hat) <- subjects
  colnames(u_hat) <- c("u1", "u2")

  # Return value ----
  list(
    df_long = df_long,
    est = est,
    u_hat = u_hat,
    logLik = -opt2$value
  )
}


fit_model_with_retries <- function(
    df_long,
    fixed_pars,
    lower_bounds,
    upper_bounds,
    max_tries = 20,
    return_first_sucess = FALSE,
    verbose = TRUE
) {

  set.seed(1234)

  best_fit <- NULL
  best_logLik <- -Inf

  for (attempt in seq_len(max_tries)) {

    if (verbose) {
      message("Attempt ", attempt, " / ", max_tries)
    }

    init <- generate_random_init_from_bounds(lower_bounds, upper_bounds)

    fit <- tryCatch(
      fit_model(
        df_long,
        fixed_pars = fixed_pars,
        default_init = init
      ),
      error = function(e) {
        if (verbose) {
          message("Failed: ", conditionMessage(e))
        }
        NULL
      }
    )

    if (!is.null(fit)) {
      if (verbose) {
        message("Success | logLik = ", round(fit$logLik, 3))
      }

      if (is.finite(fit$logLik) && fit$logLik > best_logLik) {
        best_logLik <- fit$logLik
        best_fit <- fit

        if (verbose) {
          message("Updating with new best fit")
        }
      }
      if(return_first_sucess) {
        return(best_fit)
      }
    }
  }

  if (is.null(best_fit)) {
    stop("All ", max_tries, " attempts failed. Model did not converge.")
  }

  if (verbose) {
    message("Best logLik across tries: ", round(best_logLik, 3))
  }

  best_fit
}

# Iterate over drug list ----

estimated_param_list <- list()

fixed_pars <- list(
  # eta_early = 1,
  # a1 = 1,
  # beta0_2 = 0
  # gamma_early = 0,
  # gamma_late = 0
)

# fixed_pars_per_drug <- list(
#   "Anti IL-1 agents" = list(eta_early = 0),
#   "Corticosteroids" = list(),
#   "Colchicine" = list()
# )


# Loop over each drug class # "NSAIDs",
drug_classes_of_interest <-  c("Anti IL-1 agents",
                               "Corticosteroids", "Colchicine")
for (drug_class in drug_classes_of_interest) {

  message("Estimating model for drug class: ", drug_class)

  est_list <- fit_model_with_retries(
    df_long = get_df_long(drug_df_list[[drug_class]]),
    fixed_pars      = fixed_pars,
    lower_bounds    = lower_bounds_example,
    upper_bounds    = upper_bounds_example
  )


  print(plot_drug_probabilities(est_list, title = drug_class))

  # Store the parameter estimates
  estimated_param_list[[drug_class]] <- est_list
}

# Display Results ----

# Check results
estimated_param_list

plot_all_drug_classes_facet(estimated_param_list)


for (nm in names(estimated_param_list)) {
  cat("Model:", nm, "\n")
  cat("  logLik:", estimated_param_list[[nm]]$logLik, "\n")
  cat("  est:\n")
  print(estimated_param_list[[nm]]$est)
  cat("\n")
}


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



est_long <- lapply(names(estimated_param_list), function(nm) {
  est <- estimated_param_list[[nm]]$est
  data.frame(
    model = nm,
    parameter = names(est),
    estimate = unname(est),
    stringsAsFactors = FALSE
  )
}) |>
  bind_rows()

# Wide: parameters = rows, models = columns
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

est_wide |>
  gt(rowname_col = "parameter") |>
  fmt_number(
    columns = everything(),
    decimals = 2
  ) |>
  tab_header(
    title = "Estimated Model Parameters"
  ) |>
  opt_table_outline()

