# Model ----

fit_model_lite <- function(df_long,
                           nGH = 40,
                           fixed_pars = list(),
                           full_param_names = lite_param_names,
                           default_init = default_init_example) {


  # Gauss–Hermite quadrature ----
  gh <- gauss.quad.prob(nGH, dist = "normal")
  nodes <- gh$nodes
  weights <- gh$weights

  # Distinguish parameters to optimize vs fixed
  train_pars <- setdiff(full_param_names, names(fixed_pars))
  init_pars  <- default_init[train_pars]

  # Compute negative log liklihood for optimization
  negLogLik <- function(pars, data, nodes, weights, fixed_pars = list(), eps = 1e-300) {

    # Combine fixed and optimized parameters
    full_pars <- pars
    for (p in names(fixed_pars)) {
      full_pars[[p]] <- fixed_pars[[p]]
    }

    last_pars <<- full_pars

    beta0_1 <- full_pars["beta0_1"]
    beta0_2 <- full_pars["beta0_2"]
    a1       <- full_pars["a1"]
    a2       <- full_pars["a2"]

    # Use exp so that these values are always positive
    sigma_u  <- exp(full_pars["log_sigma"])
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

    conditional_odds <- function(t, u_i) {
      Omega1 <- exp(beta0_1 + a1 * u_i)
      Omega2 <- exp(beta0_2 + a2 * u_i)

      Lambda1 <- Lambda_1(t)
      Lambda2 <- Lambda_2(t)

      Omega1 * Lambda_1(t) + Omega2 * Lambda_2(t)
    }

    subjects <- unique(data$Subject_ID)
    total_loglik <- 0

    for (i in subjects) {
      subdat <- data[data$Subject_ID == i, ]
      y_i <- subdat$Binary_outcome
      t_i <- subdat$Time

      # Compute log-likelihood for each quadrature node
      logL_k <- sapply(nodes, function(k) {
        u_k <- sigma_u * k
        odds <- conditional_odds(t_i, u_k)

        pi_k <- odds / (1 + odds)
        pi_k <- pmin(pmax(pi_k, 1e-12), 1 - 1e-12)
        sum(y_i * log(pi_k) + (1 - y_i) * log(1 - pi_k))
      })

      maxlog <- max(logL_k)
      Li <- exp(maxlog) * sum(weights * exp(logL_k - maxlog))
      total_loglik <- total_loglik + log(Li + eps)
    }

    return(-total_loglik)
  }


  # Optimization ----
  # BFGS optimizer often fails, so making a variable to capture the last set of params
  # right before it fails
  last_pars <- NULL

  cat("Starting L-BFGS-B refinement\n")
  opt2 <- tryCatch(
    optim(
      init_pars, negLogLik,
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

  # Exponentiate all "log_" variables to guarantee positive
  log_params <- grep("^log_", names(est), value = TRUE)

  for (lp in log_params) {
    plain_name <- sub("^log_", "", lp)      # remove "log_" prefix
    est[plain_name] <- exp(est[lp])        # exponentiate
    est <- est[!names(est) %in% lp]        # remove the original log_ entry
  }

  # Empirical Bayes to find random effects ----
  # sigma_u <- exp(est["log_sigma"]))
  # t_half_early <- exp(est["log_t_half_early"])
  # t_half_late  <- exp(est["log_t_half_late"])

  conditional_odds_fun <- function(t, u_i) {
    exp(est["beta0_1"] + est["a1"] * u_i) *
      get_early_phase(t, est["t_half_early"], est["eta_early"], est["gamma_early"]) +
      exp(est["beta0_2"] + est["a2"] * u_i) *
      get_late_phase(t, est["t_half_late"], est["eta_late"], est["gamma_late"])
  }

  subjects <- unique(df_long$Subject_ID)

  estimate_u_i <- function(y, t) {
    obj <- function(u) {
      pi <- conditional_odds_fun(t, u) /
        (1 + conditional_odds_fun(t, u))
      -sum(y * log(pi) + (1 - y) * log(1 - pi)) +
        u^2 / (2 * est["sigma"]^2)
    }
    optim(0, obj, method = "BFGS")$par
  }

  u_hat <- sapply(subjects, function(i) {
    sub <- subset(df_long, Subject_ID == i)
    estimate_u_i(sub$Binary_outcome, sub$Time)
  })
  names(u_hat) <- subjects


  # Return value ----
  multi_mix_model_lite <- list(
    df_long = df_long,
    est = est,
    u_hat = u_hat,
    logLik = -opt2$value
  )

  class(multi_mix_model_lite) <- "multi_mix_model_lite"
  return(multi_mix_model_lite)

}

fit_model_with_retries_lite <- function(
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
      fit_model_lite(
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
