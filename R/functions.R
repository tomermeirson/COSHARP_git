# R/functions.R

# --- Helper Functions for Survival Simulation ---

#' Bootstrap Confidence Interval for Proportion of Lost Significance
bootstrap.sig.prop = function(df) {
  sig <- as.integer(df > 0.05)

  boot_stat <- function(data, indices) {
    mean(data[indices])
  }

  set.seed(11111)
  res <- boot::boot(
    data = sig,
    statistic = boot_stat,
    R = 1000
  )

  if (sd(res$t) == 0) {
    ci <- boot::boot.ci(res, type = "perc")
    ci$bca[4] <- NA
    ci$bca[5] <- NA
  } else {
    ci = boot::boot.ci(res, type = c("perc", "bca"))
  }

  results <- tibble::tibble(
    mean = res$t0,
    sd = sd(res$t),
    ci_low_perc = ci$perc[4],
    ci_high_perc = ci$perc[5],
    ci_low_bca = ci$bca[4],
    ci_high_bca = ci$bca[5]
  )
  return(results)
}

#' Log-cumulative hazard function for counterfactual simulation
logcumhaz <- function(t, x, betas, fit, HR_subseq) {
  prog_time <- x[["prog_time"]]

  h1_df <- summary(fit, type = "cumhaz", t = t, tidy = TRUE)
  H1_val <- as.numeric(h1_df$est[1])

  if (t <= prog_time) {
    return(log(H1_val))
  } else {
    H1_prog_df <- summary(fit, type = "cumhaz", t = prog_time, tidy = TRUE)
    H1_prog <- as.numeric(H1_prog_df$est[1])
    H_combined <- max(H1_prog + HR_subseq * (H1_val - H1_prog), 1e-10)
    return(log(H_combined))
  }
}

#' Match OS and PFS endpoints using Hungarian Algorithm
match.endpoints = function(dat0, dat.sur0) {
  dat0$id.os = 1:nrow(dat0)
  dat.sur0$id.prog = 1:nrow(dat.sur0)

  A = dat0
  B = dat.sur0 %>%
    dplyr::rename(prog_time = time, prog_time_status = status) %>%
    dplyr::select(-'arm')

  penalty <- 1e6
  cost_mat <- outer(B$prog_time, A$time, FUN = function(tb, ta) {
    ifelse(ta > tb, ta - tb, penalty + abs(ta - tb))
  })

  assignment <- clue::solve_LSAP(cost_mat)
  assignment <- as.integer(assignment)

  inv_assignment <- integer(nrow(A))
  inv_assignment[assignment] <- seq_len(nrow(A))

  final_matched <- cbind(A, B[inv_assignment, ])
  return(final_matched)
}

#' Sample the crossover deficit patients
sample.cross.deficit = function(dat.match, cross.prop) {
  no.of.pd = floor(cross.prop * nrow(dat.match))

  if (sum(dat.match$prog_time_status == 1) < no.of.pd) {
    dat0.cross.i.pd = dplyr::filter(dat.match, prog_time_status == 1) %>%
      dplyr::pull(id.os)
    no.of.non.pd = no.of.pd - length(dat0.cross.i.pd)
    dat0.cross.i.non.pd = dplyr::filter(dat.match, prog_time_status == 0) %>%
      dplyr::pull(id.os) %>%
      sample(no.of.non.pd)
    dat0.cross.i = c(dat0.cross.i.pd, dat0.cross.i.non.pd)
  } else {
    dat0.cross.i = dplyr::filter(dat.match, prog_time_status == 1) %>%
      dplyr::pull(id.os) %>%
      sample(no.of.pd)
  }

  dat0.cross <- dat.match[dat.match$id.os %in% dat0.cross.i, ]
  # Add negligible jitter
  dat0.cross$time <- dat0.cross$time +
    runif(nrow(dat0.cross), min = 1e-6, max = 1e-5)
  return(dat0.cross)
}

#' Find best fit K for flexsurvspline
find.best.fit.k = function(dat.tmp) {
  k_vals <- 1:5
  fits <- lapply(k_vals, function(k) {
    flexsurv::flexsurvspline(Surv(time, status) ~ 1, data = dat.tmp, k = k)
  })
  aics <- sapply(fits, AIC)
  k.best <- k_vals[which.min(aics)]
  cat("Optimal k by AIC is", k.best, "\n")
  return(k.best)
}

#' Simulate counterfactual crossover
sim.crossover = function(dat0.cross, HR_subseq_treatment, k.best) {
  k_value <- k.best

  # Iteratively reduce k if fitting fails
  while (k_value > 0) {
    fit_attempt <- tryCatch(
      {
        flexsurv::flexsurvspline(
          Surv(time, status) ~ 1,
          data = dat0.cross,
          k = k_value
        )
      },
      warning = function(w) return(NULL),
      error = function(e) return(NULL)
    )
    if (!is.null(fit_attempt)) {
      break
    } else {
      k_value <- k_value - 1
    }
  }

  fit.1L <- flexsurv::flexsurvspline(
    Surv(time, status) ~ 1,
    data = dat0.cross,
    k = k_value
  )

  try({
    sim = simsurv::simsurv(
      logcumhazard = logcumhaz,
      x = dat0.cross,
      fit = fit.1L,
      HR_subseq = HR_subseq_treatment,
      betas = NULL,
      maxt = max(dat0.cross$time),
      interval = c(1E-20, 1E20)
    )
  })

  sim %<>% dplyr::rename(time = eventtime)
  sim$arm = 0
  return(sim)
}

#' Combine simulated data with original data
update.dat0.combine.dat1 = function(dat.match, sim, dat0.cross, dat1) {
  dat0.non.cross.tmp = dat.match[!dat.match$id.os %in% dat0.cross$id.os, ]
  dat0.non.cross.tmp %<>% dplyr::select(arm, time, status)

  sim.cross = sim %>% dplyr::select(time, status, arm)
  dat0.update = dplyr::bind_rows(dat0.non.cross.tmp, sim.cross)
  dat.sim = dplyr::bind_rows(dat0.update, dat1)
  return(dat.sim)
}

#' Extract Survival Statistics (HR, RMST)
extract.survival.data <- function(dat.temp, NCT.number) {
  cox_model <- survival::coxph(Surv(time, status) ~ arm, data = dat.temp)

  hr <- exp(coef(cox_model))
  hr.lo = summary(cox_model)$conf.int[3]
  hr.hi = summary(cox_model)$conf.int[4]
  p_value_cox <- summary(cox_model)$coefficients[, "Pr(>|z|)"]

  rmst = survRM2::rmst2(
    dat.temp$time,
    dat.temp$status,
    dat.temp$arm,
    tau = NULL
  )
  rmstd = rmst$unadjusted.result[1, 1]
  rmstd.p = rmst$unadjusted.result[1, 4]

  return(data.frame(
    hr = hr,
    hr.lo = hr.lo,
    hr.hi = hr.hi,
    hr.p = p_value_cox,
    rmstd = rmstd,
    rmstd.p = rmstd.p
  ))
}
