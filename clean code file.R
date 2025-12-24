library(magrittr)
library(dplyr)
library(tibble)
library(purrr)
library(survival)
library(flexsurv)
library(simsurv)
library(survRM2)
library(pbapply)
library(boot)
library(clue)
library(rstudioapi)
library(tidyr)

path.dir = dirname(getSourceEditorContext()$path) %>% dirname # get current directory
dat <- read.csv(paste0(path.dir, '/MONALEESA7 OS.csv'))
dat.sur <- read.csv(paste0(path.dir, '/MONALEESA7 PFS.csv'))
# --- Initial Data Jittering ---
# Adding a small random jitter to event times to break ties,
# which can prevent numerical issues in survival model fitting.
set.seed(11111)

dat$time <- dat$time + runif(nrow(dat), min = 1e-6, max = 1e-5)
dat.sur$time <- dat.sur$time + runif(nrow(dat.sur), min = 1e-6, max = 1e-5)

###################################################
### Function for CI % proportion
###################################################

# Bootstrap Confidence Interval for Proportion of Lost Significance
#
#  Calculates the mean and bootstrap confidence intervals (percentile
#  and BCa) for the proportion of simulation replicates where the p-value
#  for the Hazard Ratio (HR) is > 0.05 (i.e., where significance is "lost").
#
#  df is a vector of p-values for the Hazard Ratio from all simulation replicates.
#
#  return a tibble with the mean proportion, standard deviation, and
#  95% confidence interval bounds (percentile and BCa).
#
bootstrap.sig.prop = function(df) {
  # Convert p-values into a binary vector (1 = not significant, 0 = significant)
  sig <- as.integer(df > 0.05)

  # 1) Define the statistic function for the bootstrap
  # The statistic is the mean of the sampled binary vector (the proportion of 1s)
  boot_stat <- function(data, indices) {
    # data[indices] picks a bootstrap sample of the sig vector
    mean(data[indices])
  }

  # 2) Run the bootstrap simulation
  set.seed(11111) # for reproducibility
  res <- boot(
    data = sig, # the 0/1 indicator vector
    statistic = boot_stat,
    R = 1000 # number of bootstrap replications
  )

  # 3) Inspect results
  res$t0 # original estimate (should equal mean(sig))
  sd(res$t) # bootstrap estimate of the SE

  # 4) Get bootstrap CIs (percentile & BCa)
  # Handle cases where the standard deviation is zero (e.g., if all p-values are < 0.05)
  if (sd(res$t) == 0) {
    ci <- boot.ci(res, type = "perc")
    # Manually set BCa bounds to NA if variance is zero (BCa relies on non-zero variance)
    ci$bca[4] <- NA
    ci$bca[5] <- NA
  } else {
    ci = boot.ci(res, type = c("perc", "bca"))
  }

  # 5) Consolidate results into a tibble
  results <- tibble(
    mean = res$t0, # Original estimate (mean proportion)
    sd = sd(res$t), # Bootstrap estimate of the standard error
    ci_low_perc = ci$perc[4], # Percentile 95% CI lower bound
    ci_high_perc = ci$perc[5], # Percentile 95% CI upper bound
    ci_low_bca = ci$bca[4], # BCa 95% CI lower bound
    ci_high_bca = ci$bca[5] # BCa 95% CI upper bound
  )
  return(results)
}


###Function for creating data for crossover simulation

logcumhaz <- function(t, x, betas, fit, HR_subseq) {
  #  Defines the log-cumulative hazard H(t) for a patient who
  #  progresses at 'prog_time' and would have continued on 1st-line therapy
  #  (counterfactual) instead of receiving 2nd-line treatment.
  #
  #  The cumulative hazard is piecewise:
  #  - For t <= prog_time: H(t) is based on the baseline 1st-line hazard (H1).
  #  - For t > prog_time: H(t) is H1(prog_time) + HR_subseq * [H1(t) - H1(prog_time)].
  #    This formula applies the subsequent-line HR to the *excess* hazard
  #    accumulated after progression, assuming the subsequent treatment has
  #    an effect relative to the baseline hazard of the first line.
  #
  #  t is the time point at which to calculate the cumulative hazard.
  #  x a data frame row containing the patient's progression time ('prog_time').
  #  betas not used here (for models with covariates).
  #  fit the fitted flexible parametric survival model (flexsurvspline) for
  #  the control arm, used to estimate the baseline hazard (H1).
  #  HR_subseq the assumed Hazard Ratio for the subsequent-line treatment
  # (crossover therapy), relative to the first-line control.
  #
  #  The function returns the natural logarithm of the counterfactual cumulative hazard H(t).

  prog_time <- x[["prog_time"]]

  # Retrieve cumulative hazard from the fitted baseline model (H1_val = H1(t))
  h1_df <- summary(fit, type = "cumhaz", t = t, tidy = TRUE)
  H1_val <- as.numeric(h1_df$est[1])

  if (t <= prog_time) {
    # Before progression: use the baseline cumulative hazard H1(t)
    return(log(H1_val))
  } else {
    # At or after progression:    # Retrieve cumulative hazard at progression time
    H1_prog_df <- summary(fit, type = "cumhaz", t = prog_time, tidy = TRUE)
    H1_prog <- as.numeric(H1_prog_df$est[1])

    #  piecewise formula to compute the combined cumulative hazard.
    H_combined <- max(H1_prog + HR_subseq * (H1_val - H1_prog), 1e-10)

    return(log(H_combined))
  }
}


match.endpoints = function(dat0, dat.sur0) {
  #  Matches each patient's Overall Survival (OS) record to their
  #  Progression-Free Survival (PFS) record using the Hungarian algorithm to
  #  solve the Linear Sum Assignment Problem (LSAP). This is robust to
  #  slight mismatches or non-perfect indexing.
  #
  #  The cost matrix penalizes pairings where OS time is NOT greater than PFS time.
  #
  #  dat0 The OS data frame for the control arm.
  #  dat.sur0 the PFS data frame for the control arm.
  #
  #  Returns a single merged data frame with matched OS and PFS records.

  # Create unique IDs for matching
  dat0$id.os = 1:nrow(dat0)
  dat.sur0$id.prog = 1:nrow(dat.sur0)

  A = dat0
  B = dat.sur0
  # Rename PFS columns for clarity in the merged dataset
  B %<>% rename(prog_time = time, prog_time_status = status) %>% select(-'arm')

  # We build a cost matrix (B rows x A columns) for the assignment problem.
  # A pairing is "cheap" when OS time (ta) > PFS time (tb) with cost = ta - tb.
  # A pairing is "expensive" otherwise (penalty + absolute time difference).
  penalty <- 1e6
  cost_mat <- outer(B$prog_time, A$time, FUN = function(tb, ta) {
    # If OS time > PFS time (the desired scenario)
    ifelse(
      ta > tb,
      ta - tb,
      # Otherwise, apply a large penalty
      penalty + abs(ta - tb)
    )
  })

  # Solve the assignment problem using the Hungarian algorithm (LSAP)
  # 'assignment' maps B rows to A columns (assignment[i] = matched A row for B row i)
  assignment <- solve_LSAP(cost_mat)
  assignment <- as.integer(assignment)

  # Create an inverse mapping: 'inv_assignment' maps A rows to B rows
  inv_assignment <- integer(nrow(A))
  inv_assignment[assignment] <- seq_len(nrow(A))

  # Combine matched rows: each row in A gets the corresponding matched row from B.

  final_matched <- cbind(A, B[inv_assignment, ])
  return(final_matched)
}

sample.cross.deficit = function(dat.match, cross.prop) {
  #  Samples the control arm patients whose OS will be adjusted.
  #  The number of sampled patients is determined by 'cross.prop'.
  #  Patients with observed progression (prog_time_status == 1) are prioritized
  #  before sampling from non-progressors/censored patients.
  #
  #  dat.match The matched OS and PFS data frame for the control arm.
  #  cross.prop The proportion of the control arm that crossed over (the deficit).
  #
  #  returns a subset of 'dat.match' containing the sampled crossover deficit patients.

  # Calculate the target number of patients to adjust
  no.of.pd = floor(cross.prop * nrow(dat.match))
  # Case 1: The number of observed progressors is LESS than the target deficit size
  if (sum(dat.match$prog_time_status == 1) < no.of.pd) {
    # Use all progressors
    dat0.cross.i.pd = filter(dat.match, prog_time_status == 1) %>% pull(id.os)
    # Sample the remainder from non-progressors (who are assumed to have progressed at censoring)
    no.of.non.pd = no.of.pd - length(dat0.cross.i.pd)
    dat0.cross.i.non.pd = filter(dat.match, prog_time_status == 0) %>%
      pull(id.os) %>%
      sample(no.of.non.pd)
    # Combine indices
    dat0.cross.i = c(dat0.cross.i.pd, dat0.cross.i.non.pd)
  } else {
    # if the number of cross-over defict is lower than the nubmer of events sample them
    # Case 2: The number of observed progressors is GREATER than or equal to the target deficit size
    # Sample 'no.of.pd' patients from the progressors
    dat0.cross.i = filter(dat.match, prog_time_status == 1) %>%
      pull(id.os) %>%
      sample(no.of.pd)
  }

  # Subset dataframe using the sampled indices
  dat0.cross <- dat.match[dat.match$id.os %in% dat0.cross.i, ]
  # Add negligible jitter to time to avoid bugs in flexsurv or simsurv
  dat0.cross$time <- dat0.cross$time +
    runif(nrow(dat0.cross), min = 1e-6, max = 1e-5)
  return(dat0.cross)
}


find.best.fit.k = function(dat.tmp) {
  #  Fits flexsurvspline models with k=1 to k=5 (the number of
  #  internal knots) and selects the k that minimizes the Akaike Information
  #  Criterion (AIC).
  #
  #  dat.tmp the survival data frame for the control arm.
  #
  #  returns the optimal number of internal knots (k).

  k_vals <- 1:5

  # Fit a flexsurvspline model for each k
  fits <- lapply(k_vals, function(k) {
    flexsurvspline(Surv(time, status) ~ 1, data = dat.tmp, k = k)
  })

  # Extract AIC for each fit
  aics <- sapply(fits, AIC)

  # Pick the k with lowest AIC
  k.best <- k_vals[which.min(aics)]
  cat("Optimal k by AIC is", k.best, "\n")

  return(k.best)
}

sim.crossover = function(dat0.cross, HR_subseq_treatment, k.best) {
  #  Fits a flexsurvspline model to the control arm and then simulates
  #  counterfactual OS times for the crossover deficit patients using the
  #  piecewise log-cumulative hazard function.
  #
  #  dat0.cross the sampled crossover deficit patients.
  #  HR_subseq_treatment the assumed subsequent-line HR (must be > 0).
  #  k.best the initial optimal number of knots (k) for flexsurvspline.
  #  maxt.multiplier multiplier for the maximum simulation time (to ensure
  #  a complete curve; defaults to 1).
  #
  #  returns a tibble with the simulated OS times and event status for the deficit patients.

  k_value <- k.best
  # Loop to attempt fitting with decreasing k if the model fails (e.g., non-positive definite Hessian)
  while (k_value > 0) {
    fit_attempt <- tryCatch(
      {
        flexsurv::flexsurvspline(
          Surv(time, status) ~ 1,
          data = dat0.cross,
          k = k_value
        )
      },
      warning = function(w) {
        # Return NULL to trigger the k-decrement
        return(NULL) # Continue trying lower k
      },
      error = function(e) {
        # Return NULL to trigger the k-decrement
        return(NULL) # Continue trying lower k
      }
    )
    if (!is.null(fit_attempt)) {
      break # Exit loop if successful
    } else {
      k_value <- k_value - 1
      cat("Trying k =", k_value, "\n")
    }
  }
  #Fit first-line model #
  fit.1L <- flexsurvspline(
    Surv(time, status) ~ 1,
    data = dat0.cross,
    k = k_value
  )

  # simulate counterfactual survival
  try(
    {
      sim = simsurv(
        logcumhazard = logcumhaz,
        x = dat0.cross,
        fit = fit.1L,
        HR_subseq = HR_subseq_treatment,
        betas = NULL,
        # Max time for simulation should be at least as long as the max observed time
        maxt = max(dat0.cross$time),
        interval = c(1E-20, 1E20)
      )
    }
  )

  sim %<>% rename(time = eventtime)
  sim$arm = 0

  return(sim)
}

update.dat0.combine.dat1 = function(dat.match, sim, dat0.cross, dat1) {
  #  Replaces the original OS data of the crossover deficit patients
  #  in the control arm with their newly simulated counterfactual OS data.
  #  The resulting adjusted control arm is then combined with the original
  #  experimental arm data.
  #
  #  dat.match the full matched control arm data (OS and PFS).
  #  sim the simulated counterfactual survival data for the deficit patients.
  #  dat0.cross the subset of patients whose data was simulated.
  #  dat1 the original survival data for the experimental arm.
  #
  #  returns the final combined data frame for survival analysis.

  # 1. Identify and keep the control arm patients who were NOT adjusted
  dat0.non.cross.tmp = dat.match[!dat.match$id.os %in% dat0.cross$id.os, ]
  dat0.non.cross.tmp %<>% select(arm, time, status)
  # 2. Extract the simulated data (the adjusted subset)
  sim.cross = sim %>% select(time, status, arm)
  # 3. Combine the non-adjusted control arm data with the simulated data
  dat0.update = bind_rows(dat0.non.cross.tmp, sim.cross)

  # 4. Combine the adjusted control arm (dat0.update) with the experimental arm (dat1)
  dat.sim = bind_rows(dat0.update, dat1)
}

extract.survival.data <- function(dat.temp, NCT.number) {
  #  Fits a Cox model to the combined data to calculate Hazard Ratio (HR)
  #  and performs Restricted Mean Survival Time (RMST) analysis.
  #
  #  dat.temp The combined, adjusted survival data.
  #  NCT.number The trial identifier (for record keeping).
  #
  #  returns a data frame summarizing the HR, RMST-D, and survival probabilities.

  # Extract KM curve data for the median HR simulation
  cox_model <- coxph(Surv(time, status) ~ arm, data = dat.temp)

  # Extract HR and CI
  hr <- exp(coef(cox_model))
  hr.lo = summary(cox_model)$conf.int[3]
  hr.hi = summary(cox_model)$conf.int[4]
  p_value_cox <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
  rmst = rmst2(dat.temp$time, dat.temp$status, dat.temp$arm, tau = NULL) # non-time restricted RMST
  rmstd = rmst$unadjusted.result[1, 1] # RMST-D
  rmstd.hi = rmst$unadjusted.result[1, 3] # 95% CI-high
  rmstd.lo = rmst$unadjusted.result[1, 2] # 95% CI-low
  rmstd.p = rmst$unadjusted.result[1, 4] # p-value

  return(data.frame(
    hr = hr,
    hr.lo = hr.lo,
    hr.hi = hr.hi,
    hr.p = p_value_cox,
    rmstd = rmstd,
    rmstd.p = rmstd.p
  ))
}


# --- Simulation Parameters (Inputs) ---
# NOTE: The values for these variables (cross.prop, OS.HR.subsequent.line, NCT.number, dat, dat.sur)
# must be defined and loaded before this code block is run.

# cross.prop: Proportion of control arm patients assumed to have crossed over.
# OS.HR.subsequent.line: The assumed Hazard Ratio for the subsequent line of therapy.
# NCT.number: ClinicalTrials.gov identifier for the trial.
# dat: Data frame for Overall Survival (OS).
# dat.sur: Data frame for Progression-Free Survival (PFS).

replicates = 10

cross.prop = 0.43
HR_subseq_treatment = 0.73
NCT.number = 'NCT02278120'
# --- Data Preparation ---
dat0 = dat %>% filter(arm == 0) # Control Arm OS
dat1 = dat %>% filter(arm == 1) # Experimental Arm OS
dat.sur0 = dat.sur %>% filter(arm == 0) # Control Arm PFS
# Step 1: Match OS and PFS data for the control arm
dat.match = match.endpoints(dat0, dat.sur0) # Match OS and PFS data
# Step 2: Determine the optimal 'k' for the flexsurvspline model based on the *original* control arm data
k.best = find.best.fit.k(dat0)
# --- Simulation Execution ---
set.seed(11111)
results = pblapply(1:replicates, function(i) {
  # run the simulation N times
  try({
    # 1. Sample the crossover deficit patients
    dat0.cross = sample.cross.deficit(dat.match, cross.prop)
    # 2. Simulate counterfactual survival for the sampled patients
    sim.result = sim.crossover(dat0.cross, HR_subseq_treatment, k.best)
    # 3. Replace original survival data with simulated data and combine with experimental arm
    dat.update = update.dat0.combine.dat1(
      dat.match,
      sim.result,
      dat0.cross,
      dat1
    )
    # 4. Summarize the survival outcome for this replicate
    summary = extract.survival.data(dat.update, NCT.number) # summarize the survival

    list(dat.update = dat.update, summary = summary)
  })
})

# --- Post-Simulation Analysis ---

# 1. Tidy the summary results from all replicates
summary.df = tibble(
  id = seq_along(results),
  summary = map(results, ~ pluck(.x, "summary", .default = list()))
) %>%
  unnest_wider(summary)

dat.update.list <- map(results, "dat.update")

# 2. Calculate the median adjusted HR
median.hr = median(summary.df$hr, na.rm = T)
# 3. Identify the replicate whose HR is closest to the median HR (for plotting/reporting a representative curve)
median.hr.i <- match(
  quantile(summary.df$hr, 0.5, na.rm = TRUE, type = 1),
  summary.df$hr
)
sim.with.median.hr = dat.update.list[median.hr.i] %>% .[[1]]
# 4. Calculate the proportion of replicates where significance is lost (hr.p > 0.05)
mean.hr.lost.sign <- mean(summary.df$hr.p > 0.05)
# 5. Bootstrap the proportion of lost significance
mean_data = bootstrap.sig.prop(summary.df$hr.p)
# 6. Extract full survival data for the simulation with the median HR (the representative result)
out = extract.survival.data(sim.with.median.hr, NCT.number)
# 7. Final output summary
out1 = cbind(out, mean.hr.lost.sign, NCT.number, mean_data)
# Print the final result
print(out1)
