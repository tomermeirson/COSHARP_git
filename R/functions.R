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

#Formatted p-v according to request
format_p_value <- function(p_value) {
  if (p_value < 0.001) {
    return('p < 0.001')
  } else if (p_value < 0.01) {
    return(paste0('p = ', round(p_value, 3)))
  } else {
    return(paste0('p = ', format(round(p_value, 2), nsmall = 2)))
  }
}

plot_3_arm_km <- function(
  simulated_data,
  original_data,
  plot_label,
  output_path
) {
  cox_model <- coxph(Surv(time, status) ~ arm, data = simulated_data)

  # Extract HR and CI
  hr <- exp(coef(cox_model))
  hr.lo = summary(cox_model)$conf.int[3]
  hr.hi = summary(cox_model)$conf.int[4]
  p_value_cox <- summary(cox_model)$coefficients[, "Pr(>|z|)"]

  annot.hr = paste0(
    'HR = ',
    round(hr, 2),
    ' (',
    round(hr.lo, 2),
    '-',
    round(hr.hi, 2),
    ');',
    '\n',
    format_p_value(p_value_cox)
  )

  original_data$arm <- as.numeric(as.character(original_data$arm))
  original_data <- original_data %>% mutate(group = ifelse(arm == 0, 0, 1))

  T0 = simulated_data %>% filter(arm == 0)
  T0 <- T0 %>% mutate(group = 2)

  dat.full = rbind(original_data, T0)
  dat.full$group = factor(
    dat.full$group,
    levels = c(1, 0, 2),
    labels = c("Experimental", "Control", "Control counterfactual")
  )
  # group=0, arm 0 in original data, group=1, arm 1 in original, group=2, arm 0 in crossover

  KM_full <- survfit(Surv(time, status) ~ group, data = dat.full)
  time.break = 12
  f.size = 36
  f.x = 36

  KM.comb <- KM_full %>%
    ggsurvplot(
      dat.full,
      palette = c('#FC766AFF', 'blue', '#5B84B1FF'),
      linetype = c("solid", "solid", "dashed"),
      conf.int = F,
      risk.table = TRUE,
      break.time.by = time.break,
      legend.title = '',
      # legend.title = "Outcome groups",
      legend.labs = c("Experimental", "Control", "Control counterfactual"),
      legend = 'none',
      # xlim = c(0,88),
      # legend = "left",
      ggtheme = theme_survminer(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 48)
      ),
      xlab = 'Time (Months)',
      #title = 'K-M',
      fontsize = f.size,
      font.x = f.x,
      font.y = f.x,
      size = 3,
      risk.table.fontsize = 14
    )

  KM.comb$plot <- KM.comb$plot +
    annotate(
      "text_npc",
      # x = 10, y = 0.2, # x and y coordinates of the text
      npcx = 0.05,
      npcy = 0.1, # x and y coordinates of the text
      label = annot.hr,
      size = 14
    ) +
    theme(
      axis.text.x = element_text(size = 36),
      axis.text.y = element_text(size = 36),
      axis.title.x = element_text(size = 36),
      axis.title.y = element_text(size = 36),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
    )

  KM.comb$table <- KM.comb$table +
    theme(
      text = element_text(size = 36),
      axis.text.x = element_text(size = 36),
      axis.title.x = element_text(size = 36),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
    )
  KM.comb$table$theme$axis.text.y$size <- rel(1)

  KM.comb1 = plot_grid(
    KM.comb$plot,
    KM.comb$table,
    nrow = 2,
    rel_heights = c(3, 1),
    align = 'hv'
  )
  plot(KM.comb1)
  # ggsave(paste0(output_path,'combined plot', '_', plot_label, '_noNAR.pdf'),print(KM.comb$plot), width = 14,height = 10, scale=1.37)

  # ggsave(paste0(output_path,'combined plot', '_', plot_label, '_NAR.pdf'),print(KM.comb1), width = 14,height = 12, scale=1.37)
}
