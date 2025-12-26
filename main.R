# main.R

# 1. Load Libraries -------------------------------------------------------
library(here)
library(magrittr)
library(dplyr)
library(tibble)
library(purrr)
library(survival)
library(flexsurv)
library(simsurv)
library(survRM2)
library(survminer)
library(pbapply)
library(boot)
library(clue)
library(tidyr)
library(ggpp)
library(cowplot)

# 2. Source Functions -----------------------------------------------------
# This loads all the logic from your separate file
source(here("R", "functions.R"))

# 3. Control Center (Parameters) ------------------------------------------
SIM_PARAMS <- list(
  replicates = 1, # Number of simulation runs
  cross_prop = 0.43, # Proportion of crossover
  HR_subseq = 0.73, # Hazard Ratio for subsequent line
  NCT_number = 'NCT02278120',
  seed = 11111 # Master seed
)

# 4. Load & Preprocess Data -----------------------------------------------
# Using here() makes this work on everyone's computer
dat <- read.csv(here("data", "MONALEESA7 OS.csv"))
dat.sur <- read.csv(here("data", "MONALEESA7 PFS.csv"))

# Jittering to prevent ties
set.seed(SIM_PARAMS$seed)
dat$time <- dat$time + runif(nrow(dat), min = 1e-6, max = 1e-5)
dat.sur$time <- dat.sur$time + runif(nrow(dat.sur), min = 1e-6, max = 1e-5)

dat0 <- dat %>% filter(arm == 0) # Control Arm OS
dat1 <- dat %>% filter(arm == 1) # Experimental Arm OS
dat.sur0 <- dat.sur %>% filter(arm == 0) # Control Arm PFS

# 5. Run Simulation -------------------------------------------------------

# Step A: Match Endpoints
dat.match <- match.endpoints(dat0, dat.sur0)

# Step B: Find optimal K (based on original data)
k.best <- find.best.fit.k(dat0)

# Step C: Parallel Simulation Loop
set.seed(SIM_PARAMS$seed)
results <- pblapply(1:SIM_PARAMS$replicates, function(i) {
  try({
    # 1. Sample crossover deficit
    dat0.cross <- sample.cross.deficit(dat.match, SIM_PARAMS$cross_prop)

    # 2. Simulate counterfactual survival
    sim.result <- sim.crossover(dat0.cross, SIM_PARAMS$HR_subseq, k.best)

    # 3. Combine data
    dat.update <- update.dat0.combine.dat1(
      dat.match,
      sim.result,
      dat0.cross,
      dat1
    )

    # 4. Summarize
    summary <- extract.survival.data(dat.update, SIM_PARAMS$NCT_number)

    list(dat.update = dat.update, summary = summary)
  })
})

# 6. Post-Simulation Analysis ---------------------------------------------

# Tidy results
summary.df <- tibble(
  id = seq_along(results),
  summary = map(results, ~ pluck(.x, "summary", .default = list()))
) %>%
  unnest_wider(summary)

dat.update.list <- map(results, "dat.update")

# Median statistics
median.hr <- median(summary.df$hr, na.rm = TRUE)
median.hr.i <- match(
  quantile(summary.df$hr, 0.5, na.rm = TRUE, type = 1),
  summary.df$hr
)
sim.with.median.hr <- dat.update.list[[median.hr.i]]

# Bootstrap significance
mean.hr.lost.sign <- mean(summary.df$hr.p > 0.05)
mean_data <- bootstrap.sig.prop(summary.df$hr.p)

# Final summary of the representative (median) run
out <- extract.survival.data(sim.with.median.hr, SIM_PARAMS$NCT_number)
out1 <- cbind(
  out,
  mean.hr.lost.sign,
  NCT_number = SIM_PARAMS$NCT_number,
  mean_data
)

# 7. Print Results --------------------------------------------------------
print(out1)

plot_3_arm_km(sim.with.median.hr, dat)
