# Survival Analysis: Crossover Adjustment (MONALEESA-7)

## Overview
This project performs a counterfactual simulation to estimate Overall Survival (OS) adjusted for treatment crossover, using the **MONALEESA-7** trial data (NCT02278120).

## Project Structure
* `main.R`: The primary execution script. Run this to perform the analysis.
* `R/functions.R`: Contains the custom logic (RPSFT/Piecewise constant hazard simulation).
* `data/`: **(Local Only)** Place your raw CSV files here.
* `output/`: Generated results.

## Setup Instructions

### 1. Environment
This project uses `renv` to sync package versions.
1. Open this folder in Positron or RStudio.
2. Run `renv::restore()` in the console to install the required libraries.

### 2. Data
1. Create a `data/` folder in the project root.
2. Add the following files:
   - `MONALEESA7 OS.csv`
   - `MONALEESA7 PFS.csv`

## Usage
Open `main.R`. You can modify simulation parameters (replicates, crossover proportion) in the "Control Center" at the top of the file.