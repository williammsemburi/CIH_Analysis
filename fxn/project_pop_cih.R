#' Demographic Projection with Health Interventions
#' This script provides a descriptive, tidyverse-only implementation that reproduces the outputs of project_pop.R exactly (same function signature and return object names), while using clearer names internally.
#'
#' @author Tarun Shankar Choudhary
#  email: tarun.choudhary@uib.no; tschoudhary@outlook.com
#  Based on original codes by Jan-Magnus Økland (jan-magnus.okland@uib.no)
#' @version 1.0.0
#' Dependencies: tidyverse
#'
#' Required globals (same as project_pop_fc.R):
#' - date_year (integer base year)
#' - locations (data.frame with columns iso3, location_name)
#' - wpp.pars (list of country parameters with base.pop, mx, asfr, mig, srb)
#' - wpp2024.actual (data.frame with iso3 and WPP reference data)
#' - pregnancies.stillbirths (data.frame: iso3, stillbirth.rate, birthtopreg.scale)
#' - fert.impact(iso3, c1t, scaleup, calc.fert, fert.targets)
#' - get.lt(mx, parm = NULL, dws = NULL, qual.adj = FALSE)
#' - HYs(...)
#'
#' Returns exactly the same keys as project_pop_fc.R:
#' - period.df, cohort.mat, actual.wpp, hale.mat, asfr

# ====
# Loading and if needed installing packages ####
# List all required packages
required_packages <- c("tidyverse")

# Identify any packages that are not installed
missing_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]

# Install missing packages, with error handling
if (length(missing_packages) > 0) {
  message(
    "Installing missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
  tryCatch(
    {
      install.packages(missing_packages, repos = "http://cran.us.r-project.org")
    },
    error = function(e) {
      # This code runs if install.packages() fails
      stop(
        paste(
          "Failed to install packages:",
          paste(missing_packages, collapse = ", "),
          "\nPlease install them manually and rerun the script. Error:",
          e$message
        ),
        call. = FALSE
      )
    }
  )
}

# Load all required packages, suppressing startup messages
suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
})

message("All required packages are successfully loaded.")

# ====
# Descriptive constants ####
# can be changed if needed
sex_labels <- c("Female", "Male")

age_groups <- list(
  all_ages = 0:100,
  reproductive_ages = 11:55 # matches project_pop_fc.R births computation
)

india_states <- c(
  "IN-AP",
  "IN-AR",
  "IN-AS",
  "IN-BR",
  "IN-CG",
  "IN-DL",
  "IN-GA",
  "IN-GJ",
  "IN-HP",
  "IN-HR",
  "IN-JH",
  "IN-JK",
  "IN-KA",
  "IN-KL",
  "IN-MH",
  "IN-ML",
  "IN-MN",
  "IN-MP",
  "IN-MZ",
  "IN-NL",
  "IN-OD",
  "IN-OT",
  "IN-PB",
  "IN-RJ",
  "IN-SK",
  "IN-TN",
  "IN-TR",
  "IN-TS",
  "IN-UK",
  "IN-UP",
  "IN-WB"
)

# ====
# Validation and helper functions ####

validate_inputs <- function(iso3c, mub, mua, yldb, ylda) {
  if (!is.character(iso3c) || nchar(iso3c) < 3) {
    stop("iso3c must be a valid ISO code string (e.g., 'TZA', 'ZNZ', 'IN-AP').")
  }
  for (nm in c("mub", "mua", "yldb", "ylda")) {
    x <- get(nm)
    if (!is.array(x) || length(dim(x)) != 3) {
      stop(sprintf("%s must be a 3D array [cohort, age, sex].", nm))
    }
  }
  needed <- c(
    "date_year",
    "locations",
    "wpp.pars",
    "wpp2024.actual",
    "pregnancies.stillbirths",
    "fert.impact",
    "get.lt",
    "HYs"
  )
  missing <- needed[!vapply(needed, exists, logical(1), inherits = TRUE)]
  if (length(missing)) {
    stop("Missing required globals: ", paste(missing, collapse = ", "))
  }
}

# Convert cohort -> period with descriptive return (list of two 2D matrices)
# Always pad using baseline first row (lower triangular), then extract diagonals.
cohort_to_period <- function(baseline_matrix, adjusted_matrix) {
  n_cohorts <- nrow(baseline_matrix)
  n_ages <- ncol(baseline_matrix)

  # Build lower-triangular extension from baseline first row
  baseline_first_row <- baseline_matrix[1, ]
  baseline_lower <- matrix(
    rep(baseline_first_row, each = 100),
    nrow = 100,
    ncol = n_ages,
    byrow = FALSE
  )
  baseline_lower[col(baseline_lower) > row(baseline_lower)] <- NA

  baseline_combined <- rbind(baseline_lower, baseline_matrix)
  adjusted_combined <- rbind(baseline_lower, adjusted_matrix)

  period_adjusted <- matrix(
    NA_real_,
    nrow = n_cohorts,
    ncol = n_ages
  )

  period_baseline <- matrix(
    NA_real_,
    nrow = n_cohorts,
    ncol = n_ages
  )

  for (cohort_year in 1:n_cohorts) {
    idx <- n_cohorts - cohort_year + 1
    rng <- idx:(idx + 100) # 101 ages
    period_adjusted[cohort_year, ] <- diag(adjusted_combined[rng, ])
    period_baseline[cohort_year, ] <- diag(baseline_combined[rng, ])
  }

  list(adjusted = period_adjusted, baseline = period_baseline)
}

# Convert period -> cohort with descriptive return (list of two 2D matrices)
period_to_cohort <- function(period_adjusted, period_baseline) {
  n_cohorts <- nrow(period_adjusted)
  n_ages <- ncol(period_adjusted)

  adjusted_canvas <- matrix(
    NA_real_,
    nrow = n_cohorts + n_ages - 1,
    ncol = n_ages
  )
  baseline_canvas <- matrix(
    NA_real_,
    nrow = n_cohorts + n_ages - 1,
    ncol = n_ages
  )

  for (cohort_year in 1:n_cohorts) {
    idx <- n_cohorts - cohort_year + 1
    rng <- idx:(idx + 100)
    diag(adjusted_canvas[rng, ]) <- period_adjusted[cohort_year, ]
    diag(baseline_canvas[rng, ]) <- period_baseline[cohort_year, ]
  }

  list(
    adjusted = adjusted_canvas[n_ages:(n_cohorts + n_ages - 1), ],
    baseline = baseline_canvas[n_ages:(n_cohorts + n_ages - 1), ]
  )
}

# Survival probabilities from mortality rates (clip > 1)
calculate_survival_probabilities <- function(mortality_rates) {
  rates <- mortality_rates
  rates[rates > 1] <- 1 - 1e-6
  get.lt(rates, "Sx")
}

# Convert 3D array [year, sex, age] to long data.frame
array3d_to_long <- function(
  a,
  metric,
  source,
  date_year,
  sex_labels = sex_labels
) {
  dims <- dim(a) # [ny, ns, na]
  if (length(dims) != 3) {
    stop("Expected a 3D array [year, sex, age].")
  }
  ny <- dims[1]
  ns <- dims[2]
  na <- dims[3]

  # Reorder to [age, sex, year] so that as.vector varies age fastest, then sex, then year
  a2 <- aperm(a, c(3, 2, 1)) # [age, sex, year]
  val <- as.vector(a2)

  # The grid must be created with the same varying order as the flattened vector:
  # age varies fastest, then sex, then year.
  grid <- expand.grid(
    age = 0:(na - 1L), # fastest
    sex = seq_len(ns), # middle
    year = seq_len(ny) # slowest
  )

  data.frame(
    age = grid$age,
    year = date_year + grid$year,
    sex = factor(grid$sex, labels = sex_labels),
    value = val,
    metric = metric,
    source = source,
    cause = "All",
    stringsAsFactors = FALSE
  )
}


# Convert 4D array [cause, year, sex, age] to long data.frame
array4d_to_long <- function(
  a,
  metric,
  source,
  cause_names,
  date_year,
  sex_labels = sex_labels
) {
  dims <- dim(a) # [nc, ny, ns, na]
  if (length(dims) != 4) {
    stop("Expected a 4D array [cause, year, sex, age].")
  }
  nc <- dims[1]
  ny <- dims[2]
  ns <- dims[3]
  na <- dims[4]

  # Reorder to [age, sex, year, cause] so that as.vector varies age fastest and cause slowest
  a2 <- aperm(a, c(4, 3, 2, 1)) # [na, ns, ny, nc]
  val <- as.vector(a2)

  grid <- expand.grid(
    age = 0:(na - 1L), # fastest
    sex = seq_len(ns),
    year = seq_len(ny),
    cause = seq_len(nc) # slowest
  )

  data.frame(
    age = grid$age,
    year = date_year + grid$year,
    sex = factor(grid$sex, labels = sex_labels),
    cause = as.character(factor(grid$cause, labels = cause_names)),
    value = val,
    metric = metric,
    source = source,
    stringsAsFactors = FALSE
  )
}

# Fertility impacts with ZNZ and India subnational handling
fertility_impact_for_iso <- function(
  iso3c,
  c1t,
  scaleup,
  calc.fert,
  fert.targets
) {
  if (iso3c != "ZNZ" && !(iso3c %in% india_states)) {
    fert.impact(iso3c, c1t, scaleup, calc.fert, fert.targets)
  } else if (iso3c == "ZNZ") {
    fert.impact("TZA", c1t, scaleup, calc.fert, fert.targets) %>%
      mutate(iso = "ZNZ", country = "Zanzibar")
  } else {
    fert.impact("IND", c1t, scaleup, calc.fert, fert.targets) %>%
      mutate(iso = iso3c, country = "India")
  }
}

# Project population for one year (descriptive)
project_one_year <- function(
  current_population, # [sex, age]
  mortality_rates, # [sex, age]
  survival_probabilities, # [sex, age]
  fertility_rates_by_age, # [age]
  sex_ratio_at_birth,
  migration_rates, # [sex, age]
  n_ages
) {
  new_population <- array(0, dim = dim(current_population))
  deaths <- array(0, dim = dim(current_population))
  continuous_deaths <- array(0, dim = dim(current_population))
  births <- c(0, 0) # female, male

  age_idx <- 2:n_ages
  prev_age_idx <- 1:(n_ages - 1)

  # Age progression with mortality
  deaths[, age_idx] <- current_population[, prev_age_idx] *
    mortality_rates[, age_idx]

  continuous_deaths[, age_idx] <- current_population[, prev_age_idx] *
    (1 - survival_probabilities[, age_idx])

  new_population[, age_idx] <- current_population[, prev_age_idx] *
    survival_probabilities[, age_idx]

  # Terminal age group (100+)
  deaths[, n_ages] <- deaths[, n_ages] +
    current_population[, n_ages] * mortality_rates[, n_ages]

  continuous_deaths[, n_ages] <- continuous_deaths[, n_ages] +
    current_population[, n_ages] * (1 - survival_probabilities[, n_ages])

  new_population[, n_ages] <- new_population[, n_ages] +
    current_population[, n_ages] * survival_probabilities[, n_ages]

  # Births (age 0) — use mid-year female population 10–54 (indices 11:55)
  reproductive_ages <- age_groups$reproductive_ages
  reproductive_pop_mid <- 0.5 *
    (new_population[1, reproductive_ages] +
      current_population[1, reproductive_ages])
  total_births <- sum(reproductive_pop_mid * fertility_rates_by_age)

  births[2] <- total_births * sex_ratio_at_birth / (1 + sex_ratio_at_birth) # male
  births[1] <- total_births - births[2] # female

  # Apply newborn mortality
  deaths[, 1] <- births * mortality_rates[, 1]
  continuous_deaths[, 1] <- births * (1 - survival_probabilities[, 1])
  new_population[, 1] <- births * survival_probabilities[, 1]

  # Migration
  new_population <- new_population * (1 + migration_rates)

  list(
    population = new_population,
    deaths = deaths,
    continuous_deaths = continuous_deaths,
    births = births
  )
}

# ====
# Main function ####
#' Project population with health intervention effects
#' @return list(period.df, cohort.mat, actual.wpp, hale.mat, asfr)
# New, descriptive parameters
project_pop <- function(
  iso3c = "TZA",
  mortality_baseline = NULL,
  mortality_adjusted = NULL,
  cause_deaths_baseline = NULL,
  cause_deaths_adjusted = NULL,
  scaleup = 25,
  calibration_method = NULL,
  contraceptive_target = NULL,
  calculate_fertility_impact = NULL,
  fertility_targets = NULL,
  disability_baseline = NULL,
  disability_adjusted = NULL,
  ...
) {
  # Backward-Compatibility Mapping Layer
  dots <- list(...)

  map_arg <- function(new_val, old_name, arg_list) {
    if (!is.null(arg_list[[old_name]])) {
      warning(sprintf(
        "Argument '%s' is deprecated. Please use the new descriptive name.",
        old_name
      ))
      return(arg_list[[old_name]])
    }
    return(new_val)
  }

  # Map old arguments to their internal variable names
  mub <- map_arg(mortality_baseline, "mub", dots)
  mua <- map_arg(mortality_adjusted, "mua", dots)
  deathb <- map_arg(cause_deaths_baseline, "deathb", dots)
  deatha <- map_arg(cause_deaths_adjusted, "deatha", dots)
  calib <- map_arg(calibration_method, "calib", dots)
  c1t <- map_arg(contraceptive_target, "c1t", dots)
  calc.fert <- map_arg(calculate_fertility_impact, "calc.fert", dots)
  fert.targets <- map_arg(fertility_targets, "fert.targets", dots)
  yldb <- map_arg(disability_baseline, "yldb", dots)
  ylda <- map_arg(disability_adjusted, "ylda", dots)

  # Check that essential arguments are not NULL
  if (is.null(mub) || is.null(mua) || is.null(yldb) || is.null(ylda)) {
    stop(
      "Essential arguments (mortality/disability baseline/adjusted) are missing."
    )
  }

  validate_inputs(iso3c, mub, mua, yldb, ylda)

  # Dimensions and timeline
  years <- (date_year + 1) + 0:(100 + scaleup)
  n_years <- length(years)
  ages <- age_groups$all_ages
  n_ages <- length(ages)
  sex_count <- length(sex_labels)

  # WPP parameters
  wpp_params <- wpp.pars[[iso3c]]
  wpp_mortality_target <- wpp_params$mx[1:n_years, , ]
  wpp_asfr <- wpp_params$asfr[1:n_years, ]
  wpp_migration <- wpp_params$mig[1:n_years, , ]
  wpp_srb <- wpp_params$srb[1:n_years]
  base_population <- wpp_params$base.pop

  # Convert cohort->period for all-cause mortality and YLDs (baseline vs adjusted)
  mortality_period_baseline <- array(NA_real_, dim = dim(wpp_mortality_target))
  mortality_period_adjusted <- array(NA_real_, dim = dim(wpp_mortality_target))
  disability_period_baseline <- array(NA_real_, dim = dim(wpp_mortality_target))
  disability_period_adjusted <- array(NA_real_, dim = dim(wpp_mortality_target))

  for (s in 1:sex_count) {
    mort_pair <- cohort_to_period(
      baseline_matrix = mub[,, s],
      adjusted_matrix = mua[,, s]
    )
    mortality_period_adjusted[, s, ] <- mort_pair$adjusted
    mortality_period_baseline[, s, ] <- mort_pair$baseline

    yld_pair <- cohort_to_period(
      baseline_matrix = yldb[,, s],
      adjusted_matrix = ylda[,, s]
    )
    disability_period_adjusted[, s, ] <- yld_pair$adjusted
    disability_period_baseline[, s, ] <- yld_pair$baseline
  }

  # Fertility impacts
  fertility_df <- fertility_impact_for_iso(
    iso3c,
    c1t,
    scaleup,
    calc.fert,
    fert.targets
  )

  # Initialize projection arrays (period, sex, age for population and deaths; period, sex for births)
  population_baseline <- array(NA_real_, dim = c(n_years, sex_count, n_ages))
  population_adjusted <- array(NA_real_, dim = c(n_years, sex_count, n_ages))
  deaths_baseline <- array(NA_real_, dim = c(n_years, sex_count, n_ages))
  deaths_adjusted <- array(NA_real_, dim = c(n_years, sex_count, n_ages))
  continuous_deaths_baseline <- array(
    NA_real_,
    dim = c(n_years, sex_count, n_ages)
  )
  continuous_deaths_adjusted <- array(
    NA_real_,
    dim = c(n_years, sex_count, n_ages)
  )
  births_baseline <- array(NA_real_, dim = c(n_years, sex_count))
  births_adjusted <- array(NA_real_, dim = c(n_years, sex_count))
  mortality_rates_baseline <- array(
    NA_real_,
    dim = c(n_years, sex_count, n_ages)
  )
  mortality_rates_adjusted <- array(
    NA_real_,
    dim = c(n_years, sex_count, n_ages)
  )

  population_baseline[1, , ] <- base_population
  population_adjusted[1, , ] <- base_population

  # Scale to WPP
  mortality_scale_to_wpp <- wpp_mortality_target / mortality_period_baseline

  # Initial year mortality (match FC behavior)
  if (calib == "Fertility Migration, Base Year") {
    mortality_rates_baseline[1, , ] <- mortality_period_baseline[1, , ]
    mortality_rates_adjusted[1, , ] <- mortality_period_baseline[1, , ]
  } else {
    mortality_rates_baseline[1, , ] <- mortality_scale_to_wpp[1, , ] *
      mortality_period_baseline[1, , ]
    # In FC script, mort.a[1,,] and mort.b[1,,] are identical. Replicating that.
    mortality_rates_adjusted[1, , ] <- mortality_scale_to_wpp[1, , ] *
      mortality_period_baseline[1, , ]
  }

  # ASFR output array (index 1 = Baseline ASFR; index 2 = Adjusted ASFR)
  asfr <- array(NA_real_, dim = c(2, ncol(wpp_asfr), n_years - 1))

  # Projection loop
  for (t in 2:n_years) {
    current_scale <- mortality_scale_to_wpp[t, , ]
    current_asfr <- wpp_asfr[t, ]
    current_srb <- wpp_srb[t]
    current_migration <- wpp_migration[t, , ]

    if (calib == "Fertility Migration, Base Year" && t >= 2) {
      current_scale[,] <- 1
      current_asfr <- wpp_asfr[2, ]
      current_srb <- wpp_srb[2]
      current_migration[,] <- wpp_migration[2, , ]
    }
    # if (calib == "Mortality Fertility Migration, Base Year" && t >= 2) {
    #   current_scale <- mortality_scale_to_wpp[1, , ]
    #   current_asfr <- wpp_asfr[2, ]
    #   current_srb <- wpp_srb[2]
    #   current_migration[,] <- wpp_migration[2, , ]
    # }
    if (calib == "Mortality Fertility Migration, Base Year" && t >= 2) {
      current_scale <- mortality_scale_to_wpp[t, , ]
      current_asfr <- wpp_asfr[t, ]
      current_srb <- wpp_srb[t]
      current_migration[,] <- wpp_migration[t, , ]
    }

    # Fertility adjustments
    baseline_asfr <- current_asfr
    if (t <= scaleup) {
      lambda <- fertility_df$lambda[t]
      gamma <- fertility_df$gam[t]
    } else {
      lambda <- fertility_df$lambda[scaleup + 1]
      gamma <- fertility_df$gam[scaleup + 1]
    }
    adjusted_asfr <- baseline_asfr *
      (1 - lambda) +
      baseline_asfr * lambda * (1 - gamma)

    # Mortality for this year
    current_mortality_baseline <- current_scale *
      mortality_period_baseline[t, , ]
    current_mortality_adjusted <- current_scale *
      mortality_period_adjusted[t, , ]
    current_mortality_baseline[current_mortality_baseline > 1] <- 1 - 1e-6
    current_mortality_adjusted[current_mortality_adjusted > 1] <- 1 - 1e-6

    mortality_rates_baseline[t, , ] <- current_mortality_baseline
    mortality_rates_adjusted[t, , ] <- current_mortality_adjusted

    # Survival probabilities
    survival_probs_baseline <- t(apply(
      current_mortality_baseline,
      1,
      calculate_survival_probabilities
    ))
    survival_probs_adjusted <- t(apply(
      current_mortality_adjusted,
      1,
      calculate_survival_probabilities
    ))

    # Baseline projection step
    baseline_step <- project_one_year(
      current_population = population_baseline[t - 1, , ],
      mortality_rates = current_mortality_baseline,
      survival_probabilities = survival_probs_baseline,
      fertility_rates_by_age = baseline_asfr,
      sex_ratio_at_birth = current_srb,
      migration_rates = current_migration,
      n_ages = n_ages
    )
    population_baseline[t, , ] <- baseline_step$population
    deaths_baseline[t - 1, , ] <- baseline_step$deaths
    continuous_deaths_baseline[t - 1, , ] <- baseline_step$continuous_deaths
    births_baseline[t - 1, ] <- baseline_step$births

    # Adjusted projection step
    adjusted_step <- project_one_year(
      current_population = population_adjusted[t - 1, , ],
      mortality_rates = current_mortality_adjusted,
      survival_probabilities = survival_probs_adjusted,
      fertility_rates_by_age = adjusted_asfr,
      sex_ratio_at_birth = current_srb,
      migration_rates = current_migration,
      n_ages = n_ages
    )
    population_adjusted[t, , ] <- adjusted_step$population
    deaths_adjusted[t - 1, , ] <- adjusted_step$deaths
    continuous_deaths_adjusted[t - 1, , ] <- adjusted_step$continuous_deaths
    births_adjusted[t - 1, ] <- adjusted_step$births

    # Store ASFR (index 1 = baseline, 2 = adjusted)
    asfr[1, , t - 1] <- baseline_asfr
    asfr[2, , t - 1] <- adjusted_asfr
  }

  # Cause-specific deaths (if provided)
  n_causes <- length(deathb)
  if (n_causes > 0) {
    cause_mortality_period_baseline <- array(
      NA_real_,
      dim = c(n_causes, dim(wpp_mortality_target))
    )
    cause_mortality_period_adjusted <- array(
      NA_real_,
      dim = c(n_causes, dim(wpp_mortality_target))
    )

    for (c in 1:n_causes) {
      for (s in 1:sex_count) {
        # Note: FC script calls get.period(deathb, deatha), so mata=deathb, matb=deatha
        cp <- cohort_to_period(
          baseline_matrix = deathb[[c]][,, s],
          adjusted_matrix = deatha[[c]][,, s]
        )
        # FC script assigns period.mat[1] (from deathb) to mxbase.cse
        # and period.mat[2] (from deatha) to mxadj.cse. This is correct.
        cause_mortality_period_baseline[c, , s, ] <- cp$baseline
        cause_mortality_period_adjusted[c, , s, ] <- cp$adjusted
      }
    }

    deaths_cause_specific_baseline <- array(
      0,
      dim = dim(cause_mortality_period_baseline)
    )
    deaths_cause_specific_adjusted <- array(
      0,
      dim = dim(cause_mortality_period_adjusted)
    )

    for (c in 1:n_causes) {
      ratio_baseline <- cause_mortality_period_baseline[c, , , ] /
        mortality_rates_baseline
      ratio_adjusted <- cause_mortality_period_adjusted[c, , , ] /
        mortality_rates_adjusted
      deaths_cause_specific_baseline[c, , , ] <- ratio_baseline *
        deaths_baseline
      deaths_cause_specific_adjusted[c, , , ] <- ratio_adjusted *
        deaths_adjusted
    }
  }

  # ====
  # Period outputs (long) and fertility/contraception metrics
  # Flows (deaths, continuous deaths, births) exist for n_years - 1 intervals
  valid_years <- 1:(n_years - 1)

  # ====

  pop_adjusted_df <- array3d_to_long(
    population_adjusted,
    "Population",
    "Adjusted",
    date_year = date_year,
    sex_labels = sex_labels
  )
  pop_baseline_df <- array3d_to_long(
    population_baseline,
    "Population",
    "Baseline",
    date_year = date_year,
    sex_labels = sex_labels
  )
  deaths_adjusted_df <- array3d_to_long(
    deaths_adjusted[valid_years, , , drop = FALSE],
    "Deaths",
    "Adjusted",
    date_year = date_year,
    sex_labels = sex_labels
  )
  deaths_baseline_df <- array3d_to_long(
    deaths_baseline[valid_years, , , drop = FALSE],
    "Deaths",
    "Baseline",
    date_year = date_year,
    sex_labels = sex_labels
  )
  mx_adjusted_df <- array3d_to_long(
    mortality_rates_adjusted,
    "mx",
    "Adjusted",
    date_year = date_year,
    sex_labels = sex_labels
  )
  mx_baseline_df <- array3d_to_long(
    mortality_rates_baseline,
    "mx",
    "Baseline",
    date_year = date_year,
    sex_labels = sex_labels
  )
  mx.in_adjusted_df <- array3d_to_long(
    mortality_period_adjusted,
    "mx",
    "Adjusted",
    date_year = date_year,
    sex_labels = sex_labels
  )
  mx.in_baseline_df <- array3d_to_long(
    mortality_period_baseline,
    "mx",
    "Baseline",
    date_year = date_year,
    sex_labels = sex_labels
  )


  # Births (2D arrays -> long)
  births_baseline_df <- {
    ny <- n_years - 1
    ns <- ncol(births_baseline)
    grid <- expand.grid(sex = seq_len(ns), year = seq_len(ny))
    val <- as.vector(births_baseline[valid_years, , drop = FALSE]) # [year, sex], year fastest
    data.frame(
      age = 0L,
      year = date_year + grid$year,
      sex = factor(grid$sex, labels = sex_labels),
      value = val,
      metric = "Births",
      source = "Baseline",
      cause = "All",
      stringsAsFactors = FALSE
    )
  }

  births_adjusted_df <- {
    ny <- n_years - 1
    ns <- ncol(births_adjusted)
    grid <- expand.grid(sex = seq_len(ns), year = seq_len(ny))
    val <- as.vector(births_adjusted[valid_years, , drop = FALSE]) # [year, sex], year fastest
    data.frame(
      age = 0L,
      year = date_year + grid$year,
      sex = factor(grid$sex, labels = sex_labels),
      value = val,
      metric = "Births",
      source = "Adjusted",
      cause = "All",
      stringsAsFactors = FALSE
    )
  }

  period_core1 <- bind_rows(
    pop_adjusted_df,
    pop_baseline_df,
    deaths_adjusted_df,
    deaths_baseline_df,
    mx_adjusted_df,
    mx_baseline_df,
    births_baseline_df,
    births_adjusted_df
  ) %>%
    mutate(calibration = calib) %>%
    select(age, year, sex, value, metric, source, cause, calibration)

  period_core2 <- bind_rows(
    mx.in_adjusted_df,
    mx.in_baseline_df
  ) %>%
    mutate(calibration = "Unscaled") %>%
    select(age, year, sex, value, metric, source, cause, calibration)
  
  period_core <- rbind(period_core1, period_core2)
  
  # Cause-specific period long (if any)
  if (n_causes > 0) {
    cause_names <- names(deathb)
    deaths_a_cause_df <- array4d_to_long(
      deaths_cause_specific_adjusted[, valid_years, , , drop = FALSE],
      "Deaths",
      "Adjusted",
      cause_names,
      date_year,
      sex_labels
    )
    deaths_b_cause_df <- array4d_to_long(
      deaths_cause_specific_baseline[, valid_years, , , drop = FALSE],
      "Deaths",
      "Baseline",
      cause_names,
      date_year,
      sex_labels
    )
    period_core <- bind_rows(
      period_core,
      deaths_a_cause_df %>% mutate(calibration = calib),
      deaths_b_cause_df %>% mutate(calibration = calib)
    )
  }

  # Fertility: pregnancies and stillbirths (match FC logic)
  preg_still <- pregnancies.stillbirths
  if (iso3c == "ZNZ") {
    preg_still <- preg_still %>%
      filter(iso3 == "TZA") %>%
      mutate(iso3 = "ZNZ", iso = "ZNZ", country = "Zanzibar")
  }
  if (iso3c %in% india_states) {
    preg_still <- pregnancies.stillbirths %>%
      filter(iso3 == "IND") %>%
      mutate(iso3 = iso3c, iso = iso3c, country = "India")
  }

  births_summary <- bind_rows(births_baseline_df, births_adjusted_df) %>%
    group_by(year, source) %>%
    summarise(births = sum(value), .groups = "drop")

  births_preg_still_df <- births_summary %>%
    bind_cols(preg_still %>% filter(iso3 == iso3c) %>% slice(1)) %>%
    mutate(
      stillbirths = round(births * stillbirth.rate),
      pregnancies = round(births * birthtopreg.scale),
      births = round(births),
      calibration = calib
    ) %>%
    select(
      iso3,
      year,
      source,
      calibration,
      births,
      stillbirths,
      pregnancies
    ) %>%
    pivot_longer(
      cols = c(births, stillbirths, pregnancies),
      names_to = "metric",
      values_to = "value"
    ) %>%
    rename(group = source)

  # Contraception metrics (match FC construction)
  female_15_49 <- period_core %>%
    filter(sex == "Female", age > 14, age < 50, metric == "Population") %>%
    group_by(year, source) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = source, values_from = value)

  contraception_metrics_df <- fertility_df %>%
    select(-any_of(c("rho", "gam"))) %>%
    left_join(female_15_49, by = "year") %>%
    mutate(
      Baseline = Baseline * lambda * c0,
      Adjusted = Adjusted * lambda * c1,
      lambda = NULL,
      c0 = NULL,
      c1 = NULL,
      country = NULL
    ) %>%
    rename(iso3 = iso) %>%
    pivot_longer(
      cols = -c(iso3, year, Adjusted, Baseline),
      names_to = "metric",
      values_to = "proportion"
    ) %>%
    mutate(
      Baseline = round(Baseline * proportion),
      Adjusted = round(Adjusted * proportion)
    ) %>%
    select(-proportion) %>%
    pivot_longer(
      cols = c(Baseline, Adjusted),
      names_to = "group",
      values_to = "val"
    ) %>%
    arrange(group, metric, year) %>%
    select(iso3, metric, group, year, val) %>%
    mutate(calibration = calib)

  contraception_metrics_df <- contraception_metrics_df %>%
    bind_rows(
      births_preg_still_df %>%
        filter(year %in% unique(fertility_df$year)) %>%
        select(iso3, metric, group, year, val = value, calibration)
    ) %>%
    arrange(group, metric, year) %>%
    rename(value = val, source = group) %>%
    mutate(age = NA_integer_, sex = "Female", cause = "pregnancy") %>%
    select(age, year, sex, value, metric, source, cause, calibration)

  # Final period.df (FC-compatible)
  period.df <- bind_rows(period_core, contraception_metrics_df) %>%
    mutate(age = as.numeric(age)) %>%
    arrange(year, age, sex) %>%
    tibble::as_tibble() # Ensure output is a tibble

  # Actual WPP
  actual.wpp <- wpp2024.actual %>% filter(iso3 == iso3c)

  # Return with the exact same keys as project_pop_fc.R plus new mortality outputs
  list(
    # Existing outputs
    period.df = period.df,
    actual.wpp = actual.wpp,
    asfr = asfr
  )
}
