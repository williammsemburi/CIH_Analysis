# =============================================================================
# 1_Run_CIH.R
#
# CIH analysis pipeline — runs in two stages:
#
#   Stage 1 (parallel):  For each country and intervention ID, runs the CIH
#                        model and writes per-ID output files to:
#                          CIH/output/fullq70/ids/
#                          CIH/output/fullcost/ids/
#
#   Stage 2 (sequential): For each country, runs the full-package model
#                         (all IDs together, coverage target = 1) and writes
#                         country-level output files to:
#                           CIH/output/fulldf/
#                           CIH/output/fullq70/
#                           CIH/output/fullcost/
#
# Timing results are saved to CIH/output/timing_results.fst.
#
# Dependencies:
#   - FairChoices-v3/ repo must be a sibling directory of CIH/
#   - CIH/input/ must contain country_region_mapping.csv
#   - FairChoices-v3/data/db/ must contain offline_<ISO3>.rda files
# =============================================================================

rm(list = ls())

# Set working directory to FairChoices-v3 (sibling of CIH/)
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("FairChoices-v3")

library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(data.table)
library(fst)

# Load shared parameters: iso3s, no_vmc_countries, region_colors
# iso3s has already had no_gbd_data countries removed
source("../CIH/config.R")

# -----------------------------------------------------------------------------
# Compute resources
# -----------------------------------------------------------------------------
core_data <- "data"
num_cores  <- 0.5 * parallel::detectCores()
Sys.setenv(OMP_NUM_THREADS  = 1)
Sys.setenv(BLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS  = 1)

# -----------------------------------------------------------------------------
# Intervention ID list
# Base set of intervention IDs run for every country.
# For no_vmc_countries, ID 43 (VMC) is excluded — see loop below.
# -----------------------------------------------------------------------------
idlist_base <- c(
   27L,  28L,  29L,  30L,  31L,  32L,  33L,  34L,  36L,  37L,  38L,  39L,
   40L,  41L,  42L,  43L,  44L,  45L,  46L,  47L,  54L,  56L,  57L,  58L,  59L,
   60L,  61L,  62L,  63L,  64L,  65L,  66L,  67L,  68L,  69L,  70L,  71L,  72L,
   73L,  74L,  82L,  83L,  84L,  85L,  86L,  87L,  88L,  89L,  90L,  91L,  92L,
   93L,  94L,  95L,  96L,  97L,  98L,  99L, 100L, 101L, 102L, 103L, 104L,
  105L, 106L, 107L, 108L, 109L, 110L, 111L, 112L, 113L, 114L, 115L,
  116L, 117L, 118L, 119L, 120L, 121L, 122L, 123L, 124L, 125L, 126L,
  127L, 128L, 129L, 130L, 131L, 132L, 133L, 134L, 135L, 136L, 137L,
  138L, 139L, 140L, 141L, 142L, 143L, 144L, 145L, 146L, 147L, 148L,
  149L, 150L, 151L, 152L, 153L, 154L, 156L, 157L, 158L, 159L, 160L,
  161L, 163L, 164L, 165L, 168L, 169L, 178L, 179L, 180L, 190L, 200L,
  206L, 224L, 226L, 227L, 228L, 229L
)

message(paste("Total countries to process:", length(iso3s)))
message(paste("IDs per country:",            length(idlist_base)))
message(paste("Using",                       num_cores, "cores for parallel processing"))

# =============================================================================
# STAGE 1 — Parallel per-ID run
# =============================================================================

# Helper to start (or restart) the parallel cluster and load required
# packages and scripts on all workers
start_cluster <- function(n_cores) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    library(fst)
    source("config/init.R")
    source("scripts/load_model.R")
    source("scripts/get_pins.R")
    load("data/taxonomy/taxGroups.rda")
    load("data/wb/2024/country2024.rda")
  })
  cl
}

cl             <- start_cluster(num_cores)
timing_results <- list()

for (w in 1:length(iso3s)) {
  cat("\n")
  print("##################################################")

  start_time <- Sys.time()

  # Restart the cluster every 10 countries to free memory
  if (w %% 10 == 1 && w > 1) {
    message(paste("\n=== RESTARTING CLUSTER after", w-1, "countries ==="))
    message("=== 3 minute cooling break ===\n")

    stopCluster(cl)
    Sys.sleep(180)  # 3 min cooling

    cl <- start_cluster(num_cores)
  }

  iso3 <- iso3s[w]

  # Drop the VMC intervention (ID 43) for countries where it is not applicable
  if (iso3 %in% no_vmc_countries) {
    idlist <- setdiff(idlist_base, 43L)
  } else {
    idlist <- idlist_base
  }

  message(paste0("ISO3C: ", iso3, " - ", w, " of ", length(iso3s),
                 " (", round(100*w/length(iso3s), 2), "%)"))

  # Pre-load shared offline data for this country
  message("  Pre-loading shared data...")
  offline_file <- paste0(core_data, "/db/offline_", iso3, ".rda")
  load(offline_file)  # loads: coverages_chosen, cln

  coverages_chosen_full <- coverages_chosen %>%
    filter(!id == 223) %>%
    mutate(target = baseline + 0.1, target = ifelse(target > 1, 1, target))
  cln_full <- cln
  rm(coverages_chosen, cln, offline_file)

  message(paste("  Processing", length(idlist), "IDs in parallel..."))

  clusterExport(cl, c("core_data", "iso3", "idlist",
                      "coverages_chosen_full", "cln_full"),
                envir = environment())

  res <- foreach(
    theid = idlist,
    .packages = c("data.table", "dplyr", "fst"),
    .errorhandling = "stop"
  ) %dopar% {

    coverages_chosen <- coverages_chosen_full %>% filter(id %in% theid)
    cln              <- cln_full

    assign("theid",            theid,            envir = .GlobalEnv)
    assign("coverages_chosen", coverages_chosen, envir = .GlobalEnv)
    assign("cln",              cln,              envir = .GlobalEnv)
    cat("Worker", Sys.getpid(), "- ISO", iso3, "- ID", theid, "\n")

    suppressWarnings(suppressMessages(
      capture.output(
        source("../CIH/fxn/run_cih_funcs.R"),
        file = nullfile()
      )
    ))

    write_fst(full.q70,   paste0("../CIH/output/fullq70/ids/full.q70_",    iso3, "_", theid, ".fst"))
    write_fst(full.costs, paste0("../CIH/output/fullcost/ids/full.costs_", iso3, "_", theid, ".fst"))

    rm(full.q70, full.costs, coverages_chosen, cln)
    gc(verbose = FALSE)

    list(iso3 = iso3, id = theid, status = "done")
  }

  end_time     <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  message(paste("✓ Finished:", iso3, "| Time taken:", round(elapsed_time, 2), "mins"))

  timing_results[[iso3]] <- data.frame(
    iso3          = iso3,
    start_time    = start_time,
    end_time      = end_time,
    duration_mins = round(elapsed_time, 2),
    condition     = "Complete",
    stringsAsFactors = FALSE
  )

  # Cleanup before next iteration
  rm(res, coverages_chosen_full, cln_full, start_time, end_time, elapsed_time)
  rm(list = ls(pattern = "^full\\."))
  gc(verbose = FALSE)
}

stopCluster(cl)

# =============================================================================
# STAGE 2 — Sequential full-package run (all IDs, coverage target = 1)
# =============================================================================

library(fst)
source("config/init.R")
source("scripts/load_model.R")
source("scripts/get_pins.R")
load("data/taxonomy/taxGroups.rda")
load("data/wb/2024/country2024.rda")

# Project population baseline (required before the per-country loop)
suppressMessages(
  source("../CIH/fxn/project_pop_cih.R")
)

# Re-declare compute settings (environment may have changed after Stage 1)
core_data <- "data"
num_cores  <- 0.5 * parallel::detectCores()
Sys.setenv(OMP_NUM_THREADS  = 1)
Sys.setenv(BLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS  = 1)

# Re-read country list (iso3s may have been modified in Stage 1 environment)
iso3s       <- fread("../CIH/input/country_region_mapping.csv") |> pull(ISO3)
no_gbd_data <- c("COK", "NIU", "TKL", "PRK", "TWN", "VEN", "CPV", "PSE")

no_vmc_countries <- c(
  "BHS", "BRN", "CIV", "COD", "COG", "CPV", "CZE", "FSM",
  "GMB", "IRN", "KGZ", "KOR", "LAO", "LCA", "MKD",
  "PRK", "PSE", "STP", "SVK", "SWZ", "SYR", "TUR", "TWN",
  "VCT", "VEN", "VIR", "YEM"
)

iso3s <- iso3s[!iso3s %in% no_gbd_data]

idlist_base <- c(
   27L,  28L,  29L,  30L,  31L,  32L,  33L,  34L,  36L,  37L,  38L,  39L,
   40L,  41L,  42L,  43L,  44L,  45L,  46L,  47L,  54L,  56L,  57L,  58L,  59L,
   60L,  61L,  62L,  63L,  64L,  65L,  66L,  67L,  68L,  69L,  70L,  71L,  72L,
   73L,  74L,  82L,  83L,  84L,  85L,  86L,  87L,  88L,  89L,  90L,  91L,  92L,
   93L,  94L,  95L,  96L,  97L,  98L,  99L, 100L, 101L, 102L, 103L, 104L,
  105L, 106L, 107L, 108L, 109L, 110L, 111L, 112L, 113L, 114L, 115L,
  116L, 117L, 118L, 119L, 120L, 121L, 122L, 123L, 124L, 125L, 126L,
  127L, 128L, 129L, 130L, 131L, 132L, 133L, 134L, 135L, 136L, 137L,
  138L, 139L, 140L, 141L, 142L, 143L, 144L, 145L, 146L, 147L, 148L,
  149L, 150L, 151L, 152L, 153L, 154L, 156L, 157L, 158L, 159L, 160L,
  161L, 163L, 164L, 165L, 168L, 169L, 178L, 179L, 180L, 190L, 200L,
  206L, 224L, 226L, 227L, 228L, 229L
)

print("##################################################")
message(paste("  Processing all countries for full package..."))

timing_results <- list()

for (w in 1:length(iso3s)) {
  cat("\n")
  print("##################################################")

  start_time <- Sys.time()
  iso3 <- iso3s[w]

  message(paste0("Full ISO3C: ", iso3, " - ", w, " of ", length(iso3s),
                 " (", round(100*w/length(iso3s), 2), "%)"))

  # Drop VMC intervention for countries where it is not applicable
  if (iso3 %in% no_vmc_countries) {
    idlist <- setdiff(idlist_base, 43L)
  } else {
    idlist <- idlist_base
  }

  theid <- idlist

  # Load and prepare data (target = 1 for full-package scenario)
  offline_file <- paste0(core_data, "/db/offline_", iso3, ".rda")
  load(offline_file)  # loads: coverages_chosen, cln

  coverages_chosen <- coverages_chosen %>%
    filter(!id == 223) %>%
    mutate(target = 1) %>%
    filter(id %in% theid)

  suppressMessages(
    source("../CIH/fxn/run_cih_funcs.R")
  )

  # Subset mortality/population output and tag with country code
  full.df <- period.df |>
    dplyr::filter(metric %in% c("mx", "Population") & year <= 2050) |>
    dplyr::mutate(iso3c = iso3) |>
    dplyr::arrange(source, metric, cause, year, sex, age)

  write_fst(full.df,    paste0("../CIH/output/fulldf/full.df_",      iso3, ".fst"))
  write_fst(full.q70,   paste0("../CIH/output/fullq70/full.q70_",   iso3, ".fst"))
  write_fst(full.costs, paste0("../CIH/output/fullcost/full.costs_", iso3, ".fst"))

  # Cleanup
  rm(full.df, full.q70, full.costs, period.df, coverages_chosen, cln)
  end_time     <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  message(paste("✓ Finished:", iso3, "| Time taken:", round(elapsed_time, 2), "mins"))

  timing_results[[paste(iso3, 2)]] <- data.frame(
    iso3          = paste(iso3, 2),
    start_time    = start_time,
    end_time      = end_time,
    duration_mins = round(elapsed_time, 2),
    condition     = "Complete",
    stringsAsFactors = FALSE
  )

  rm(start_time, end_time, elapsed_time)
}

# =============================================================================
# Save timing summary
# =============================================================================

timing_results_df <- data.table::rbindlist(timing_results, fill = TRUE)
print(timing_results_df)
write_fst(timing_results_df, "../CIH/output/timing_results.fst")

message("\n========== ALL COUNTRIES PROCESSED ==========")
message(paste("Total time:",         round(sum(timing_results_df$duration_mins),  2), "mins"))
message(paste("Average per country:", round(mean(timing_results_df$duration_mins), 2), "mins"))
