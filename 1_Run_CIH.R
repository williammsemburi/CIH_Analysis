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

# Load shared parameters: iso3s, no_vmc_countries, idlist_base,
# core_data, num_cores, thread settings
source("../CIH/config.R")

message(paste("Total countries to process:", length(iso3s)))
message(paste("IDs per country:",            length(idlist_base)))
message(paste("Using",                       num_cores, "cores"))

# =============================================================================
# STAGE 1 — Parallel per-ID run
# =============================================================================

# Helper to initialise (or reinitialise) the parallel cluster
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

for (w in seq_along(iso3s)) {

  cat("\n")
  message("##################################################")
  start_time <- Sys.time()

  # Restart the cluster every 10 countries to free memory
  if (w %% 10 == 1 && w > 1) {
    message(paste("=== Restarting cluster after", w - 1, "countries ==="))
    message("=== 3-minute cooling break ===")
    stopCluster(cl)
    Sys.sleep(180)
    cl <- start_cluster(num_cores)
  }

  iso3 <- iso3s[w]

  # Drop the VMC intervention for countries where it is not applicable
  idlist <- if (iso3 %in% no_vmc_countries) setdiff(idlist_base, 43L) else idlist_base

  message(paste0("ISO3C: ", iso3, " — ", w, " of ", length(iso3s),
                 " (", round(100 * w / length(iso3s), 2), "%)"))

  # Pre-load shared offline data for this country
  message("  Pre-loading shared data...")
  offline_file <- paste0(core_data, "/db/offline_", iso3, ".rda")
  load(offline_file)  # loads: coverages_chosen, cln

  coverages_chosen_full <- coverages_chosen |>
    filter(!id == 223) |>
    mutate(target = baseline + 0.1, target = ifelse(target > 1, 1, target))
  cln_full <- cln
  rm(coverages_chosen, cln, offline_file)

  # Export country-level data to workers
  clusterExport(cl, c("core_data", "iso3", "idlist",
                      "coverages_chosen_full", "cln_full"),
                envir = environment())

  message(paste("  Processing", length(idlist), "IDs in parallel..."))

  res <- foreach(
    theid        = idlist,
    .packages    = c("data.table", "dplyr", "fst"),
    .errorhandling = "stop"
  ) %dopar% {

    coverages_chosen <- coverages_chosen_full |> filter(id %in% theid)
    cln              <- cln_full

    # Assign to global env so sourced functions can find them
    assign("theid",            theid,            envir = .GlobalEnv)
    assign("coverages_chosen", coverages_chosen, envir = .GlobalEnv)
    assign("cln",              cln,              envir = .GlobalEnv)

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
  message(paste("✓ Finished:", iso3, "| Time:", round(elapsed_time, 2), "mins"))

  timing_results[[iso3]] <- data.frame(
    iso3          = iso3,
    start_time    = start_time,
    end_time      = end_time,
    duration_mins = round(elapsed_time, 2),
    condition     = "Complete",
    stringsAsFactors = FALSE
  )

  rm(res, coverages_chosen_full, cln_full, start_time, end_time, elapsed_time)
  rm(list = ls(pattern = "^full\\."))
  gc(verbose = FALSE)
}

stopCluster(cl)

# =============================================================================
# STAGE 2 — Sequential full-package run (all IDs, coverage target = 1)
# =============================================================================

source("config/init.R")
source("scripts/load_model.R")
source("scripts/get_pins.R")
load("data/taxonomy/taxGroups.rda")
load("data/wb/2024/country2024.rda")

# Project population baseline (required before the per-country loop)
suppressMessages(source("../CIH/fxn/project_pop_cih.R"))

message("##################################################")
message("Processing all countries — full package...")

for (w in seq_along(iso3s)) {

  cat("\n")
  message("##################################################")
  start_time <- Sys.time()

  iso3   <- iso3s[w]
  idlist <- if (iso3 %in% no_vmc_countries) setdiff(idlist_base, 43L) else idlist_base
  theid  <- idlist

  message(paste0("Full ISO3C: ", iso3, " — ", w, " of ", length(iso3s),
                 " (", round(100 * w / length(iso3s), 2), "%)"))

  # Load and prepare data (target = 1 for full-package scenario)
  offline_file <- paste0(core_data, "/db/offline_", iso3, ".rda")
  load(offline_file)  # loads: coverages_chosen, cln

  coverages_chosen <- coverages_chosen |>
    filter(!id == 223) |>
    mutate(target = 1) |>
    filter(id %in% theid)

  suppressMessages(source("../CIH/fxn/run_cih_funcs.R"))

  # Subset and annotate mortality/population output
  full.df <- period.df |>
    dplyr::filter(metric %in% c("mx", "Population") & year <= 2050) |>
    dplyr::mutate(iso3c = iso3) |>
    dplyr::arrange(source, metric, cause, year, sex, age)

  write_fst(full.df,    paste0("../CIH/output/fulldf/full.df_",      iso3, ".fst"))
  write_fst(full.q70,   paste0("../CIH/output/fullq70/full.q70_",   iso3, ".fst"))
  write_fst(full.costs, paste0("../CIH/output/fullcost/full.costs_", iso3, ".fst"))

  rm(full.df, full.q70, full.costs, period.df, coverages_chosen, cln)
  end_time     <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  message(paste("✓ Finished:", iso3, "| Time:", round(elapsed_time, 2), "mins"))

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
message(paste("Total time:          ", round(sum(timing_results_df$duration_mins), 2), "mins"))
message(paste("Average per country: ", round(mean(timing_results_df$duration_mins), 2), "mins"))
