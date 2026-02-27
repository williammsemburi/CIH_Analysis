# =============================================================================
# config.R
# Shared parameters for the CIH analysis pipeline.
# Sourced by 1_Run_CIH.R and 2_CIH_Paper_Plots.R.
#
# NOTE: This script assumes the working directory has already been set to
# FairChoices-v3/ (as done at the top of each calling script).
#
# NOTE: idlist_base is NOT defined here because the two scripts derive it
# differently:
#   - 1_Run_CIH.R         uses a hardcoded integer vector
#   - 2_CIH_Paper_Plots.R reads it from Intervention_map.csv (Core only)
# Each script defines its own idlist_base after sourcing this file.
# =============================================================================

library(data.table)
library(dplyr)

# -----------------------------------------------------------------------------
# Country lists
# -----------------------------------------------------------------------------

# Full ISO3 list from the region mapping input
iso3s <- fread("../CIH/input/country_region_mapping.csv") |> pull(ISO3)

# Countries excluded due to missing GBD data
no_gbd_data <- c("COK", "NIU", "TKL", "PRK", "TWN", "VEN", "CPV", "PSE")

# Countries excluded from the VMC (voluntary medical male circumcision)
# intervention — ID 43 is dropped for these countries
no_vmc_countries <- c(
  "BHS", "BRN", "CIV", "COD", "COG", "CPV", "CZE", "FSM",
  "GMB", "IRN", "KGZ", "KOR", "LAO", "LCA", "MKD",
  "PRK", "PSE", "STP", "SVK", "SWZ", "SYR", "TUR", "TWN",
  "VCT", "VEN", "VIR", "YEM"
)

# Remove countries with no GBD data from the processing list
iso3s <- iso3s[!iso3s %in% no_gbd_data]

# -----------------------------------------------------------------------------
# Region colors
# Used by 2_CIH_Paper_Plots.R for time-series and Figure Y charts.
# Note: the regional overview map (intro.map) uses a separate color palette
# defined within that function.
# -----------------------------------------------------------------------------
region_colors <- list(
  "Sub-Saharan Africa"                 = "#6FA8DC",
  "India"                              = "#6AA84F",
  "World"                              = "#000000",
  "China"                              = "#BF9000",
  "United States"                      = "#F1C232",
  "North Atlantic"                     = "#6DD8E8",
  "Central and Eastern Europe"         = "#E06666",
  "Central Asia"                       = "#8E7CC3",
  "Latin America and Caribbean"        = "#F6B26B",
  "Middle East and North Africa"       = "#93C47D",
  "Western Pacific and Southeast Asia" = "#C27BA0"
)
