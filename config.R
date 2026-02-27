# =============================================================================
# config.R
# Shared parameters for the CIH analysis pipeline.
# Sourced by 1_Run_CIH.R and 2_CIH_Paper_Plots.R.
#
# NOTE: This script assumes the working directory has already been set to
# FairChoices-v3/ (as done at the top of each calling script).
# =============================================================================

library(data.table)
library(dplyr)

# -----------------------------------------------------------------------------
# Compute resources
# -----------------------------------------------------------------------------
core_data <- "data"
num_cores  <- 0.5 * parallel::detectCores()

# Restrict multithreading to avoid conflicts with foreach/doParallel
Sys.setenv(OMP_NUM_THREADS  = 1)
Sys.setenv(BLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS  = 1)

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
# Intervention ID list
# Base set of intervention IDs run for every country.
# For no_vmc_countries, ID 43 is excluded (see scripts).
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

# -----------------------------------------------------------------------------
# Region colors (used by 2_CIH_Paper_Plots.R for all time-series charts)
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
