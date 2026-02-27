rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##########################################################
# Install required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyr, dplyr, data.table, countrycode, devtools, readxl, stringr)

# # # Creating the WPP file
#  if(!require("wpp2024")){
#    options(timeout = 600)
#    install_github("PPgp/wpp2024")
#  }
#  library(wpp2024)
# # ##########################################################
# #
# # ##########################################################
# # # Baseline demography from WPP2024
# # ##########################################################
# #
# # # # Population numbers
#  data(popAge1dt)
#  data(popprojAge1dt)
# # #
# # # # Mortality rates
#  data(mx1dt)
# # #
# # # # Fertility rates and sex-ratio-at-birth
#  data(tfr1dt)
#  data(tfrproj1dt)
#  data(percentASFR1dt)
#  data(sexRatio1dt)
# # #
# # # # location names
#  data(UNlocations)
# # #
# # # # Other data
#  data(misc1dt)
#  data(miscproj1dt)
# # #
# wpp2024.all <- list(popprojAge1dt = rbind(popAge1dt, popprojAge1dt |> select(names(popAge1dt))),
#                       mx1dt = mx1dt,
#                       tfrproj1dt = rbind(tfr1dt, tfrproj1dt |> select(names(tfr1dt))),
#                       percentASFR1dt = percentASFR1dt,
#                       sexRatio1dt = sexRatio1dt,
#                       UNlocations = UNlocations,
#                       miscproj1dt = rbind(misc1dt, miscproj1dt |> select(names(misc1dt))))
# save(wpp2024.all, file = "wpp2024.all.Rda")

load("wpp2024.all.Rda")

###########################################################################################
# Short list of locations
locations <- suppressWarnings(wpp2024.all$UNlocations %>%
  mutate(iso3 = countrycode(name, "country.name", "iso3c")) %>%
  rename(location_name = name) %>%
  filter(!is.na(iso3) &
           location_name != "Less developed regions, excluding China" &
           country_code %in% unique(wpp2024.all$popprojAge1dt$country_code)) %>%
  select(location_name, country_code, iso3))

###########################################################################################
# Combine the projection and historical datasets
# Begin projection from the year 2025

iso3s            <- fread("country_region_mapping.csv") |> pull(ISO3)
no_gbd_data      <- c("COK","NIU","TKL","PRK","TWN","VEN","CPV","PSE")
iso3s<- iso3s[! iso3s %in% no_gbd_data]

pop.dt <- wpp2024.all$popprojAge1dt %>%
  right_join(locations, by = join_by(country_code)) %>%
  filter(year %in% 1950:2050 & !is.na(age) & !country_code == 2093) %>%
  rename(Female = popF, Male = popM) %>%
  select(-c(pop, name, country_code)) %>%
  gather(sex, Nx, -location_name, -iso3, -year, -age) %>%
  mutate(Nx = Nx*1e3, year = as.numeric(year), location_name = NULL) %>%
  arrange(iso3, year, sex, age)

# Mortality rates are full series to 2100
mort.dt <- wpp2024.all$mx1dt %>%
  right_join(locations, by = join_by(country_code)) %>%
  filter(year %in% 1950:2050 & !is.na(age) & !country_code == 2093) %>%
  rename(Female = mxF, Male = mxM) %>%
  select(-c(mxB, name, country_code)) %>%
  gather(sex, mx, -location_name, -iso3, -year, -age) %>%
  mutate(year = as.numeric(year), location_name = NULL) %>%
  arrange(iso3, year, sex, age)

retro.df <- left_join(mort.dt, pop.dt,
                      by = join_by(year, age, iso3, sex)) |>
  mutate(Nx = ifelse(Nx == 0, 1e-6, Nx)) |>
  filter(year %in% 1950:2050) |>
  rename(Population = Nx, iso3c = iso3) |>
  select(age, year, sex, mx, Population, iso3c) |>
  gather(metric, value, -c(age, year, sex, iso3c)) |>
  mutate(source = "WPP2024", cause = "All", calibration = source) |>
  select("age", "year", "sex" , "value",
         "metric", "source", "cause", "iso3c",
         "calibration")

save(retro.df, file = "retro.df.Rda")
