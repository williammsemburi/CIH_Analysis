# =============================================================================
# 2_CIH_Paper_Plots.R
#
# Generates all plots and data objects for the CIH paper.
# Outputs are saved to CIH/report/plots.Rda for use by Results_Report.Rmd.
#
# Required objects saved to plots.Rda:
#   df_labeled, region_colors, intro.map.hc, cih.p1, cih.p2,
#   deaths.costs.agg, map1, df.is, locs, plot.q70.trend, plot_ts_metric
#
# Dependencies:
#   - Stage 2 of 1_Run_CIH.R must have been run first
#   - CIH/output/ must contain fulldf/, fullq70/, fullcost/ outputs
#   - CIH/input/ must contain country_region_mapping.csv, gni.csv,
#     Intervention_map.csv, intervention_map.csv, cih_modules.csv
#   - FairChoices-v3/data/ must contain retro.df.Rda
# =============================================================================

rm(list = ls())

# Set working directory to FairChoices-v3 (sibling of CIH/)
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
setwd("FairChoices-v3")

library(fst)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)
library(data.table)
library(countrycode)
library(readxl)
library(highcharter)
library(tidyr)
library(stringr)
library(purrr)

# Load shared parameters: iso3s, no_vmc_countries, idlist_base, region_colors
source("../CIH/config.R")

# =============================================================================
# Load model output data
# =============================================================================

region_mapping <- fread("../CIH/input/country_region_mapping.csv") |>
  mutate(UN.region = countrycode(ISO3, "iso3c", "region"))

# Helper: read an .fst file and attach the ISO3 code parsed from the filename
read_with_iso3 <- function(path) {
  iso <- str_match(basename(path), "_([A-Z]{3})\\.fst$")[, 2]
  df  <- read_fst(path, as.data.table = TRUE)
  df[, iso3c := iso]
  return(df)
}

file_list1 <- list.files("../CIH/output/fulldf/",
                         pattern = "^full.df_.*\\.fst$",
                         full.names = TRUE)
full.dfs <- rbindlist(lapply(file_list1, read_with_iso3), use.names = TRUE, fill = TRUE)

load("../CIH/input/retro.df.Rda")

# =============================================================================
# Life table function
# Computes standard life table columns from a vector of age-specific mx.
# Pass parm = "ex" or parm = "lx" to extract a single column.
# Pass qual.adj = TRUE with dws to compute HLE (health-adjusted life expectancy).
# =============================================================================
get.lt <- function(mx, parm = NULL, dws = NULL, qual.adj = NULL) {

  mx    <- mx[!is.na(mx)]
  n_age <- length(mx)
  nx    <- rep(1, n_age)
  qx    <- 1 - exp(-nx * mx); qx[n_age] <- 1
  ax    <- nx + 1 / mx - nx / qx
  lx    <- c(1, cumprod(1 - qx)[1:(n_age - 1)])
  dx    <- c(rev(diff(rev(lx))), lx[1] - sum(rev(diff(rev(lx)))))
  nLx   <- nx * lx - (nx - ax) * dx
  Tx    <- rev(cumsum(rev(nLx)))
  ex    <- Tx / lx

  if (is.null(qual.adj)) {
    Sx        <- nLx / c(1, nLx)[1:n_age]
    Sx[n_age] <- nLx[n_age] / (nLx[n_age - 1] + nLx[n_age])
    df.lt     <- data.table(age = 0:(n_age - 1), ax, mx, lx, qx, dx, nLx, Tx, ex, Sx)
    if (is.null(parm)) df.lt else df.lt |> pull(parm)
  } else {
    qual <- (1 - dws)[!is.na(1 - dws)]
    nLxh <- nLx * qual
    Txh  <- rev(cumsum(rev(nLxh)))
    exh  <- Txh / lx
    exh[1]
  }
}

# =============================================================================
# Country-level mortality indicators (e0, 45q15, 70q0, 5q0, 1q0)
# series = "Retro" uses WPP historical data; "Model" uses CIH model output
# =============================================================================
get.q70.is <- function(series = "Retro") {

  if (series == "Retro") {

    region.q70 <- rbind(
      retro.df |> filter(year %in% 2000:2050 & iso3c %in% iso3s) |> mutate(source = "Baseline"),
      retro.df |> filter(year %in% 2000:2050 & iso3c %in% iso3s) |> mutate(source = "Adjusted")
    ) |>
      left_join(fread("../CIH/input/country_region_mapping.csv") |> rename(iso3c = ISO3)) |>
      mutate(calibration = "WPP2024") |>
      spread(metric, value) |>
      mutate(
        mx         = ifelse(mx > 1, 1 - 1e-9, mx),
        Population = ifelse(is.na(Population) | Population == 0, 1e-6, Population),
        Deaths     = mx * Population
      ) |>
      gather(metric, value, -c(age, year, sex, source, cause, iso3c, calibration, Region)) |>
      mutate(calibration = "WPP2024")

  } else {

    region.q70 <- full.dfs |>
      filter(cause == "All") |>
      left_join(fread("../CIH/input/country_region_mapping.csv") |> rename(iso3c = ISO3)) |>
      mutate(
        metric      = ifelse(calibration == "Unscaled", "mx.in", metric),
        calibration = ifelse(calibration == "Unscaled",
                             "Mortality Fertility Migration, Base Year", calibration)
      ) |>
      spread(metric, value) |>
      mutate(Population.in = ifelse(year == 2027, Population, NA)) |>
      group_by(iso3c, calibration, source, age, sex) |>
      mutate(Population.in = mean(Population.in, na.rm = TRUE)) |>
      ungroup() |>
      mutate(
        mx            = ifelse(mx > 1, 1 - 1e-9, mx),
        Population    = ifelse(is.na(Population) | Population == 0, 1e-6, Population),
        mx.in         = ifelse(mx.in > 1, 1 - 1e-9, mx.in),
        Population.in = ifelse(is.na(Population.in) | Population.in == 0, 1e-6, Population.in),
        Deaths        = mx * Population,
        Deaths.in     = mx.in * Population.in
      ) |>
      gather(metric, value, -c(age, year, sex, source, cause, iso3c, calibration, Region)) |>
      mutate(
        calibration = ifelse(metric %in% c("Population.in", "Deaths.in", "mx.in"),
                             "Unscaled", calibration),
        metric = ifelse(calibration == "Unscaled" & metric == "Population.in", "Population", metric),
        metric = ifelse(calibration == "Unscaled" & metric == "Deaths.in",     "Deaths",     metric),
        metric = ifelse(calibration == "Unscaled" & metric == "mx.in",         "mx",         metric)
      )
  }

  DT <- as.data.table(region.q70)

  DT_rates <- DT[
    metric %in% c("Deaths", "Population") & cause == "All" & year <= 2050,
    .(value = sum(value, na.rm = TRUE)),
    by = .(age, year, source, metric, calibration, iso3c)
  ]

  DT_rates <- dcast(
    DT_rates,
    age + year + source + iso3c + calibration ~ metric,
    value.var = "value"
  )[, mx := Deaths / Population][order(iso3c, source, calibration, year, age)]

  df_70q0 <- DT_rates[, {
    exs <- get.lt(mx, "ex")
    lxs <- get.lt(mx, "lx")
    list(
      e0  = round(exs[1],                          1),
      q45 = round(100  * (1 - lxs[61] / lxs[16]), 1),
      q70 = round(100  * (1 - lxs[71] / lxs[1]),  1),
      q5  = round(1000 * (1 - lxs[6]  / lxs[1]),  1),
      q1  = round(1000 * (1 - lxs[2]  / lxs[1]),  1)
    )
  }, by = .(iso3c, source, calibration, year)]

  df_70q0 <- melt(df_70q0,
    id.vars      = c("iso3c", "source", "calibration", "year"),
    measure.vars = c("e0", "q45", "q70", "q5", "q1"),
    variable.name = "metric", value.name = "value"
  )

  dcast(df_70q0, iso3c + calibration + year + metric ~ source, value.var = "value") |>
    arrange(iso3c, calibration, metric, year) |>
    as.data.frame()
}

# =============================================================================
# Regional mortality indicators (e0, 45q15, 70q0, 5q0, 1q0)
# series = "Retro" uses WPP; "Model" uses CIH model output
# =============================================================================
get.q70 <- function(series = "Retro") {

  if (series == "Retro") { yl <- 2050; cse <- "Alls" } else { yl <- 2025; cse <- "All" }

  region.q70 <- rbind(
    retro.df |> filter(year %in% 2000:2050 & iso3c %in% iso3s) |> mutate(source = "Baseline"),
    retro.df |> filter(year %in% 2000:2050 & iso3c %in% iso3s) |> mutate(source = "Adjusted"),
    full.dfs  |> filter(cause == cse)
  ) |>
    left_join(fread("../CIH/input/country_region_mapping.csv") |> rename(iso3c = ISO3)) |>
    spread(metric, value) |>
    mutate(
      mx         = ifelse(mx > 1, 1 - 1e-9, mx),
      Population = ifelse(is.na(Population) | Population == 0, 1e-6, Population),
      Deaths     = mx * Population
    ) |>
    gather(metric, value, -c(age, year, sex, source, cause, Region, iso3c, calibration))

  DT <- as.data.table(region.q70)

  DT_rates <- DT[
    metric %in% c("Deaths", "Population") & cause == "All" & year <= 2050,
    .(value = sum(value, na.rm = TRUE)),
    by = .(age, year, source, metric, calibration, Region)
  ]

  DT_rates2 <- DT_rates[
    metric %in% c("Deaths", "Population"),
    .(value = sum(value, na.rm = TRUE)),
    by = .(age, year, source, metric, calibration)
  ]

  DT_rates <- rbind(DT_rates, DT_rates2 |> mutate(Region = "World"))

  DT_rates <- dcast(
    DT_rates,
    age + year + source + Region + calibration ~ metric,
    value.var = "value"
  )[, mx := Deaths / Population][order(calibration, Region, source, year, age)]

  # get.lt is redefined locally here to avoid data.table scoping issues
  # with the outer environment. The logic is identical to the top-level get.lt.
  get.lt <- function(mx, parm = NULL, dws = NULL, qual.adj = NULL) {
    mx    <- mx[!is.na(mx)]
    n_age <- length(mx)
    nx    <- rep(1, n_age)
    qx    <- 1 - exp(-nx * mx); qx[n_age] <- 1
    ax    <- nx + 1 / mx - nx / qx
    lx    <- c(1, cumprod(1 - qx)[1:(n_age - 1)])
    dx    <- c(rev(diff(rev(lx))), lx[1] - sum(rev(diff(rev(lx)))))
    nLx   <- nx * lx - (nx - ax) * dx
    Tx    <- rev(cumsum(rev(nLx)))
    ex    <- Tx / lx
    if (is.null(qual.adj)) {
      Sx        <- nLx / c(1, nLx)[1:n_age]
      Sx[n_age] <- nLx[n_age] / (nLx[n_age - 1] + nLx[n_age])
      df.lt     <- data.table(age = 0:(n_age - 1), ax, mx, lx, qx, dx, nLx, Tx, ex, Sx)
      if (is.null(parm)) df.lt else df.lt |> pull(parm)
    } else {
      qual <- (1 - dws)[!is.na(1 - dws)]
      nLxh <- nLx * qual; Txh <- rev(cumsum(rev(nLxh))); exh <- Txh / lx; exh[1]
    }
  }

  df_70q0 <- DT_rates[, {
    exs <- get.lt(mx, "ex")
    lxs <- get.lt(mx, "lx")
    list(
      e0  = round(exs[1],                          1),
      q45 = round(100  * (1 - lxs[61] / lxs[16]), 1),
      q70 = round(100  * (1 - lxs[71] / lxs[1]),  1),
      q5  = round(1000 * (1 - lxs[6]  / lxs[1]),  1),
      q1  = round(1000 * (1 - lxs[2]  / lxs[1]),  1)
    )
  }, by = .(calibration, Region, source, year)]

  df_70q0 <- melt(df_70q0,
    id.vars       = c("calibration", "Region", "source", "year"),
    measure.vars  = c("e0", "q45", "q70", "q5", "q1"),
    variable.name = "metric", value.name = "value"
  )

  dcast(df_70q0, calibration + Region + year + metric ~ source, value.var = "value") |>
    arrange(calibration, Region, metric, year) |>
    as.data.frame()
}

# =============================================================================
# Build indicator datasets
# =============================================================================

df.wpp    <- get.q70("Retro")
df.cih    <- get.q70("Model")
df.wpp.is <- get.q70.is("Retro")
df.cih.is <- get.q70.is("Model")

# Regional time series (wide → long, recode metric labels)
df.wpp.cih <- rbind(df.wpp, df.cih) |>
  gather(grouping, value, -c(calibration, Region, year, metric)) |>
  arrange(grouping, metric, calibration, Region, year) |>
  select(grouping, metric, calibration, Region, year, value)

# Country-level time series
df.is <- rbind(df.cih.is, df.wpp.is) |>
  gather(grouping, value, -c(iso3c, year, calibration, metric)) |>
  arrange(grouping, calibration, metric, iso3c, year) |>
  select(grouping, calibration, metric, iso3c, year, value) |>
  mutate(
    location    = countrycode(iso3c, "iso3c", "country.name"),
    metric      = case_when(
      metric == "e0"  ~ "e0",
      metric == "q45" ~ "45q15",
      metric == "q70" ~ "70q0",
      metric == "q5"  ~ "5q0",
      metric == "q1"  ~ "1q0"
    ),
    calibration = ifelse(calibration == "Mortality Fertility Migration, Base Year",
                         "Scaled to WPP", calibration)
  )

# Labeled regional dataset for the interactive time-series chart
df_labeled <- df.wpp.cih |>
  group_by(grouping, metric, Region, calibration) |>
  arrange(year, .by_group = TRUE) |>
  mutate(
    calibration = ifelse(calibration == "Mortality Fertility Migration, Base Year",
                         "Scaled to WPP", calibration),
    dataLabel = case_when(
      year %in% c(2000, 2019, 2025, 2050) ~ as.character(value),
      TRUE ~ NA_character_
    ),
    metric = case_when(
      metric == "e0"  ~ "e0",
      metric == "q45" ~ "45q15",
      metric == "q70" ~ "70q0",
      metric == "q5"  ~ "5q0",
      metric == "q1"  ~ "1q0",
      TRUE            ~ metric
    )
  ) |>
  ungroup() |>
  arrange(grouping, metric, calibration, Region, year)

# =============================================================================
# Interactive time-series chart function (used by Shiny in Results_Report.Rmd)
# =============================================================================
plot_ts_metric <- function(met, grp, dfin = df_labeled,
                           reg_cols = region_colors, calib) {

  dfs  <- dfin |> filter(grouping == grp, metric == met)
  dfs1 <- dfs  |> filter(calibration == "WPP2024" & year < 2027)
  dfs2 <- dfs  |> filter(calibration == calib)
  dfs  <- rbind(dfs1, dfs2)

  if (nrow(dfs) == 0)
    stop(sprintf("No data for grouping = '%s', metric = '%s'.", grp, met))

  sublab <- switch(grp,
    "Adjusted" = "Mortality from WPP2024 (2000-2025) then FC impact (2026-2050)",
    "Baseline" = "Mortality from WPP2024 (2000-2025) then held constant (2026-2050)",
    "Measure calculated using WPP2024 mortality estimates for all years"
  )

  mets <- switch(met,
    "1q0"   = "1q0 (per 1,000 live births)",
    "5q0"   = "5q0 (per 1,000 live births)",
    "e0"    = "e0 (Life expectancy at birth, years)",
    "45q15" = "45q15 (Adult mortality, %)",
    "70q0"  = "70q0 (Probability of premature death, %)",
    met
  )

  value_suffix <- if (met %in% c("45q15", "70q0")) "%" else ""

  hc <- highchart() |>
    hc_xAxis(
      title = list(text = "Year"),
      min = min(dfs$year, na.rm = TRUE), max = max(dfs$year, na.rm = TRUE),
      tickInterval = 5, gridLineWidth = 1, gridLineColor = "#f0f0f0"
    ) |>
    hc_yAxis(
      title = list(text = mets),
      labels = list(format = "{value}"),
      gridLineWidth = 1, gridLineColor = "#e6e6e6"
    ) |>
    hc_title(text = paste0("Trends in ", met, " from 2000 to 2050 (", grp, ")"),
             style = list(fontSize = "16px", fontWeight = "bold")) |>
    hc_subtitle(text = sublab) |>
    hc_tooltip(
      shared = TRUE, crosshairs = TRUE,
      valueDecimals = 1, valueSuffix = value_suffix,
      backgroundColor = "rgba(255, 255, 255, 0.95)", borderWidth = 1
    ) |>
    hc_plotOptions(
      line = list(marker = list(enabled = FALSE), lineWidth = 2.5,
                  states = list(hover = list(lineWidthPlus = 1)))
    )

  for (region in unique(dfs$Region)) {
    region_data <- dfs |> filter(Region == region) |> arrange(year)
    data_points <- lapply(seq_len(nrow(region_data)), function(i) {
      point <- list(x = region_data$year[i], y = region_data$value[i])
      if (!is.na(region_data$dataLabel[i])) {
        point$dataLabels <- list(
          enabled = TRUE,
          formatter = JS(sprintf("function(){ return '%s'; }", region_data$dataLabel[i])),
          style = list(fontSize = "11px", fontWeight = "bold",
                       color = "black", textOutline = "2px white")
        )
      }
      point
    })

    hc <- hc |>
      hc_add_series(
        name   = region, type = "line",
        color  = reg_cols[[region]],
        data   = data_points,
        marker = list(enabled = FALSE)
      )
  }

  hc |>
    hc_legend(enabled = TRUE, align = "center", verticalAlign = "bottom",
              layout = "horizontal", itemStyle = list(fontSize = "11px")) |>
    hc_add_theme(hc_theme_smpl()) |>
    hc_exporting(
      enabled = TRUE, sourceWidth = 700, sourceHeight = 600,
      buttons = list(contextButton = list(
        menuItems = c("downloadPNG", "downloadJPEG", "downloadPDF")
      ))
    )
}

# =============================================================================
# Country-level 70q0 trend chart (used by Shiny in Results_Report.Rmd)
# =============================================================================
plot.q70.trend <- function(df.is, loc) {

  plot.df   <- df.is |> filter(location == loc, metric == "70q0") |>
    arrange(year) |> spread(grouping, value)

  p.df.wpp  <- plot.df |> filter(calibration == "WPP2024")
  p.df.uns  <- plot.df |> filter(calibration == "Unscaled")
  p.df.sca  <- plot.df |> filter(calibration == "Scaled to WPP")

  # Exponential path to halve 70q0 from its 2019 level by 2050
  v1        <- p.df.wpp$Baseline[20]
  r         <- (1 / 31) * log(0.5)
  v.thresh  <- round(v1 * exp(seq_len(31) * r), 2)

  p.df.tar <- data.frame(
    calibration = "Target", metric = "70q0", iso3c = p.df.wpp$iso3c[1],
    year = 2020:2050, location = loc, Adjusted = v.thresh, Baseline = v.thresh
  )

  highchart() |>
    hc_chart(type = "spline") |>
    hc_add_series(list_parse2(data.frame(x = p.df.sca$year, y = p.df.sca$Baseline)),
                  name = "Baseline (Scaled)",   type = "spline",
                  color = "#ff7f0e", marker = list(enabled = FALSE)) |>
    hc_add_series(list_parse2(data.frame(x = p.df.sca$year, y = p.df.sca$Adjusted)),
                  name = "Adjusted (Scaled)",   type = "spline",
                  color = "#1f77b4", marker = list(enabled = FALSE)) |>
    hc_add_series(list_parse2(data.frame(x = p.df.uns$year, y = p.df.uns$Baseline)),
                  name = "Baseline (Unscaled)", type = "spline",
                  color = "#ff7f0e", dashStyle = "Dash", marker = list(enabled = FALSE)) |>
    hc_add_series(list_parse2(data.frame(x = p.df.uns$year, y = p.df.uns$Adjusted)),
                  name = "Adjusted (Unscaled)", type = "spline",
                  color = "#1f77b4", dashStyle = "Dash", marker = list(enabled = FALSE)) |>
    hc_add_series(list_parse2(data.frame(x = p.df.tar$year, y = p.df.tar$Adjusted)),
                  name = "Target", type = "spline",
                  color = "#d62728", dashStyle = "Dash", marker = list(enabled = FALSE)) |>
    hc_add_series(list_parse2(data.frame(x = p.df.wpp$year, y = p.df.wpp$Adjusted)),
                  name = "WPP2024", type = "spline",
                  color = "grey", marker = list(enabled = FALSE)) |>
    hc_xAxis(title = list(text = "Year")) |>
    hc_yAxis(softMin = 0, title = list(text = "70q0 (PPD%)")) |>
    hc_tooltip(shared = TRUE, valueSuffix = "%") |>
    hc_title(text = paste0("Projected 70q0 for ", loc)) |>
    hc_legend(enabled = TRUE) |>
    hc_add_theme(hc_theme_smpl()) |>
    hc_exporting(enabled = TRUE)
}

locs <- sort(unique(df.is$location))

# Quick check
plot.q70.trend(df.is, "Tanzania")

# =============================================================================
# Choropleth map — baseline 70q0 by country (2019)
# =============================================================================

map.data.list <- df_labeled |>
  filter(calibration == "WPP2024" & grouping == "Baseline" &
         year == 2019 & metric == "70q0") |>
  mutate(code3 = countrycode(Region, "region", "iso3c")) |>
  rename(value = value)

stops <- data.frame(
  q = seq(0, 1, length.out = 7),
  c = c("#d4f4dd", "#a8ddb5", "#fef0bc", "#fdcc8a", "#fc8d59", "#d73027", "#7f0000"),
  stringsAsFactors = FALSE
) |> list_parse2()

map1 <- highchart(type = "map") |>
  hc_add_series_map(
    worldgeojson, map.data.list,
    joinBy = c("iso3", "code3"),
    value = "value", name = "Baseline 70q0 (%)",
    borderColor = "#000000", borderWidth = 0.2
  ) |>
  hc_title(text = "Baseline Estimates of 70q0 by Country (2019)",
           style = list(fontSize = "16px", fontWeight = "bold")) |>
  hc_colorAxis(
    stops = stops, min = 0, max = 100, layout = "horizontal",
    labels = list(format = "{value}")
  ) |>
  hc_legend(align = "center", verticalAlign = "bottom",
            layout = "horizontal", width = 600) |>
  hc_mapNavigation(enabled = TRUE) |>
  hc_add_theme(hc_theme_smpl()) |>
  hc_exporting(enabled = TRUE)

# =============================================================================
# Regional overview map (color by region, labels with POP / PPD / GNI)
# =============================================================================

ppd_2023 <- df_labeled |>
  filter(calibration == "WPP2024" & grouping == "Baseline" &
         year == 2023 & metric == "70q0") |>
  select(Region, value) |>
  rename(PPD = value)

pop_2023 <- retro.df |>
  filter(year == 2023 & metric == "Population") |>
  group_by(iso3c) |>
  summarise(Population = sum(value), .groups = "drop")

message("Loading GNI and population data...")

gni_data <- fread("../CIH/input/gni.csv") |>
  rename(name = "Country or Area", gni = Value) |>
  filter(name != "Kosovo") |>
  mutate(iso3c = countrycode(name, "country.name", "iso3c")) |>
  select(iso3c, gni) |>
  group_by(iso3c) |>
  summarise(gni = sum(gni), .groups = "drop") |>
  right_join(pop_2023, by = join_by(iso3c)) |>
  filter(iso3c %in% iso3s) |>
  left_join(region_mapping |> select(ISO3, Region) |> rename(iso3c = ISO3),
            by = join_by(iso3c)) |>
  rename(population = Population)

regional_summary <- rbind(gni_data, gni_data |> mutate(Region = "World")) |>
  group_by(Region) |>
  summarise(
    GNI                  = sum(gni, na.rm = TRUE) / sum(population, na.rm = TRUE),
    population_millions  = sum(population, na.rm = TRUE) / 1e6,
    n_countries          = n(),
    .groups = "drop"
  ) |>
  left_join(ppd_2023, by = "Region") |>
  mutate(GNI = round(GNI, 0))

message("\nRegional Summary:")
print(regional_summary, n = Inf)

map_data <- region_mapping |>
  select(ISO3, Region) |>
  rename(iso3c = ISO3) |>
  left_join(regional_summary |> filter(Region != "World"), by = "Region") |>
  mutate(
    pop_display = round(population_millions, 0),
    gni_display = format(GNI, big.mark = " "),
    ppd_display = round(PPD, 1)
  ) |>
  select(iso3c, Region, population_millions, PPD, GNI,
         pop_display, gni_display, ppd_display)

intro.map <- function() {

  # Distinct color palette for the regional overview map
  intro_region_colors <- list(
    "Central and Eastern Europe"         = "#3B3B8F",
    "China"                              = "#8B9B4B",
    "Central Asia"                       = "#D98B9B",
    "India"                              = "#2B6B3B",
    "Latin America and Caribbean"        = "#4B9B9B",
    "Middle East and North Africa"       = "#9B5B9B",
    "North Atlantic"                     = "#7BB5D5",
    "Sub-Saharan Africa"                 = "#6B9BC5",
    "United States"                      = "#C5B56B",
    "Western Pacific and Southeast Asia" = "#7B3B7B"
  )

  map_data_colored <- map_data |>
    mutate(color = unlist(intro_region_colors[Region]))

  hc <- hcmap(
    map = "custom/world-highres3",
    data = map_data_colored, value = "PPD",
    joinBy = c("iso-a3", "iso3c"), name = "PPD",
    dataLabels = list(enabled = FALSE),
    borderColor = "#FFFFFF", borderWidth = 0.5, nullColor = "#E0E0E0"
  ) |>
    hc_chart(marginTop = 10, marginBottom = 80) |>
    hc_colorAxis(
      stops   = color_stops(n = length(intro_region_colors),
                            colors = unlist(intro_region_colors)),
      visible = FALSE
    ) |>
    hc_title(text = NULL) |>
    hc_tooltip(
      useHTML  = TRUE,
      formatter = JS("function() {
        if (this.point.Region) {
          return '<b>' + this.point['iso3c'] + '</b><br>' +
                 'Region: ' + this.point.Region + '<br>' +
                 'PPD: '    + this.point.ppd_display + '%<br>' +
                 'GNI: $'   + this.point.gni_display + '<br>' +
                 'Population: ' + this.point.pop_display + ' M';
        } else {
          return '<b>' + this.point['iso3c'] + '</b><br>No data';
        }
      }")
    ) |>
    hc_mapNavigation(enabled = TRUE) |>
    hc_legend(
      enabled = TRUE, layout = "horizontal",
      align = "center", verticalAlign = "bottom",
      title = list(text = ""), itemStyle = list(fontSize = "9px"),
      symbolHeight = 10, symbolWidth = 10,
      itemMarginTop = 1, itemMarginBottom = 1
    ) |>
    hc_caption(
      useHTML = TRUE,
      text = paste0(
        "<div style='text-align:left; font-size:12px; margin-top:10px;'>",
        "<b>World</b><br>",
        "POP: ", round(regional_summary$population_millions[regional_summary$Region == "World"], 0), " M  |  ",
        "PPD: ", round(regional_summary$PPD[regional_summary$Region == "World"], 1), " %  |  ",
        "GNI: $", format(round(regional_summary$GNI[regional_summary$Region == "World"], 0), big.mark = " "),
        "<br><span style='font-size:10px;'>Population (millions, 2023) | ",
        "Probability of premature death before age 70 (2023) | ",
        "Gross national income per person, current US$ (2023)</span></div>"
      )
    ) |>
    hc_credits(enabled = FALSE)

  # Regional label positions (placed in ocean / empty space)
  label_positions <- data.frame(
    Region = c("United States", "Latin America and Caribbean", "North Atlantic",
               "Sub-Saharan Africa", "Middle East and North Africa",
               "Central and Eastern Europe", "Central Asia", "China",
               "India", "Western Pacific and Southeast Asia"),
    lat    = c( 36, -30,  58, -10,  26,  60, 42,  35,  18, -10),
    lon    = c(-125, -85,   5,   8,  40,  70, 60, 118,  70, 140)
  )

  label_data <- regional_summary |>
    filter(Region != "World") |>
    left_join(label_positions, by = "Region") |>
    mutate(label_text = paste0(
      "<b>", Region, "</b><br>",
      "POP: ", round(population_millions, 0), " M<br>",
      "PPD: ", round(PPD, 1), " %<br>",
      "GNI: $", format(round(GNI, 0), big.mark = " ")
    ))

  hc <- hc |>
    hc_plotOptions(series = list(dataLabels = list(enabled = FALSE)))

  for (i in seq_len(nrow(label_data))) {
    hc <- hc |>
      hc_add_series(
        type = "mappoint",
        data = list(list(lat = label_data$lat[i], lon = label_data$lon[i],
                         name = label_data$Region[i])),
        dataLabels = list(
          enabled = TRUE, useHTML = TRUE, format = label_data$label_text[i],
          style = list(fontSize = "8px", fontWeight = "normal", textOutline = "none"),
          backgroundColor = "rgba(255,255,255,0.85)",
          borderWidth = 0.3, borderColor = "#cccccc", padding = 2, borderRadius = 1
        ),
        marker = list(enabled = FALSE),
        enableMouseTracking = FALSE, showInLegend = FALSE
      )
  }

  for (i in seq_along(intro_region_colors)) {
    hc <- hc |>
      hc_add_series(
        type = "mappoint", name = names(intro_region_colors)[i], data = list(),
        color = intro_region_colors[[i]],
        marker = list(symbol = "circle", radius = 5),
        showInLegend = TRUE, enableMouseTracking = FALSE
      )
  }

  hc |>
    hc_exporting(
      enabled = TRUE, sourceWidth = 700, sourceHeight = 600,
      buttons = list(contextButton = list(
        menuItems = c("downloadPNG", "downloadJPEG", "downloadPDF")
      ))
    )
}

intro.map.hc <- intro.map()

# =============================================================================
# Figure X — Intervention cost-effectiveness by condition (ICER range plot)
# =============================================================================

cause_mapping <- fread("../CIH/input/intervention_map.csv")
cih_mapping   <- fread("../CIH/input/cih_modules.csv")

q70_path  <- "../CIH/output/fullq70/ids"
cost_path <- "../CIH/output/fullcost/ids"

# Compute ICER for each intervention ID × country combination
get_icers <- function(q70_path, cost_path, iso3) {

  idlist <- if (iso3 %in% no_vmc_countries) setdiff(idlist_base, 43L) else idlist_base

  icer.list <- list()
  for (theid in idlist) {

    q70_id <- read_fst(paste0(q70_path, "/full.q70_",    iso3, "_", theid, ".fst")) |>
      mutate(intervention_id = theid) |>
      rename(q70.base = Baseline, q70.adj = Adjusted)

    cost_id <- read_fst(paste0(cost_path, "/full.costs_", iso3, "_", theid, ".fst")) |>
      mutate(iso3c = iso3) |>
      group_by(iso3c, year, source, intervention_id) |>
      summarise(cost = sum(cost, na.rm = TRUE), .groups = "drop") |>
      spread(source, cost) |>
      rename(cost.base = Baseline, cost.adj = Adjusted)

    df.id <- left_join(q70_id, cost_id, by = join_by(year, iso3c, intervention_id)) |>
      filter(year == 2050)

    delta.cost <- 1e-6 * (sum(df.id$cost.adj) - sum(df.id$cost.base))
    delta.q70  <- df.id |> pull(q70.base) - df.id |> pull(q70.adj)

    icer.list[[paste(theid)]] <- data.frame(
      iso3c = iso3, id = theid,
      delta.q70 = delta.q70, delta.cost = delta.cost,
      icer = delta.q70 / delta.cost
    )
  }

  rbindlist(icer.list) |>
    mutate(icer = ifelse(is.na(icer), 0, icer)) |>
    arrange(-icer)
}

# Collect ICERs for all countries
all.icers.l <- list()
for (iso in iso3s) {
  message(iso)
  all.icers.l[[iso]] <- get_icers(q70_path, cost_path, iso)
}

all.icers.df <- rbindlist(all.icers.l) |>
  left_join(pop_2023, by = join_by(iso3c)) |>
  filter(delta.cost != 0) |>
  rename(pop = Population) |>
  mutate(wicer = pop * icer)

# Aggregate by CIH module
loop_input_1 <- cause_mapping |>
  filter(package_status == "Core") |>
  select(id, cih_module, category)

cih_modules <- unique(loop_input_1$cih_module)

loop_1_l <- list()
for (mods in cih_modules) {
  ids   <- loop_input_1 |> filter(cih_module == mods) |> pull(id) |> unique()
  d.sum <- all.icers.df |>
    filter(id %in% ids) |>
    group_by(id) |>
    summarise(wicer = sum(wicer, na.rm = TRUE), pop = sum(pop, na.rm = TRUE), .groups = "drop")

  loop_1_l[[mods]] <- cbind(
    d.sum |> mutate(wicer = wicer / pop) |>
      summarise(min.icer = min(wicer, na.rm = TRUE), max.icer = max(wicer, na.rm = TRUE)),
    d.sum |> summarise(wicer = sum(wicer, na.rm = TRUE), pop = sum(pop, na.rm = TRUE)) |>
      mutate(med.icer = wicer / pop, wicer = NULL, pop = NULL, cih_module = mods)
  )
}

loop_1_df <- rbindlist(loop_1_l) |>
  left_join(loop_input_1 |> select(cih_module, category) |> distinct())

# Aggregate by CIH condition
cih_conditions <- unique(cih_mapping$condition_name)

loop_2_l <- list()
for (conds in cih_conditions) {
  message(conds)
  cih_mods_for_cond <- cih_mapping |> filter(condition_name == conds) |> pull(cih_module)

  d.sum.l <- list()
  for (mods in cih_mods_for_cond) {
    ids <- loop_input_1 |> filter(cih_module == mods) |> pull(id) |> unique()
    d.sum.l[[paste(conds, mods)]] <- all.icers.df |>
      filter(id %in% ids) |>
      group_by(id) |>
      summarise(wicer = sum(wicer, na.rm = TRUE), pop = sum(pop, na.rm = TRUE), .groups = "drop") |>
      mutate(condition_name = conds)
  }
  d.sum <- rbindlist(d.sum.l)

  loop_2_l[[conds]] <- cbind(
    d.sum |> mutate(wicer = wicer / pop) |>
      summarise(min.icer = min(wicer, na.rm = TRUE), max.icer = max(wicer, na.rm = TRUE)),
    d.sum |> summarise(wicer = sum(wicer, na.rm = TRUE), pop = sum(pop, na.rm = TRUE)) |>
      mutate(med.icer = wicer / pop, wicer = NULL, pop = NULL, condition_name = conds)
  )
}

groups <- left_join(loop_input_1, cih_mapping,
                    relationship = "many-to-many", by = join_by(cih_module)) |>
  select(condition_name, category) |> unique()

loop_2_df <- rbindlist(loop_2_l) |> left_join(groups, by = join_by(condition_name))

# Build the ICER range chart
get.dots <- function(loop_2_df) {

  stopifnot(all(c("min.icer", "max.icer", "med.icer", "condition_name", "category") %in%
                  names(loop_2_df)))

  df <- loop_2_df |>
    mutate(condition_name = as.character(condition_name),
           category       = as.character(category)) |>
    arrange(category, condition_name) |>
    mutate(y = row_number() - 1L)

  y_map <- df |> select(y, condition_name)
  cats  <- unique(df$category)
  stopifnot(length(cats) == 2)

  df1    <- df |> filter(category == cats[1])
  df2    <- df |> filter(category == cats[2])
  sep_y  <- max(df1$y) + 0.5
  cap_h  <- 0.2
  x_max  <- max(df$max.icer, df$med.icer, na.rm = TRUE) * 1.05

  make_segments <- function(d) {
    unlist(lapply(seq_len(nrow(d)), function(i)
      list(list(x = d$min.icer[i], y = d$y[i]),
           list(x = d$max.icer[i], y = d$y[i]),
           list(x = NA, y = NA))), recursive = FALSE)
  }

  make_caps <- function(d, which = c("min", "max"), cap_h = 0.2) {
    which <- match.arg(which)
    unlist(lapply(seq_len(nrow(d)), function(i) {
      x <- if (which == "min") d$min.icer[i] else d$max.icer[i]
      list(list(x = x, y = d$y[i] - cap_h),
           list(x = x, y = d$y[i] + cap_h),
           list(x = NA, y = NA))
    }), recursive = FALSE)
  }

  make_medians <- function(d) {
    lapply(seq_len(nrow(d)), function(i)
      list(x = d$med.icer[i], y = d$y[i], name = d$condition_name[i],
           min_icer = d$min.icer[i], max_icer = d$max.icer[i], med_icer = d$med.icer[i]))
  }

  label_formatter <- JS(sprintf(
    "function() { var m = %s; return m[this.value] === undefined ? '' : m[this.value]; }",
    jsonlite::toJSON(setNames(y_map$condition_name, y_map$y), auto_unbox = TRUE)
  ))

  col1 <- "#4C78A8"  # blue  — infections / maternal
  col2 <- "#F58518"  # orange — NCDs / injuries

  add_range_series <- function(hc, segs, caps_l, caps_r, meds, col, cat_name) {
    hc |>
      hc_add_series(name = cat_name, type = "line", data = segs, color = col,
                    lineWidth = 2.5, marker = list(enabled = FALSE),
                    enableMouseTracking = FALSE, showInLegend = FALSE) |>
      hc_add_series(name = NULL, type = "line", data = caps_l, color = col,
                    lineWidth = 2.5, marker = list(enabled = FALSE),
                    enableMouseTracking = FALSE, showInLegend = FALSE) |>
      hc_add_series(name = NULL, type = "line", data = caps_r, color = col,
                    lineWidth = 2.5, marker = list(enabled = FALSE),
                    enableMouseTracking = FALSE, showInLegend = FALSE) |>
      hc_add_series(name = cat_name, type = "scatter", data = meds, color = col,
                    marker = list(symbol = "circle", radius = 6,
                                  lineWidth = 1, lineColor = "#FFFFFF"))
  }

  highchart() |>
    hc_title(text = "Intervention cost-effectiveness by condition",
             style = list(fontSize = "18px", fontWeight = "bold")) |>
    hc_subtitle(text = "ICER = Δ70q0 / ΔCost (in $M), for a standardised 10% coverage increase",
                style = list(fontSize = "12px", color = "#666666")) |>
    hc_xAxis(title = list(text = "ICER (Δ70q0 per $1M)"), min = 0, max = x_max,
             gridLineWidth = 1, gridLineColor = "#e6e6e6") |>
    hc_yAxis(
      title = list(text = "Condition"),
      min = min(df$y) - 0.5, max = max(df$y) + 0.5,
      tickInterval = 1, reversed = TRUE, gridLineWidth = 0,
      labels = list(formatter = label_formatter, style = list(fontSize = "11px")),
      plotLines = list(list(value = sep_y, width = 1.5, dashStyle = "Dash",
                            color = "#999999", zIndex = 5))
    ) |>
    hc_chart(type = "line", marginLeft = 240, marginRight = 40) |>
    hc_legend(enabled = TRUE, align = "center", verticalAlign = "bottom",
              itemStyle = list(fontSize = "11px")) |>
    hc_tooltip(
      shared = FALSE, useHTML = TRUE,
      formatter = JS("function() {
        if (this.point.name) {
          return '<b>' + this.point.name + '</b><br/>' +
                 'Min ICER: '    + this.point.min_icer.toFixed(4) + '<br/>' +
                 'Median ICER: ' + this.point.med_icer.toFixed(4) + '<br/>' +
                 'Max ICER: '    + this.point.max_icer.toFixed(4);
        }
        return false;
      }")
    ) |>
    add_range_series(make_segments(df1), make_caps(df1, "min", cap_h),
                     make_caps(df1, "max", cap_h), make_medians(df1), col1, cats[1]) |>
    add_range_series(make_segments(df2), make_caps(df2, "min", cap_h),
                     make_caps(df2, "max", cap_h), make_medians(df2), col2, cats[2])
}

cih.p1 <- get.dots(loop_2_df)

# =============================================================================
# Figure Y — Regional 70q0 trend relative to 2019 (50×50 target chart)
# =============================================================================

df.dq70 <- rbind(
  df.wpp |> filter(year >= 2019, year <= 2026),
  df.cih |> filter(calibration == "Mortality Fertility Migration, Base Year", year > 2026)
) |>
  arrange(Region, year) |>
  filter(metric == "q70") |>
  mutate(denom = ifelse(year == 2019, Baseline, NA)) |>
  group_by(Region) |>
  mutate(denom = mean(denom, na.rm = TRUE)) |>
  ungroup() |>
  mutate(value = Adjusted / denom,
         Baseline = NULL, Adjusted = NULL, denom = NULL, metric = NULL) |>
  rename(region = Region)

cih.p2 <- highchart() |>
  hc_chart(reflow = TRUE) |>
  hc_title(text  = "Projected reduction in premature mortality (70q0)",
           style = list(fontSize = "16px", fontWeight = "bold")) |>
  hc_subtitle(text = "Values shown relative to 2019 baseline (Adjusted / 2019 level)") |>
  hc_xAxis(title = list(text = "Year"), min = 2019, max = 2050, tickInterval = 5,
           gridLineWidth = 1, gridLineColor = "#f0f0f0") |>
  hc_yAxis(
    title = list(text = "70q0 relative to 2019 level"),
    min = 0.4, max = 1.4, tickInterval = 0.1,
    gridLineWidth = 1, gridLineColor = "#e6e6e6",
    plotLines = list(list(
      value = 0.5, color = "#666666", dashStyle = "Dash", width = 2, zIndex = 5,
      label = list(text = "50% of 2019 level", align = "",
                   style = list(color = "#666666", fontSize = "10px"))
    ))
  ) |>
  hc_tooltip(shared = TRUE, valueDecimals = 2, valueSuffix = " × baseline") |>
  hc_legend(enabled = TRUE, align = "center", verticalAlign = "bottom",
            layout = "horizontal", itemStyle = list(fontSize = "9px"), maxColumns = 4) |>
  hc_plotOptions(
    series = list(label = list(enabled = FALSE), turboThreshold = 0),
    line   = list(marker = list(enabled = FALSE), lineWidth = 2.5)
  ) |>
  hc_add_series_list(
    df.dq70 |>
      group_by(region) |>
      group_map(~ list(
        name  = .y$region,
        type  = "line",
        color = region_colors[[.y$region]],
        data  = lapply(seq_len(nrow(.x)), function(i) list(.x$year[i], .x$value[i]))
      )) |>
      setNames(NULL)
  ) |>
  hc_credits(enabled = FALSE) |>
  hc_exporting(
    enabled = TRUE, sourceWidth = 700, sourceHeight = 600,
    buttons = list(contextButton = list(
      menuItems = c("downloadPNG", "downloadJPEG", "downloadPDF")
    ))
  ) |>
  hc_boost(enabled = FALSE) |>
  hc_add_theme(hc_theme_smpl())

# =============================================================================
# Summary table — deaths averted and costs by region
# =============================================================================

gni_tots <- rbind(gni_data, gni_data |> mutate(Region = "World")) |>
  group_by(Region) |>
  summarise(GNI = sum(gni, na.rm = TRUE),
            population_millions = sum(population, na.rm = TRUE) / 1e6,
            n_countries = n(), .groups = "drop")

delta.q70 <- rbind(
  df.wpp |> filter(year == 2019, metric == "q70"),
  df.cih |> filter(calibration == "Mortality Fertility Migration, Base Year",
                   year == 2050, metric == "q70")
) |>
  arrange(Region, year) |>
  select(Region, year, Adjusted) |>
  spread(year, Adjusted) |>
  rename(q70.2019 = "2019", q70.2050 = "2050") |>
  mutate(dq70 = round(100 * (q70.2050 - q70.2019) / q70.2019, 1))

get.deaths.agg <- function() {
  deaths.agg <- full.dfs |>
    filter(cause == "All" & calibration == "Mortality Fertility Migration, Base Year") |>
    left_join(fread("../CIH/input/country_region_mapping.csv") |> rename(iso3c = ISO3)) |>
    spread(metric, value) |>
    mutate(
      mx         = ifelse(mx > 1, 1 - 1e-9, mx),
      Population = ifelse(is.na(Population) | Population == 0, 1e-6, Population),
      Deaths     = mx * Population
    ) |>
    gather(metric, value, -c(age, year, sex, source, cause, Region, iso3c, calibration)) |>
    as.data.table()

  deaths.agg <- rbind(
    deaths.agg[metric %in% c("Deaths", "Population") & cause == "All" & year <= 2050,
               .(value = sum(value, na.rm = TRUE)),
               by = .(age, year, source, metric, Region)],
    deaths.agg[metric %in% c("Deaths", "Population"),
               .(value = sum(value, na.rm = TRUE)),
               by = .(age, year, source, metric)] |>
      mutate(Region = "World")
  ) |>
    filter(metric == "Deaths") |>
    group_by(source, Region) |>
    summarise(Deaths = sum(value), .groups = "drop") |>
    spread(source, Deaths) |>
    mutate(deaths.averted = Baseline - Adjusted, Baseline = NULL, Adjusted = NULL)

  deaths.agg
}

deaths.agg <- get.deaths.agg()

file_list_costs <- list.files("../CIH/output/fullcost/",
                               pattern = "^full.costs_.*\\.fst$", full.names = TRUE)
full.costs <- rbindlist(lapply(file_list_costs, read_with_iso3), use.names = TRUE, fill = TRUE)

get.costs.agg <- function() {
  costs.agg <- full.costs |>
    left_join(fread("../CIH/input/country_region_mapping.csv") |> rename(iso3c = ISO3)) |>
    group_by(Region, source) |>
    summarise(cost = sum(cost, na.rm = TRUE), .groups = "drop")

  rbind(costs.agg,
        costs.agg |> group_by(source) |>
          summarise(cost = sum(cost), .groups = "drop") |>
          mutate(Region = "World")) |>
    spread(source, cost) |>
    mutate(incremental.cost = Adjusted - Baseline, Adjusted = NULL, Baseline = NULL)
}

costs.agg <- get.costs.agg()

deaths.costs.agg <- gni_tots |>
  left_join(delta.q70)  |>
  left_join(deaths.agg) |>
  left_join(costs.agg)  |>
  mutate(
    ave.cost.per.gni = 100 * (incremental.cost / length(2027:2050)) / GNI,
    incremental.cost = incremental.cost * 1e-9,
    deaths.averted   = deaths.averted   * 1e-6
  ) |>
  select(Region, q70.2019, q70.2050, dq70, deaths.averted, incremental.cost, ave.cost.per.gni)

# =============================================================================
# Validate and save all required objects to report/plots.Rda
# =============================================================================

required <- c(
  "df_labeled", "region_colors", "intro.map.hc", "cih.p1", "cih.p2",
  "deaths.costs.agg", "map1", "df.is", "locs", "plot.q70.trend", "plot_ts_metric"
)

missing <- setdiff(required, ls(envir = globalenv()))
if (length(missing) > 0) {
  stop("Cannot save — missing objects: ", paste(missing, collapse = ", "))
}

save(list = required, file = "../CIH/report/plots.Rda")
message("Saved report/plots.Rda with: ", paste(required, collapse = ", "))
