get.q70 <- function(period.df, iso3) {
  DT <- as.data.table(period.df)

  # Aggregate and compute mx in one go
  DT_rates <- DT[
    metric %in% c("Deaths", "Population") &
      cause == "All" &
      year <= 2050,
    .(value = sum(value)),
    by = .(age, year, source, metric)
  ]

  # Wide format and mx
  DT_rates <- dcast(
    DT_rates,
    age + year + source ~ metric,
    value.var = "value"
  )[
    , mx := Deaths / Population
  ][
    order(source, year, age)
  ]

  # Life table per (source, year)
  df_70q0 <- DT_rates[
    ,
    {
      lxs <- get.lt(mx, "lx")
      q70 <- if (length(lxs) >= 71) 100 * (1 - lxs[71] / lxs[1]) else NA_real_
      .(q70 = round(q70, 1))
    },
    by = .(source, year)
  ]

  # Cast to wide and add iso
  df_70q0 <- dcast(df_70q0, year ~ source, value.var = "q70")
  df_70q0[, iso3c := iso3]

  as.data.frame(df_70q0)
}


iso3_ <- iso3
loc_ <- country[which(country$iso3 == iso3), "country"]
generate_epi_matrices(root_folder_src)

if(!223 %in% coverages_chosen$id){
  eecc_baseline<<-eecc_target<<-0.05
} else{
  eecc_baseline<<- as.numeric(coverages_chosen[which(coverages_chosen$id==223),"baseline"])
  eecc_target  <<- as.numeric(coverages_chosen[which(coverages_chosen$id==223),"target"])
}
if(!eecc_flag) eecc_baseline<-eecc_target<-1

infos_out<-infos(coverages_chosen=coverages_chosen,iso3c = iso3_)
model_out <- model_cascade(coverages_chosen=coverages_chosen,
                           iso3=iso3,tax_groups,einfo = infos_out$einfo,
                           cinfo = infos_out$cinfo,eecc_baseline=eecc_baseline,
                           eecc_target=eecc_target)

period.df <- model_out$proj_out$period.df

all.costs.in <- get.costs(period.df, model_out$cinfo, withdis)

pins <- all.costs.in$costs.by.row    %>%
  group_by(id) %>%
  group_by(year, source, intervention_id,
           id, sex, cause, treated.fraction,
           unit.cost) %>%
  summarise(
    total = sum(tcount, na.rm = TRUE),
    pin = sum(pin, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(year >= year0 & year <= yearn) %>%
  group_by(id, year, source) %>%
  summarise(pin = sum(pin, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

cinfo <- left_join(model_out$cinfo[, setdiff(colnames(model_out$cinfo), "pin")],
                   pins %>% filter(year == year0, source == "Baseline"), by = "id")

costs.in <- all.costs.in$costs.by.row

cids <- cinfo[which(cinfo$intervention_id %in% coverages_chosen$id), "intervention_id"]

coverage.in.b <- data.frame(intervention_id = cids,
                            source = rep("Baseline", length(cids)))

coverage.in.a <- data.frame(intervention_id = cids,
                            source = rep("Adjusted", length(cids)))
c0 <- left_join(coverage.in.b,
                coverages_chosen,
                by = c("intervention_id" = "id"))
c0 <- c0[, c("intervention_id", "source", "abaseline", "abaseline")]

colnames(c0) <- c("intervention_id", "source", "c0", "ct")
ct <- left_join(coverage.in.a,
                coverages_chosen,
                by = c("intervention_id" = "id"))
ct <- ct[, c("intervention_id", "source", "abaseline", "target")]
colnames(ct) <- c("intervention_id", "source", "c0", "ct")
coverage.in <- unique(rbind(c0, ct))

discount=1+3/100

full.costs <- left_join(costs.in, coverage.in, by = join_by(source, intervention_id)) %>%
  mutate(
    ci = ifelse(
      year > year0 + scaleup - 1,
      c0,
      c0 + (year + 1 - year0) * (ct - c0) / scaleup / (discount)^(1:scaleup)
    ),
    #
    pcov = pin * ci,
    pcov.cost = pin.cost * ci
  ) %>%
  group_by(year,
           #
           source,
           intervention_id,
           id,
           sex,
           cause,
           treated.fraction,
           unit.cost,
           c0,
           ci,
           ct) %>%
  summarise(
    total = sum(tcount),
    pin = sum(pin),
    pcov = sum(pcov),
    cost = sum(pcov.cost),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(year >= year0 & year <= yearn) %>%
  suppressWarnings()

full.q70   <- get.q70(period.df, iso3)

