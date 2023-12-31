winter_run_fitness <- function(
    known_adults,
    seeds,
    params,
    surv_adult_enroute_int,
    surv_juv_rear_int,
    surv_juv_rear_contact_points,
    surv_juv_rear_prop_diversions,
    surv_juv_rear_total_diversions,
    surv_juv_bypass_int,
    surv_juv_delta_int,
    surv_juv_delta_contact_points,
    surv_juv_delta_total_diverted,
    surv_juv_outmigration_sj_int,
    ocean_entry_success_int
) {

  params_init <- params

  params_init$..surv_adult_enroute_int = surv_adult_enroute_int
  params_init$..surv_juv_rear_int = rep(surv_juv_rear_int, 31)
  params_init$..surv_juv_rear_contact_points = surv_juv_rear_contact_points
  params_init$..surv_juv_rear_prop_diversions = surv_juv_rear_prop_diversions
  params_init$..surv_juv_rear_total_diversions = surv_juv_rear_total_diversions
  params_init$..surv_juv_bypass_int = surv_juv_bypass_int
  params_init$..surv_juv_delta_int = surv_juv_delta_int
  params_init$..surv_juv_delta_contact_points = surv_juv_delta_contact_points
  params_init$..surv_juv_delta_total_diverted = surv_juv_delta_total_diverted
  params_init$..surv_juv_outmigration_sj_int = surv_juv_outmigration_sj_int
  params_init$..ocean_entry_success_int = rep(ocean_entry_success_int, 31)

  keep <- c(1, 3)

  tryCatch({
    preds <- winter_run_model(mode = "calibrate",
                              seeds = seeds,
                              stochastic = FALSE,
                              ..params = params_init)

    known_nats <- known_adults[keep, 6:19, drop = FALSE] * (1 - winterRunDSM::params$proportion_hatchery[keep])
    mean_escapent <-rowMeans(known_nats, na.rm = TRUE)

    sse <- sum(((preds[keep,, drop = FALSE] - known_nats)^2)/mean_escapent, na.rm = TRUE)

    return(sse)
  },
  error = function(e) return(1e12),
  warning = function(w) return(1e12)
  )
}

#
# x <- sample(seq(.1, 3, by = .01), size = 11, replace = TRUE)
#
# winter_run_fitness(
#   known_adults = DSMCalibrationData::grandtab_observed$winter,
#   seeds = DSMCalibrationData::grandtab_imputed$winter,
#   params = winterRunDSM::params,
#   x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10],
#   x[11]
# )
