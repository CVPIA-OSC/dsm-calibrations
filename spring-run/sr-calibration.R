library(springRunDSM)
library(GA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

source("spring-run/sr-fitness.R")
source("spring-run/sr-update-params.R")

params <- DSMCalibrationData::set_synth_years(springRunDSM::params_2022)
params$prey_density <- rep("max")
best_previous_solution <- readr::read_rds("calibration/calibration-results-2023.rds")@solution


lo_bounds <- rep(-10, 27)
lo_bounds[1] <- 2.5
lo_bounds[13:15] <- 0
hi_bounds <- rep(10, 27)
# Perform calibration --------------------
res <- ga(type = "real-valued",
          fitness =
            function(x) -spring_run_fitness(
              known_adults = DSMCalibrationData::grandtab_observed$spring,
              seeds = DSMCalibrationData::grandtab_imputed$spring,
              params = params,
              x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10],
              x[11], x[12], x[13], x[14], x[15], x[16], x[17], x[18], x[19],
              x[20], x[21], x[22], x[23], x[24], x[25], x[26], x[27]
            ),
          lower = lo_bounds,
          upper = hi_bounds,
          popSize = 150,
          maxiter = 10000,
          run = 50,
          parallel = TRUE,
          pmutation = .4) # <- remove this argument if wanting to start from zero
                                                # its a good idea to start from scractch when doing "annual" calibrations
                                                # if you make a small tweak to the fitness function and want to re-run calibration you
                                                # pass best previous solution to start from there and see if changes make a difference

readr::write_rds(res, paste0("spring-run/result-last-week-max-prey-1-", format(Sys.time(), "%Y-%m-%d_%H-%M"), ".rds"))

# res <- readr::read_rds("winter-run/result-max-cor-556-2023-08-11_133255.rds")

# Evaluate Results ------------------------------------

keep <- c(2, 3, 6, 7, 10, 12, 19, 20)
r1_solution <- res@solution[1, ]

r1_params <- update_params(x = r1_solution, springRunDSM::params_2022)
r1_params <- DSMCalibrationData::set_synth_years(r1_params)
r1_params$prey_density <- rep("max", 31)
r1_sim <- spring_run_model(seeds = DSMCalibrationData::grandtab_imputed$spring, mode = "calibrate",
                           ..params = r1_params,
                           stochastic = FALSE)


r1_nat_spawners <- as_tibble(r1_sim[keep, ,drop = F]) %>%
  mutate(watershed = DSMscenario::watershed_labels[keep]) %>%
  gather(year, spawners, -watershed) %>%
  mutate(type = "simulated",
         year = readr::parse_number(year) + 5)



r1_observed <- as_tibble((1 - springRunDSM::params$proportion_hatchery[keep]) * DSMCalibrationData::grandtab_observed$spring[keep,, drop=F]) %>%
  mutate(watershed = DSMscenario::watershed_labels[keep]) %>%
  gather(year, spawners, -watershed) %>%
  mutate(type = "observed", year = as.numeric(year) - 1997) %>%
  filter(!is.na(spawners),
         year > 5)



r1_eval_df <- bind_rows(r1_nat_spawners, r1_observed)


r1_eval_df %>%
  ggplot(aes(year, spawners, color = type)) + geom_line() + facet_wrap(~watershed, scales = "free_y")

r1_eval_df %>%
  spread(type, spawners) %>%
  ggplot(aes(observed, simulated)) + geom_point() +
  geom_abline(intercept = 0, slope = 1)

x <- r1_eval_df %>%
  spread(type, spawners) %>%
  filter(!is.na(observed)) %>%
  group_by(watershed) %>%
  summarise(
    r = cor(observed, simulated, use = "pairwise.complete.obs")
  )


r1_eval_df %>%
  spread(type, spawners) %>%
  filter(!is.na(observed)) %>%
  summarise(
    r = cor(observed, simulated, use = "pairwise.complete.obs")
  )

