# remotes::install_github("cvpia-osc/DSMflow")
# remotes::install_github("cvpia-osc/DSMtemperature")
# remotes::install_github("cvpia-osc/DSMhabitat")
# remotes::install_github("cvpia-osc/DSMscenario")
# remotes::install_github("cvpia-osc/winterRunDSM")

library(winterRunDSM)
library(GA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

source("winter-run/wr-fitness.R")
source("winter-run/wr-update-params.R")

params <- DSMCalibrationData::set_synth_years(winterRunDSM::params_2022)
params$prey_density <- rep("max", 31)


wr_lo_bounds <- rep(-10, 11)
wr_lo_bounds[c(0, 3:5)] <- 0

wr_hi_bounds <- rep(15, 11)

# Perform calibration --------------------
res <- ga(type = "real-valued",
          fitness =
            function(x) -winter_run_fitness(
              known_adults = DSMCalibrationData::grandtab_observed$winter,
              seeds = DSMCalibrationData::grandtab_imputed$winter,
              params = params,
              x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10],
              x[11]
            ),
          lower = wr_lo_bounds,
          upper = wr_hi_bounds,
          popSize = 150,
          maxiter = 10000,
          run = 50,
          parallel = TRUE,
          pmutation = .4)

readr::write_rds(res, paste0("winter-run/result-max-last-week", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".rds"))

# Evaluate Results ------------------------------------
res <- read_rds("winter-run/result-max-cor-556-2023-08-11_133255.rds")
r1_solution <- res@solution[1, ]
keep <- c(1, 3)

new_params <- update_params(x = r1_solution, winterRunDSM::params_2022)
calib_params <- DSMCalibrationData::set_synth_years(new_params)
calib_params$prey_density <- rep("max", 31)
calib_sim <- winter_run_model(seeds = DSMCalibrationData::grandtab_imputed$winter, mode = "calibrate",
                           ..params = calib_params,
                           stochastic = FALSE,
                           calib_return_all = T)


nat_spawners <- as_tibble(calib_sim[keep, ,drop = F]) %>%
  mutate(watershed = DSMscenario::watershed_labels[keep]) %>%
  gather(year, spawners, -watershed) %>%
  mutate(type = "simulated",
         year = readr::parse_number(year) + 2002)

observed <- as_tibble((1 - winterRunDSM::params$proportion_hatchery[c(1, 3)]) * DSMCalibrationData::grandtab_observed$winter[keep, , drop = FALSE]) %>%
  mutate(watershed = DSMscenario::watershed_labels[keep]) %>%
  gather(year, spawners, -watershed) %>%
  mutate(type = "observed", year = as.numeric(year)) %>%
  filter(!is.na(spawners),
         year >= 2005)



eval_df <- bind_rows(nat_spawners, observed) %>%
  filter(!(year %in% 2005:2006))


eval_df %>%
  ggplot(aes(year, spawners, color = type)) +
  geom_line() +
  facet_wrap(~watershed, scales = "free_y") +
  scale_x_continuous(breaks = seq(2002, 2019, by = 2))

eval_df %>%
  spread(type, spawners) %>%
  ggplot(aes(observed, simulated)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Observed vs Predicted updated",
       x = "Observed Natural Spawners",
       y = "Predicted Natural Spawners") +
  xlim(0, 20000) +
  ylim(0, 20000)

eval_df %>%
  spread(type, spawners) %>%
  filter(!is.na(observed)) %>%
  group_by(watershed) %>%
  summarise(
    r = cor(observed, simulated, use = "pairwise.complete.obs")
  ) %>% arrange(desc(abs(r)))

eval_df %>%
  spread(type, spawners) %>%
  filter(!is.na(observed)) %>%
  summarise(
    r = cor(observed, simulated, use = "pairwise.complete.obs")
  )


