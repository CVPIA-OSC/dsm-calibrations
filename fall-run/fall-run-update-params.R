fr_update_params <- function(x, params) {

  surv_adult_enroute = x[1]
  surv_adult_prespawn = x[2]
  surv_egg_to_fry = x[3]
  bypass_surv_juv = x[4]
  upsac_surv_juv = x[5]
  butte_surv_juv = x[6]
  clear_surv_juv = x[7]
  deer_surv_juv = x[8]
  mill_surv_juv = x[9]
  sac_surv_juv = x[10]
  feather_and_bear_surv_juv = x[11]
  yuba_surv_juv = x[12]
  american_surv_juv = x[13]
  deltatribs_surv_juv = x[14]
  moke_surv_juv = x[15]
  merced_surv_juv = x[16]
  stan_surv_juv = x[17]
  tuol_surv_juv = x[18]
  sj_surv_juv = x[19]
  surv_juv_rear_contact_points = x[20]
  surv_juv_rear_prop_diversions = x[21]
  surv_juv_rear_total_diversions = x[22]
  surv_juv_bypass_int = x[23]
  surv_juv_delta_int = x[24]
  surv_juv_delta_contact_points = x[25]
  surv_juv_delta_total_diverted = x[26]
  surv_juv_outmigration_sj_int = x[27]
  default_ocean_entry_surv = x[28]
  upsac_ocean_entry_surv = x[29]
  butte_ocean_entry_surv = x[30]
  deer_ocean_entry_surv = x[31]
  mill_ocean_entry_surv = x[32]
  bear_feather_ocean_entry = x[33]
  yuba_ocean_entry_surv = x[34]
  american_ocean_entry_surv = x[35]
  deltatribs_ocean_entry_surv = x[36]
  moke_ocean_entry_surv = x[37]
  merced_ocean_entry_surv = x[38]
  stan_ocean_entry_surv = x[39]
  tuol_ocean_entry_surv = x[40]

  params$..surv_adult_enroute_int = surv_adult_enroute
  params$..surv_adult_prespawn_int = surv_adult_prespawn
  params$..surv_egg_to_fry_int = surv_egg_to_fry
  params$..surv_juv_rear_int = c(`Upper Sacramento River` = upsac_surv_juv,
                                 `Antelope Creek` = deer_surv_juv,
                                 `Battle Creek` = deer_surv_juv,
                                 `Bear Creek` = deer_surv_juv,
                                 `Big Chico Creek` = deer_surv_juv,
                                 `Butte Creek` = butte_surv_juv,
                                 `Clear Creek` = clear_surv_juv,
                                 `Cottonwood Creek` = deer_surv_juv,
                                 `Cow Creek` = deer_surv_juv,
                                 `Deer Creek` = deer_surv_juv,
                                 `Elder Creek` = deer_surv_juv,
                                 `Mill Creek` = mill_surv_juv,
                                 `Paynes Creek` = deer_surv_juv,
                                 `Stony Creek` = deer_surv_juv,
                                 `Thomes Creek` = deer_surv_juv,
                                 `Upper-mid Sacramento River` = sac_surv_juv,
                                 `Sutter Bypass` = bypass_surv_juv,
                                 `Bear River` = feather_and_bear_surv_juv,
                                 `Feather River` = feather_and_bear_surv_juv,
                                 `Yuba River` = yuba_surv_juv,
                                 `Lower-mid Sacramento River` = sac_surv_juv,
                                 `Yolo Bypass` = bypass_surv_juv,
                                 `American River` = american_surv_juv,
                                 `Lower Sacramento River` = sac_surv_juv,
                                 `Calaveras River` = deltatribs_surv_juv,
                                 `Cosumnes River` = deltatribs_surv_juv,
                                 `Mokelumne River` = moke_surv_juv,
                                 `Merced River` = merced_surv_juv,
                                 `Stanislaus River` = stan_surv_juv,
                                 `Tuolumne River` = tuol_surv_juv,
                                 `San Joaquin River` = sj_surv_juv)
  params$..surv_juv_rear_contact_points = surv_juv_rear_contact_points
  params$..surv_juv_rear_prop_diversions = surv_juv_rear_prop_diversions
  params$..surv_juv_rear_total_diversions = surv_juv_rear_total_diversions
  params$..surv_juv_bypass_int = surv_juv_bypass_int
  params$..surv_juv_delta_int = surv_juv_delta_int
  params$..surv_juv_delta_contact_points = surv_juv_delta_contact_points
  params$..surv_juv_delta_total_diverted = surv_juv_delta_total_diverted
  params$..surv_juv_outmigration_sj_int = surv_juv_outmigration_sj_int
  params$..ocean_entry_success_int = c(
    `Upper Sacramento River` = upsac_ocean_entry_surv,
    `Antelope Creek` = default_ocean_entry_surv,
    `Battle Creek` = default_ocean_entry_surv,
    `Bear Creek` = default_ocean_entry_surv,
    `Big Chico Creek` = default_ocean_entry_surv,
    `Butte Creek` = butte_ocean_entry_surv,
    `Clear Creek` = default_ocean_entry_surv,
    `Cottonwood Creek` = default_ocean_entry_surv,
    `Cow Creek` = default_ocean_entry_surv,
    `Deer Creek` = deer_ocean_entry_surv,
    `Elder Creek` = default_ocean_entry_surv,
    `Mill Creek` = mill_ocean_entry_surv,
    `Paynes Creek` = default_ocean_entry_surv,
    `Stony Creek` = default_ocean_entry_surv,
    `Thomes Creek` = default_ocean_entry_surv,
    `Upper-mid Sacramento River` = default_ocean_entry_surv,
    `Sutter Bypass` = default_ocean_entry_surv,
    `Bear River` = bear_feather_ocean_entry,
    `Feather River` = bear_feather_ocean_entry,
    `Yuba River` = yuba_ocean_entry_surv,
    `Lower-mid Sacramento River` = default_ocean_entry_surv,
    `Yolo Bypass` = default_ocean_entry_surv,
    `American River` = american_ocean_entry_surv,
    `Lower Sacramento River` = default_ocean_entry_surv,
    `Calaveras River` = deltatribs_ocean_entry_surv,
    `Cosumnes River` = deltatribs_ocean_entry_surv,
    `Mokelumne River` = moke_ocean_entry_surv,
    `Merced River` = merced_ocean_entry_surv,
    `Stanislaus River` = stan_ocean_entry_surv,
    `Tuolumne River` = tuol_ocean_entry_surv,
    `San Joaquin River` = default_ocean_entry_surv
    )

  return(params)

}



fr_update_params_multi_route <- function(x, params) {

  surv_adult_enroute = x[1]
  surv_adult_prespawn = x[2]
  surv_egg_to_fry = x[3]
  bypass_surv_juv = x[4]
  upsac_surv_juv = x[5]
  butte_surv_juv = x[6]
  clear_surv_juv = x[7]
  deer_surv_juv = x[8]
  mill_surv_juv = x[9]
  sac_surv_juv = x[10]
  feather_and_bear_surv_juv = x[11]
  yuba_surv_juv = x[12]
  american_surv_juv = x[13]
  deltatribs_surv_juv = x[14]
  moke_surv_juv = x[15]
  merced_surv_juv = x[16]
  stan_surv_juv = x[17]
  tuol_surv_juv = x[18]
  sj_surv_juv = x[19]
  surv_juv_rear_contact_points = x[20]
  surv_juv_rear_prop_diversions = x[21]
  surv_juv_rear_total_diversions = x[22]
  surv_juv_bypass_int = x[23]
  surv_juv_delta_int = x[24]
  surv_juv_delta_contact_points = x[25]
  surv_juv_delta_total_diverted = x[26]
  surv_juv_outmigration_sj_int = x[27]
  default_ocean_entry_surv = x[28]
  upsac_ocean_entry_surv = x[29]
  butte_ocean_entry_surv = x[30]
  deer_ocean_entry_surv = x[31]
  mill_ocean_entry_surv = x[32]
  bear_feather_ocean_entry = x[33]
  yuba_ocean_entry_surv = x[34]
  american_ocean_entry_surv = x[35]
  deltatribs_ocean_entry_surv = x[36]
  moke_ocean_entry_surv = x[37]
  merced_ocean_entry_surv = x[38]
  stan_ocean_entry_surv = x[39]
  tuol_ocean_entry_surv = x[40]
  habitat_cap = x[41]
  floodplain_habitat_cat = x[42]

  params$..surv_adult_enroute_int = surv_adult_enroute
  params$..surv_adult_prespawn_int = surv_adult_prespawn
  params$..surv_egg_to_fry_int = surv_egg_to_fry
  params$..surv_juv_rear_int = c(`Upper Sacramento River` = upsac_surv_juv,
                                 `Antelope Creek` = deer_surv_juv,
                                 `Battle Creek` = deer_surv_juv,
                                 `Bear Creek` = deer_surv_juv,
                                 `Big Chico Creek` = deer_surv_juv,
                                 `Butte Creek` = butte_surv_juv,
                                 `Clear Creek` = clear_surv_juv,
                                 `Cottonwood Creek` = deer_surv_juv,
                                 `Cow Creek` = deer_surv_juv,
                                 `Deer Creek` = deer_surv_juv,
                                 `Elder Creek` = deer_surv_juv,
                                 `Mill Creek` = mill_surv_juv,
                                 `Paynes Creek` = deer_surv_juv,
                                 `Stony Creek` = deer_surv_juv,
                                 `Thomes Creek` = deer_surv_juv,
                                 `Upper-mid Sacramento River` = sac_surv_juv,
                                 `Sutter Bypass` = bypass_surv_juv,
                                 `Bear River` = feather_and_bear_surv_juv,
                                 `Feather River` = feather_and_bear_surv_juv,
                                 `Yuba River` = yuba_surv_juv,
                                 `Lower-mid Sacramento River` = sac_surv_juv,
                                 `Yolo Bypass` = bypass_surv_juv,
                                 `American River` = american_surv_juv,
                                 `Lower Sacramento River` = sac_surv_juv,
                                 `Calaveras River` = deltatribs_surv_juv,
                                 `Cosumnes River` = deltatribs_surv_juv,
                                 `Mokelumne River` = moke_surv_juv,
                                 `Merced River` = merced_surv_juv,
                                 `Stanislaus River` = stan_surv_juv,
                                 `Tuolumne River` = tuol_surv_juv,
                                 `San Joaquin River` = sj_surv_juv)
  params$..surv_juv_rear_contact_points = surv_juv_rear_contact_points
  params$..surv_juv_rear_prop_diversions = surv_juv_rear_prop_diversions
  params$..surv_juv_rear_total_diversions = surv_juv_rear_total_diversions
  params$..surv_juv_bypass_int = surv_juv_bypass_int
  params$..surv_juv_delta_int = surv_juv_delta_int
  params$..surv_juv_delta_contact_points = surv_juv_delta_contact_points
  params$..surv_juv_delta_total_diverted = surv_juv_delta_total_diverted
  params$..surv_juv_outmigration_sj_int = surv_juv_outmigration_sj_int
  params$..ocean_entry_success_int = c(
    `Upper Sacramento River` = upsac_ocean_entry_surv,
    `Antelope Creek` = default_ocean_entry_surv,
    `Battle Creek` = default_ocean_entry_surv,
    `Bear Creek` = default_ocean_entry_surv,
    `Big Chico Creek` = default_ocean_entry_surv,
    `Butte Creek` = butte_ocean_entry_surv,
    `Clear Creek` = default_ocean_entry_surv,
    `Cottonwood Creek` = default_ocean_entry_surv,
    `Cow Creek` = default_ocean_entry_surv,
    `Deer Creek` = deer_ocean_entry_surv,
    `Elder Creek` = default_ocean_entry_surv,
    `Mill Creek` = mill_ocean_entry_surv,
    `Paynes Creek` = default_ocean_entry_surv,
    `Stony Creek` = default_ocean_entry_surv,
    `Thomes Creek` = default_ocean_entry_surv,
    `Upper-mid Sacramento River` = default_ocean_entry_surv,
    `Sutter Bypass` = default_ocean_entry_surv,
    `Bear River` = bear_feather_ocean_entry,
    `Feather River` = bear_feather_ocean_entry,
    `Yuba River` = yuba_ocean_entry_surv,
    `Lower-mid Sacramento River` = default_ocean_entry_surv,
    `Yolo Bypass` = default_ocean_entry_surv,
    `American River` = american_ocean_entry_surv,
    `Lower Sacramento River` = default_ocean_entry_surv,
    `Calaveras River` = deltatribs_ocean_entry_surv,
    `Cosumnes River` = deltatribs_ocean_entry_surv,
    `Mokelumne River` = moke_ocean_entry_surv,
    `Merced River` = merced_ocean_entry_surv,
    `Stanislaus River` = stan_ocean_entry_surv,
    `Tuolumne River` = tuol_ocean_entry_surv,
    `San Joaquin River` = default_ocean_entry_surv
  )

  params$..habitat_capacity = habitat_cap
  params$..floodplain_capacity = floodplain_habitat_cat

  return(params)

}





