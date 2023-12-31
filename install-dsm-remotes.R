# install dep packages
remotes::install_github("cvpia-osc/DSMhabitat")
remotes::install_github("cvpia-osc/DSMflow")
remotes::install_github("cvpia-osc/DSMtemperature")
remotes::install_github("cvpia-osc/DSMCalibrationData")
remotes::install_github("cvpia-osc/DSMscenario")

# install dsm packages
remotes::install_github("cvpia-osc/fallRunDSM@movement-hypothesis")
remotes::install_github("cvpia-osc/winterRunDSM")
remotes::install_github("cvpia-osc/springRunDSM@calibration")

# install genetic alog
install.packages("GA")
