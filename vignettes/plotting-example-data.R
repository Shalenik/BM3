## ----setup, include=FALSE-----------------------------------------------------
library(BM3)
library(ggplot2)


## ----load-data----------------------------------------------------------------
climate_models <- readRDS(system.file("data", "climate_models.rds", package = "BM3"))
hadcrut5_annual <- readRDS(system.file("data", "hadcrut5_annual.rds", package = "BM3"))


## ----plot-climate-models------------------------------------------------------
if (is.data.frame(climate_models)) {
  ggplot(climate_models, aes(x = .data[[1]], y = .data[[2]])) +
    geom_line() +
    labs(title = "Climate Models Data", x = names(climate_models)[1], y = names(climate_models)[2])
} else {
  plot(climate_models, main = "Climate Models Data")
}


## ----plot-hadcrut5------------------------------------------------------------
if (is.data.frame(hadcrut5_annual)) {
  ggplot(hadcrut5_annual, aes(x = .data[[1]], y = .data[[2]])) +
    geom_line() +
    labs(title = "HadCRUT5 Annual Data", x = names(hadcrut5_annual)[1], y = names(hadcrut5_annual)[2])
} else {
  plot(hadcrut5_annual, main = "HadCRUT5 Annual Data")
}

