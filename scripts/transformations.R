

vars <- c('ele')
for (var in vars){
  m <-rast(paste0('/Users/dennistpw/Projects/AsGARD/data/environmental_modelling/circuitscape_data/',var,'/',var,'.processed.untransformed.tif.asc'))
  r_mean <- global(m, fun = "mean", na.rm = TRUE)[1, 1]
  r_sd <- global(m, fun = "sd", na.rm = TRUE)[1, 1]
  r_z <- (m - r_mean) / r_sd
  r_z_shifted <- r_z +10
  fn <- paste0('/Users/dennistpw/Projects/AsGARD/data/environmental_modelling/circuitscape_data/',var,'/',var,'.processed.z.untransformed.asc')
  writeRaster(r_z_shifted, filename=fn, overwrite=TRUE)
}
var<-'mdr'
v<-paste0('/Users/dennistpw/Projects/AsGARD/data/environmental_modelling/circuitscape_data/',var,'/',var,'.processed.untransformed.tif.asc')
m <- rast(v)
r_mean <- global(m, fun = "mean", na.rm = TRUE)[1, 1]
r_sd <- global(m, fun = "sd", na.rm = TRUE)[1, 1]
r_z <- (m - r_mean) / r_sd
r_z_shifted <- r_z +10
plot(r_z_shifted)

library(terra)
f<- rast('/Users/dennistpw/Projects/AsGARD/data/environmental_modelling/circuitscape_data/mdr/mdr.processed.z.untransformed.asc')
hist(f)

f_half <- aggregate(f, fact = 2, fun = mean)  # Aggregate by a factor of 2 using the mean function
#writeRaster(f_half, '/Users/dennistpw/Projects/AsGARD/data/environmental_modelling/circuitscape_data/ele/ele.processed.untransformed.tif.asc', overwrite=TRUE)


d<-terra::rast('~/Projects/AsGARD/data/environmental_modelling/circuitscape_data/ele/ele.asc')
f_half <- terra::aggregate(d, fact = 2, fun = mean)  # Aggregate by a factor of 2 using the mean function

terra::writeRaster(f_half,'~/Projects/AsGARD/data/environmental_modelling/circuitscape_data/ele/ele.dscale.asc')
