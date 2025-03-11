library(terra)


# Read MTM data
mtm_locs <- read.csv('~/Projects/AsGARD/metadata/mtm_invasive_vector_atlas_locs.csv')

# Read CEASE locs
site_status <- fread('~/Projects/AsGARD/data/feems_20240920/status_table.csv')


# Read rasters
evi<-terra::rast('~/Projects/AsGARD/data/environmental_data_rasters/evi/modis_evi.rawmerged.tif')
friction <- terra::rast('~/Projects/AsGARD/data/environmental_data_rasters/transit_friction/201501_Global_Travel_Speed_Friction_Surface_2015.tif')
bio1 <- terra::rast('~/Projects/AsGARD/data/environmental_data_rasters/climate/CHELSA_bio1_1981-2010_V.2.1.tif')


# Function to transform temperature into mosquito development rate based on the model in Villena et al
temp_to_mdr <- function(raster){
  # Load, crop raster

  # Now apply the transformation
  c_est <- 7.138852e-05 # Constant that determines the curvature of the Briére distribution, estimated from the model of MDR in the paper above
  # Calculate tmin and tmax from the raster
  tmin <- min(values(raster, na.rm=TRUE))
  tmax <- max(values(raster, na.rm=TRUE))
  
  mdr_transformation <- function(temp) { # Equation from Vallena, 2022
    c_est * temp * ((temp - tmin) * sqrt((tmax - temp)) * (tmax > temp))
  }
  # Apply the function to the raster
  rast_transformed <- app(raster, mdr_transformation)
  return(rast_transformed)
}

mdr <-temp_to_mdr(bio1)

# Now let's transform to be between 0 and 1 again, and crop to the desired extent

# Define extent
extent <- ext(1,68,-19,38)

# Force same resolution
mdr_t <- resample(mdr, bio1_t, method = "bilinear")
fric_t <- resample(friction, bio1_t, method = "bilinear")
evi_t <- resample(evi, bio1_t, method = "bilinear")

# Now crop
bio1_t <- crop(bio1_t, extent)
fric_t <- crop(fric_t, extent)
evi_t <- crop(evi_t, extent)

#now transform

# Function to transform and crop raster
normalise_raster <- function(raster) {
  #Normalise MDR to be between 0 and 1
  r_min <- min(values(raster), na.rm = TRUE)
  r_max <- max(values(raster), na.rm = TRUE)
  r_normalized <- (raster - r_min) / (r_max - r_min)
  return(r_normalized)
}


# Normalise
fric_t <- 1-normalise_raster(fric_t)
mdr_t <- normalise_raster(mdr_t)
evi_t <- normalise_raster(evi_t)

# Cumulative raster
cumulative_permissivity <- fric_t+bio1_t+evi_t

# Now wrangle WHO and CEASE data, plot map with borders and water bodies etc






# Now we need to force same grid and resolution
#Choose a reference (bio1 has good resolution)










ext(fric_t)
ext(evi_t)
ext(mdr_t)
