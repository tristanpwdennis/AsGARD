#####
#script to process all the different rasters and types and
#produce, for each, a single raster, cropped to the analysis area
#bounded between 0 and 1, with any necessary transformation applied
#Tristan Dennis 2024-12-13
#####

# Import packages
pkg <- c('sf','gdistance','rnaturalearth','terra','parallel')
# install.packages(pkg) # Uncomment if packages are not installed # nolint
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)

# Define some utility vars - projection and bbox
crs.geo <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
crs.geo.terra <- "+init=epsg:4326"
analysis_region <- c(xmin=27, xmax=47, ymin=6, ymax=21)
analysis_region.terra <- ext(analysis_region)

# Start off by loading the ocean from rnaturalearth
ocean <- ne_download(scale=10, type='ocean', category = 'physical')
ocean <- st_transform(ocean, crs=crs.geo)
ocean_vect <- vect(ocean) #transform to spatial vector
ocean_crop <- crop(ocean_vect, ext(analysis_region)) # crop to size

#Metadata
df_samples <- read.csv("cease_combinedmetadata_qcpass.20240914.csv")
points <- df_samples %>%
    dplyr::filter(country %in% c('Djibouti', 'Sudan', 'Ethiopia')) %>%
    dplyr::select(longitude, latitude, location) %>%
    unique()

pC <- as.matrix(points[c("longitude","latitude")])
# Some utility functions

# Function to crop and mask a raster to remove ocean cells, invert if asked
crop_mask_raster <- function(raster, analysis_region, invert) {
  #raster <- rast(file) # Load
  crs(raster) <- crs.geo #Define CRS
  #basename <- gsub('.tif','',basename(file)) # Get basename for output
  r_cropped <- crop(raster, ext(analysis_region))  # Crop to size
  r_cropped <- clamp(r_cropped, lower = 0)  # Rescale the raster to set negative values to zero

  # Invert the raster if 'invert' is TRUE
  if (invert) {
    r_cropped <- max(r_cropped[], na.rm = TRUE) - r_cropped
  }
    r_masked <- mask(r_cropped, ocean_crop, inverse = TRUE)  # Apply ocean mask
    return(r_masked)
}

# Function to transform temperature into mosquito development rate based on the model in Villena et al
temp_to_mdr <- function(file, invert, outfilename){
    # Load, crop raster
    raster <- crop_mask_raster(file, analysis_region = analysis_region, invert=invert)

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

    #Normalise MDR to be between 0 and 1
    r_min <- min(values(rast_transformed), na.rm = TRUE)
    r_max <- max(values(rast_transformed), na.rm = TRUE)
    #r_normalized <- (rast_transformed - r_min) / (r_max - r_min)
    # WO
    basename <- gsub('.tif','',basename(file))
    writeRaster(r_normalized, paste0('/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/processed/',outfilename,'.processed.tif'), overwrite = TRUE) # Write cropped, masked raster
    return(rast_transformed)
}

# Function to normalise a raster to bound the values between zero and 1 (e.g. will account for extreme values)
normalise_raster <- function(file, invert, outfilename){
    # Load, crop raster
    raster <- crop_mask_raster(file, analysis_region = analysis_region, invert=invert)
    # Normalize the raster to 0-1 range
    r_min <- min(values(raster), na.rm = TRUE)
    r_max <- max(values(raster), na.rm = TRUE)
    r_normalized <- (raster - r_min) / (r_max - r_min)
    # WO
   # basename <- gsub('.tif','',basename(file))
    writeRaster(raster, paste0('/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/processed/',outfilename,'.processed.tif'), overwrite = TRUE) # Write cropped, masked raster
    return(raster)
    }

# Function to create transition layer for a raster, and writeout
create_transitionlayer <- function(raster, varname) {
    print(paste0('calculating transition matrices and cost distances for: ', varname))
    transition_layer <- gdistance::transition(raster::raster(raster), mean, directions = 8)
    #transition_layer_c <- gdistance::geoCorrection(transition_layer, type='c', multpl=FALSE, scl=TRUE)
    transition_layer_r <- gdistance::geoCorrection(transition_layer, type='r', multpl=FALSE, scl=TRUE)

    #cosDist <- costDistance(transition_layer_c, pC)
    #resdist <- commuteDistance(transition_layer_r, pC)
    #write.table(cosDist, paste0('/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/paths_output/',varname,'.cosdist.txt')) # Write cropped, masked transitionlayer as RDS)
    #write.table(resdist, paste0('/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/paths_output/',varname,'.resdist.txt')) # Write cropped, masked transitionlayer as RDS)

    #saveRDS(transition_layer_c, paste0('/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/transitionlayers/',varname,'.transitionlayer.c.Rds')) # Write cropped, masked transitionlayer as RDS
    saveRDS(transition_layer_r, paste0('/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/transitionlayers/',varname,'.transitionlayer.r.Rds')) # Write cropped, masked transitionlayer as RDS

}

#Function to merge multiple raster tiles and reproject into WGS84
merge_reproject_raster_tiles <- function(directory, outfilename, outdir){
    raster_files <- list.files(directory, pattern = "\\.tif$", full.names = TRUE)
    raster_collection <- sprc(lapply(raster_files, rast))
    final_raster <- merge(raster_collection)
    #final_raster <- hdfs_merge[[2]]
    final_raster <- project(final_raster, "+init=epsg:4326")
    writeRaster(final_raster, paste0(outdir, outfilename))
}

#start with temperature - for this we need to convert temperature to mosquito development rate according to Vallena 2022
#bio1_file <- '/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/chelsa/CHELSA_bio1_1981-2010_V.2.1.tif'
#bio1_transformed <- temp_to_mdr(bio1_file, invert=FALSE, outfilename='mdr')
#bio1_normalised <- normalise_raster(bio1_file, invert=FALSE, outfilename='bio1')

#create_transitionlayer(bio1_transformed, varname='mdr')

#relative humidity
#hurs_file <- '/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/chelsa/CHELSA_hurs_mean_1981-2010_V.2.1.tif'
#hurs_transformed <- normalise_raster(hurs_file, invert=FALSE, outfilename='hurs')
#create_transitionlayer(hurs_transformed, varname='hurs')

#friction
#friction_file <- '/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/malaria_ap/201501_Global_Travel_Speed_Friction_Surface_2015.tif'
#friction_transformed <- normalise_raster(friction_file, invert=TRUE, outfilename='friction')
#create_transitionlayer(friction_transformed, varname='friction')

#population - this comes in tiles of Africa and Asia, so we need to merge these and then process
#pop_file <- '/home/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/WorldPop/ppp_2016_1km_Aggregated.tif'
#merge_reproject_raster_tiles(directory = pop_dir, outfilename = 'pop_raster_merged.tif', outdir = pop_dir)
#pop_transformed <- normalise_raster(pop_file, invert=FALSE, outfilename='pop')
#create_transitionlayer(pop_transformed, varname='pop')

#enhanced vegetation index - gives better resolution when veg is sparse than NDVI - need to merge multiple tiles and then normalise
#evi_dir <- '~/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/MODIS/modis_evi/'
#outfile <- 'evi_raster_aster.merged.tif'
#evi_file <- paste0(evi_dir,outfile)
#merge_reproject_raster_tiles(directory = evi_dir, outfilename = outfile, outdir = evi_dir)
#evi_transformed <- normalise_raster(evi_file, invert=FALSE, outfilename='evi')
#create_transitionlayer(evi_transformed, varname='evi')

#elevation - need to merge multiple tiles and then normalise
#ele_dir <- '~/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/ASTER/'
#outfile <- 'ele_raster_aster.merged.tif'
#merge_reproject_raster_tiles(directory = ele_dir, outfilename = outfile, outdir = ele_dir)
#ele_file <- paste0(ele_dir,outfile )
#ele_transformed <- normalise_raster(ele_file, invert=TRUE,outfilename='ele')
#create_transitionlayer(ele_transformed, varname='ele')

#slope (e.g. how hilly the landscape is)
outfile <- 'slope_gdem.merged.tif'
slope_dir <- '~/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/GDEM2010/'
#merge_reproject_raster_tiles(directory = slope_dir, outfilename = outfile, outdir = slope_dir, outfilename='slope')
slope_file <- paste0(slope_dir,outfile)
slope_r <- rast(slope_file)
#slope_transformed <- normalise_raster(slope_file, invert=TRUE, outfilename = outfile) #this is high resolution - 3 arc seconds - need to convert to 30 arc seconds (approx 1km at the equator)
slope_transformed_downscaled <- terra::aggregate(slope_r, fact = 10, fun = mean) # Use 'mean' or another function depending on your data
cropslope <- crop_mask_raster(slope_transformed_downscaled, analysis_region, invert=TRUE)
writeRaster(cropslope, '/mnt/user_shares/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/processed/slope.processed.tif.processed.untransformed.asc', overwrite=TRUE)

#writeRaster(slope_transformed_downscaled, '/mnt/user_shares/dennist/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/processed/slope.processed.tif', overwrite=TRUE)
#create_transitionlayer(slope_transformed_downscaled, varname='slope')
#rast('~/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/processed/slope_gdem.merged.tif.processed.tif')

#cattle density (stephensi like cattle) - we need to upscale to match the resolution of the other rasters (eg 1km/30s, 0.008 degrees)
#cattlefile <- '~/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/GLW4-2020.D-DA.CTL.tif'
#cattle_transformed <- normalise_raster(cattlefile, invert=FALSE, outfilename='cattle')
#cattle_upscaled <- disagg(cattle_transformed, fact=10)
#aligned_raster <- resample(cattle_upscaled, bio1_transformed)
#create_transitionlayer(cattle_transformed, varname='cattle')

#now we need to upscale the raster so that lcp works
#irrigation density
#irrigationfile <- '~/lstm_scratch/network_scratch/cease_workspace/environmental_analysis_202412/rasters/preprocessed/gmia.tif'
#irrigation_transformed <- normalise_raster(irrigationfile, invert=FALSE, outfilename='irrigation')
#create_transitionlayer(irrigation_transformed, varname='irrigation')
