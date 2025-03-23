# Import packages
pkg <- c('sf', 'tidyverse','rnaturalearth','terra')
#install.packages(pkg) # Uncomment if packages are not installed
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)

# Load fst data
fst <- read.csv('~/Projects/AsGARD/data/fst_dist.csv')
# Sample metadata
df_samples <- read.csv('~/Projects/AsGARD/metadata/cease_combinedmetadata_qcpass.20240914.csv')
points <- df_samples %>%
  dplyr::filter(country %in% c('Djibouti', 'Sudan', 'Ethiopia')) %>%
  dplyr::select(longitude, latitude, location, country) %>%
  unique()


### start by loading example processed raster, and defining the shape around which no mirgation occurs
# here we will use the western border between Ethiopia & Eritrea and Sudan

#load country admin borders
border_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_0_sovereignty/ne_10m_admin_0_sovereignty.shp")
ethshape <- border_shape[border_shape$ADM0_A3 == 'ETH']

world <- ne_countries(scale = "large", returnclass = "sf")
ethiopia_eritrea <- world[world$name %in% c("Ethiopia", "Eritrea", "Djibouti","Somaliland","Somalia"), ]

#load an example raster
example_r <-  rast('~/Projects/AsGARD/data/environmental_modelling/circuitscape_data/mdr/mdr.asc')
head(as.data.frame(example_r))

points_sf <- st_as_sf(points, coords = c("longitude", "latitude"), crs = 4326)


#Cropping the line
# Now we want to set the notional barrier to migration to be between sudan and not sudan
# This roughly corresponds to the mountains that we hypothesise are a major barrier to stephensi migration
# This is the northernmost point of eritrea, and the southernmost point of the line where it bisects the lower limit
# of the raster

#merge ethiopia and eritrea polygons
combined_polygon <- st_union(ethiopia_eritrea)

poly_coords <- matrix(c(
  35, 3,  # Point 1 (longitude, latitude)
  38.585112812853204, 18.028833819929034,  # Point 2
  50, 18.028833819929034,   # Point 3
  50, 3,   # Point 4, 
  35, 3    # Closing the polygon (back to Point 1)
), ncol = 2, byrow = TRUE)

polygon <- st_sfc(st_polygon(list(poly_coords)), crs = 4326)
crop_poly <- as(polygon, "SpatVector")
plot(crop_poly)


combined_polygon_terra <- as(combined_polygon, "SpatVector")
plot(combined_polygon_terra)
lines(crop_poly, add=TRUE)




# Plot raster with country borders overlaid
ggplot() +
  # Raster plot
  geom_raster(data = as.data.frame(example_r, xy = TRUE), aes(x = x, y = y, fill = lyr.1)) +
  scale_fill_viridis_c() + 
  # Borders plot
  geom_sf(data = combined_polygon, fill = NA, color = "red", size = 1) +
  theme_minimal() +
  geom_sf(data = points_sf, aes(color = "black"), size = 3) +
  labs(title = "Raster with Borders of Ethiopia and Eritrea")

#now let;s set raster cells touching the line to 
vectpoly<-vect(combined_polygon)
line_cells <- extract(example_r, vectpoly, buffer = 1)  # 'buffer' can help include cells that touch the line
line_cells <- unique(unlist(line_cells))

example_r[line_cells] <- NA


#rasterise exclusion zone and output
# Define a resolution for the output raster (e.g., 1000 x 1000 cells)

example_r
rasterized

res <- 0.008333333  # Resolution (adjust as needed)
# Create an empty raster with the same extent as the SpatVector
ext_vector <- ext(e)
r <- rast(ext_vector, resolution = res)
rasterized <- rasterize(e, r, field = 0) #create raster of eclusion zone with -1 values
rasterized
plot(rasterized)
writeRaster(rasterized,'~/Projects/AsGARD/data/environmental_modelling/circuitscape_data/exclude_region.asc', overwrite=TRUE)



########
#raster processing
########

#For each raster we need to load it, crop it to the right size, mask areas of ocean

#define some utility vars
crs.geo <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
analysis_region <- c(xmin=27, xmax=47, ymin=6, ymax=21)

vars <- c('friction','ele', 'evi', 'hurs', 'mdr', 'pop', 'slope')

#vars <- c('ele', 'evi', 'hurs', 'mdr', 'pop', 'slope')

# Load the first variable to get the shared location data
init.pdat <- readRDS(paste0('~/Projects/AsGARD/data/environmental_modelling/least_cost_path_data/lcp_means_data/', vars[1], '.pathdata.Rds'))

# Create the base data frame with shared location data
shared_data <- data.frame(
  location_x = sapply(init.pdat[1:253], function(x) x$p1_loc),
  location_y = sapply(init.pdat[1:253], function(x) x$p2_loc),
  latitude_x = sapply(init.pdat[1:253], function(x) x$p1_latlon[1]),
  longitude_x = sapply(init.pdat[1:253], function(x) x$p1_latlon[2]),
  latitude_y = sapply(init.pdat[1:253], function(x) x$p2_latlon[1]),
  longitude_y = sapply(init.pdat[1:253], function(x) x$p2_latlon[2])
)

init.pdat
# Loop over variables and add columns to the shared data frame
for (var in vars) {
  # Load the data for the current variable
  init.pdat <- readRDS(paste0('~/Projects/AsGARD/data/environmental_modelling/least_cost_path_data/lcp_means_data/', var, '.pathdata.Rds'))
  
  # Extract and calculate the desired values
  shared_data[[paste0("lcp.mean.", var)]] <- sapply(init.pdat[1:253], function(x) mean(x$lcp_values[2, ], na.rm=TRUE))
  shared_data[[paste0("lcp.total.", var)]] <- sapply(init.pdat[1:253], function(x) sum(x$lcp_values[2, ],na.rm=TRUE))
  shared_data[[paste0("line.mean.", var)]] <- sapply(init.pdat[1:253], function(x) mean(x$line_values[2, ], na.rm=TRUE))
  shared_data[[paste0("line.total.", var)]] <- sapply(init.pdat[1:253], function(x) sum(x$line_values[2, ],na.rm=TRUE))
}

#calculate distance
#calcuklate geographic distance between points
shared_data$pdist <- mapply(function(lat1, lon1, lat2, lon2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2))
}, shared_data$latitude_x, shared_data$longitude_x, shared_data$latitude_y, shared_data$longitude_y)
shared_data$pdist <- shared_data$pdist / 1000


# Load wind
mean.wind <- read.table('~/Projects/AsGARD/data/environmental_modelling/least_cost_path_data/wind_time_means.txt')
mean.wind[lower.tri(mean.wind)] <- NA #Mask lower triangle
colnames(mean.wind) <- points$location
mean.wind <- as.data.frame(mean.wind)
mean.wind$site_x <- points$location
wind_lcp_longform <- mean.wind %>% pivot_longer(cols=1:23, names_to = 'site_y', values_to = 'wind.travel.time')
# Get rid of NA

wind_lcp_longform <- wind_lcp_longform[!is.na(wind_lcp_longform$wind.travel.time),]
colnames(wind_lcp_longform) <- c('location_x','location_y','wind.travel.time')
shared_data <- left_join(shared_data, wind_lcp_longform)


#Define common loc var
fst$loc_fac <- paste0(fst$location_x, fst$location_y)
shared_data$loc_fac <- paste0(shared_data$location_x, shared_data$location_y)



result <- fst %>%
  full_join(shared_data, by = c("location_x" = "location_x", "location_y" = "location_y")) %>% 
  full_join(shared_data, by = c("location_x" = "location_y", "location_y" = "location_x"))



# Create a combined 'location_pair' column in both dataframes
fst_df <- fst %>%
  mutate(location_pair = pmin(location_x, location_y) %>%
           paste(pmax(location_x, location_y), sep = "-"))

winddf <- shared_data %>%
  mutate(location_pair = pmin(location_x, location_y) %>%
           paste(pmax(location_x, location_y), sep = "-"))

# Perform the join using the 'location_pair' column
result <- full_join(fst_df, winddf, by = "location_pair")


location_df <- unique(df_samples[,c('location','country','latitude','longitude')])


result_metadata <- left_join(result, location_df, by=c('popa' = 'location')) %>% 
  left_join(.,location_df, by=c('popb' = 'location'))

result_metadata$countryfact <- paste0(result_metadata$country.x, result_metadata$country.y)

lcp_plot_df=result_metadata %>% filter(country.x!='Yemen' & country.y != 'Yemen') #%>% filter(countryfact == 'EthiopiaDjibouti' | countryfact ==  'EthiopiaEthiopia' | countryfact== 'SudanSudan')

lcp_plot_df<-lcp_plot_df %>% filter(countryfact == 'SudanSudan' | countryfact == 'EthiopaEthiopia')
pdist <- ggplot(lcp_plot_df,aes(x=pdist, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'Geographic distance (km)')+
  theme(legend.position = "none")

ele<- ggplot(lcp_plot_df,aes(x=lcp.mean.ele, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'elevation')+
  theme(legend.position = "none")



wind <- ggplot(lcp_plot_df,aes(x=wind.travel.time, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'mean wind travel time')+
  theme(legend.position = "none")



hurs<-ggplot(lcp_plot_df,aes(x=lcp.mean.hurs, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'relative humidity')+
  theme(legend.position = "none")



evi <- ggplot(lcp_plot_df,aes(x=lcp.mean.evi, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'enhanced vegetation index')+
  theme(legend.position = "none")



pop <- ggplot(lcp_plot_df,aes(x=lcp.mean.pop, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'human population')+
  theme(legend.position = "none")



mdr <- ggplot(lcp_plot_df,aes(x=lcp.mean.mdr, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'mosquito development rate')+
  theme(legend.position = "none")

friction <- ggplot(lcp_plot_df,aes(x=lcp.mean.friction, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'friction')+
  theme(legend.position = "none")

slope_leg <- ggplot(lcp_plot_df,aes(x=lcp.mean.slope, y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  labs(title = 'slope')


slopeleg <- cowplot::get_legend(slope_leg)

cowplot::plot_grid(pdist, wind, ele, hurs, evi, pop, friction,mdr, slope_leg+theme(legend.position = "none"), slopeleg)

#let's have a wee look at the matrices we generated to see if they match


mat <- as.matrix(readRDS('~/Projects/AsGARD/data/environmental_modelling/least_cost_path_data/cost_distance_matrices/ele.lcp.Rds'))
colnames(mat) <- points$location
mat[lower.tri(mat)] <- NA #Mask lower triangle
mat <- as.data.frame(mat)
mat$location_x <- points$location
mat_long <- mat %>% pivot_longer(cols=1:23, names_to = 'location_y', values_to = 'lcp.ele')
mat_long <- mat_long[!is.na(mat_long$lcp.ele),]

location_df <- unique(df_samples[,c('location','country','latitude','longitude')])

result_metadata <- left_join(mat_long, location_df, by=c('location_x' = 'location')) %>% 
  left_join(.,location_df, by=c('location_y' = 'location'))


mat_df <- result_metadata %>%
  mutate(location_pair = pmin(location_x, location_y) %>%
           paste(pmax(location_x, location_y), sep = "-"))

# Perform the join using the 'location_pair' column
result <- full_join(fst_df, mat_df, by = "location_pair")

ggplot(result, aes(x=fst, y=lcp.ele))+
  geom_point()





vars <- c('mdr','ele','pop','cattle','hurs','evi','friction','slope')

varcols <- list()
for (i in seq_along(vars)){
  col <- paste0(vars[[i]],'.circuit')
  mdrmat <- fread(paste0('~/Projects/AsGARD/data/environmental_modelling/circuitscape_data/',vars[[i]],'/',vars[[i]],'_resistances.out'))
  mdrmat[lower.tri(mdrmat)] <- NA #Mask lower triangle
  mdrmat<- as.data.frame(mdrmat[-1,-1])
  colnames(mdrmat) <- points$location
  mdrmat$location_x <- points$location
  mat_long <- mdrmat %>% pivot_longer(cols=1:23, names_to = 'location_y', values_to = col)
  mat_long <- mat_long[!is.na(mat_long[,3]),]
  varcols[[i]] <- mat_long
  }

#add wind
wind <- fread('~/Projects/AsGARD/data/environmental_modelling/least_cost_path_data/wind_time_means.txt')[,-1]
wind[lower.tri(wind)] <- NA #Mask lower triangle
colnames(wind) <- points$location
wind$location_x <- points$location
wind_long <- wind %>% pivot_longer(cols=1:23, names_to = 'location_y', values_to = 'mean.wind.traveltime')
varcols[[9]] <- wind_long

merged_df <- Reduce(function(x, y) merge(x, y, by = c("location_x", "location_y")), varcols)

# add wind


#wind_long <- wind_long[!is.na(wind_long[,3]),]

location_df <- unique(df_samples[,c('location','country','latitude','longitude')])


result <- fst %>%
  full_join(merged_df, by = c("location_x" = "location_x", "location_y" = "location_y")) %>% 
  full_join(merged_df, by = c("location_x" = "location_y", "location_y" = "location_x"))



# Create a combined 'location_pair' column in both dataframes
fst_df <- fst %>%
  mutate(location_pair = pmin(location_x, location_y) %>%
           paste(pmax(location_x, location_y), sep = "-"))

merged_df <- merged_df %>%
  mutate(location_pair = pmin(location_x, location_y) %>%
           paste(pmax(location_x, location_y), sep = "-"))

# Perform the join using the 'location_pair' column
result <- full_join(fst_df, merged_df, by = "location_pair")


#location_df <- unique(df_samples[,c('location','country','latitude','longitude')])


result_metadata <- left_join(result, location_df, by=c('popa' = 'location')) %>% 
  left_join(.,location_df, by=c('popb' = 'location'))

result_metadata$countryfact <- paste0(result_metadata$country.x, result_metadata$country.y)

plot_df=result_metadata %>% filter(country.x!='Yemen' & country.y != 'Yemen') %>% 
  filter(countryfact == 'SudanSudan' | countryfact == 'EthiopiaEthiopia' | countryfact == 'EthiopiaDjibouti')

plot_df$invasion <- ifelse(plot_df$countryfact %in% c("EthiopiaEthiopia", "EthiopiaDjibouti", "DjiboutiEthiopia"), "Ethiopia", "Sudan")




#plot_df=result %>% filter(country.x!='Yemen' & country.y != 'Yemen') %>% filter(countryfact == 'EthiopiaDjibouti' | countryfact ==  'EthiopiaEthiopia' | countryfact== 'SudanSudan')
pdist.circuit <- ggplot(plot_df, aes(x=dist,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")


mdr.circuit <- ggplot(plot_df, aes(x=mdr.circuit,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")


hurs.circuit <- ggplot(plot_df, aes(x=hurs.circuit,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")


ele.circuit <- ggplot(plot_df, aes(x=ele.circuit,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")


evi.circuit <- ggplot(plot_df, aes(x=evi.circuit,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")


pop.circuit <- ggplot(plot_df, aes(x=pop.circuit,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")


slope.circuit <- ggplot(plot_df, aes(x=slope.circuit,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()

cattle.circuit <- ggplot(plot_df, aes(x=sqrt(cattle.circuit),y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")


friction.circuit <- ggplot(plot_df, aes(x=friction.circuit,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")

wind.lcp <- ggplot(plot_df, aes(x=mean.wind.traveltime,y=fst, colour=countryfact))+
  geom_point()+
  theme_classic()+
  theme(legend.position = "none")

legend <- cowplot::get_legend(slope.circuit)

cowplot::plot_grid(pdist.circuit, wind.lcp, ele.circuit, hurs.circuit, evi.circuit, pop.circuit, friction.circuit,mdr.circuit, cattle.circuit, slope.circuit+theme(legend.position = "none"), legend)






library(corrplot)
M <- cor(plot_df[,c('dist','mean.wind.traveltime','ele.circuit','hurs.circuit','evi.circuit','pop.circuit','friction.circuit','mdr.circuit','slope.circuit', 'cattle.circuit')])
corrplot(M, method = 'number') # colorful number

#function for plotting residuals
plotresids <- function(model){
  modsimout <- simulateResiduals(fittedModel = model, plot = F)
  residuals(modsimout, quantileFunction = qnorm, outlierValues = c(0,1))
  plot(modsimout)
}


plot_df$fst_lin <- plot_df$fst / (1-plot_df$fst) #Linearise Fst
model_df <- plot_df

# function to rescale all variables


variable_df <- model_df[c('mean.wind.traveltime','dist','mdr.circuit','ele.circuit','pop.circuit','cattle.circuit','hurs.circuit','evi.circuit','friction.circuit','slope.circuit')]
variable_df$cattle.circuit[variable_df$cattle.circuit == -1] <- NA


col_list <- list()
for (i in seq(1:ncol(variable_df))){
  col <- variable_df[,i]
  scaled_col <- col / max(col, na.rm = TRUE)
  col_list[[i]] <- scaled_col
}
rescaled_df <- data.frame(do.call(cbind,col_list))
colnames(rescaled_df) <- c('mean.wind.traveltime','dist','mdr.circuit','ele.circuit','pop.circuit','cattle.circuit','hurs.circuit','evi.circuit','friction.circuit','slope.circuit')
rescaled_model_df <- cbind(model_df[,c('invasion','fst_lin')], rescaled_df)

#model_df$wind.scaled <- model_df$mean.wind.traveltime / max(model_df$mean.wind.traveltime)
#model_df$dist.scaled <- model_df$dist / max(model_df$dist)

library(lme4)
library(brms)

#without country
model0 <- glm(fst_lin~dist, data=plot_df)
model1 <- glm(fst_lin ~ dist + cattle.circuit +pop.circuit + friction.circuit + evi.circuit + mean.wind.traveltime, data=plot_df) # model with everythiung
model2 <- glm(fst_lin ~ cattle.circuit +pop.circuit + friction.circuit + evi.circuit + mean.wind.traveltime, data=plot_df) #minus dist
model3 <- glm(fst_lin ~ dist + cattle.circuit +pop.circuit + friction.circuit + evi.circuit + mean.wind.traveltime, data=plot_df, family=gaussian(link="log")) #log gaussian dist
model4 <- glm(fst_lin ~ cattle.circuit +pop.circuit + friction.circuit + evi.circuit + mean.wind.traveltime, data=plot_df, family=gaussian(link="log")) #as above but mionus distance
model5 <- glm(fst_lin ~ cattle.circuit + friction.circuit + evi.circuit + mean.wind.traveltime, data=plot_df, family=gaussian(link="log")) #as above but mionus population
model6 <- glm(fst_lin ~ cattle.circuit + friction.circuit + evi.circuit, data=plot_df, family=gaussian(link="log")) #as above but minus wind

#
c.model0 <- glm(fst_lin~dist*countryfact, data=plot_df)
c.model1 <- glmer(fst_lin ~ dist + cattle.circuit +pop.circuit + friction.circuit + evi.circuit + mean.wind.traveltime + (1|invasion), data=plot_df, family=gaussian(link="log")) # model with everythiung
summary(c.model1)


library("rstanarm")
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


stan1 <- rstanarm::stan_glmer(
  fst_lin ~ dist + cattle.circuit +pop.circuit + friction.circuit + evi.circuit + ele.circuit + slope.circuit + mean.wind.traveltime + hurs.circuit + mdr.circuit + (1|invasion), 
  data = rescaled_model_df,
  family = gaussian(link = "log"),
  chains = 4, iter = 1000 )  

summary(stan1)
plot(stan1)
plot(stan1, plotfun = "trace")
prior_summary(stan1)
pp_check(stan1, plotfun = "stat", stat = "mean")


stan2 <- rstanarm::stan_glm(
  fst_lin ~ mean.wind.traveltime*invasion, 
  data = rescaled_model_df,
  family = gaussian(link = "log"),
  chains = 4, iter = 1000 )  



summary(stan2)
plot(stan2)
plot(stan2, plotfun = "trace")
prior_summary(stan2)
pp_check(stan2, plotfun = "stat", stat = "mean")

predict(stan2, new_data = data.frame(mean.wind.traveltime = seq(0,1,20), invasion = 'Ethiopia'))

write.csv(rescaled_model_df, '~/Desktop/rescaled_model_df.csv')

dim(rescaled_model_df)

fst_pred_eth = exp(stan2$coefficients[1] + stan2$coefficients[2]*seq(0,1,0.02) + stan2$coefficients[3]*0 +  stan2$coefficients[4]*0 )

fst_pred_sud = stan2$coefficients[1] + stan2$coefficients[2]*seq(0,1,0.02) + stan2$coefficients[3]*1 +  stan2$coefficients[4]*1*seq(0,1,0.02) 


ggplot(rescaled_model_df, aes(x=mean.wind.traveltime,y=fst_lin, colour=invasion))+
  geom_point()+
  theme_classic()

plot(rescaled_model_df$mean.wind.traveltime, rescaled_model_df$fst, ylim=c(0,0.15), xlim=c(0,1))


freq_model_1 <- glmer(fst_lin ~ dist + cattle.circuit +pop.circuit + friction.circuit + evi.circuit + ele.circuit + slope.circuit + mean.wind.traveltime + hurs.circuit + mdr.circuit + (1|invasion), data=rescaled_model_df) 
summary(freq_model_1)



fst_pred = posterior_predict(stan2, newdat = expand.grid(
  invasion = c("Sudan"),
  dist.scaled =  seq(min(model_df$dist.scaled),max(model_df$dist.scaled),length=6),
  wind.scaled = seq(min(model_df$wind.scaled),max(model_df$wind.scaled),length=6),
  friction.scaled = seq(min(model_df$friction.scaled),max(model_df$friction.scaled),length=6)))
  
  
newdat = expand.grid(
  invasion = c("Sudan"),
  dist.scaled =  seq(min(model_df$dist.scaled),max(model_df$dist.scaled),length=6),
  wind.scaled = seq(min(model_df$wind.scaled),max(model_df$wind.scaled),length=6),
  friction.scaled = seq(min(model_df$friction.circuit, na.rm=TRUE),max(model_df$friction.circuit),length=6))

distance.xaxis <- seq(min(model_df$dist.scaled),max(model_df$dist.scaled),length=6)
y1pred.sudan <- fst_pred[1,1:21]
#y1pred.sudan <- fst_pred[1,22:42]

#plot(x=distance.xaxis,y=y1pred.ethiopia)
plot(x=distance.xaxis,y=y1pred.sudan)


#plot model predictions for sudan
plot(x=distance.xaxis,y=y1pred.sudan)
ss = sample(1:2000,100,replace=FALSE)
for(i in 1:length(ss)){
  lines(fst_pred[ss[i],1:21] ~ distance.xaxis)
}


points(model_df$fst_lin ~ model_df$dist.scaled, col="blue")
hist(model_df$fst_lin)

#try logit normal dist
#linearise fst
#include temp
#maybe a beta distribition

