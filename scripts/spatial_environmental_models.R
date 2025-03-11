#-------------------------------------------------------------------------------#
# This script contains R and Rstan code required to reproduce figure 2          #
# and the associated supplementary material for the AsGARD MS                   #
# Modelling spatial and environmental genomic data for An. stephensi            #
#-------------------------------------------------------------------------------#
# Tristan Dennis, 12.2024

# Import libs
libs <- c('tidyverse','data.table','ggthemes','ggquiver','bayesplot','ggeffects','cowplot','maps','terra','sf','geosphere','rstanarm','loo','projpred',
  'tidybayes','corrplot','posterior','lme4','MuMIn','DHARMa','akima','scales','glmmTMB','dggridR')
lapply(libs, library, character.only = TRUE)

# Read metadata
df_samples = read.csv('~/Projects/AsGARD/metadata/cease_combinedmetadata_noqc.20241228.csv')

# Define some helper variables

# Define the coordinate reference system. Getting this wrong really screws everything up so I've been probably a bit overzealous in specifying the CRS
# Pretty much everywhere where there is the option to - even when many objects in spatial analyses inherit it as a property
crs.geo <- 4326

# Region extent
invasive_coords <- c(xmin=27, ymin=6, ymax=30, xmax=52)

# Sampling locations
points <- df_samples %>% filter(country %in% c('Sudan','Ethiopia','Djibouti')) %>%  dplyr::select(country, location,pop_code, latitude, longitude) %>% unique()

# Let's start by plotting the fEEMS output and sampling data on a map

# Load border data
border_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_0_sovereignty/ne_10m_admin_0_sovereignty.shp")

# Generate a dggs specifying an intercell spacing of ~150 miles
dggs <- dgconstruct(spacing=75, metric=TRUE, resround='nearest', aperture = 4, topology = "TRIANGLE",  projection = "ISEA")

# Make square box covering sampled region in KSA/HoA/Yemen
invasivebbox = sf::st_bbox(invasive_coords, crs = crs.geo)
invasiveshape <- st_as_sfc(invasivebbox)

# WO
write.csv(invasiveshape[[1]][[1]], '/Users/dennistpw/Projects/AsGARD/data/feems_20240920/hoa_ksa_outline.csv', quote = FALSE, row.names = FALSE)
st_write(invasiveshape, '/Users/dennistpw/Projects/AsGARD/data/feems_20240920/hoa_ksa_outline.shp', append=FALSE)

# Make cropped grid
su_grid <- dgshptogrid(dggs, '/Users/dennistpw/Projects/AsGARD/data/feems_20240920/hoa_ksa_outline.shp')
st_write(su_grid, '/Users/dennistpw/Projects/AsGARD/data/feems_20240920/hoa_ksa_outline.TRI.75k.shp', append=FALSE)

# These grids can be used for feems

#######################
# Plot fEEMS output on a map
#######################

# Function for reading feems output
prepare_data <- function(edge_file, node_file){
  edges <- fread(edge_file, col.names = c("from_id", "to_id", "edge_weight"))
  nodes <- fread(node_file, col.names = c("Longitude", "Latitude","N")) %>% mutate(V1 = row_number() - 1)
  
  # Convert necessary columns to integer
  edges$from_id <- as.integer(edges$from_id)
  edges$to_id <- as.integer(edges$to_id)
  nodes$V1 <- as.integer(nodes$V1)
  
  # Join edges and nodes data to get the start and end points of each edge
  edges <- edges %>%
    left_join(nodes, by = c("from_id" = "V1")) %>%
    left_join(nodes, by = c("to_id" = "V1"), suffix = c(".from", ".to")) %>%
    mutate(weight = log10(edge_weight)-mean(log10(edge_weight)))
  
  # Create a list of linestrings, each defined by a pair of points
  edges$geometry <- mapply(function(lon_from, lat_from, lon_to, lat_to) {
    st_linestring(rbind(c(lon_from, lat_from), c(lon_to, lat_to)))
  }, edges$Longitude.from, edges$Latitude.from, edges$Longitude.to, edges$Latitude.to, SIMPLIFY = FALSE)
  
  # Convert edges to an sf object
  edges_sf <- st_as_sf(edges)
  
  # Set the CRS
  st_crs(edges_sf) <- crs.geo
  
  # Convert nodes data.table to an sf object
  nodes_sf <- st_as_sf(nodes, coords = c("Longitude", "Latitude"), crs = crs.geo)
  
  list(edges_sf = edges_sf, nodes_sf = nodes_sf)
}

# Load feems output
edge_file = '/Users/dennistpw/Projects/AsGARD/data/feems_20240920/incksaedgew.csv'
node_file = '/Users/dennistpw/Projects/AsGARD/data/feems_20240920/incksanodepos.csv'

edges <- fread(edge_file, col.names = c("from_id", "to_id", "edge_weight"))
nodes <- fread(node_file, col.names = c("Longitude", "Latitude","N")) %>% mutate(V1 = row_number() - 1)

# Convert necessary columns to integer
edges$from_id <- as.integer(edges$from_id)
edges$to_id <- as.integer(edges$to_id)
nodes$V1 <- as.integer(nodes$V1)

# Join edges and nodes data to get the start and end points of each edge
edges <- edges %>%
  left_join(nodes, by = c("from_id" = "V1")) %>%
  left_join(nodes, by = c("to_id" = "V1"), suffix = c(".from", ".to")) %>%
  mutate(weight = log10(edge_weight)-mean(log10(edge_weight)))

# Create a list of linestrings, each defined by a pair of points
edges$geometry <- mapply(function(lon_from, lat_from, lon_to, lat_to) {
  st_linestring(rbind(c(lon_from, lat_from), c(lon_to, lat_to)))
}, edges$Longitude.from, edges$Latitude.from, edges$Longitude.to, edges$Latitude.to, SIMPLIFY = FALSE)

# Convert edges to an sf object
edges_sf <- st_as_sf(edges)

# Set the CRS
st_crs(edges_sf) <- 4326

# Convert nodes data.table to an sf object
nodes_sf <- st_as_sf(nodes, coords = c("Longitude", "Latitude"), crs = 4326)

# Load shapefiles gor plotting
admin_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")
border_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_0_sovereignty/ne_10m_admin_0_sovereignty.shp")
river_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp")
lake = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster//ne_50m_lakes/ne_50m_lakes.shp")
ocean_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_50m_ocean/ne_50m_ocean.shp")

# Make new lat and long for more concise plotting
df_samples$newlong <- trunc(df_samples$longitude * 100) / 100
df_samples$newlat <- trunc(df_samples$latitude * 100) / 100

# Load collection site status data from WP1a
site_status <- fread('~/Projects/AsGARD/data/feems_20240920/status_table.csv')

# Plot with feems weighted graph grid on top
world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))

# Plotting below

# Now let's analyse some spatial data
# Start by loading the data on diversity and differentiation

fst <- fread('~/Projects/AsGARD/data/fst_dist.csv')
roh <- fread('~/Projects/AsGARD/data/roh_portdist.csv')
pi <- fread('Projects/AsGARD/data/pi_portdist.csv')

#now let's do a test for isolation by distance (a mantel test is often appropriate but these data were generated using two triangular rather than two square matrices so lm prob ok too?)
#fst %>% filter(country.x != 'AdenCity' & country.y != 'AdenCity')
#subset df_samples to include only samples from HoA
df_samples_hoa <- df_samples[df_samples$country %in% c('Djibouti','Sudan','Ethiopia'),]

#make 'invasive pop' factor
df_samples$InvasivePop <- ifelse(df_samples$country == 'DjiboutiCity' | df_samples$country == 'Ethiopia', 'Ethiopia & Djibouti','Sudan')

#quick look at structuring of fst by country
locs <- df_samples %>% dplyr::select(location, country, InvasivePop) %>% unique()

fst <- left_join(fst, locs, by=c('popa'='location')) %>% left_join(., locs, by=c('popb'='location'))
fst$countrycomp <- paste0(fst$country.x, ':',fst$country.y)

fst$isdiffcountry = ifelse(fst$country.x == fst$country.y,'samecountry','diffcountry')
fst$putativeinvasion = ifelse(fst$InvasivePop.x == fst$InvasivePop.y,'Same','Different')

# Transform Fst to unbound and normalise (Slatkin 1997 I think)
fst$fst_lin <- fst$fst / (1-fst$fst)


fst_plot<- fst %>% filter(country.x != 'Yemen' & country.y != 'Yemen') %>% 
  ggplot(., aes(x=dist,y=fst_lin, colour=countrycomp))+
  geom_point(size=4,pch=20, alpha=0.6)+
  theme_classic()+
  scale_color_brewer(palette = 'Paired')
fst_plot

fst %>% filter(country.x != 'Yemen' & country.y != 'Yemen') %>% 
  ggplot(., aes(y=fst, fill=as.factor(countrycomp)))+
  geom_histogram()+
  facet_grid(~countrycomp)


############
# Fit linear models of:
# Isolation by distance
# Genetic diversity as a function from distance to main port
###########

#function for plotting residuals
plotresids <- function(model){
  modsimout <- simulateResiduals(fittedModel = model, plot = F)
  residuals(modsimout, quantileFunction = qnorm, outlierValues = c(0,1))
  plot(modsimout)
}

# Model for Fst as function of pairwise pop dist
fst_model <- glm(fst_lin ~ dist, data=fst)
summary(fst_model)
plotresids(fst_model) #this is fine. prob a little heteroscedasticity 

# Model for Fst as function of dist with ineraction for invasive pop
fst_model_int <- glm(fst_lin ~ dist*putativeinvasion, data=fst)
summary(fst_model_int)
plotresids(fst_model_int) 

# Model for Fst as function of dist with fixed effect for invasive pop
fst_model_fix <- glm(fst_lin ~ dist+putativeinvasion, data=fst)
summary(fst_model_fix)
plotresids(fst_model_fix) 

# Check AIC
AIC(fst_model, fst_model_fix, fst_model_int)
#Fst model fix has lower AIC, ANOVA rejects interac

#Get R-Squared
r.squaredGLMM(fst_model)
r.squaredGLMM(fst_model_fix)
r.squaredGLMM(fst_model_int)

# Plot below

#very nice, let's take a look at roh and pi

#firts calculate the distance between points and djibouti
# Apply distHaversine to each row of the dataframe
roh$dist_to_nearest_port <- mapply(function(lat1, lon1, lat2, lon2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2))
}, roh$latitude, roh$longitude, roh$lat_nearest_port, roh$long_nearest_port)
roh$dist_to_nearest_port <- roh$dist_to_nearest_port / 1000

#do for pi
# Apply distHaversine to each row of the dataframe
pi$dist_to_nearest_port <- mapply(function(lat1, lon1, lat2, lon2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2))
}, pi$latitude, pi$longitude, pi$lat_nearest_port, pi$long_nearest_port)
pi$dist_to_nearest_port <- pi$dist_to_nearest_port / 1000

#now run some models
#coerce to factor
pi$InvasionPop <- as.factor(pi$InvasionPop)

#now let's do pi
hist(pi$pi)
#pi with gam dist
mod_pi_gam = glm(formula = pi ~ dist_to_nearest_port, data = pi, family=Gamma(link = "log"))
mod_pi_fix_gam = glm(formula = pi ~ dist_to_nearest_port+InvasionPop, data = pi, family=Gamma(link = "log"))
mod_pi_int_gam = glm(formula = pi ~ dist_to_nearest_port*InvasionPop, data = pi, family=Gamma(link = "log"))
AIC(mod_pi_gam, mod_pi_fix_gam, mod_pi_int_gam)
anova(mod_pi_gam, mod_pi_fix_gam, mod_pi_int_gam) #fixed effect model is adequate for roh

####roh
#model fitting

#filter djib
#roh <- dplyr::filter(roh, country != 'Djibouti')
hist(roh$fROH) #normal looking dist, bit of right skew?
roh$InvasionPop = as.factor(roh$InvasionPop)
roh_baseinroh <- roh$fROH * (93706023+88747589+22713616) #coerce to number of bases in the genome so that can be modelled as gaussian
mod_roh_gaus <- glm(formula = roh_baseinroh ~ dist_to_nearest_port, data = roh, family="gaussian")
mod_roh_fix_gaus <- glm(formula = roh_baseinroh ~ dist_to_nearest_port+InvasionPop, data = roh, family="gaussian")
plotresids(mod_roh_fix_gaus) #looks not terrible but not perfect. let's try some different distributions. roh data are between 0 and 1 and slightly skewed right so let's try binom logit link
mod_roh_int_gaus <- glm(formula = roh_baseinroh ~ dist_to_nearest_port*InvasionPop, data = roh, family="gaussian")
plotresids(mod_roh_int_gaus) #llooks best so far
AIC(mod_roh_gaus, mod_roh_fix_gaus, mod_roh_int_gaus) #fixed effect model is adequate for roh
anova(mod_roh_gaus, mod_roh_fix_gaus, mod_roh_int_gaus) #fixed effect model is adequate for roh, but very slightly, and residuals look better for interaction model

r.squaredGLMM(mod_roh_int_gaus) # better rsquared too


################
#plot model outputs
#################
#roh_predictions = predict_response(mod_roh_int_gaus, terms=~dist_to_nearest_port*InvasionPop)
#pi_predictions = predict_response(mod_pi_int_gam, terms=~dist_to_nearest_port*InvasionPop)

# For ROH
# Generate new data for predictions

# Generate new data for predictions
roh_pred <- expand.grid(
  dist_to_nearest_port = seq(min(roh$dist_to_nearest_port, na.rm = TRUE), 
                             max(roh$dist_to_nearest_port, na.rm = TRUE), 
                             length.out = 100),
  InvasionPop = unique(roh$InvasionPop)
)

# Get model predictions with 95% confidence intervals
preds <- predict(mod_roh_int_gaus, roh_pred, type = "response", se.fit = TRUE)
roh_pred$fit <- preds$fit
roh_pred$lwr <- preds$fit - 1.96 * preds$se.fit  # Lower 95% CI
roh_pred$upr <- preds$fit + 1.96 * preds$se.fit  # Upper 95% CI

# Ensure InvasionPop is a factor (if categorical)
roh_pred$InvasionPop <- as.factor(roh_pred$InvasionPop)
roh$InvasionPop <- as.factor(roh$InvasionPop)



# Generate new data for predictions - PI
pi_pred <- expand.grid(
  dist_to_nearest_port = seq(min(pi$dist_to_nearest_port, na.rm = TRUE), 
                             max(pi$dist_to_nearest_port, na.rm = TRUE), 
                             length.out = 100),
  InvasionPop = unique(pi$InvasionPop)
)

# Get model predictions with 95% confidence intervals
preds <- predict(mod_pi_int_gam, roh_pred, type = "response", se.fit = TRUE)
pi_pred$fit <- preds$fit
pi_pred$lwr <- preds$fit - 1.96 * preds$se.fit  # Lower 95% CI
pi_pred$upr <- preds$fit + 1.96 * preds$se.fit  # Upper 95% CI

# Ensure InvasionPop is a factor (if categorical)
pi_pred$InvasionPop <- as.factor(pi_pred$InvasionPop)
pi$InvasionPop <- as.factor(pi$InvasionPop)

# Plot below
# Plot Fst vs distance
fst_modelplot<- ggplot(fst, aes(x = dist/1000, y = fst, color = putativeinvasion)) +
  geom_point(alpha=0.6) +
  theme_classic()+
  geom_smooth(method = "lm", se = TRUE, ) +
  scale_y_continuous(breaks=c(0, 0.1, 0.2))+
  labs(colour= 'Invasion',
       x = "Distance (km)",
       y = "Fst/(1-Fst)") +
  #scale_x_continuous(breaks = extended_range_breaks()(fst$dist/1000)) +
  #scale_y_continuous(breaks = extended_range_breaks()(fst$fst),labels = number_format(accuracy =0.01)) +
  scale_color_brewer(palette = "Set1")+
  theme(axis.text = element_text(size = 10, family = "Arial"),
        axis.title = element_text(size=12, family = "Arial"),
        plot.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 10, family='Arial'),  # Increase legend text size
        legend.title = element_text(size = 12, family='Arial'))

# Plot roh
rohpred <- ggplot() +
  # Original data points
  geom_point(data = roh, aes(x = dist_to_nearest_port, y = roh_baseinroh, color = InvasionPop), alpha = 0.6) +
  # Model fit lines
  geom_line(data = roh_pred, aes(x = dist_to_nearest_port, y = fit, color = InvasionPop), size = 1) +
  # Confidence interval ribbon
  geom_ribbon(data = roh_pred, aes(x = dist_to_nearest_port, ymin = lwr, ymax = upr, fill = InvasionPop), 
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(x = "Distance to Nearest Port", y = "bp in ROH")+
  scale_color_manual(values = c("darkblue", "lightblue")) +  # Customize colors
  scale_fill_manual(values = c("darkblue", "lightblue")) +
  #theme(plot.margin = margin(2, 1.5, 1.5, 1.5))+
  theme(axis.text = element_text(size = 10, family = "Arial"),
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        axis.title = element_text(size=12, family = "Arial"),
        legend.text = element_text(size = 10, family='Arial'),  # Increase legend text size
        legend.title = element_text(size = 12, family='Arial'))+
  theme_classic()





# Plot pi
pipred <- ggplot() +
  # Original data points
  geom_point(data = pi, aes(x = dist_to_nearest_port, y = pi, color = InvasionPop), alpha = 0.6) +
  # Model fit lines
  geom_line(data = pi_pred, aes(x = dist_to_nearest_port, y = fit, color = InvasionPop), size = 1) +
  # Confidence interval ribbon
  geom_ribbon(data = pi_pred, aes(x = dist_to_nearest_port, ymin = lwr, ymax = upr, fill = InvasionPop), 
              alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(x = "Distance to Nearest Port", y = "π")+
  scale_y_continuous(breaks=c(0.007, 0.008,0.009, 0.01), labels = scientific_format())+
  scale_color_manual(values = c("darkblue", "lightblue")) +  # Customize colors
  scale_fill_manual(values = c("darkblue", "lightblue")) +
  labs(y='π',x='Dist. to main port (km)', colour="Invasion\nPopulation", title="")+
  #theme(plot.margin = margin(2, 1.5, 1.5, 1.5))+
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        axis.title = element_text(size=12),
        legend.text = element_text(size = 10),  # Increase legend text size
        legend.title = element_text(size = 12))+
  theme_classic()

# Load ggrastr for good plotting
# Set dpi
options("ggrastr.default.dpi" = 500) 
library(ggrastr)
sample_map <- ggplot()+
  rasterize(geom_sf(data = border_shape, fill = '#f7f7f7',col =gray(0.1))) + 
  geom_sf(data=st_geometry(lake),colour = '#4a80f5', fill='#9bbff4')+
  geom_sf(data=st_geometry(river_shape),colour = '#4a80f5', fill='#9bbff4')+
  rasterize(geom_sf(data=edges_sf, aes(colour=log10(edges_sf$edge_weight)), alpha=0.7))+
  scale_color_gradient2(low = "#ad3c07",mid='white', high = "#4dc1ff") +       # Color gradient for weight
  geom_point(data=site_status, aes(x=as.numeric(long), y=as.numeric(lat), fill=Status), pch=21, size=4)+
  scale_fill_manual(values = c('#808080','#2e8fff','#fcba03'))+
  # coord_sf(xlim=c(28, 52),ylim=c(4,28), expand=FALSE)+
  labs(colour="log10(w)", x='Long.', y='Lat.') +
  #theme_void() +
  coord_sf(xlim = c(28, 52), ylim = c(4, 28), expand = FALSE)+
theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        #legend.position = 'none',
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "top", 
        axis.line = element_blank(),
        axis.text = element_text(size = 10),  # Increase legend text size
        #axis.title = element_blank(),
        legend.text = element_text(size = 10),  # Increase legend text size
        legend.title = element_text(size = 12)) +  
  theme(legend.key=element_blank())

ggsave(filename = '~/Projects/AsGARD/figures/feems_map.svg', sample_map, width=7, height=7)



# Align right-hand plots properly
modelplots <- cowplot::plot_grid(fst_modelplot, pipred, rohpred, 
                                 ncol = 1, 
                                 align = "v", 
                                 axis = "tb",
                                 rel_heights = c(1, 1.15, 1)) # Ensures equal height
ggsave(filename = '~/Projects/AsGARD/figures/modelplots.svg', modelplots, width=4, height=7)

# Fix the overall layout
final_plot <- cowplot::plot_grid(sample_map, modelplots, 
                                 ncol = 2, 
                                 rel_widths = c(2, 1),  # Adjust widths if needed
                                 align = "hv",  # Ensures vertical & horizontal alignment
                                 axis = "tblr")


print(final_plot)



r.squaredGLMM(mod_pi_int_gam)

r.squaredGLMM(mod_roh_int_gaus)


