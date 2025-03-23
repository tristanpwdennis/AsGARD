library(terra)
library(windscape)

# Load wind data downloaded from ERA5
wind<-rast("~/Projects/AsGARD/data/environmental_data_rasters/wind_temp/data_stream-oper_stepType-instant.nc")

# Create wind series
series <- wind_series(wind, order = "uuvv")

# Create wind rose, this takes a fair bit
rose <- wind_rose(series, trans = 1)

# Plot wind rose
plot(rose)

# Get points
df_samples = read.csv('~/Projects/AsGARD/metadata/cease_combinedmetadata_qcpass.20240914.csv')
points <- df_samples %>% 
            filter(country  %in% c("Sudan","Ethiopia","Djibouti")) %>% 
            dplyr::select(country, location,analysis_pop, latitude, longitude) %>% 
            unique()
sites <- cbind(x=points$longitude,
               y=points$latitude)

# Check cell overlap and distance
check_cell_distance(rose, sites)
# Looks like we have only one site pair in the same grid cell (prob Khartoum)
# These aren't huge errors but they'll add some noise to our wind estimates
# Let's downscale by a factor of 5 and see what we get in terms of error rates
rose <- downscale(rose, 5)
check_cell_distance(rose, sites) # Looks better

# Now we build the wind graph, and calculate the LCP between sites
downwind <- wind_graph(rose, direction = "downwind")
upwind <- wind_graph(rose, direction = "upwind")
wind_time <- least_cost_distance(downwind, sites)

# Get estimates of pairwise mean tind travel time
wind_conn <- pairwise_means(wind_time)

#wo
write.csv(wind_conn,'~/Projects/AsGARD/data/environmental_modelling/least_cost_path_data/wind_time_means.txt')

# Now let's try and plot wind average data
# Separate the stack into u and v components
n_layers <- nlyr(wind)
u_stack <- wind[[1:(n_layers / 2)]]
v_stack <- wind[[(n_layers / 2 + 1):n_layers]]
u_stack

# Compute layer-wise averages
u_avg <- mean(u_stack)
v_avg <- mean(v_stack)

# Create a data frame of values and coordinates
u_df <- as.data.frame(u_avg, xy = TRUE, na.rm = TRUE)
v_df <- as.data.frame(v_avg, xy = TRUE, na.rm = TRUE)

# Merge u and v data
wind_df <- merge(u_df, v_df, by = c("x", "y"), suffixes = c("_u", "_v"))
writeRaster(rast(wind_df), '~/Projects/AsGARD/data/environmental_data_rasters/wind_mean.tif')
# Normalize Wind Vectors for Consistent Arrow Length
#wind_df <- wind_df %>%
#  mutate(
#    magnitude = sqrt(mean_u^2 + mean_v^2),  # Compute wind speed
#    norm_u = mean_u / magnitude,            # Normalize U
#    norm_v = mean_v / magnitude             # Normalize V
 # )

wind_df <- wind_df %>%
  mutate(
    magnitude = sqrt(mean_u^2 + mean_v^2),  # Compute wind speed
    scaled_u = mean_u / max(magnitude) * magnitude,  # Scale U by intensity
    scaled_v = mean_v / max(magnitude) * magnitude   # Scale V by intensity
  )

#subsample for ease
set.seed(42)  # For reproducibility
wind_df_sampled <- wind_df %>%
  sample_frac(0.5)  # Sample 10% of the points

ggplot(wind_df_sampled, aes(x = x, y = y)) +
  geom_segment(
    aes(
      xend = x + scaled_u, 
      yend = y + scaled_v, 
      size = magnitude  # Thickness proportional to wind magnitude
    ),
    arrow = arrow(length = unit(0.05, "cm")),  # Adjust arrowhead size
    color = "darkgrey", alpha = 0.2           # Dark grey with transparency
  ) +
  scale_size_continuous(range = c(0.2, 1)) +  # Adjust thickness range
  labs(title = "Wind Quiver Plot", x = "Longitude", y = "Latitude") +
  theme_minimal()

