#######
#plot nice map and table for stephensi collections
#tristan dennis 06.11.24
#######

library(tidyverse)
library(sf)
library(RColorBrewer)
library("ggrepel")
library("ggspatial")
library("ggsci")
library("data.table")
library(kableExtra)

#read sample metadata
df_samples = read.csv('~/Projects/AsGARD/metadata/cease_combinedmetadata_noqc.20250212.csv')
#load shapefiles
admin_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")
border_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_10m_admin_0_sovereignty/ne_10m_admin_0_sovereignty.shp")
river_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp")
lake = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster//ne_50m_lakes/ne_50m_lakes.shp")
ocean_shape = st_read("~/Projects/cease/cease_wp1c/analysis_data/raster/ne_50m_ocean/ne_50m_ocean.shp")

#subset border shape to bounding box to avoid mental SVG sizes
sf::sf_use_s2(FALSE)

#define plotting bbox and clip other spatial data accordingly.
bounding_box <- st_bbox(c(xmin = 25, xmax = 85, ymin = 2, ymax = 37), crs = st_crs(border_shape))
border_shape <- st_crop(border_shape, bounding_box)
ocean_shape <- st_crop(ocean_shape, bounding_box)
lakeshape <- st_crop(lake, bounding_box)

#reformat lat and long
df_samples$newlong <- trunc(df_samples$longitude * 100) / 100
df_samples$newlat <- trunc(df_samples$latitude * 100) / 100

#------------------#
# Make nice table #
#------------------#


loc_df <- df_samples %>% select(country,location, admin1_name, latitude, longitude ) %>% unique()
write.csv(loc_df, '~/Projects/AsGARD/locations.csv')

sampdf <- df_samples %>%
  dplyr::select(country, location, pop_code, latitude, longitude) %>%
  unique()

counted_df <- df_samples %>%
  group_by(country,pop_code ) %>%
  count()

# Define colnames and do some ordering for the table
colnames(counted_df) <- c('Country','Cohort', 'Count')

countries_ordered <- c('India','Pakistan','Afghanistan','Iran','SaudiArabia','Yemen','Djibouti','Sudan','Ethiopia')
counted_df$Country <- factor(counted_df$Country, levels=countries_ordered)
counted_df <- counted_df[order(counted_df$Country), ]

# Define colors
cols_table <- c('#96172e', '#f03e5e','#ff7f00','#ff7f00', '#c57fc9','#c27a88', 
                '#6a3d9a', '#cab2d6','#CC7722','#507d2a', '#fccf86','#007272', 
                '#33a02c', '#a6cee3')

counted_df$`Cohort` <- mapply(function(val, color) {
  cell_spec(val, background = color)
}, counted_df$`Cohort`, cols_table)
# Apply cell_spec only to the 'Pop. Code' column
#counted_df$`Pop. Code` <- cell_spec(counted_df$`Pop. Code`, 
#                                    align = 'c', 
#                                    background = cols_table, 
#                                    color = 'white')
# Generate the table
counted_df %>%
  #arrange(Country) %>% 
  kbl(escape = FALSE) %>%
  collapse_rows(columns = 1, valign = "top") %>%
  kable_styling(font_size = 18,  # Increase font size
                full_width = F) %>%  # Optional: to prevent table from stretching to full width
  column_spec(1, width = "2cm,", background = "white") %>%  # Adjust column widths
  column_spec(2, width = "2cm") %>%  #%>%  # White background for 'Country'
  #column_spec(3, width = "2.5cm",) %>%  # White background for 'Country'
 save_kable(file = '~/Projects/AsGARD/figures/poptable.png', zoom = 10) 


##########
#plot nice map of sample locations
##########

#plot range line
steprange <- sf::read_sf('~/Projects/AsGARD/data/EO_Asia_Pacific/stephensi.shp')

pop_cols <- c('#ff7f00','#507d2a','#007272','#33a02c','#a6cee3','#96172e','#f03e5e','#c57fc9','#c27a88','#6a3d9a','#cab2d6','#fccf86','#CC7722')
#get countries included in analysis 
rangecountries = c('SDN', 'ETH', 'IND', 'DJI', 'YEM', 'PAK', 'AFG', 'SAU', 'IRN','IRQ','LKA','KWT','OMN','SOM','SOL','ERI','KEN','MMR','BGD','QAT','ARE')
sample_map <- ggplot()+
  geom_sf(data = border_shape, fill='#e3e3e3')+#fill = ifelse(border_shape$ADM0_A3 %in% rangecountries, "#efe3c3",'white'),col =gray(0.7)) + # geom_sf(data=
#  geom_sf(data=st_geometry(ocean_shape), colour = '#4a80f5', fill='#9bbff4')+
  geom_sf(data=st_geometry(lake),colour = '#4a80f5', fill='#9bbff4')+
  geom_sf(data = steprange, fill="#918e8e",linewidth=0, alpha=0.4) + # geom_sf(data=st_geometry(river_shape), colour = '#4a80f5', fill='#9bbff4', alpha=0.7)+
  geom_point(data=sampdf, aes(x=as.numeric(longitude), y=as.numeric(latitude), colour=pop_code),  alpha = 0.7, size=4)+#colour=query))+#size=sequenced_count))+
  coord_sf(xlim=c(25, 83),ylim=c(3,35.5), expand=FALSE)+
  labs(x="Longitude", y="Latitude", colour='Admin1')+
  annotation_scale(location = "br", width_hint = 0.5) +
  labs(x = "Longitude", y = "Latitude", colour="Population") +
  #geom_text(data=sampdf, aes(x=longitude, y=latitude, label = location))+
  scale_color_manual(values = pop_cols)+
  theme_void() +
  theme(panel.background = element_rect(fill = "#9bbff4"),
        legend.position = 'none',
        #legend.position = "bottom", 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold")) +  theme(legend.key=element_blank())
sample_map

ggsave('~/Projects/AsGARD/figures/samplemap.svg', plot=sample_map, width = 2000, height=1000, units = 'px')








