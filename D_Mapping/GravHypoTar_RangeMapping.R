# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GRAVESII AND HYPOLEUCOIDES RANGE MAPS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is used to make a map of the ranges of two of the species in
# the study (Q. gravesii and hypoleucoides). Overlaid on this map are points
# corresponding to the samples included in the study. Two versions of the map
# are generated: one using hillshade, and one using coloring to represent the DEM.
# Both of these are are saved as PNG files to the disk.

# Load necessary libraries
pacman::p_load(sf, terra, rgbif, ggplot2, dplyr, geojsonsf, rnaturalearth,
               rnaturalearthdata, elevatr, ggnewscale, patchwork, cowplot)
# Specify and set filepath
QUTA_mappingDir <- '/home/akoontz/Documents/QUTA_TxRedOaks/Code/D_Mapping/'
setwd(QUTA_mappingDir)
# Specify image output directory
imageOut <- 
  '/home/akoontz/Documents/QUTA_TxRedOaks/Documentation/Images/GravHypoTar-RangeMaps_20251223/'
# Specify width and length dimensions for images generated
imageWidth <- 1300
imageLength <- 700

# %%% DOWNLOAD SHAPEFILES FOR GRAVESII/HYPOLEUCODIES RANGES ----
# This code call commands for downloading geojson files corresponding to the
# ranges of Quercus gravesii and hypoleucoides. These shape files were sourced
# from this GitHub repo: https://github.com/wpetry/USTreeAtlas?tab=readme-ov-file
# These commands only need to be run once, hence the file.exists
if(!file.exists("QUGR.geojson")){
  download.file(
    url = "https://github.com/wpetry/USTreeAtlas/blob/main/geojson/quergrav.geojson",
    destfile = "QUGR.geojson"
  )
}
if(!file.exists("QUHY.geojson")){
  download.file(
    url = "https://github.com/wpetry/USTreeAtlas/blob/main/geojson/querhypo.geojson",
    destfile = "QUHY.geojson"
  )
}

# %%% LOAD RANGE POLYGONS AND POINTS ----
# Read in polygons (st_read), confirming consistent CRS (st_transform)
QUGR_range <- st_transform(st_read("QUGR.geojson"), 4326)
QUHY_range <- st_transform(st_read("QUHY.geojson"), 4326)
# Read in CSV of sample occurrences. This is spreadsheet containing the 
# occurrences of gravesii, hypoleucoides, and tardifolia individuals. 

# NOTE: this spreadsheet does NOT include known gravesii and hypoleucoides
# hybrids, or gravesii individuals found on limestone or "Langtry's" specimens.
# Additionally, samples with coordinates of NA have been removed.
samples <- read.csv("GravHypoTar_RangeMapSamples.csv")
# Convert data.frame into sf object
samples_sf <- st_as_sf(
  samples,
  coords = c("Long", "Lat"),
  crs = 4326
)
# Specify order of taxa
samples_sf$Taxon <- factor(
  samples_sf$Taxon,
  levels = c("gravesii", "hypoleucoides", "hybrid")
)
# Combine ranges into a single object
ranges <- bind_rows(
  QUGR_range %>% mutate(species = "Quercus gravesii"),
  QUHY_range %>% mutate(species = "Quercus hypoleucoides")
)

# %%% GENERATE BASEMAP ----
# Specify countries, for political boundaries and DEM limits
countries <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) |>
  dplyr::filter(iso_a3 %in% c("USA", "MEX")) |>
  st_transform(4326)
# Check if DEM has already been downloaded and processed; if it has, just read
# in the existin .Rdata file. Rasterizing the DEM is time-consuming, hence these commands.
if(file.exists('DEM.Rdata')){
  dem <- readRDS('DEM.Rdata')
} else{
  # Download DEM
  dem <- get_elev_raster(
    locations = countries,
    z = 6,            # ~1 km resolution
    clip = "locations"
  )
  # Rasterize and crop DEM, then convert to data.frame
  dem <- rast(dem)
  dem <- crop(dem, ext(-115, -95, 25, 38))
  # Save original DEM to disk
  saveRDS(dem, file='DEM.Rdata')
}
# Generate object for political boundaries
countries <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_transform(4326)

# %%% CROP TO FOCUS AREA ----
# Specify a bounding box, with an optional buffer around it
sf::sf_use_s2(FALSE)
bbox <- st_bbox(samples_sf)
bbox <- bbox + c(-2.4, -2.4, 2.4, 2.4)
bbox_sfc <- st_as_sfc(bbox)
# Crop all spatial layers
countries_crop <- st_crop(countries, bbox)
ranges_crop    <- st_crop(ranges, bbox)
dem_crop <- crop(
  dem,
  ext(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])
)

# %%% BUILD HILLSHADE LAYER OFF OF DEM ----
# Ensure DEM is projected (hillshade requires projected CRS)
dem_proj <- project(dem, "EPSG:5070")  # NAD83 / Conus Albers
# Derive slope and aspect
slope  <- terrain(dem_proj, v = "slope", unit = "radians")
aspect <- terrain(dem_proj, v = "aspect", unit = "radians")
# Hillshade
hillshade <- shade(slope, aspect, angle = 45, direction = 315)
# Reproject back to WGS84 for plotting
hillshade <- project(hillshade, "EPSG:4326")
# Crop to study extent
hillshade <- crop(
  hillshade,
  ext(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])
)
# Convert to data.frame
hillshade_df <- as.data.frame(hillshade, xy = TRUE)
colnames(hillshade_df) <- c("lon", "lat", "hillshade")
# Specify DEM background
dem_df <- as.data.frame(dem_crop, xy = TRUE)
colnames(dem_df) <- c("lon", "lat", "elevation")

# %%% PLOTTING, VERSION A: WITH HILLSHADE %%% ----
# Specify colors for sampling points
sample_colors <- c(
  gravesii = "darkred",
  hypoleucoides = "darkblue",
  hybrid = "green"
)

# Primary plotting call
main_map <- ggplot() +
  # Elevation background
  geom_raster(
    data = hillshade_df,
    aes(lon, lat, fill = hillshade),
    alpha = 0.6
  ) +
  scale_fill_gradient(
    low = "white",
    high = "grey20",
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  # Country boundaries
  geom_sf(
    data = countries_crop,
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  # Species ranges
  geom_sf(
    data = ranges_crop,
    aes(fill = species),
    color = "black",
    alpha = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "Quercus gravesii" = "firebrick3",
      "Quercus hypoleucoides" = "steelblue3"
    ),
    name = "Species range"
  ) +
  # Sample points (change shapes AND colors); specify non-hybrids 
  # first, to make sure hybrids are visible
  geom_sf(
    data = dplyr::filter(samples_sf, Taxon != "hybrid"),
    aes(shape = Taxon, color = Taxon),
    size = 2.75
  ) +
  # Specify hybrids, with slightly larger size
  geom_sf(
    data = dplyr::filter(samples_sf, Taxon == "hybrid"),
    aes(shape = Taxon, color = Taxon),
    size = 3.25
  ) +
  # Use specified colors
  scale_color_manual(
    values = sample_colors,
    name = "Taxon"
  ) +
  # Customize plotting characters for different taxa
  scale_shape_manual(
    values = c(
      gravesii = 16,        # filled circle
      hypoleucoides = 17,   # filled triangle
      hybrid = 12            # square
    ),
    name = "Taxon"
  ) +
  # Axes labels and theme
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  # Adjust spacing of the legends, bringing them closer to the main map
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(5, "pt") 
  ) +
  # Increase sizes of points in the legend
  guides(
    shape = guide_legend(override.aes = list(size = 5))
  )

# Adjust text sizes, in different parts of the figure
main_map <- main_map +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13)
  )

# Inset map
inset_map <- ggplot() +
  geom_sf(
    data = countries,
    fill = "grey90",
    color = "black",
    linewidth = 0.3
  ) +
  geom_sf(
    data = bbox_sfc,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  coord_sf(xlim = c(-125, -85), ylim = c(15, 50)) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

# Specify margins around inset
inset_map <- inset_map +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(t = -135, r = 5, b = 7, l = 5)
  )

# Extract the legends from the main map
legend_grob <- get_legend(main_map)
main_map_nolegend <- main_map +
  theme(legend.position = "none")
# Specify margins around main map
main_map_nolegend <- main_map_nolegend +
  theme(
    plot.margin = margin(t = 4, r = 0, b = 5, l = -20)
  )

# Stack the legends and the inset (to build a legend column). First value in rel_heights adjusts inset height
legend_column <- plot_grid(
  legend_grob,
  inset_map,
  ncol = 1,
  rel_heights = c(2.6, 1)
)
# Adjusting spacing
legend_column <- legend_column +
  theme(
    plot.margin = margin(t = -120, r = 5, b = 5, l = -120)
  )
# Combined primary map and inset map (with primary above)
final_plot <- plot_grid(
  main_map_nolegend,
  legend_column,
  ncol = 2,
  rel_widths = c(5.2, 0.8)
)

final_plot

# Generate the final plot, as PDF
png(file = paste0(imageOut, "GravHypoTar_Hillshade_RangeMap.png"), 
    width = imageWidth, height = imageLength)
final_plot
dev.off()

# %%% PLOTTING, VERSION B: WITHOUT HILLSHADE (DEM) %%% ----
# Specify colors for sampling points
sample_colors <- c(
  gravesii = "darkred",
  hypoleucoides = "darkblue",
  hybrid = "green"
)

# Primary plotting call
main_map <- ggplot() +
  geom_raster(
  data = dem_df,
  aes(lon, lat, fill = elevation),
  alpha = 0.8
  ) +
  # Specify elevation gradient coloration
  scale_fill_gradientn(
    colors = c("grey95", "grey80", "grey60", "grey40"),
    name = "Elevation (m)"
  ) +
  ggnewscale::new_scale_fill() +
  # Country boundaries
  geom_sf(
    data = countries_crop,
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  # Species ranges
  geom_sf(
    data = ranges_crop,
    aes(fill = species),
    color = "black",
    alpha = 0.35
  ) +
  scale_fill_manual(
    values = c(
      "Quercus gravesii" = "firebrick3",
      "Quercus hypoleucoides" = "steelblue3"
    ),
    name = "Species range"
  ) +
  # Sample points (change shapes AND colors); specify non-hybrids 
  # first, to make sure hybrids are visible
  geom_sf(
    data = dplyr::filter(samples_sf, Taxon != "hybrid"),
    aes(shape = Taxon, color = Taxon),
    size = 2.75
  ) +
  # Specify hybrids, with slightly larger size
  geom_sf(
    data = dplyr::filter(samples_sf, Taxon == "hybrid"),
    aes(shape = Taxon, color = Taxon),
    size = 3.25
  ) +
  # Use specified colors
  scale_color_manual(
    values = sample_colors,
    name = "Taxon"
  ) +
  # Customize plotting characters for different taxa
  scale_shape_manual(
    values = c(
      gravesii = 16,        # filled circle
      hypoleucoides = 17,   # filled triangle
      hybrid = 12            # square
    ),
    name = "Taxon"
  ) +
  # Axes labels and theme
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  # Adjust spacing of the legends, bringing them closer to the main map
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(5, "pt") 
  ) +
  # Increase sizes of points in the legend
  guides(
    shape = guide_legend(override.aes = list(size = 5))
  )

# Adjust text sizes, in different parts of the figure
main_map <- main_map +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 13)
  )

# Inset map
inset_map <- ggplot() +
  geom_sf(
    data = countries,
    fill = "grey90",
    color = "black",
    linewidth = 0.3
  ) +
  geom_sf(
    data = bbox_sfc,
    fill = NA,
    color = "black",
    linewidth = 0.8
  ) +
  coord_sf(xlim = c(-125, -85), ylim = c(15, 50)) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

# Specify margins around inset
inset_map <- inset_map +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(t = -135, r = 5, b = 7, l = 5)
  )

# Extract the legends from the main map
legend_grob <- get_legend(main_map)
main_map_nolegend <- main_map +
  theme(legend.position = "none")
# Specify margins around main map
main_map_nolegend <- main_map_nolegend +
  theme(
    plot.margin = margin(t = 4, r = 0, b = 5, l = -20)
  )

# Stack the legends and the inset (to build a legend column). First value in rel_heights adjusts inset height
legend_column <- plot_grid(
  legend_grob,
  inset_map,
  ncol = 1,
  rel_heights = c(2.6, 1)
)
# Adjusting spacing
legend_column <- legend_column +
  theme(
    plot.margin = margin(t = -120, r = 5, b = 5, l = -120)
  )
# Combined primary map and inset map (with primary above)
final_plot <- plot_grid(
  main_map_nolegend,
  legend_column,
  ncol = 2,
  rel_widths = c(5.2, 0.8)
)

final_plot

# Generate the final plot, as PDF
png(file = paste0(imageOut, "GravHypoTar_DEM_RangeMap.png"), 
    width = imageWidth, height = imageLength)
final_plot
dev.off()
