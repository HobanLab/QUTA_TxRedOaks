# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% QUERCUS TARDIFOLIA MAPPING SCRIPT %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is used to make a map of the samples used in the Q. tardifolia Texas red oak analyses.
# It reads in a spreadsheet that contains latitude and longitude values, and then uses leaflet to build
# a map based on those values. Points are colored based on their taxonomic designation (which, note, has
# been significantly cleaned up from original designations...)

# Load necessary libraries
pacman::p_load(dplyr, rnaturalearth, rnaturalearthdata, readr, leaflet, htmlwidgets)

# Specify filepath to sample coordinates, and read in
QUTA_coordinateDir <- '/home/akoontz/Documents/QUTA_TxRedOaks/Code/D_Mapping/QUTA_TxRedOak_MappingList.csv'
QUTA_coordinateData <- read_csv(QUTA_coordinateDir)
# Filter out rows with NA Latitude or Longitude
QUTA_coordinateData <- QUTA_coordinateData %>% filter(!is.na(LatitudeUnprotected) & !is.na(LongitudeUnprotected))

# # Define a color palette based only on the remaining species
species_present <- unique(QUTA_coordinateData$TaxaUpdated)
# Specify a vector of color names. The grays are meant to represent taxa which are less frequent/central to analyses
QUTA_colors <- c('coral1', 'antiquewhite', 'blue4', 'chartreuse4', 'chartreuse', 'cadetblue1', 'darkorchid', 'bisque',
                 'darksalmon', 'darkorange4', 'darkolivegreen', 'darkgoldenrod2', 'gray')
# Create a named color vector (species -> color mapping)
species_colors <- setNames(QUTA_colors, species_present)
# Define a color function for Leaflet
species_pal <- colorFactor(palette = species_colors, domain = species_present)

# Create the Leaflet map
QUTA_map <- leaflet(QUTA_coordinateData) %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%  # Use a nice basemap
  addCircleMarkers(
    ~LongitudeUnprotected, ~LatitudeUnprotected,
    color = ~species_pal(TaxaUpdated),
    radius = 6,
    stroke = TRUE,
    fillOpacity = 0.8,
    popup = ~paste("<b>Sample:</b>", SampleName, "<br><b>Species:</b>", TaxaUpdated),
    group = ~TaxaUpdated  # Grouping for layer control
  ) %>%
  addLegend(
    position = "topright",  # Initially place legend at top right
    pal = species_pal,
    values = species_present,  # Only include present species
    title = "Taxa",
    opacity = 1,
    layerId = "legend"
  ) %>%
  addLayersControl(
    overlayGroups = species_present,  # Only include present species
    options = layersControlOptions(collapsed = FALSE),
    position = "topright"  # Place layers control at top right
  )

# Use JavaScript to adjust the position of the legend
QUTA_map <- onRender(QUTA_map, "
  function(el, x) {
    var legend = document.querySelector('.leaflet-top.leaflet-right .leaflet-control');
    var layersControl = document.querySelector('.leaflet-control-layers');
    if (legend && layersControl) {
      legend.style.marginTop = layersControl.clientHeight + 'px';
    }
  }
")

# Display the map
QUTA_map
