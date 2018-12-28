# load necessary R packages 
library(rhdf5)
library(rgdal)
library(raster)
library(ggplot2)
library(tidyr)
library(sf)
library(dplyr)
library(data.table) 

# set working directory
setwd("~/github/jubilant-waffle/code/")

# load any local functions in external files 
source("supporting_functions.R")

# code for NEON site 
site_code <- 'NIWO'

# directory with hyperspectral .h5 files
h5_dir <- paste0('../data/', site_code, '/hyperspectral/')

# directory with shapefiles (tree stem locations and crown polygons)
shapefile_dir <- paste0('../data/', site_code, '/shapefiles/')

# define the output directory. If it doesn't exist already, create it.
check_create_dir('../output/') # create top level "output" directory
out_dir <- paste0('../output/', site_code, '/')
check_create_dir(out_dir) # create output folder for site

# define the "bad bands" wavelength ranges in nanometers, where atmospheric 
# absorption creates unreliable reflectance values. 
bad_band_window_1 <- c(1340, 1445)
bad_band_window_2 <- c(1790, 1955)

# read tree polygons file 
tree_polygons <- rgdal::readOGR(dsn = shapefile_dir,
                                layer = "polygons_checked_overlap")

# read the tree stem locations file
tree_points <- rgdal::readOGR(dsn = shapefile_dir,
                              layer = "mapped_stems_final")

# convert polygons and tree locations to SF objects
tree_polygons_sf <- sf::st_as_sf(tree_polygons)
tree_points_sf <- sf::st_as_sf(tree_points) 

# isolate the tree location coordinates 
tree_coords <- tree_points_sf %>% 
  sf::st_coordinates() %>% 
  as.data.frame()

# add new columns for the tree location coordinates 
tree_points_sf$X <- tree_coords$X
tree_points_sf$Y <- tree_coords$Y

# add empty columns for the min and max coordinates for each polygon
tree_polygons_sf$xmin <- NA 
tree_polygons_sf$xmax <- NA 
tree_polygons_sf$ymin <- NA 
tree_polygons_sf$ymax <- NA 

# add the min, max X and Y values to each polygon for filtering 
for (i in 1:nrow(tree_polygons_sf)) {
  tree_polygons_sf$xmin[i] <- as.numeric(sf::st_bbox(tree_polygons_sf$geometry[i])[1])
  tree_polygons_sf$ymin[i] <- as.numeric(sf::st_bbox(tree_polygons_sf$geometry[i])[2])
  tree_polygons_sf$xmax[i] <- as.numeric(sf::st_bbox(tree_polygons_sf$geometry[i])[3])
  tree_polygons_sf$ymax[i] <- as.numeric(sf::st_bbox(tree_polygons_sf$geometry[i])[4])
}

# merge the polygons with tree locations;
# rename the geometry columns to be more descriptive 
tree_polygons_points <- merge(as.data.frame(tree_polygons_sf),
                              as.data.frame(tree_points_sf)
                              [,c("indvdID", "X","Y","geometry")],
                              by="indvdID") %>% 
  dplyr::rename(geometry.polygon = geometry.x, 
                geometry.point = geometry.y)





# loop through h5 files ---------------------------------------------------

# get the names of all HDF5 files to iterate through
h5_list <- list.files(path = h5_dir, full.names = TRUE)
h5_list <- h5_list[grepl("*.h5", h5_list)]

# loop through h5 files 
for (h5 in h5_list) {
  
  print(h5)
  
  # list the contents of HDF5 file
  h5_struct <- rhdf5::h5ls(h5, all=T)
  
  # construct the string using "/Reflectance/Metadata/Coordinate_System",
  # without explicitly using a site code 
  crs_tag <- h5_struct$group[grepl("/Reflectance/Metadata/Coordinate_System", 
                                   h5_struct$group)][1] 
  
  # read coordinate reference system data
  crs_info <- rhdf5::h5read(h5, crs_tag)
  
  # convert "UTM" to lowercase "utm" for proper usage later
  crs_info$Proj4 <- CRS(chartr("UTM", "utm", crs_info$Proj4))
  
  # get attributes for the Reflectance dataset.
  # construct the string using "/Reflectance/Reflectance_Data"" 
  refl_tag <- paste0(h5_struct$group[grepl("/Reflectance", 
                                           h5_struct$group)][1],
                     "/Reflectance_Data")
  
  # read the reflectance metadata
  refl_info <- rhdf5::h5readAttributes(h5,refl_tag)
  
  # get the dimensions of the reflectance data
  n_rows <- refl_info$Dimensions[1]
  n_cols <- refl_info$Dimensions[2]
  n_bands <- refl_info$Dimensions[3]
  
  # print dimensions 
  print(paste0("# Rows: ", as.character(n_rows)))
  print(paste0("# Columns: ", as.character(n_cols)))
  print(paste0("# Bands: ", as.character(n_bands)))
  
  # read the wavelengths of the hyperspectral image bands
  wavelength_tag <- paste0(h5_struct$group[grepl("/Reflectance/Metadata/Spectral_Data", 
                                                 h5_struct$group)][1],
                           "/Wavelength")
  wavelengths <- rhdf5::h5read(h5,
                               wavelength_tag)
  
  # define spatial extent: extract resolution and origin coordinates
  map_info <- unlist(strsplit(crs_info$Map_Info, 
                              split = ", "))
  res_x <- as.numeric(map_info[6])
  res_y <- as.numeric(map_info[7])
  x_min <- as.numeric(map_info[4])
  y_max <- as.numeric(map_info[5])
  
  # calculate the maximum X and minimum Y values 
  x_max <- (x_min + (n_cols * res_x))
  y_min <- (y_max - (n_rows * res_y))
  tile_extent <- raster::extent(x_min, x_max, y_min, y_max)
  print("tile extent")
  print(tile_extent)
  
  # figure out which trees are within the current tile 
  polygons_in <- tree_polygons_points %>% 
    dplyr::filter(xmin >= tile_extent@xmin & 
                    xmax < tile_extent@xmax & 
                    ymin >= tile_extent@ymin & 
                    ymax < tile_extent@ymax)
  
  print(paste0(as.character(nrow(polygons_in))," polygons in current tile"))
  
  # if no polygons are within the current tile, skip to the next one
  if (nrow(polygons_in)==0){
    print("no trees located within current tile... skipping to next tile")
    next
  }
  
  # convert from SF obect to SpatialPolygons object for clipping
  polygons_in_sp <- sf::as_Spatial(polygons_in$geometry.polygon,
                                   IDs = as.character(polygons_in$indvdID))
  points_in_sp <- sf::as_Spatial(polygons_in$geometry.point,
                                 IDs = as.character(polygons_in$indvdID))
  
  
  
  # read reflectance data for all bands
  refl <- rhdf5::h5read(h5, refl_tag,
                        index = list(1:n_bands, 1:n_cols, 1:n_rows))
  
  # view and apply scale factor to convert integer values to reflectance [0,1]
  # and data ignore value
  scale_factor <- refl_info$Scale_Factor
  data_ignore <- refl_info$Data_Ignore_Value
  refl[refl == data_ignore] <- NA 
  refl_scaled <- refl / scale_factor
  
  # create georeferenced raster using band 1 
  r1 <- (refl_scaled[1,,]) # convert first band to matrix
  r1 <- raster::raster(r1, crs = crs_info$Proj4)
  extent(r1) <- tile_extent
  
  # start the raster stack with first band 
  s <- raster::stack(r1)
  
  # loop through bands and create a giant rasterstack with 426 (n_bands) bands
  for(b in 2:n_bands){
    print(b)
    
    # create raster with current band
    r <- (refl_scaled[b,,]) # convert to matrix
    r <- raster::raster(r, crs = crs_info$Proj4)
    extent(r) <- tile_extent
    
    # add additional band to the stack with the addLayer function
    s <- raster::addLayer(s, r)
    
  }
  
  # adjust the names for each layer in raster stack to correspond to wavelength
  names(s) <- round(wavelengths)
  
  # clip the hyperspectral raster stack with the polygons within current tile.
  # the returned objects are data frames, each row corresponds to a pixel in the
  # hyperspectral imagery. The ID number refers to which polygon or point that 
  # the pixel belongs to. A large polygon will lead to many extracted pixels
  # (many rows in the output data frame)
  extracted_point_spectra <- raster::extract(s, points_in_sp, df = TRUE)
  extracted_polygon_spectra <- raster::extract(s, polygons_in_sp, df = TRUE)
  
  # create polygon metadata data frame 
  polygon_metadata <- data.frame(individualID = polygons_in$indvdID,
                                 scientificName = polygons_in$scntfcN,
                                 taxonID = polygons_in$taxonID,
                                 maxCrownDiameter = polygons_in$crownDm,
                                 height = polygons_in$height,
                                 X = polygons_in$X,
                                 Y = polygons_in$Y,
                                 # create ID column to pair the polygon metadata with  spectra
                                 ID =  1:nrow(polygons_in))
  
  # create a list of increasing integer counts to keep track of how many rows 
  # (pixels or spectra) belong to each tree 
  for (i in unique(extracted_polygon_spectra$ID)){
    if(i==1){
      counts = 1:sum(extracted_polygon_spectra$ID==i)
    }
    else{
      counts = append(counts, 1:sum(extracted_polygon_spectra$ID==i))
    }
  }
  
  # combine the additional data with each spectrum for writing to file
  spectra_write <- merge(polygon_metadata,
                         as.data.frame(extracted_polygon_spectra),
                         by="ID") %>% 
    mutate(spectra_count = counts)%>% 
    select(ID, spectra_count, everything())
  
  # take a look at the first rows of the spectra data to write 
  head(spectra_write[,1:10] %>% select(-c(scientificName,X,Y)))
  
  
  # write the spectral data to file for future analysis 
  write.csv(spectra_write, file = paste0(out_dir,
                                         "spectral_reflectance_",
                                         as.character(tile_extent@xmin),"_",
                                         as.character(tile_extent@ymin),".csv")) 

}
  
  
  

  
  

# read and plot the spectra .csv files  -----------------------------------

# get a list of the .csv file per tile containing woody veg stemse
csvs <- list.files(out_dir, full.names = TRUE)
csvs <- csvs[grepl("*000.csv", csvs)]

# combine all .csv data into a single data frame 
for (i in 1:length(csvs)){
  print(csvs[i])
  csv <- read.csv(csvs[i])
  
  if(i==1){
    spectra.all <- csv
  } else {
    spectra.all <- rbind(spectra.all, csv)
  }
}



# write ALL the spectra to a single file 
write.csv(spectra.all,
          file=paste0(out.dir.spectra,site, "_spectral_reflectance_ALL.csv"))

# write the exact wavelengths to file for future use 
write.table(data.frame(wavelengths = wavelengths),
            paste(out.dir.spectra,"wavelengths.txt"),
            sep="\n",
            row.names=FALSE)

# read wavelengths if not previously created
#wavelengths = read.table(paste(out.dir.spectra,"wavelengths.txt"),
#                         sep="\n",
#                         skip = 1,
#                         col.names = 'wavelength')

# Plot spectra and color line by species ----------------------------------

