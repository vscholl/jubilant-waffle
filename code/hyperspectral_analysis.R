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
  
  # create a georeferenced rasterstack using the current hyperspectral tile
  s <- stack_hyperspectral(h5)
  
  # figure out which trees are within the current tile by 
  polygons_in <- tree_polygons_points %>% 
    dplyr::filter(xmin >= extent(s)[1] & 
                  xmax < extent(s)[2] & 
                  ymin >= extent(s)[3] & 
                  ymax < extent(s)[4])
  
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
  
  ### clip the hyperspectral raster stack with the polygons within current tile.
  # the returned objects are data frames, each row corresponds to a pixel in the
  # hyperspectral imagery. The ID number refers to which polygon or point that 
  # the pixel belongs to. A large polygon will lead to many extracted pixels
  # (many rows in the output data frame)
  
  # stem point locations 
  extracted_point_spectra <- raster::extract(s, points_in_sp, df = TRUE)
  
  # checked_overlap polygons generated using the neon_veg workflow 
  extracted_polygon_spectra <- raster::extract(s, polygons_in_sp, df = TRUE)
  
  # the buffer parameter can be used to include cells around each point of a
  # given size. the buffer parameter can be specified as a vector of the length
  # of the number of points. 
  
  # maxCrownDiameter (buffer of (maxCrownDiameter / 2))
  buffers_mxDm <- tree_polygons_points$crownDm / 2
  extracted_spectra_buffer_mxDm <- raster::extract(s, 
                                                   points_in_sp,
                                                   buffer = buffers_mxDm,
                                                   df = TRUE)
  
  # 50% max crown diameter (buffer of (maxCrownDiameter / 4))
  #buffers_50percent <- tree_polygons_points$crownDm / 4
  #extracted_spectra_buffer_50percentDm <- raster::extract(s, 
  #                                                 points_in_sp,
  #                                                 buffer = buffers_50percent,
  #                                                 df = TRUE)
  # GETTING WEIRD ERROR: 
  #Error in (function (..., deparse.level = 1)  : 
  #number of columns of matrices must match (see arg 6)
  
  # Try substituting this 50% buffer parameter instead with the polygons
  # created using neon_veg workflow with a 50% smaller diameter?????? 
  
  
  
  
  ### write spectra to file 
  
  # stem point locations 
  write_spectra_to_file(spectra = extracted_point_spectra,
                        polygons_in = polygons_in,
                        filename_out = paste0(out_dir,
                                              "spectral_reflectance_",
                                              as.character(s@extent[1]),"_",
                                              as.character(s@extent[3]),"_",
                                              "stem_points",
                                              ".csv"))
  
  # checked_overlap polygons generated using the neon_veg workflow 
  write_spectra_to_file(spectra = as.data.frame(extracted_polygon_spectra),
                        polygons_in = polygons_in,
                        filename_out = paste0(out_dir,
                                              "spectral_reflectance_",
                                              as.character(s@extent[1]),"_",
                                              as.character(s@extent[3]),"_",
                                              "polygons_checked_overlap_max_diameter",
                                              ".csv"))
  
  # maxCrownDiameter (buffer of (maxCrownDiameter / 2))
  write_spectra_to_file(spectra = as.data.frame(extracted_spectra_buffer_mxDm),
                        polygons_in = polygons_in,
                        filename_out = paste0(out_dir,
                                              "spectral_reflectance_",
                                              as.character(s@extent[1]),"_",
                                              as.character(s@extent[3]),"_",
                                              "buffer_max_diameter",
                                              ".csv"))
  
}
  



# read and plot the spectra .csv files  -----------------------------------

# get a list of the .csv file per tile containing woody veg stemse
csvs <- list.files(out_dir, full.names = TRUE)

# csvs <- csvs[grepl("*000.csv", csvs)] # from back when there was only one collection of csvs

# specify a description that the different shapefile iterations are named by
out_description <- "stem_points" # stem point locations 
out_description <- "polygons_checked_overlap_max_diameter" # checked_overlap polygons
out_description <- "buffer_max_diameter" # buffer of (maxCrownDiameter / 2)

# refine the output csv selection 
csvs <- csvs[grepl(paste0("*000_", out_description, ".csv"), csvs)]

# combine all .csv data into a single data frame 
for (i in 1:length(csvs)){
  print(csvs[i])
  csv <- read.csv(csvs[i])
  
  if(i==1){
    spectra_all <- csv
  } else {
    # add a bias value to the ID column, so in the end
    # the ID values will range from 1 to n_trees
    csv$ID <- csv$ID + max(spectra_all$ID)
    spectra_all <- rbind(spectra_all, csv)
  }
}

# remove the unneccessary column "X.1"
spectra_all <- spectra_all %>% dplyr::select(-X.1)

# write ALL the spectra to a single file 
write.csv(spectra_all,
          file=paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                      out_description,".csv"))

# write the exact wavelengths to file for future use 
write.table(data.frame(wavelengths = wavelengths),
            paste(out_dir,"wavelengths.txt"),
            sep="\n",
            row.names=FALSE)

# read wavelengths if not previously created
#wavelengths = read.table(paste(out.dir.spectra,"wavelengths.txt"),
#                         sep="\n",
#                         skip = 1,
#                         col.names = 'wavelength')









# Testing: crop HS data using single polygon ------------------------------

# read wavelengths
wavelengths = as.numeric(unlist(read.table(paste0(out_dir,"wavelengths.txt"),
                                           sep="\n",
                                           skip = 1,
                                           col.names = 'wavelength')))


p <- polygons_in_sp[1]
h5_cropped <- raster::crop(s, p)
dim(h5_cropped)

# get indices of RGB bands within the raster stack 
wl_r <- which.min(abs(wavelengths - 669)) # R 
wl_g <- which.min(abs(wavelengths - 549)) # G 
wl_b <- which.min(abs(wavelengths - 474)) # B 

h5_rgb <- raster::subset(h5_cropped, subset = c(wl_r, wl_g, wl_b)) 


plotRGB(h5_rgb,
        r = 1, g = 2, b = 3,
        stretch = "hist",
        axes = TRUE,
        main="RGB Composite",
        xlab="Easting (m)",
        ylab="Northing (m)",
        cex.main=2)

# write cropped RGB RasterBrick to a tif
writeRaster(h5_rgb,
            paste0(out_dir,"myStack.tif"), 
            format="GTiff",
            overwrite=TRUE)

