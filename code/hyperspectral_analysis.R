# load necessary R packages 
library(rhdf5)
library(rgdal)
library(raster)
library(ggplot2)
library(tidyr)
library(sf)
library(dplyr)
library(data.table) 
library(stringr) # str_split function
library(randomForest) # randomForest function
library(ggbiplot) # for PCA visualization with a biplot

# set working directory
setwd("~/github/jubilant-waffle/code/")

# load any local functions in external files 
source("supporting_functions.R")

# code for NEON site 
site_code <- 'NIWO'

# define the "bad bands" wavelength ranges in nanometers, where atmospheric 
# absorption creates unreliable reflectance values. 
bad_band_window_1 <- c(1340, 1445)
bad_band_window_2 <- c(1790, 1955)

# taxon ID's to predict
taxonList <- c("ABLAL","PICOL","PIEN","PIFL2")

# define the output directory. If it doesn't exist already, create it.
check_create_dir('../output/') # create top level "output" directory
out_dir <- paste0('../output/', site_code, '/')
check_create_dir(out_dir) # create output folder for site


# shapefile sets to test --------------------------------------------------

# shapefiles were generated using the neon_veg workflow
# where each entry (point or polygon) corresponds to a tree in the woody vegetation 
# in-situ field NEON data set.

# directory with shapefiles to test (tree stem locations and crown polygons)
shapefile_dir <- paste0('../data/', site_code, '/shapefiles/')

# define the subdirectory (destination or dsn when reading shapefile)
# and the layer name (the filename before the .shp extension) 
# in a vector for each of the shapefile scenarios. 

# all stem points at NIWO that also have crown diameter measurements. 
# this includes multi-bole entries, which result in duplicated spectra. 
allStems_layer <- c("allStems",
                    "shapefiles_maxDiameter/",
                    "mapped_stems_with_crown_diameter")

# polygons with max crown diameter, one generated for each stem point at NIWO.
# just as for the "allStems" layer above, the shapfile used here
# includes multi-bole entries, since this is how the data come when
# downloaded straight from the data portal. 
allPolygons_maxDiameter_layer <- c("allPolygons_maxDiameter",
                                   "shapefiles_maxDiameter/",
                                   "polygons_all")

# polygons with half max diameter, one generated for each stem point at NIWO. 
allPolygons_halfDiameter_layer <- c("allPolygons_halfDiameter",
                                    "shapefiles_halfDiameter/",
                                    "polygons_all")

# neon_veg workflow polygons generated with max diameter 
neonvegPolygons_maxDiameter_layer <- c("neonvegPolygons_maxDiameter",
                                       "shapefiles_maxDiameter/",
                                       "polygons_clipped_overlap")

# neon_veg workflow polygons generated with half the max diameter  
neonvegPolygons_halfDiameter_layer <- c("neonvegPolygons_halfDiameter",
                                        "shapefiles_halfDiameter/",
                                        "polygons_clipped_overlap")

# neon_veg workflow stems, corresponding to polygons generated with max diameter
neonvegStems_maxDiameter_layer <- c("neonvegStems_maxDiameter",
                                    "shapefiles_maxDiameter/",
                                    "mapped_stems_final")

# create a data frame containing the description, directory, and shapefile name
# for each of the shapefile scenarios to be tested 
shapefileLayerNames <- as.data.frame(rbind(allStems_layer,
                                           allPolygons_halfDiameter_layer,
                                           allPolygons_maxDiameter_layer,
                                           neonvegStems_maxDiameter_layer,
                                           neonvegPolygons_halfDiameter_layer,
                                           neonvegPolygons_maxDiameter_layer),
                                     stringsAsFactors = FALSE) %>% 
          `colnames<-`(c("description", "dsn", "layer")) 
rownames(shapefileLayerNames) <- 1:nrow(shapefileLayerNames)

print("shapefile layers to test:")
shapefileLayerNames



# input AOP airborne remote sensing data --------------------------------------

# assemble all of the NEON AOP remote sensing data layers into a giant data cube.
# a set of spectral and structural features from this data cube will be extracted
# for each sample (pixel) that intersects with the shapefile entries.  

# specify the paths for each data directory
# hyperspectral .h5 and .rds
h5_dir <- paste0('../data/', site_code, '/hyperspectral/') 
# CHM geotiffs 
chm_dir <- paste0('../data/', site_code, '/chm/')      
# slope geotiffs
slope_dir <- paste0('../data/', site_code, '/slope/')   
# aspect geotiffs
aspect_dir <- paste0('../data/', site_code, '/aspect/') 
# rgb image geotiffs
rgb_dir <- paste0('../data/', site_code, '/rgb/')
# vegetation index .tifs 
vegIndices_dir <- paste0('../data/', site_code, '/vegIndices/') 

# list the files in each data directory; filter results based on file type
# hyperspectral data - list the .h5 files 
h5_list <- list.files(path = h5_dir, full.names = TRUE)
h5_list <- h5_list[grepl("*.h5", h5_list)] 
# canopy height model data - list the .tif files 
chm_list <- list.files(path = chm_dir, full.names = TRUE)
chm_list <- chm_list[grepl("*CHM.tif$", chm_list)]
# lidar-derived slope data - list the .tif files 
slope_list <- list.files(path = slope_dir, full.names = TRUE)
slope_list <- slope_list[grepl("*slope.tif$", slope_list)]
# lidar-derived aspect data - list the .tif files
aspect_list <- list.files(path = aspect_dir, full.names = TRUE)
aspect_list <- aspect_list[grepl("*aspect.tif$", aspect_list)]
# digital camera rgb data - list the .tif files 
rgb_list <- list.files(path = rgb_dir, full.names = TRUE)
rgb_list<- rgb_list[grepl("*image.tif$", rgb_list)]
# hyperspectral-derived veg indices - list the subdirectories for each index 
vegIndices_list <- list.dirs(path = vegIndices_dir, full.names = TRUE)
vegIndices_names <- c("ARVI","EVI","NDLI","NDNI","NDVI","PRI","SAVI")


# read wavelengths from text file 
# created in supporting_functions/stack_hyperspectral
wavelengths = as.numeric(unlist(read.table(paste0(out_dir,"wavelengths.txt"),
                                           sep="\n",
                                           skip = 1,
                                           col.names = 'wavelength')))



# create data cubes with AOP-derived features ---------------------------

start_time <- Sys.time() # start the timer 

# loop through the tiles; build up a data cube for each tile 

for (h5 in h5_list) {
  
  # print current hyperspectral filename 
  print(h5)
  
  # each current hyperspectral tile must be read and stacked into a 
  # georeferenced rasterstack object (so it can be clipped with point / polygon
  # shapefiles). The process of creating a rasterstack takes a while for 
  # each tile, so after creating each rasterstack once, each object gets 
  # written to a file. 
  
  # Build up the rasterstack filename by parsing out the easting/northing
  # coordinates from the current h5 filename.
  
  # parse the UTM easting and northing values from the current h5 filename
  easting <- stringr::str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][5]
  northing <- stringr::str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][6]
  # combine them with an underscore; use this to find corresponding tiles 
  # of various remote sensing data
  east_north_string <- paste0(easting,"_",northing)
  
  
  # hyperspectral and lidar features ----------------------------------------
  
  # Build up the h5 rasterstack filename
  rasterstack_filename <- paste0(h5_dir, "rasterstack_",
                                 east_north_string, ".rds")
  
  print(paste("rasterstack filename: ", rasterstack_filename))
  
  # check to see if a .rds file already exists for the current tile
  # within the hyperspectral directory - just HS reflectance bands so far 
  if (file.exists(rasterstack_filename)){
    
    # if it exists, read that instead of re-generating the same rasterstack.
    message("reading rasterstack (already created for current tile)...")
    # restore / read the rasterstack from file
    s <- readRDS(file = rasterstack_filename)
    
  } else{
    
    # if it doesn't exist, generate the rasterstack. 
    message("creating rasterstack for current tile...")
    # create a georeferenced rasterstack using the current hyperspectral tile
    s <- stack_hyperspectral(h5)
    # save the rasterstack to file 
    saveRDS(s, file = rasterstack_filename)
  }
  
  # read the corresponding remote sensing data layers for current tile
  print("Reading lidar-derived layers for the training data cube...")
  chm <- raster::raster(grep(east_north_string, chm_list, value=TRUE))
  slope <- raster::raster(grep(east_north_string, slope_list, value=TRUE))
  aspect <- raster::raster(grep(east_north_string, aspect_list, value=TRUE))
  
  # for the vegetation indices, go into the corresponding folder for current tile
  # and get a list of all the vegIndex geotiffs. then read all of those geotiffs 
  # into a single raster stack.
  print("Reading vegetation indices for the training data cube...")
  vegIndices <- raster::stack(list.files(grep(east_north_string, 
                                              vegIndices_list, 
                                              value=TRUE), 
                                         pattern="tif$", full.names=TRUE))
  
  # set the raster name for each layer to be simply the name of the data 
  # (i.e. "aspect") as opposed to the full filename 
  # (i.e. ""NEON_D13_NIWO_DP3_452000_4431000_aspect")
  names(chm) <- "chm"
  names(slope) <- "slope"
  names(aspect) <- "aspect"
  # name each of the vegetation index layers based on the last piece of each 
  # respective filename, e.g. "NDVI" and 
  names(vegIndices) <- sapply(stringr::str_split(names(vegIndices),"_"),tail,1)
  
  # RGB features ------------------------------------------------------------
  
  # The RGB data tile has 10,000 x 10,000 pixels, 10cm spatial resolution. 
  # All other layers tiles have 1,000 x 1,000 pixels, 1 meter spatial resolution. 
  # Aggregate red, green, blue intensity within each coarser grid cell using
  # statistics such as mean and standard deviation. 
  
  rgb_features_filename <- paste0(rgb_dir, "rgb_features_",
                                  east_north_string, ".rds")
  
  print(paste("rgb_features: ", rgb_features_filename))
  
  # check if a .rds file already exists for the current feature data cube
  if (file.exists(rgb_features_filename)){
    
    # if it exists, read that instead of re-generating the same rasterstack.
    message("reading rgb_features (already created for current tile)...")
    
    # restore / read the rasterstack from file
    rgb_features <- readRDS(file = rgb_features_filename)
    
  } else{
    
    # if it doesn't exist, create the features from the aop data to file 
    message("creating rgb_features for current tile...")
    
    # rgb data has 3 bands. read each one individually 
    rgb_red <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), band = 1) 
    rgb_green <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), band = 2) 
    rgb_blue <- raster::raster(grep(east_north_string, rgb_list, value = TRUE), band = 3) 
    
    # The "fact" parameter of the raster::aggregate function is the number of cells
    # in each direction (horizontal and vertically) to aggregate across.
    # Since the RGB data has a spatial resolution that is 1/10th of the 
    # other data layers (10cm compared to 1m), fact should be 10 to produce 
    # an output raster with 1000 x 1000 pixels. 
    
    # mean intensity per 1m x 1m grid cell 
    rgb_meanR <- raster::aggregate(rgb_red, fact = 10, fun = mean)
    rgb_meanG <- raster::aggregate(rgb_green, fact = 10, fun = mean)
    rgb_meanB <- raster::aggregate(rgb_blue, fact = 10, fun = mean)
    
    # standard deviation of intensity per 1m x 1m grid cell 
    rgb_sdR <- raster::aggregate(rgb_red, fact = 10, fun = sd)
    rgb_sdG <- raster::aggregate(rgb_green, fact = 10, fun = sd)
    rgb_sdB <- raster::aggregate(rgb_blue, fact = 10, fun = sd)
    
    # (mean + SD) gives some idea about relative position on the spectrum
    # and variation rather than just variation (standard deviation on its own)
    rgb_mean_sd_R <- (rgb_meanR + rgb_sdR)
    rgb_mean_sd_G <- (rgb_meanG + rgb_sdG)
    rgb_mean_sd_B <- (rgb_meanB + rgb_sdB)
    
    
    # future work: would be awesome to calculate texture features: 
    # https://cran.r-project.org/web/packages/radiomics/vignettes/TextureAnalysis.html
    
    
    # to confirm the order of the red, green, and blue intensities,
    # stack all 3 bands into a RasterStack and use the plotRGB function.
    # the colors appear natural, so r=1, g=2, b=3
    #rgb_stack <- raster::stack(rgb_red,rgb_green, rgb_blue)
    #raster::plotRGB(rgb_stack, r = 1, g = 2, b = 3)
    
    # set the names of each layer to reflect the metric it contains
    names(rgb_meanR) <- "rgb_meanR" # mean intensity 
    names(rgb_meanG) <- "rgb_meanG"
    names(rgb_meanB) <- "rgb_meanB"
    names(rgb_sdR) <- "rgb_sdR"     # standard deviation of intensity 
    names(rgb_sdG) <- "rgb_sdG"
    names(rgb_sdB) <- "rgb_sdB"
    names(rgb_mean_sd_R) <- "rgb_mean_sd_R" # mean plus standard deviation
    names(rgb_mean_sd_G) <- "rgb_mean_sd_G"
    names(rgb_mean_sd_B) <- "rgb_mean_sd_B"
    
    # stack up all the RGB features
    rgb_features <- raster::stack(rgb_meanR, rgb_meanG, rgb_meanB,
                                  rgb_sdR, rgb_sdG, rgb_sdB,
                                  rgb_mean_sd_R, rgb_mean_sd_G, rgb_mean_sd_B)
    
    # save the stacked features from the aop data to file 
    saveRDS(rgb_features, file = rgb_features_filename)
  }
  
  # Create the pixel number grid as a layer to add to the data cube. 
  # this one keeps track of individual pixel ID's
  # to avoid duplicate spectra being extracted. Basically, assign an integer ID
  # to each pixel in the 1000x1000 raster. This raster needs to have the same 
  # dimensions, extent, crs as the other layers so they can be stacked together. 
  # create a vector of IDs from 1 to the number of pixels in one band (#rows x #cols)
  pixelID <- 1:(nrow(s) * ncol(s))
  # add tile east, north coordinates - how to do this if raster values must be numeric? 
  #pixelID <- paste(pixelID, east_north_string, sep="_") 
  # reshape this 1D vector into a 2D matrix 
  dim(pixelID) <- c(nrow(s),ncol(s))
  # create a raster layer of pixel numbers 
  pixelNumbers <- raster::raster(pixelID, crs = crs(s))
  extent(pixelNumbers) <- extent(s)
  names(pixelNumbers) <- "pixelNumber"
  
  
  # now, all of the hyperspectral data files have been read in for the current
  # tile. add each one to the hyperspectral data stack along with the 
  # layer to keep track pixel number within the tile. 
  stacked_aop_data <- raster::addLayer(s, chm, slope, aspect, vegIndices, 
                                       rgb_features, pixelNumbers)
  print("Stacked AOP data for current tile. ")
  
  
  # clip data cube - extract features ----------------------------------------
  
  # loop through shapefile sets 
  for(i in 1:nrow(shapefileLayerNames)){ 
    print(paste0("Currently extracting features for tree points / polygons in:  ", 
                 shapefileLayerNames$description[i]))
    
    # read the shapefile layer 
    shp <- rgdal::readOGR(dsn = paste0(shapefile_dir,shapefileLayerNames$dsn[i]),
                          layer = shapefileLayerNames$layer[i])
    
    # convert to SF object
    shp_sf <- sf::st_as_sf(shp)
    
    # add columns for the center location of each tree 
    shp_coords <- shp_sf %>% 
      sf::st_centroid() %>%  # get the centroids first for polygon geometries 
      sf::st_coordinates() %>% 
      as.data.frame()
    
    # add new columns for the tree location coordinates 
    shp_sf$X <- shp_coords$X
    shp_sf$Y <- shp_coords$Y
    
    # figure out which trees are within the current tile by comparing each
    # X,Y coordinate to the extent of the current tile 
    trees_in <- shp_sf %>% 
      dplyr::filter(X >= extent(s)[1] & 
                      X < extent(s)[2] & 
                      Y >= extent(s)[3] & 
                      Y  < extent(s)[4])
    
    print(paste0(as.character(nrow(trees_in))," trees in current tile"))
    
    # if no polygons are within the current tile, skip to the next one
    if (nrow(trees_in)==0){
      print("no trees located within current tile... skipping to next shapefile")
      next
    }
    
    # convert from SF obect to Spatial object for clipping
    trees_in_sp <- sf::as_Spatial(trees_in,
                                  IDs = as.character(trees_in$indvdID))
    
    # clip the hyperspectral raster stack with the polygons within current tile.
    # the returned objects are data frames, each row corresponds to a pixel in the
    # hyperspectral imagery. The ID number refers to which tree that the 
    # the pixel belongs to. A large polygon will lead to many extracted pixels
    # (many rows in the output data frame), whereas tree stem points will
    # lead to a single extracted pixel per tree. 
    print("Extracting features for each tree from the data cube... ")
    extracted_spectra <- raster::extract(stacked_aop_data, 
                                         trees_in_sp, 
                                         df = TRUE)
    
    # TO DO: 
    # adjust this extract step to only get pixels WITHIN each tree polygon,
    # also try calculating the percentage that each pixel is within a polygon
    # and keep only pixels with > 50% overlap 
    
    # merge the extracted spectra and other data values with the tree info 
    tree_metadata <- data.frame(trees_in) %>% 
      mutate(ID = 1:nrow(trees_in))
    
    # create a list of increasing integer counts to keep track of how many rows 
    # (pixels or spectra) belong to each tree 
    for (j in unique(extracted_spectra$ID)){
      if(j==1){
        counts = 1:sum(extracted_spectra$ID==j)
      }
      else{
        counts = append(counts, 1:sum(extracted_spectra$ID==j))
      }
    }
    
    # combine the additional data with each spectrum for writing to file.
    # remove the geometry column to avoid issues when writing to csv later 
    spectra_write <- merge(tree_metadata,
                           extracted_spectra,
                           by="ID") %>% 
      mutate(spectra_count = counts)%>% 
      select(ID, spectra_count, everything()) %>% 
      select(-geometry)
    
    # write extracted spectra and other remote sensing data values to file 
    write.csv(spectra_write, 
              file = paste0(out_dir,
                            "extracted_features_",
                            east_north_string, "_",
                            shapefileLayerNames$description[i], ".csv")) 
  }
  
}

end_time <- Sys.time()
elapsed <- end_time - start_time
print("Elapsed time: ")
print(elapsed)



# combine all extracted features into a single 
# .csv for each of the shapefile sets.
for(i in 1:nrow(shapefileLayerNames)){
  
  csvs <- list.files(out_dir, full.names = TRUE)
  
  # refine the output csv selection 
  csvs <- csvs[grepl(paste0("*000_", shapefileLayerNames$description[i], ".csv"), csvs)]
  
  # combine all .csv data into a single data frame 
  for (c in 1:length(csvs)){
    print(csvs[c])
    csv <- read.csv(csvs[c])
    
    if(c==1){
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
  
  # write ALL the spectra to a single .csv file 
  write.csv(spectra_all,
            file=paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                        shapefileLayerNames$description[i],".csv"))
  
}



# ribbon plots ------------------------------------------------------------

# read wavelengths if not previously created
wavelengths = as.numeric(unlist(read.table(paste0(out_dir,"wavelengths.txt"),
                                           sep="\n",
                                           skip = 1,
                                           col.names = 'wavelength')))

# loop through each shapefile name, read the .csv containing spectral reflectance 
# for all trees within each data set, generate a ribbon plot and write to image file
# in the figures/ output directory 
for(i in 1:nrow(shapefileLayerNames)){
  extracted_features_filename <- paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                                        shapefileLayerNames$description[i],".csv")
  
  createRibbonPlot(wavelengths, extracted_features_filename)
 
}



# calculate JM distance - spectral separability ---------------------------

for(i in 1:nrow(shapefileLayerNames)){
  extracted_features_filename <- paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                                        shapefileLayerNames$description[i],".csv")
  
  distances <- jmDist(extracted_features_filename)
  print(distances)
  
}



# spectral data with reflectance per wavelength (column), with each
# row being a different pixel
x <- features %>% dplyr::select(c(wavelength_lut$xwavelength)) %>% as.matrix()

# class labels for each pixel. taxonIDs are originally factors.
# reclassify them as numbers. 
classes <- features$taxonID
classNumbers <- matrix(nrow = length(classes))
classNumbers[classes == "ABLAL"] <- 1
classNumbers[classes == "PICOL"] <- 2
classNumbers[classes == "PIEN"] <- 3
classNumbers[classes == "PIFL2"] <- 4
classes <- classNumbers

# spectral separability ---------------------------------------------------

# https://www.rdocumentation.org/packages/fpc/versions/2.1-11.1/topics/bhattacharyya.dist 
#https://www.rdocumentation.org/packages/varSel/versions/0.1/topics/JMdist
library(varSel)

# create a data frame where each row is a reflectance spectrum and each 
# coloumn is a wavelength
X <- features %>% dplyr::select(c(wavelength_lut$xwavelength))
X <- X[,1:100]
  
# a column vector of the lables. length(g) is equal to nrow(X).
classes <- features$taxonID
classNumbers <- matrix(nrow = length(classes))
classNumbers[classes == "ABLAL"] <- 1
classNumbers[classes == "PICOL"] <- 2
classNumbers[classes == "PIEN"] <- 3
classNumbers[classes == "PIFL2"] <- 4
g <- classNumbers

sep <- varSel::JMdist(g,X)


# Jeffries-Matusita distance ----------------------------------------------

# https://stats.stackexchange.com/questions/106325/jeffries-matusita-distance-for-14-variables

# Compute the Mahalanobis distance between two vectors.
mahalanobis <- function(m1, m2, sigma) {m <- m1 - m2; m %*% solve(sigma, m)}

# Compute the Bhattacharyya distance between two multivariate normal distributions
# given by their means and covariance matrices.
bhattacharyya <- function(m1, s1, m2, s2) {
  d <- function(u) determinant(u, logarithm=TRUE)$modulus # Log determinant of matrix u
  s <- (s1 + s2)/2                                        # mean covariance matrix
  mahalanobis(m1, m2, s)/8 + (d(s) - d(s1)/2 - d(s2)/2)/2
}

# Re-express the Bhattacharyya distance as the Jeffries-Matusita distance.
# Values range from 0 (poor separability) to 2 (great separability)
jeffries.matusita <- function(...) 2*(1-exp(-bhattacharyya(...)))

#------------------------------------------------------------------------------------#
# Generate sets of sample data for d bands aggregated into classes.
d <- 14                          # Number of bands
n <- rep(1000, 5)                # Class pixel counts
n.class <- length(n)             # Number of classes
require(MASS)                    # For generating multivariate normals
set.seed(17)                     # Allows reproducible results

# Create random mean vectors for the classes
mu <- round(matrix(rnorm(d*n.class, 128, 1), ncol=n.class, byrow=TRUE), 0)
# Initialize the data structure {x, classes}
x <- matrix(double(), ncol=d, nrow=0)
classes <- integer()
# Generate the random data
for (i in 1:n.class) {
  # Create a random valid covariance matrix for this class
  f <- svd(matrix(rnorm(d^2), ncol=d))
  sigma <- t(f$v) %*% diag(rep(10, d)) %*% f$v
  # Generate d-variate normals
  x <- rbind(x, mvrnorm(n[i], mu[, i], sigma))
  classes <- c(classes, rep(i, n[i]))
}



# Given a set of bands (as the columns of x) and a grouping variable in `class`,
# compute the class means and covariance matrices (the "signatures").
classes.stats <- by(x, classes, function(y) list(mean=apply(y, 2, mean), cov=cov(y)))

# Compute the J-M distances between the classes.
distances <- matrix(0.0, n.class, n.class)
for (i in 2:n.class) {
  m1 <- classes.stats[[i]]$mean; s1 <- classes.stats[[i]]$cov
  for (j in 1:(i-1)) {
    m2 <- classes.stats[[j]]$mean; s2 <- classes.stats[[j]]$cov
    distances[i,j] <- distances[j,i] <- jeffries.matusita(m1,s1,m2,s2)
  }
}
print(distances)
#------------------------------------------------------------------------------------#



# Jeffries-Matusita distance using NIWO spectra ---------------------------

# spectral data with reflectance per wavelength (column), with each
# row being a different pixel
x <- features %>% dplyr::select(c(wavelength_lut$xwavelength)) %>% as.matrix()

# class labels for each pixel. taxonIDs are originally factors.
# reclassify them as numbers. 
classes <- features$taxonID
classNumbers <- matrix(nrow = length(classes))
classNumbers[classes == "ABLAL"] <- 1
classNumbers[classes == "PICOL"] <- 2
classNumbers[classes == "PIEN"] <- 3
classNumbers[classes == "PIFL2"] <- 4
classes <- classNumbers

# number of unique classes 
n.class <- length(unique(classes))

# calculate the mean and covariance for all spectra within each class 
classes.stats <- by(x, classes, function(y) list(mean=apply(y, 2, mean), cov=cov(y)))

# Compute the J-M distances between the classes.
distances <- matrix(0.0, n.class, n.class)
for (i in 2:n.class) {
  m1 <- classes.stats[[i]]$mean; s1 <- classes.stats[[i]]$cov
  for (j in 1:(i-1)) {
    m2 <- classes.stats[[j]]$mean; s2 <- classes.stats[[j]]$cov
    distances[i,j] <- distances[j,i] <- jeffries.matusita(m1,s1,m2,s2)
  }
}
print(distances)










# Interspecies comparison Boxplots of variables per species -------------------

gg_alpha <- 1

# ASPECT 
ggplot(data = features, aes(x = taxonID, y = aspect, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: Aspect")

ggsave(paste0(out_dir, outDescription, "boxplot_ASPECT_",shapefileLayerNames$description[i],".png"))

# SLOPE 
ggplot(data = features, aes(x = taxonID, y = slope, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: Slope")

ggsave(paste0(out_dir, outDescription, "boxplot_SLOPE_",shapefileLayerNames$description[i],".png"))


# CHM 
ggplot(data = features, aes(x = taxonID, y = chm, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: CHM-derived height")

ggsave(paste0(out_dir, outDescription, "boxplot_CHM_",shapefileLayerNames$description[i],".png"))


# NDVI 
ggplot(data = features, aes(x = taxonID, y = NDVI, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: NDVI")

ggsave(paste0(out_dir, outDescription, "boxplot_NDVI_",shapefileLayerNames$description[i],".png"))


# ARVI
ggplot(data = features, aes(x = taxonID, y = ARVI, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: ARVI")

ggsave(paste0(out_dir, outDescription, "boxplot_ARVI_",shapefileLayerNames$description[i],".png"))


# PRI 
ggplot(data = features, aes(x = taxonID, y = PRI, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: PRI")

ggsave(paste0(out_dir, outDescription, "boxplot_PRI_",shapefileLayerNames$description[i],".png"))


# NDLI
ggplot(data = features, aes(x = taxonID, y = NDLI, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: NDLI")

ggsave(paste0(out_dir, outDescription, "boxplot_NDLI_",shapefileLayerNames$description[i],".png"))


# NDNI
ggplot(data = features, aes(x = taxonID, y = NDNI, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: NDNI")

ggsave(paste0(out_dir, outDescription, "boxplot_NDNI_",shapefileLayerNames$description[i],".png"))


# rgb_mean_sd_B
ggplot(data = features, aes(x = taxonID, y = rgb_mean_sd_B, colour = taxonID)) +
  geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5) + 
  scale_color_brewer(palette="Spectral") + 
  ggtitle("Interspecies boxplot comparison: rgb_mean_sd_B")

ggsave(paste0(out_dir, outDescription, "boxplot_rgb_mean_sd_B_",shapefileLayerNames$description[i],".png"))








# count number of polygons and pixels per shapefile scenario ------------------
for(i in 1:nrow(shapefileLayerNames)){
  print(i)
  extracted_features_filename <- paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                                        shapefileLayerNames$description[i],".csv")
  print("Currently counting the number of individual IDs and extracted spectra in:")
  print(extracted_features_filename)
  
  # read the values extracted from the data cube
  df_orig <- read.csv(extracted_features_filename)
  
  if(i ==1){
    count <- c(length(unique(df_orig$indvdID)), length(unique(df_orig$pixelNumber)))
  } else{
    count <- rbind(count, c(length(unique(df_orig$indvdID)), length(unique(df_orig$pixelNumber))))
  }
}
countDF <- data.frame(count,
                   row.names = shapefileLayerNames$description)
colnames(countDF) <- c("nIndvdID", "nPixelNumbers")



# Random Forest Classification --------------------------------------------

# At this point, all of the shapefile scenarios have been used to extract
# features from the giant remote sensing data cube.

# specific string to name a directory to hold classification output files
#outDescription <- "rf_allSamplesPerClass_ntree500_validationSet/"
outDescription <- "rf_allSamplesPerClass_ntree500/"
outDescription <- "rf_allSamplesPerClass_ntree1000/"
outDescription <- "rf_allSamplesPerClass_ntree500_pcaInsteadOfWavelengths/"
outDescription <- "rf_allSamplesPerClass_ntree500_pca4InsteadOfWavelengths_keep10MostImpVar/"
outDescription <- "rf_allSamplesPerClass_ntree2000_pca2InsteadOfWavelengths_keep10MostImpVar/"
outDescription <- "rf_allSamplesPerClass_ntree2000_pca2InsteadOfWavelengths/"
outDescription <- "rf_allSamplesPerClass_ntree5000_pca2InsteadOfWavelengths/" 
outDescription <- "rf_allSamplesPerClass_ntree500_pca2InsteadOfWavelengths/" 


check_create_dir(paste0(out_dir,outDescription))

# RF tuning parameter, number of trees to grow. deafualt value 500
ntree <- 500
#ntree <- 2000 # using PCA the run time is significantly cut down 
#ntree <- 5000

# boolean variable. if TRUE, select random minSamples per class to reduce bias
randomMinSamples <- FALSE

# boolean variable. if TRUE, keep separate set for validation
independentValidationSet <- FALSE 
# randomly select this amount of data for training, use the rest for validation
percentTrain <- 0.8 

pcaInsteadOfWavelengths <- TRUE
nPCs <- 2 # number of PCAs to keep 

# keep most important variables and run RF again with reduced feature set 
keepMostImpVar <- FALSE

# create boxplots and write to file 
createBoxplots <- FALSE

# open a text file to record the output results 
rf_output_file <- file(paste0(out_dir,outDescription,
                                  "rf_model_summaries.txt"), "w")

# create an empty list to summarise the model performances
# --> how to access the OOD error rate from the randomForest output object???!?
#modelComparison <- data.frame(matrix(ncol = 3, nrow = nrow(shapefileLayerNames)))
#x <- c("OOB-error", "OveralAccuracy", "rankOA")
#colnames(modelComparison) <- x
rfAccuracies <- data.frame(matrix(ncol = 5, nrow = nrow(shapefileLayerNames)))
x <- c("shapefileDescription", "OA", "UA","PA","K")
colnames(rfAccuracies) <- x
rfAccuracies$shapefileDescription <- shapefileLayerNames$description


# start the timer
start_time <- Sys.time()

# loop through all shapefile sets 
for(i in 1:nrow(shapefileLayerNames)){ 
  
  # filename of current .csv file to read 
  extracted_features_filename <- paste0(out_dir, site_code, 
                                        "_spectral_reflectance_ALL_",
                                        shapefileLayerNames$description[i],
                                        ".csv")
  print("Currently training the random forest with features extracted using:")
  print(extracted_features_filename)
  
  # read the values extracted from the data cube
  df_orig <- read.csv(extracted_features_filename)
  
  # Remove any spectra that have a height == 0
  print(paste0(as.character(sum(df_orig$chm==0)), 
               " pixels have a height of 0 in the CHM"))
  print("Removing these rows from the training set ... ")
  
  # also reset the factor levels (in case there are dropped taxonID levels)
  df <- df_orig %>% filter(chm>0) %>% droplevels()
  
  # remove the bad bands from the list of wavelengths 
  remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                                 wavelengths < bad_band_window_1[2]) | 
                                (wavelengths > bad_band_window_2[1] & 
                                   wavelengths < bad_band_window_2[2])]
  
  # create a LUT that matches actual wavelength values with the column names,
  # X followed by the rounded wavelength values. Remove the rows that are 
  # within thebad band ranges. 
  wavelength_lut <- data.frame(wavelength = wavelengths,
                               xwavelength = paste0("X", 
                                                    as.character(round(wavelengths))),
                               stringsAsFactors = FALSE) %>% 
    filter(!wavelength %in% remove_bands)
  
  # features to use in the RF models.
  # this list is used to filter the columns of the data frame,
  # to remove the ones containing other metadata per tree from the model. 
  featureNames <- c("taxonID", 
                    wavelength_lut$xwavelength,
                    "chm", 
                    "slope", 
                    "aspect",
                    "ARVI",
                    "EVI",
                    "NDLI",
                    "NDNI",
                    "NDVI",
                    "PRI",
                    "SAVI",
                    "rgb_meanR",
                    "rgb_meanG",
                    "rgb_meanB",
                    "rgb_mean_sd_R",
                    "rgb_mean_sd_G",
                    "rgb_mean_sd_B")
  
  # filter the data to contain only the features of interest 
  features <- df %>% 
    dplyr::select(featureNames)
  
  # testing whether PCA yields better accuracy than individual wavelength reflectance data
  if(pcaInsteadOfWavelengths == TRUE){
    
    # remove the individual spectral reflectance bands from the training data
    features <- features %>% dplyr::select(-c(wavelength_lut$xwavelength))
    
    # PCA: calculate Principal Components 
    hs <- df %>% dplyr::select(c(wavelength_lut$xwavelength)) %>% as.matrix()
    hs_pca <- stats::prcomp(hs, center = TRUE, scale. = TRUE)
    summary(hs_pca)
    features <- cbind(features, hs_pca$x[,1:nPCs]) # add first n PCs to features data frame
    # visualize where each sample falls on a plot with PC2 vs PC1 
    ggbiplot::ggbiplot(hs_pca, 
                       obs.scale = 1, var.scale = 1, # scale the observations and variables 
                       #var.axes=FALSE, # remove arrows 
                       groups = df$taxonID, # color the points by species
                       ellipse = TRUE, # draw ellipse around each group
                       circle = TRUE # draw circle around center of data set
                        )   + 
      ggtitle("PCA biplot, PC1 and PC2") + 
      scale_color_brewer(palette="Spectral")
    # save to file 
    ggsave(paste0(out_dir, outDescription, "pcaPlot_",shapefileLayerNames$description[i],".png"))
    
  }
  
  print("Features used in current RF model: ")
  print(colnames(features))
  
  
  # count the number of samples per species 
  featureSummary <- features %>%
    dplyr::group_by(as.character(taxonID)) %>%
    dplyr::summarize(total = n()) 
  
  print("number of samples per species class")
  print(featureSummary)
  
  # randomly select <percentTrain> of data for training,
  # use the remaining samples for validation
  if(independentValidationSet){
    train <- sample(nrow(features), percentTrain*nrow(features), replace = FALSE)
    trainSet <- features[train,]
    validationSet <- features[-train,]
    summary(trainSet)
    summary(validationSet)
    features <- trainSet
  }
  
  if(randomMinSamples){
    # reduce number of samples per species to avoid classifier bias
    
    # count the minimum number of samples for a single class
    minSamples <- min(featureSummary$total)  
    print(paste0("Randomly selecting ",
                 as.character(minSamples),
                 " samples per species class to avoid classifier bias"))
    
    # isolate the samples per species
    taxon1 <- features[features$taxonID==taxonList[1],]
    taxon2 <- features[features$taxonID==taxonList[2],]
    taxon3 <- features[features$taxonID==taxonList[3],]
    taxon4 <- features[features$taxonID==taxonList[4],]
    
    # keep random minSamples of each species; merge
    taxon1 <- taxon1[sample(nrow(taxon1), minSamples), ]
    taxon2 <- taxon2[sample(nrow(taxon2), minSamples), ]
    taxon3 <- taxon3[sample(nrow(taxon3), minSamples), ]
    taxon4 <- taxon4[sample(nrow(taxon4), minSamples), ]
    
    features <- rbind(taxon1, taxon2, taxon3, taxon4)
    
  } else{
    print("Using all samples per class")
  }
  
  
  # train the RF model using the training set
  rf_startTime <- Sys.time()
  set.seed(104)
  rf_model <- randomForest::randomForest(as.factor(features$taxonID) ~ .,
                                         data=features, 
                                         importance=TRUE, 
                                         ntree=ntree) # ntree is number of trees to grow
  print("randomForest time elapsed for model training: ")
  print(Sys.time()-rf_startTime)
  
  print(rf_model)
  
  
  # use the rfUtilities package to calculate OA, PA, UA accuracies 
  library(rfUtilities)
  accuracy <- rfUtilities::accuracy(rf_model$predicted,rf_model$y)
  
  # record each accuracy metric in the table for a final comparison.
  # round each value to the nearest decimal place 
  rfAccuracies$OA[i] <- round(accuracy$PCC, 1) # Overall Accuracy
  rfAccuracies$UA[i] <-  round(accuracy$users.accuracy, 1) #User's Accuracy
  rfAccuracies$PA[i] <- round(accuracy$producers.accuracy, 1) #Producer's Accuracy
  rfAccuracies$K[i] <- round(accuracy$kappa, 3) #Cohen's Kappa 
  
  # Parallel randomForest
  #library(foreach)
  #rf_startTime <- Sys.time()
  #test <- foreach(ntree=rep(500, 6), .combine=randomForest::combine,
  #              .multicombine=TRUE, .packages='randomForest') %dopar% {
  #                randomForest(as.factor(features$taxonID) ~ .,
  #                data=features, importance=TRUE, ntree=ntree)
  #              }
  #print("randomForest time elapsed: ")
  #print(Sys.time()-rf_startTime)
  
  
  # plot error as a function of ntrees
  #plot(rf_model)
  #legend("topright", colnames(rf_model$err.rate),col=1:5,cex=0.8,fill=1:5)
  
  
  # What variables were important? --> Consult the variable importance plot. 
  
  # save varImpPlot to image file. create a filename and png object.  
  varImpFilename <- paste0(out_dir, outDescription,"varImpPlot_",shapefileLayerNames$description[i],".png")
  png(filename = varImpFilename)
  # make varImpPlot
  randomForest::varImpPlot(rf_model,
                           main = shapefileLayerNames$description[i])
  dev.off() # saves the plot to the image filename  
  
  # save RF model to file 
  save(rf_model, file = paste0(out_dir, outDescription,"rf_model_",
                               shapefileLayerNames$description[i],".RData"))
  
  # use the model to predict the species for the training set 
  predTrain <- predict(rf_model, features, type = "class")
  
  # Checking classification accuracy
  confusionTable <- table(predTrain, features$taxonID)
  print(confusionTable)
  
  # print overall accuracy
  # TO DO: modify this for the independent validation set
  print(paste0("overall accuracy predicting train set: ",
               as.character(mean(predTrain == features$taxonID))))
  
  # predict species ID for validation set 
  if(independentValidationSet){
    predValidation <- predict(rf_model, validationSet, type = "class")
    confusionTable <- table(predValidation, validationSet$taxonID)
    print(confusionTable)
    print(paste0("overall accuracy predicting validation set: ",
                 as.character(mean(predValidation == validationSet$taxonID))))
  }
  
  
  # write all relevant information to the textfile: 
  # shapefile name
  write(shapefileLayerNames$description[i], rf_output_file, append=TRUE)
  write("\n", rf_output_file, append=TRUE) #newline
  
  # number of samples per class
  featureSummary <- data.frame(featureSummary)
  colnames(featureSummary) <- c("taxonID","numberOfSamples")
  capture.output(featureSummary, file = rf_output_file, append=TRUE)
  
  # features used to describe each sample (pixel)
  write("\ndescriptive features used to train this model: ", rf_output_file, append=TRUE) #newline
  write(colnames(features), rf_output_file, append=TRUE)
  
  # RF model summary, OOB error rate 
  capture.output(rf_model, file = rf_output_file, append=TRUE)
  
  # variable importance, ordered from highest MDGini to lowest
  write("\n20 most important variables, ranked by Mean Decrease Gini: \n", 
        rf_output_file, append=TRUE)
  varImp <- as.data.frame(rf_model$importance[order(rf_model$importance[,"MeanDecreaseGini"],decreasing = TRUE),])
  colnames(varImp)[colnames(varImp) == 'MeanDecreaseAccuracy'] <- 'MDAcc'
  colnames(varImp)[colnames(varImp) == 'MeanDecreaseGini'] <- 'MDGini'
  capture.output(varImp[1:20,],
              file = rf_output_file,
              #sep = "\t",
              #row.names=TRUE,
              #col.names = TRUE,
              append=TRUE)
  
  write("\n", rf_output_file, append=TRUE)
  
  # variable importance, ordered from highest MDA to lowest
  write("\n20 most important variables, ranked by MDA: \n", 
        rf_output_file, append=TRUE)
  varImp <- varImp[order(varImp$MDAcc, decreasing=TRUE),]
  capture.output(varImp[1:20,],
                 file = rf_output_file,
                 #sep = "\t",
                 #row.names=TRUE,
                 #col.names = TRUE,
                 append=TRUE)
  
  print("Top 10 most important variables ranked by MDA")
  print(rownames(varImp[1:10,])) 
  
  # # TO DO: keep the n most important variables and run the classification again?
  # print("KEEPING TOP 5 VARIABLES AND RUNNING RF AGAIN ")
  # if (keepMostImpVar == TRUE) {
  #   mostImpVar <- rownames(varImp[1:5,])
  #   features2 <- features %>% dplyr::select(taxonID, c(mostImpVar))
  #   # run RF model again, this time with reduced features
  #   set.seed(104)
  #   rf_model2 <- randomForest::randomForest(as.factor(features2$taxonID) ~ .,
  #                                          data=features2,
  #                                          importance=TRUE,
  #                                          ntree=ntree) # ntree is number of trees to grow
  #   print(rf_model2)
  # }
  # initial test shows that this does not significantly improve the accuracy
  
  
  
  if (createBoxplots == TRUE){
    print("Creating boxplots to compare variable distributions across species...")
    gg_alpha <- 1
    
    # ASPECT 
    ggplot(data = features, aes(x = taxonID, y = aspect, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: Aspect")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_ASPECT_",shapefileLayerNames$description[i],".png"))
    
    # SLOPE 
    ggplot(data = features, aes(x = taxonID, y = slope, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: Slope")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_SLOPE_",shapefileLayerNames$description[i],".png"))
    
    
    # CHM 
    ggplot(data = features, aes(x = taxonID, y = chm, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: CHM-derived height")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_CHM_",shapefileLayerNames$description[i],".png"))
    
    
    # NDVI 
    ggplot(data = features, aes(x = taxonID, y = NDVI, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: NDVI")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_NDVI_",shapefileLayerNames$description[i],".png"))
    
    
    # ARVI
    ggplot(data = features, aes(x = taxonID, y = ARVI, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: ARVI")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_ARVI_",shapefileLayerNames$description[i],".png"))
    
    
    # PRI 
    ggplot(data = features, aes(x = taxonID, y = PRI, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: PRI")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_PRI_",shapefileLayerNames$description[i],".png"))
    
    
    # NDLI
    ggplot(data = features, aes(x = taxonID, y = NDLI, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: NDLI")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_NDLI_",shapefileLayerNames$description[i],".png"))
    
    
    # NDNI
    ggplot(data = features, aes(x = taxonID, y = NDNI, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: NDNI")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_NDNI_",shapefileLayerNames$description[i],".png"))
    
    # PC 1 
    ggplot(data = features, aes(x = taxonID, y = PC1, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: PC1")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_PC1_",shapefileLayerNames$description[i],".png"))
    
    
    # rgb_mean_sd_B
    ggplot(data = features, aes(x = taxonID, y = rgb_mean_sd_B, colour = taxonID)) +
      geom_boxplot(alpha = gg_alpha, notch = TRUE, varwidth = TRUE, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 0.5) + 
      scale_color_brewer(palette="Spectral") + 
      ggtitle("Interspecies boxplot comparison: rgb_mean_sd_B")
    
    ggsave(paste0(out_dir, outDescription, "boxplot_rgb_mean_sd_B_",shapefileLayerNames$description[i],".png"))
    
    
    
  }
  
  
  write("\n\n------------------------------\n\n", rf_output_file, append=TRUE)
  
}

# close the text file
close(rf_output_file)

# TO DO: write the features to a text tfile 

end_time <- Sys.time()
print("Elapsed time: ")
print(end_time-start_time)

# create a nice tabl to summarize the model accuracies 
#https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html 
library(kableExtra)
rfAccuracies %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover","condensed", 
                                      full_width=F, align = "center"))

#https://www.littlemissdata.com/blog/prettytables
library(formattable)
colnames(rfAccuracies) <- c("Shapefile Description",
                            "OA [%]", "UA [%]", "PA [%]", "Kappa")
formattable(rfAccuracies, 
            align =c("l","c","c","c","c"), 
            list(`shapefileDescription` = formatter(
              "span", style = ~ style(color = "grey",font.weight = "bold")) 
            ))
  