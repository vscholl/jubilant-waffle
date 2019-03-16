# supporting R functions for the neon-veg workflow


# locate_woody_veg ---------------------------------------------------------

locate_woody_veg <- function(df){
  # Calculate precise geolocations in UTM for each mapped stem
  # as described on page 8 of the NEON veg structure user guide:
  # http://data.neonscience.org/api/v0/documents/NEON_vegStructure_userGuide_vA
  #
  # Args: 
  #   df
  #     data frame containing NEON woody_utm vegetation data from the 
  #     vst_mappingandtagging table
  #
  # Returns:
  #   woody_utm
  #     input data frame with two additional columns,
  #     easting and northing coordinates of mapped stems
  
  message("\nCalculating UTM coordinates of tree stems...")
  
  # remove all rows without stem distance & azimuth data
  woody_utm <- df[complete.cases(df$stemAzimuth) & 
                    complete.cases(df$stemDistance),]
  
  
  # use the geoNEON R package to pull geolocation data from the NEON API
  # get location information for each woody_utm veg entry. 
  # concatenate fields for namedLocation and pointID
  woody_utm$namedLocationPointID <- paste(woody_utm$namedLocation, woody_utm$pointID, sep=".")
  woody_utm_loc <- def.extr.geo.os(woody_utm, 'namedLocationPointID')
  
  
  # get easting/northing of reference point ID
  ref_east <- as.numeric(woody_utm_loc$api.easting)
  ref_north <- as.numeric(woody_utm_loc$api.northing)
  theta <- (woody_utm_loc$stemAzimuth * pi) / 180
  
  # calculate easting and northing for each plant
  # add new columns to the woody_utm veg data frame 
  woody_utm$easting <- ref_east + 
    woody_utm_loc$stemDistance * sin(theta)
  woody_utm$northing <- ref_north + 
    woody_utm_loc$stemDistance * cos(theta)
  
  return(woody_utm)
  
}


# woody_df_to_shp ---------------------------------------------------------

woody_df_to_shp = function(df, coord_ref, shrink, num_sides, shp_filename){
  # Creates a circular polygon geometry for each veg entry (row) 
  # within the specified data frame of NEON vegetation structures. 
  #
  # Args: 
  #   df: data frame containing NEON veg structure entries that each
  #         include the columns: maxCanopyDiameter, easting,
  #         northing, scientificname, taxonid, individualid,
  #         stemheight.
  #   coord_ref: Class "CRS" of Coordinate Reference System Arguments.
  #         Can be obtained from existing R layer with existing
  #         "coord .ref" field by using layer@crs
  #   shrink: numeric integer or decimal; factor by which to divide 
  #         polygon radii to reduce the size of crown polygons. 
  #   num_sides: integer; number of sides in each circular polygon
  #   shp_filename: string; output filename to write shapefile 
  # 
  # Returns:
  #   spdfs:
  #         spatialPolygonsDataFrame containing crown polygons
  
  message("\nCreating circular polygons based on stem location and crown diameter...")
  
  # for each tree 
  for (i in 1:length(df$easting)){
    
    # calculate radius of circular polygon [units of km] 
    r = round(df$maxCrownDiameter[i]) / 2
    
    # reduce radius by shrink factor to reduce polygon size 
    r = signif( (r / as.numeric(shrink)), 2)
    
    # create circular polygon coordinates with center location, 
    # number of sides, units of distance, polygon calculation
    c <- circle.polygon(df$easting[i], 
                        df$northing[i], 
                        r, sides=num_sides, 
                        by.length=FALSE, 
                        units="km", 
                        poly.type = "cartesian")
    
    p = Polygon(c)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    
    # set coordinate reference system
    proj4string(sps) = coord_ref
    
    # create data to associate with the polygon coordinates 
    d = data.frame(scientificName = as.character(df$scientificName[i]),
                   taxonID = as.character(df$taxonID[i]),
                   individualID = df$individualID[i],
                   crownDiam = df$maxCrownDiameter[i],
                   height = df$height[i])
    
    # turn the SpatialPolygons object into a SpatialPolygonsDataFrame
    spdf_out = SpatialPolygonsDataFrame(sps,d)
    
    # check for the first loop iteration
    if (i> 1){
      # If it's a later loop iteration, merge polygons together 
      spdfs<-rbind(spdf_out, spdfs, makeUniqueIDs = TRUE)
    } else { 
      # If it's the first loop iteration, assign the SpatialPolygonsDataFrame to a variable
      spdfs <- spdf_out
    }
  }
  
  # write polygon(s) to shapefile  
  suppressWarnings(
    writeOGR(spdfs, 
             getwd(),
             shp_filename, 
             driver="ESRI Shapefile", 
             overwrite_layer = TRUE))
  
  return(spdfs)
}


# merge_vst_tables --------------------------------------------------------

merge_vst_tables <- function(vst_mapping, vst_individual){
  # Merges species data with vegetation structure data
  # based on matching individual ID's. 
  # 
  # inputs
  #   vst_mapping: 
  #     data frame containing scientific names of mapped plants (vst_mappingandtagging)
  #   vst_individual:
  #     data frame containing vegetation structure (vst_apparentindividual)
  #
  # output
  #   merged:
  #     input data frame with additional columns of information
  
  message("\nMerging vegetation structure tables...")  
  
  # find matching individual ID's between the tables
  i <- match(as.character(vst_mapping$individualID),
             as.character(vst_individual$individualID))
  
  # check for duplicates and keep most recent measurement? 
  
  # add height and crown diameter columns to mapping table
  merged <- vst_mapping
  merged$height <- vst_individual$height[i]
  merged$maxCrownDiameter <- vst_individual$maxCrownDiameter[i]
  # replace the date column with apparentindividual 
  # since the location of the plant will not change but the height 
  # may between measurements
  merged$date <- vst_individual$date[i]
  
  # keep only entries that have height and crown diameter values
  merged <- merged[complete.cases(merged$height) & 
                     complete.cases(merged$maxCrownDiameter),]
  
  
  
  return(merged)
}


# get_vst_crs -----------------------------------------------------------

get_vst_crs <- function(woody_path){
  
  print("Retrieving CRS...")
  
  # define path to vst_plotperyear table,
  # which contains geodetic datum and UTM zone 
  vst_path <- paste(woody_path, 
                    list.files(path = woody_path, 
                               pattern = "plotperyear"), 
                    sep="/")
  
  vst_data <- read.csv(vst_path)
  
  datum <- as.character(vst_data$geodeticDatum[1])
  zone <- gsub("[^0-9\\.]", "", vst_data$utmZone[1])
  coord_ref <- CRS(paste("+proj=utm +zone=",zone, 
                         " +datum=",datum," +units=m",sep=""))
  return(coord_ref)
}


# list_tiles_with_plants --------------------------------------------------

list_tiles_with_plants <- function(woody,out_dir){
  # generate a list of 1km x 1km tiles containing field data
  # write list to a text file in the output directory provided
  
  print("Generating list of tiles containing stems...")
  
  # get easting, northing coordinates of all plants
  e <- NA
  n <- NA
  for (i in 1:dim(woody)[1]){
    easting <- floor(as.numeric(woody$easting[i])/1000)*1000
    northing <- floor(as.numeric(woody$northing[i])/1000)*1000
    e[i] <- easting
    n[i] <- northing
  }
  
  # find unique rows repesenting tiles
  easting_northings <- data.frame(as.character(e),as.character(n))
  colnames(easting_northings) <- c('e','n')
  tiles <- unique(easting_northings[,c('e','n')])
  
  # order by ascending tile coordinates 
  tiles <- tiles %>%
    arrange(e)
  
  # write to text file 
  tile_names <- paste(tiles$e, tiles$n, sep="_")
  tiles_file <- file(paste(out_dir,"list_tiles.txt", sep=""))
  writeLines(tile_names, tiles_file)
  close(tiles_file)
  
  return(tiles)
}


# apply_area_threshold ----------------------------------------------------

apply_area_threshold <- function(df, nPix){ 
  # Applies an area threshold to remove polygons
  # with an area smaller than the number of [1m^2] pixels
  # specified by numPix. 
  #
  # Args
  #   df
  #     data frame containing woody veg entries, including
  #     the maxCrownDiameter measurement
  #
  #   nPix
  #     number of pixels describing the area required to 
  #     keep polygons. 
  
  print("Removing polygons below area threshold...")
  
  # area of image pixels [cm^2] for thresholding sub-pixel plants 
  px_area_rgb <- 25 * 25 #[cm^2]
  # gridded LiDAR products and HS pixels are 1m x 1m
  px_area_hs <- 16 * px_area_rgb
  
  # multiply area of 1 pixel by the numPix input parameter
  thresh <- px_area_hs * nPix
  
  # calculate approximate area [cm^2] of each plant based on diameter 
  # keep only values > 0
  diam_cm <- (df$maxCrownDiameter) * 100
  area_cm <- pi * ((diam_cm / 2)^2)
  df$area_cm <- area_cm
  
  # filter crowns with area < thresh
  df <- df %>%
    filter(area_cm > thresh) 
  
  return(df)
  
}


# polygon_overlap ---------------------------------------------------------

polygon_overlap <- function(df, nPix, shp_filename){
  
  message("\nChecking for overlap between polygons...")
  
  # multiply area of 1 pixel by the numPix input parameter
  # rgb pixels are 25cm, 16 of them in a single HS pixel
  thresh <- 25 * 25 * 16 * nPix # [cm^2]
  
  # order polygons from largest to smallest crown diameter
  polys_ordered <- df[order(df$crownDiam, 
                            decreasing = TRUE),]
  
  # create a polygon data frame to update with filtered/deleted/merged entries
  polys_filtered <- polys_ordered
  
  # create an empty list of polygon pairs that have been compared 
  compared_pairs <- list();
  c <- 1 # index for appending to list
  
  for (individualID in polys_ordered@data$individualID){
    
    # if this polygon was removed from the polys_filtered
    # data frame in a previous iteration, skip it 
    if(sum(polys_filtered$individualID==individualID) == 0){
      next
    }
    
    # extract current vs. all other polygon from data frame
    current_poly <- get_poly(polys_filtered, index_type = 'id', number = individualID)
    other_polys <- polys_filtered[polys_filtered$individualID!=current_poly$individualID,]
    
    # check for overlap between current polygon and all polygons
    overlap <- raster::intersect(current_poly, other_polys)
    n_overlap <- length(overlap)
    
    if(n_overlap>0){ 
      for (o in 1:n_overlap){
        
        # if current polygon ID is not in filtered set
        if(sum(polys_filtered$individualID==individualID) == 0){
          break
        }
        
        # get area and height of current and test polygons that overlap
        current_poly.area <- round(current_poly@polygons[[1]]@area, 2)
        current_poly.height <- current_poly@data$height
        test_poly <- get_poly(polys_filtered, 
                              index_type = 'id', 
                              number = overlap@data$individualID.2[[o]])
        test_poly.area <- test_poly@polygons[[1]]@area
        test_poly.height <- test_poly@data$height
        
        # combine the ID's of the current and test polygons
        # to keep track of which pairs have been compared 
        id.pair <- paste(current_poly@data$individualID,
                         test_poly@data$individualID,
                         sep = " ")
        
        # if polygon pair was already compared, skip to the next pair
        if (id.pair %in% compared_pairs) {
          
          next
          
        } else { 
          
          # add to the list of polygon pairs that have been compared 
          compared_pairs[[c]] <- id.pair
          
          c <- c + 1 # increment index 
          
          # add opposite combination of polygons
          compared_pairs[[c]] <- paste(test_poly@data$individualID,
                                       current_poly@data$individualID,
                                       sep = " ")
          # increment index 
          c <- c + 1
          
        }
        
        # if test polygon is not in filtered set, skip to the next overlapping polygon
        if(sum(polys_filtered$individualID==test_poly$individualID) == 0){
          next
        }
        
        # get area of the overlap between current polygon and test polygon
        overlap.area <- round(raster::intersect(current_poly, test_poly)@polygons[[1]]@area, 2)
        clip <- gIntersection(current_poly, test_poly, byid = TRUE, drop_lower_td = TRUE)
        
        # if total overlap
        if(current_poly.area == overlap.area | 
           test_poly.area == overlap.area | 
           (current_poly.area - overlap.area) < 1){
          if(current_poly.height > test_poly.height){ 
            # then remove test polygon from filtered list
            polys_filtered <- polys_filtered[polys_filtered$individualID!=test_poly$individualID,]
          } else{ 
            # remove current polygon
            polys_filtered <- other_polys
            break
          }
          
          
        } else { # partial overlap
          
          # test polygon is taller
          if(test_poly.height > current_poly.height){
            clipped <- raster::erase(current_poly,
                                     raster::crop(current_poly, test_poly))
            
          } else { 
            # current polygon is taller 
            clipped <- raster::erase(test_poly,
                                     raster::crop(test_poly, current_poly))
          }
          
          # if there is no clipped area, skip to the next overlap polygon
          if(length(clipped) == 0){
            
            next
            
          } else{ 
            
            # if clipped area exceeds area threshold
            if(clipped@polygons[[1]]@area * 10000 > thresh){
              
              # replace current ploygon with clipped polygon
              j <- which(polys_filtered@data$individualID == clipped$individualID)
              
              if(j==1){
                polys_filtered <-  rbind(clipped, 
                                         polys_filtered[(j+1):nrow(polys_filtered),])
                
              } else if(j==nrow(polys_filtered)) { # j is the last index
                
                polys_filtered <- rbind(polys_filtered[1:(j-1),],
                                        clipped)
                
              } else{
                polys_filtered <- rbind(polys_filtered[1:(j-1),],
                                        clipped, 
                                        polys_filtered[(j+1):nrow(polys_filtered),])
              }
              
              
            } else{
              # else, remove current polygon since clipped area is too small
              polys_filtered <- polys_filtered[polys_filtered$individualID!=clipped$individualID,]
            }
          }
        }
      }
    }
  }
  
  print(polys_filtered)
  
  # write final polygons to file after checking for overlap
  writeOGR(polys_filtered, getwd(),
           paste(shp_filename),
           driver="ESRI Shapefile", overwrite_layer = TRUE)
  
  return(polys_filtered)
  
}


# get_poly ---------------------------------------------------------------

get_poly = function(spdf, index_type, number){
  # this fuction extracts a single polygon inside a SpatialPolygonsDataFrame 
  # object based on its individual ID OR index in the data frame for testing
  # index_type should be either "id" or "index"
  # if index_type == id, then search for matching entry with that individualID
  # if index_type == index, then it use number as an index into data frame
  
  if(index_type=='id'){ # use number as individual ID 
    i <- which(spdf$individualID == number)
  } else { # use number as index
    i <- number
  }
  
  coords = spdf@polygons[[i]]@Polygons[[1]]@coords
  extra_data = as.data.frame(spdf@data[spdf@data$individualID == spdf$individualID[i],], 
                             row.names = as.character(spdf$individualID[i]))
  
  # create SpatialPolygons
  P1 = Polygon(coords)
  Ps1 = SpatialPolygons(list(Polygons(list(P1), ID = spdf$individualID[i])), 
                        proj4string=spdf@proj4string)
  
  # create SpatialPolygonsDataFrame
  Ps1 = SpatialPolygonsDataFrame(Ps1, 
                                 data = extra_data, match.ID = TRUE)
  return(Ps1)
  
}


# check_create_dir --------------------------------------------------------

check_create_dir <- function(new_dir){
  # check if directory exists. If it doesn't create it. 
  if (!dir.exists(new_dir)){
    dir.create(new_dir)
  }
}


# make_species_table ------------------------------------------------------

make_species_table <- function(df){
  # Counts the number of entries per species in the data frame. 
  # Creates a table with count, taxon ID, scientific name. 
  #
  # Args:
  #   df 
  #     data frame with the following attributes per entry: 
  #     scientificName, taxonID
  #
  # Returns:
  #   species_table 
  #     tibble containing a summary of taxon ID, scientific name, 
  #     and count information
  
  species_table  <- df %>% 
    mutate(sciNameShort = word(scientificName, 1, 2)) %>% 
    group_by(taxonID, sciNameShort) %>% 
    summarise(total = n()) %>%
    dplyr::select(total, taxonID, sciNameShort)
  
  return(species_table)
}


# df_to_shp_points --------------------------------------------------------

df_to_shp_points = function(df, coord_ref, shp_filename){
  # Creates a point shapefile for veg entries (rows) 
  # within the specified data frame of NEON vegetation structures. 
  #
  # Args: 
  #   df: data frame containing NEON veg structure entries that each
  #         include the columns: easting, northing, scientificname,
  #         taxonid, individualid
  #   coord_ref: Class "CRS" of Coordinate Reference System Arguments.
  #         Can be obtained from existing R layer with existing
  #         "coord .ref" field by using layer@crs
  #   shp_filename: string; output filename to write shapefile 
  # 
  # Returns:
  #   spdfs:
  #         spatialPointsDataFrame containing mapped stems as points
  
  message("\nCreating points based on mapped stem locations...")
  
  # select columns of interest
  if("height" %in% colnames(df)){
    stem_locations <- df %>%
      dplyr::select(easting, northing, individualID, scientificName, taxonID, height,maxCrownDiameter)
    
  } else{
    stem_locations <- df %>%
      dplyr::select(easting, northing, individualID, scientificName, taxonID)
    
  }
  
  # assign UTM coordinates to create SpatialPointsDataFrame
  coordinates(stem_locations) <- ~easting+northing
  proj4string(stem_locations) <- coord_ref
  
  # write polygon(s) to shapefile  
  suppressWarnings(
    writeOGR(stem_locations, 
             getwd(),
             shp_filename, 
             driver="ESRI Shapefile", 
             overwrite_layer = TRUE))
  
  return(stem_locations)
  
}


# apply_height_threshold --------------------------------------------------

apply_height_threshold <- function(df, ht){ 
  # Applies a height threshold to remove polygons
  # with a height smaller than ht (meters)
  #
  # Args
  #   df
  #     data frame containing woody veg entries, including
  #     the height measurement
  #
  #   ht
  #     height (m) describing the height required to 
  #     keep polygons. 
  
  print("Removing polygons below height threshold...")
  
  # filter crowns with area < thresh
  df <- df %>%
    filter(height > ht) 
  
  return(df)
  
}


# allometry_height_diam ---------------------------------------------------

allometry_height_diam <- function(df){
  # calculates linear regression models for each taxonID 
  # within the input data frame. 
  # Dependent variable: crown diameter (m) 
  # Independent variable: height (m)
  
  # libraries: ggplot2, tidyr, purrr, dplyr, broom 
  
  
  # nest the data by taxon ID: 
  # create separate list of data for each species
  by_taxonID <- stems_final %>% 
    group_by(taxonID) %>% 
    filter(n() > 1) %>% 
    nest()
  
  # function to calculate linear model of crownDiam 
  # as a function of height 
  taxon_model <- function(df){
    lm(crownDiam ~ height, data = df)
  }
  
  # for each species data set, calculate linear model
  models <- by_taxonID %>% 
    mutate(
      model = data %>% map(taxon_model)
    )
  
  # compute summary, r-squared, parameters, and observation statistics 
  # for each linear model using the "broom" library 
  models <- models %>% 
    mutate( 
      glance  = map(model, broom::glance),
      rsq     = glance %>% map_dbl("r.squared"),
      tidy    = map(model, broom::tidy),
      augment = map(model, broom::augment)
    )
  
  # plot data points and fitted models using ggplot 
  g <- models %>%
    unnest(data) %>%
    ggplot(aes(height, crownDiam)) +
    geom_point(shape=1) + 
    facet_wrap(~taxonID) +    
    geom_smooth(method=lm) +
    labs(x = "height (m)", y = "crown diameter (m)")
  print(g)
  
  return(models)
  
}

remove_duplicates <- function(df){
  # This function checks for duplicate individualID 
  # values and keeps only the most recent entry
  
  df_no_duplicates <- df %>% 
    group_by(individualID) %>%
    slice(which.max(as.Date(date)))
  
  return(df_no_duplicates)
  
}



# plt_hs_rgb -------------------------------------------------------------

plot_hs_rgb <- function(refl, wavelengths, rgb.bands, proj4, ext, plt=FALSE){
  
  # extract the red, green, and blue bands from a hyperspectral 
  # data frame and plot the RGB image
  # 
  # Args: 
  #   refl
  #     array containing hyperspectral reflectance data 
  #
  #   wavelengths 
  #     1D array of wavelengths for all hyperspectral bands
  # 
  #   rgb.bands
  #     numeric vector of wavelengths for R, G, and B bands.
  #     for example: rgb.bands = c(620, 555, 450)
  #
  #   proj4
  #     character string describing the projection information for the image 
  #
  #   ext
  #     Extent object containing the xmin, xmax, ymin, and ymax of the image
  #
  #   plt
  #     boolean variable to plot (TRUE) or not plot (FALSE) the RGB composite
  
  # RED 
  r <- get_hs_band(refl = refl, 
                   wavelengths = wavelengths,
                   wl = rgb.bands[1],
                   proj4 = proj4,
                   ext = ext,
                   plt = FALSE)
  
  wl <- 549
  g <- get_hs_band(refl = refl, 
                   wavelengths = wavelengths,
                   wl = rgb.bands[2],
                   proj4,
                   ext,
                   plt = FALSE)
  
  wl <- 474
  b <- get_hs_band(refl = refl, 
                   wavelengths = wavelengths,
                   wl = rgb.bands[3],
                   proj4,
                   ext,
                   plt = FALSE)
  
  rgb_stack <- stack(r,g,b)
  rgb_brick <- brick(rgb_stack)
  
  if (plt == TRUE){
    plotRGB(rgb_brick,
            r = 1, g = 2, b = 3,
            stretch = "lin",
            axes = TRUE,
            main="RGB Composite",
            xlab="Easting (m)",
            ylab="Northing (m)",
            cex.main=2)
  }
  
  return(rgb_brick)
  
  
}


# get_hs_band -------------------------------------------------------------

get_hs_band <- function(refl, wavelengths, wl, proj4, ext, plt=FALSE){
  # extract a HS band with the specified wavelength (wl)
  # convert band to a raster, assign CRS and geographic extent, 
  # transpose for proper plotting orientation
  #
  # Args: 
  #   refl - array with reflectance data for all wavelengths
  #   wavelengths - vector containing all wavelengths
  #   wl - wavelength(integer [nm]) to select
  #   proj4 - projection string describing CRS for spatial object
  #   ext - extent of the image to assign georeferences coordinates
  #   plt - boolean variable to plot (TRUE) or not plot (FALSE) the logged refl
  #
  # Returns:
  #   refl.t 
  #     transposed specified band of reflectance data with CRS and extent;
  #     ready to plot. 
  
  
  # get index of closest element in the wavelengths list
  wl.idx <- which.min(abs(wavelengths - wl))
  
  # extract reflectance data for the wavelength. convert array to matrix 
  refl.plot <- refl[wl.idx,,]
  
  # transpose x and y values for proper orientation in plot 
  refl.t <- t(refl.plot)
  
  # create raster and assign CRS
  refl.ras <- raster(refl.t,
                     crs = proj4)
  
  # assign UTM extent to raster
  extent(refl.ras) <- ext
  
  # plot 
  if (plt == TRUE){
    image(log(refl.ras), 
          main= paste("Band",
                      as.character(wl),
                      "nm",
                      sep= " "))
  }
  
  return(refl.ras)
  
  
}



# Stack_hyperspectral -----------------------------------------------------

stack_hyperspectral <- function(h5, out_dir){
  # This function creates a rasterstack object for the specified HDF5 
  # filename. 
  #
  # Args: 
  # h5
  #   character string filename of HDF5 file 
  # out_dir
  #   directory for output files where the wavelengths will be written to
  #   a text file for further analysis
  #
  # Returns: 
  # s
  #   RasterStack (collection of Raster layers with the same spatial extent
  #   and resolution) containing nrows x ncols x nbands, based on the 
  #   resolution and number of bands in the HDF5 file. This RasterStack
  #   can then be clipped using Spatial vector layers. 
  #
  
  
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
  # transpose the image pixels for proper orientation to match
  # the other layers. create a raster for this band and assign
  # the CRS.
  print("Transposing reflectance data for proper orientation")
  r1 <- raster::t(raster::raster(r1, crs = crs_info$Proj4))
  extent(r1) <- tile_extent
  
  # start the raster stack with first band 
  s <- raster::stack(r1)
  
  # loop through bands and create a giant rasterstack with 426 (n_bands) bands
  for(b in 2:n_bands){
    print(b)
    
    # create raster with current band
    r <- (refl_scaled[b,,]) # convert to matrix
    r <- raster::t(raster::raster(r, crs = crs_info$Proj4))
    extent(r) <- tile_extent
    
    # add additional band to the stack with the addLayer function
    s <- raster::addLayer(s, r)
    
  }
  
  # adjust the names for each layer in raster stack to correspond to wavelength
  names(s) <- round(wavelengths)
  
  # write wavelengths to a text file 
  # write the exact wavelengths to file for future use 
  write.table(data.frame(wavelengths = wavelengths),
              paste0(out_dir,"wavelengths.txt"),
              sep="\n",
              row.names=FALSE)
  
  # return the stacked hyperspectral data to clip with vector files 
  return(s)
  
}



# write_spectra_to_file ---------------------------------------------------

write_spectra_to_file <- function(spectra, trees_in, filename_out){
  # this function takes a data frame containing extracted data.
  # column names include ID (this refers to the individual plant ID),
  # columns for each wavelength of hyperspectral reflectance ("X381","X386", 
  # ... "X2505"), as well as chm (height), slope, aspect, and pixelNumber 
  # (unique integer for every pixel in the 1km x 1km raster). 
  # these values are combined with metadata (from the shapefiles list). 
  
  # create tree metadata data frame 
  # tree_metadata <- data.frame(individualID = trees_in$indvdID,
  #                             scientificName = trees_in$scntfcN,
  #                             taxonID = trees_in$taxonID,
  #                             maxCrownDiameter = trees_in$mxCrwnD,
  #                             height = trees_in$height,
  #                             X = trees_in$X,
  #                             Y = trees_in$Y,
  #                             # create ID column to pair the tree metadata with 
  #                             # extracted spectra data 
  #                             ID =  1:nrow(trees_in))
  
  tree_metadata <- data.frame(trees_in) %>% 
    mutate(ID = 1:nrow(trees_in))
  
  # create a list of increasing integer counts to keep track of how many rows 
  # (pixels or spectra) belong to each tree 
  for (i in unique(spectra$ID)){
    if(i==1){
      counts = 1:sum(spectra$ID==i)
    }
    else{
      counts = append(counts, 1:sum(spectra$ID==i))
    }
  }
  
  # combine the additional data with each spectrum for writing to file
  spectra_write <- merge(tree_metadata,
                         spectra,
                         by="ID") %>% 
    mutate(spectra_count = counts)%>% 
    select(ID, spectra_count, pixelNumber, everything())
  
  # take a look at the first rows of the spectra data to write 
  #head(spectra_write)
  
  
  # write the spectral data to file for future analysis 
  write.csv(spectra_write, file = filename_out) 
  
  return(spectra_write)
  
}



# createRibbonPlot --------------------------------------------------------


createRibbonPlot <- function(wavelengths, reflFilename){
  
  # define the "bad bands" wavelength ranges in nanometers, where atmospheric 
  # absorption creates unreliable reflectance values. 
  bad_band_window_1 <- c(1340, 1445)
  bad_band_window_2 <- c(1790, 1955)

  taxonList <- c("ABLAL","PICOL","PIEN","PIFL2")
  
  
  
  # absolute maximum reflectance to set the same ylimit for the plots
  y_max <- 0.35    #max(refl_tidy$max_reflectance, na.rm = TRUE)
  
  # remove the bad bands 
  remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                                 wavelengths < bad_band_window_1[2]) | 
                                (wavelengths > bad_band_window_2[1] & 
                                   wavelengths < bad_band_window_2[2])]
  
  # remove columns that contain "X" in their name but are not reflectance values 
  df <- as.data.frame(read.csv(reflFilename)) %>% 
          dplyr::select(-c(X.1,X,Y)) %>% 
            dplyr::filter(taxonID %in% taxonList)
  
  # filter the columns to only keep those with spectral reflectance
  spectra_all <- df %>% select( colnames(df)[ grepl( "X", names(df))] ) 
  
  # sanity check - check the number of unique entries in the spectra set 
  print(paste("There are", as.character(length(unique(df$indvdID))),
              "unique individual IDs for the spectra being plotted"))
  
  print(paste("There are", as.character(length(unique(df$pixelNumber))),
              "unique pixelNumbers for the spectra being plotted"))
  
  # calculate mean reflectance per species
  mean_reflectance <- stats::aggregate(spectra_all, 
                                       by = list(taxonID = df$taxonID),
                                       FUN = mean) 
  min_reflectance <- stats::aggregate(spectra_all, 
                                      by = list(taxonID = df$taxonID),
                                      FUN = min) 
  max_reflectance <- stats::aggregate(spectra_all, 
                                      by = list(taxonID = df$taxonID),
                                      FUN = max) 
  sd_reflectance <- stats::aggregate(spectra_all,
                                     by = list(taxonID = df$taxonID),
                                     FUN = sd)
  
  # create a LUT that matches actual wavelength values with the column names,
  # X followed by the rounded wavelength values. 
  wavelength_lut <- data.frame(wavelength = wavelengths,
                               xwavelength = paste0("X",as.character(round(wavelengths))),
                               stringsAsFactors = FALSE)
  
  # use the gather function makes wide data longer:
  # https://uc-r.github.io/tidyr 
  # so the reflectance data can easily be grouped by species, 
  # and the mean/min/max reflectance values can be selected for a ribbon plot. 
  mean_refl_tidy <- tidyr::gather(mean_reflectance,
                                  key = xwavelength,
                                  value = "mean_reflectance",
                                  X381:X2510) %>%
    dplyr::left_join(wavelength_lut, by="xwavelength") 
  
  # add on the max, min reflectance columns with the same format 
  max_refl_tidy <- tidyr::gather(max_reflectance,
                                 key = xwavelength,
                                 value = "max_reflectance",
                                 X381:X2510)
  
  min_refl_tidy <- tidyr::gather(min_reflectance,
                                 key = xwavelength,
                                 value = "min_reflectance",
                                 X381:X2510)
  
  sd_refl_tidy <- tidyr::gather(sd_reflectance,
                                key = xwavelength,
                                value = "sd_reflectance",
                                X381:X2510)
  
  # combine the mean, min, man reflectance data into one long data frame
  refl_tidy <- merge.data.frame(mean_refl_tidy,
                                max_refl_tidy) %>% 
    merge.data.frame(min_refl_tidy) %>% 
    merge.data.frame(sd_refl_tidy) %>% 
    select(-xwavelength) %>%          # remove the Xwavelength values 
    select(wavelength, everything())  # reorder to wavelength column is first
  
  
  # remove the first reflectance value 
  refl_tidy <- refl_tidy[refl_tidy$wavelength > 385,]
  
  # remove the bad bands 
  refl_tidy$mean_reflectance[refl_tidy$wavelength %in% remove_bands] <- NA
  refl_tidy$max_reflectance[refl_tidy$wavelength %in% remove_bands] <- NA
  refl_tidy$min_reflectance[refl_tidy$wavelength %in% remove_bands] <- NA
  refl_tidy$sd_reflectance[refl_tidy$wavelength %in% remove_bands] <- NA
  
  # add and subtract one standard deviation from the mean 
  refl_tidy$mean_plus_sd <- refl_tidy$mean_reflectance + refl_tidy$sd_reflectance
  refl_tidy$mean_minus_sd <- refl_tidy$mean_reflectance - refl_tidy$sd_reflectance
  # set any negative values to zero for proper plotting 
  refl_tidy$mean_minus_sd[refl_tidy$mean_minus_sd<0] <- 0
  
  # specify the colors for the reflectance curves & shading around them 
  shading_colors <- c("#d7191c", "#fdae61", "#abdda4", "#2b83ba")
  species <- sort(unique(df$taxonID)) #alphabetical so colors match plot above
  shading_alpha <- 0.4
  
  # generate the ribbon plot
  ggplot(refl_tidy, 
         aes(x = wavelength, y = mean_reflectance, color = taxonID)) + 
    
    # ABLAL
    geom_ribbon(data = refl_tidy[refl_tidy$taxonID == species[1], ],
                aes(ymin = mean_minus_sd,
                    ymax = mean_plus_sd), # std dev shading
                colour=NA,
                alpha = shading_alpha,
                fill = shading_colors[1],
                show.legend = F) + 
    
    # PICOL
    geom_ribbon(data = refl_tidy[refl_tidy$taxonID == species[2], ],
               aes(ymin = mean_minus_sd,
                    ymax = mean_plus_sd), # std dev shading
                colour=NA,
                alpha = shading_alpha,
                fill = shading_colors[2],
                show.legend = F) + 
    
    # PIEN
    geom_ribbon(data = refl_tidy[refl_tidy$taxonID == species[3], ],
                aes(ymin = mean_minus_sd,
                    ymax = mean_plus_sd), # std dev shading
                colour=NA,
                alpha = shading_alpha,
                fill = shading_colors[3],
                show.legend = F) + 
    
    # PIFL2
    geom_ribbon(data = refl_tidy[refl_tidy$taxonID == species[4], ],
                aes(ymin = mean_minus_sd,
                    ymax = mean_plus_sd), # std dev shading
                colour=NA,
                alpha = shading_alpha,
                fill = shading_colors[4],
                show.legend = F) + 
    
    
    # mean reflectance line
    # placing this after the ribbon shading so the mean curves are visible
    geom_line(size = 0.5, alpha = 1) + 
    
    scale_color_manual(values = shading_colors) + 
    
    # hide the "alpha" legend
    guides(alpha=FALSE) + 
    
    # label X and Y axes 
    labs(x = "wavelength (nm)", y = "reflectance") + 
    
    # set the y axis range to be consistent between plots
    ylim(0,y_max) + 
    
    # main plot title  
    ggtitle(paste0("Mean Hyperspectral reflectance per species: ", 
                   tools::file_path_sans_ext(str_split(basename(reflFilename),"ALL_")[[1]][2]),
                   " \n",
                   # std dev shading
                   "(shading shows one standard deviation from mean refl range per wavelength)"))
  
    # write plot to file 
    ggsave(paste0(out_dir,"/figures/","ribbon_plot_", 
                  tools::file_path_sans_ext(basename(reflFilename)), ".png"), 
           width = 10, height = 6)

  
}


# createRibbonPlot --------------------------------------------------------


createSeparateRibbonPlots <- function(wavelengths, reflFilename){
  
  # define the "bad bands" wavelength ranges in nanometers, where atmospheric 
  # absorption creates unreliable reflectance values. 
  bad_band_window_1 <- c(1340, 1445)
  bad_band_window_2 <- c(1790, 1955)
  
  taxonList <- c("ABLAL","PICOL","PIEN","PIFL2")
  
  
  
  # absolute maximum reflectance to set the same ylimit for the plots
  y_max <- 0.35    #max(refl_tidy$max_reflectance, na.rm = TRUE)
  
  # remove the bad bands 
  remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                                 wavelengths < bad_band_window_1[2]) | 
                                (wavelengths > bad_band_window_2[1] & 
                                   wavelengths < bad_band_window_2[2])]
  
  # remove columns that contain "X" in their name but are not reflectance values 
  df <- as.data.frame(read.csv(reflFilename)) %>% 
    dplyr::select(-c(X.1,X,Y)) %>% 
    dplyr::filter(taxonID %in% taxonList)
  
  # filter the columns to only keep those with spectral reflectance
  spectra_all <- df %>% select( colnames(df)[ grepl( "X", names(df))] ) 
  
  # sanity check - check the number of unique entries in the spectra set 
  print(paste("There are", as.character(length(unique(df$indvdID))),
              "unique individual IDs for the spectra being plotted"))
  
  print(paste("There are", as.character(length(unique(df$pixelNumber))),
              "unique pixelNumbers for the spectra being plotted"))
  
  # calculate mean reflectance per species
  mean_reflectance <- stats::aggregate(spectra_all, 
                                       by = list(taxonID = df$taxonID),
                                       FUN = mean) 
  min_reflectance <- stats::aggregate(spectra_all, 
                                      by = list(taxonID = df$taxonID),
                                      FUN = min) 
  max_reflectance <- stats::aggregate(spectra_all, 
                                      by = list(taxonID = df$taxonID),
                                      FUN = max) 
  sd_reflectance <- stats::aggregate(spectra_all,
                                     by = list(taxonID = df$taxonID),
                                     FUN = sd)
  
  # create a LUT that matches actual wavelength values with the column names,
  # X followed by the rounded wavelength values. 
  wavelength_lut <- data.frame(wavelength = wavelengths,
                               xwavelength = paste0("X",as.character(round(wavelengths))),
                               stringsAsFactors = FALSE)
  
  # use the gather function makes wide data longer:
  # https://uc-r.github.io/tidyr 
  # so the reflectance data can easily be grouped by species, 
  # and the mean/min/max reflectance values can be selected for a ribbon plot. 
  mean_refl_tidy <- tidyr::gather(mean_reflectance,
                                  key = xwavelength,
                                  value = "mean_reflectance",
                                  X381:X2510) %>%
    dplyr::left_join(wavelength_lut, by="xwavelength") 
  
  # add on the max, min reflectance columns with the same format 
  max_refl_tidy <- tidyr::gather(max_reflectance,
                                 key = xwavelength,
                                 value = "max_reflectance",
                                 X381:X2510)
  
  min_refl_tidy <- tidyr::gather(min_reflectance,
                                 key = xwavelength,
                                 value = "min_reflectance",
                                 X381:X2510)
  
  sd_refl_tidy <- tidyr::gather(sd_reflectance,
                                key = xwavelength,
                                value = "sd_reflectance",
                                X381:X2510)
  
  # combine the mean, min, man reflectance data into one long data frame
  refl_tidy <- merge.data.frame(mean_refl_tidy,
                                max_refl_tidy) %>% 
    merge.data.frame(min_refl_tidy) %>% 
    merge.data.frame(sd_refl_tidy) %>% 
    select(-xwavelength) %>%          # remove the Xwavelength values 
    select(wavelength, everything())  # reorder to wavelength column is first
  
  
  # remove the first reflectance value 
  refl_tidy <- refl_tidy[refl_tidy$wavelength > 385,]
  
  # remove the bad bands 
  refl_tidy$mean_reflectance[refl_tidy$wavelength %in% remove_bands] <- NA
  refl_tidy$max_reflectance[refl_tidy$wavelength %in% remove_bands] <- NA
  refl_tidy$min_reflectance[refl_tidy$wavelength %in% remove_bands] <- NA
  refl_tidy$sd_reflectance[refl_tidy$wavelength %in% remove_bands] <- NA
  
  # add and subtract one standard deviation from the mean 
  refl_tidy$mean_plus_sd <- refl_tidy$mean_reflectance + refl_tidy$sd_reflectance
  refl_tidy$mean_minus_sd <- refl_tidy$mean_reflectance - refl_tidy$sd_reflectance
  # set any negative values to zero for proper plotting 
  refl_tidy$mean_minus_sd[refl_tidy$mean_minus_sd<0] <- 0
  
  # specify the colors for the reflectance curves & shading around them 
  shading_colors <- c("#d7191c", "#fdae61", "#abdda4", "#2b83ba")
  species <- sort(unique(df$taxonID)) #alphabetical so colors match plot above
  shading_alpha <- 0.4
  
  # FACET WRAP mean reflectance curves with ribbon (+/- 1 SD)
  ggplot(refl_tidy, aes(x = wavelength, 
                        y = mean_reflectance)) + 
    facet_wrap(~taxonID) +
    geom_ribbon(alpha = 0.5, aes(ymin = mean_minus_sd, 
                                 ymax = mean_plus_sd, 
                                 group = taxonID,
                                 fill = taxonID)) + 
    geom_line(aes(color = taxonID)) + 
    scale_color_manual(values = shading_colors) + 
    scale_fill_manual(values = shading_colors) + 
    
    # label X and Y axes 
    labs(x = "wavelength (nm)", y = "reflectance") + 
    
    # set the y axis range to be consistent between plots
    #ylim(0,as.numeric(max(na.omit(refl_tidy$mean_plus_sd)))) + 
    ylim(0,y_max) + 
    
    # main plot title  
    ggtitle(paste0("Mean Hyperspectral reflectance per species: ", 
                   tools::file_path_sans_ext(str_split(basename(reflFilename),"ALL_")[[1]][2]),
                   " \n",
                   # std dev shading
                   "(shading shows one standard deviation from mean refl range per wavelength)"))
  
  # write plot to file 
  ggsave(paste0(out_dir,"/figures/","separate_ribbon_plot_", 
                tools::file_path_sans_ext(basename(reflFilename)), ".png"), 
         width = 10, height = 6)
  
  
}


# JM distance -------------------------------------------------------------


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


jmDist <- function(reflFilename){
  # classes: list of class names. can be numeric or factors. 
  # x:  each row is a sample described by features (each column)
  # mean and covariance will be computed for each unique class
  # using the spectra within x. 
  
  
  # define the "bad bands" wavelength ranges in nanometers, where atmospheric 
  # absorption creates unreliable reflectance values. 
  bad_band_window_1 <- c(1340, 1445)
  bad_band_window_2 <- c(1790, 1955)
  taxonList <- c("ABLAL","PICOL","PIEN","PIFL2")
  # remove the bad bands 
  remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                                 wavelengths < bad_band_window_1[2]) | 
                                (wavelengths > bad_band_window_2[1] & 
                                   wavelengths < bad_band_window_2[2])]
  wavelength_lut <- data.frame(wavelength = wavelengths,
                               xwavelength = paste0("X", 
                                                    as.character(round(wavelengths))),
                               stringsAsFactors = FALSE) %>% 
    filter(!wavelength %in% remove_bands)
  
  # remove columns that contain "X" in their name but are not reflectance values 
  df <- as.data.frame(read.csv(reflFilename)) %>% 
            dplyr::select(-c(X.1,X,Y)) %>% # remove extra columns that start with X 
              dplyr::filter(taxonID %in% taxonList)# keep only the species of interest
  
  x <- df %>% dplyr::select(wavelength_lut$xwavelength) # remove bad bands 
  
  # get class labels for the spectra
  classes <- df$taxonID %>% droplevels()
  
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
  
}
