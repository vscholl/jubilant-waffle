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



# read the tree stem locations file
tree_points <- rgdal::readOGR(dsn = shapefile_dir,
                              layer = "mapped_stems_final")

# read tree polygons file - circles with diameter = max crown diameter
tree_polygons <- rgdal::readOGR(dsn = shapefile_dir,
                                layer = "polygons_clipped_overlap")


# 50% crown diameter ------------------------------------------------------

# read tree polygons file - circles with diameter = half of max crown diameter.
# note that there are a different number of final polygons / stems when the
# diameter is changed; this is also a function of the area threshold. 
tree_polygons <- rgdal::readOGR(dsn = shapefile_dir,
                                         layer = "polygons_clipped_overlap_50percent")


tree_points <- rgdal::readOGR(dsn = shapefile_dir,
                              layer = "mapped_stems_final_50percent")

out_dir <- paste0('../output/', site_code, '/diam50percent/')
check_create_dir(out_dir) # create output folder for site

# ------------------------------------------------------------------------



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
  
  # each current hyperspectral tile must be read and stacked into a 
  # georeferenced rasterstack object (so it can be clipped with point / polygon
  # shapefiles). The process of creating a rasterstack takes a while for 
  # each tile, so after creating each rasterstack once, each object gets 
  # written to a file. 
  
  # Build up the rasterstack filename by parsing out the easting/northing
  # coordinates from the current h5 filename.
  rasterstack_filename <- paste0(h5_dir, "rasterstack_",
              str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][5],"_",
              str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][6],".rds")
  
  print(paste("rasterstack filename: ", rasterstack_filename))
  
  # check to see if a .rds file already exists for the current tile.
  if (file.exists(rasterstack_filename)){
    
    # if it exists, read that instead of re-generating the same rasterstack.
    message("rasterstack was already created for current tile")
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
  
  # clipped_overlap polygons generated using the neon_veg workflow 
  extracted_polygon_spectra <- raster::extract(s, polygons_in_sp, df = TRUE)
  
  # the buffer parameter can be used to include cells around each point of a
  # given size. the buffer parameter can be specified as a vector of the length
  # of the number of points. 
  
  # maxCrownDiameter (buffer of (maxCrownDiameter / 2))
#  buffers_mxDm <- tree_polygons_points$crownDm / 2
#  extracted_spectra_buffer_mxDm <- raster::extract(s, 
#                                                   points_in_sp,
#                                                   buffer = buffers_mxDm,
#                                                   df = TRUE)
  
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
  
  # 50% max crown diameter (buffer of (maxCrownDiameter / 4))
#  extracted_polygon_halfDiam_spectra <- raster::extract(s, polygons_halfDiam_in_sp, df = TRUE)
  
  
  
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
  
  # clipped_overlap polygons generated using the neon_veg workflow 
 # write_spectra_to_file(spectra = as.data.frame(extracted_polygon_spectra),
 #                       polygons_in = polygons_in,
 #                       filename_out = paste0(out_dir,
 #                                             "spectral_reflectance_",
 #                                             as.character(s@extent[1]),"_",
 #                                             as.character(s@extent[3]),"_",
 #                                             "polygons_clipped_overlap_max_diameter",
 #                                             ".csv"))
  
  # maxCrownDiameter (buffer of (maxCrownDiameter / 2))
#  write_spectra_to_file(spectra = as.data.frame(extracted_spectra_buffer_mxDm),
#                        polygons_in = polygons_in,
#                        filename_out = paste0(out_dir,
#                                              "spectral_reflectance_",
#                                              as.character(s@extent[1]),"_",
#                                              as.character(s@extent[3]),"_",
#                                              "buffer_max_diameter",
#                                              ".csv"))

  # clipped_overlap polygons generated using the neon_veg workflow,
  # 50% maxCrownDiameter size. 
  write_spectra_to_file(spectra = as.data.frame(extracted_polygon_spectra),
                        polygons_in = polygons_in,
                        filename_out = paste0(out_dir,
                        "spectral_reflectance_",
                        as.character(s@extent[1]),"_",
                        as.character(s@extent[3]),"_",
                        "polygons_clipped_overlap",
                        ".csv"))
  
  
}
  



# read and plot the spectra .csv files  -----------------------------------

# get a list of the .csv file per tile containing woody veg stemse
csvs <- list.files(out_dir, full.names = TRUE)

# csvs <- csvs[grepl("*000.csv", csvs)] # from back when there was only one collection of csvs

# specify a description that the different shapefile iterations are named by
out_description <- "stem_points" # stem point locations 
out_description <- "polygons_clipped_overlap" # clipped_overlap polygons

#out_description <- "polygons_checked_overlap_max_diameter" # checked_overlap polygons
#out_description <- "buffer_max_diameter" # buffer of (maxCrownDiameter / 2)
#out_description <- "polygons_checked_overlap_50percent_diameter" # checked_overlap, 50% max crown diameter

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

## write the exact wavelengths to file for future use --> moved into stack_hyperspectral function
#write.table(data.frame(wavelengths = wavelengths),
#            paste0(out_dir,"wavelengths.txt"),
#            sep="\n",
#            row.names=FALSE)

# read wavelengths if not previously created
#wavelengths = read.table(paste0(out_dir,"wavelengths.txt"),
#                         sep="\n",
#                         skip = 1,
#                         col.names = 'wavelength')

# read wavelengths if not previously created
wavelengths = as.numeric(unlist(read.table(paste0(out_dir,"wavelengths.txt"),
                                           sep="\n",
                                           skip = 1,
                                           col.names = 'wavelength')))




# ribbon plots for different shapefiles ------------------------------------------

# instead put each set of spectra into its own folder
# use folder name for graph titles
shapefile_list <- c("stem_points_maxdiameter"
                    ,"polygons_clipped_overlap_max_diameter"
                    ,"polygons_clipped_overlap_50percent"
                    )

# absolute maximum reflectance to set the same ylimit for the plots
y_max <- 0.35    #max(refl_tidy$max_reflectance, na.rm = TRUE)

# remove the bad bands 
remove_bands <- wavelengths[(wavelengths > bad_band_window_1[1] & 
                               wavelengths < bad_band_window_1[2]) | 
                              (wavelengths > bad_band_window_2[1] & 
                                 wavelengths < bad_band_window_2[2])]

# figure out which individual ID's appear across both data sets
# get ID's in max diameter set 
spectra_maxDiameter <- as.data.frame(read.csv(
  paste0(out_dir, "ribbon_plots_spectral_reflectance_ALL/",
         site_code, "_spectral_reflectance_ALL_", "stem_points_maxdiameter", ".csv"))) 
# get ID's in 50 percent max diameter set 
spectra_50percent <- as.data.frame(read.csv(
  paste0(out_dir, "ribbon_plots_spectral_reflectance_ALL/",
         site_code, "_spectral_reflectance_ALL_", "stem_points_50percent", ".csv"))) 

# get the ID's that are the same between the two data sets 
stems_to_use <- intersect(as.character(spectra_maxDiameter$individualID), 
                          as.character(spectra_50percent$individualID))

# create empty data frame to contain the mean spectra for each 
# species, for each shapefile data set 
df_mean_spectra <- data.frame(shpfile = c("stem", "stem","stem","stem",
"poly_maxDiameter","poly_maxDiameter","poly_maxDiameter","poly_maxDiameter",
"poly_50percent","poly_50percent","poly_50percent","poly_50percent")
                            species = c())


for (out_description in shapefile_list){
  print(paste0("creating ribbon plot for ", out_description))
  
  fname <- paste0(out_dir, "ribbon_plots_spectral_reflectance_ALL/",
                  site_code, "_spectral_reflectance_ALL_", out_description, ".csv")
  spectra_all <- as.data.frame(read.csv(fname)) %>%  select(-X.1)
  
  
  # filter down the spectra to utilize only the 
  # trees that are present in both the max crown diameter data set 
  spectra_all <- spectra_all %>% dplyr::filter(individualID %in% stems_to_use)
  
  print(paste("There are", as.character(length(unique(spectra_all$individualID))),
               "unique spectra being plotted"))
  
  # calculate mean reflectance per species
  # VS-NOTE: ADJUST THE COLUMN SELECTION to use names instead of indices (i.e. "10")
  mean_reflectance <- stats::aggregate(spectra_all[,10:ncol(spectra_all)], 
                                       by = list(taxonID = spectra_all$taxonID),
                                       FUN = mean) 
  min_reflectance <- stats::aggregate(spectra_all[,10:ncol(spectra_all)], 
                                      by = list(taxonID = spectra_all$taxonID),
                                      FUN = min) 
  max_reflectance <- stats::aggregate(spectra_all[,10:ncol(spectra_all)], 
                                      by = list(taxonID = spectra_all$taxonID),
                                      FUN = max) 
  sd_reflectance <- stats::aggregate(spectra_all[,10:ncol(spectra_all)],
                                         by = list(taxonID = spectra_all$taxonID),
                                         FUN = sd)
  
  # add and subtract 1 standard deviation from the mean. Keep taxon ID column.
  # mean_plus_sd <- cbind(taxonID = mean_reflectance$taxonID,
  #                       (mean_reflectance[2:ncol(mean_reflectance)] + 
  #                        sd_reflectance[2:ncol(mean_reflectance)]))
  # 
  # mean_minus_sd <- cbind(taxonID = mean_reflectance$taxonID,
  #                        (mean_reflectance[2:ncol(mean_reflectance)] - 
  #                        sd_reflectance[2:ncol(mean_reflectance)]))
  
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
  
  # specify the colors for the reflectance curves & shading around them 
  shading_colors <- c("#d7191c", "#fdae61", "#abdda4", "#2b83ba")
  species <- sort(unique(spectra_all$taxonID)) #alphabetical so colors match plot above
  shading_alpha <- 0.4
  
  # generate the ribbon plot
  ggplot(refl_tidy, 
         aes(x = wavelength, y = mean_reflectance, color = taxonID)) + 
    
    # shaded ribbon from min to max for each species
    # can't get the shading colors to match the lines
    #geom_ribbon(aes(ymin = min_reflectance,
    #                ymax = max_reflectance,
    #                alpha = 0.1,
    #                fill = taxonID)) + 
    
    # ABLAL
    geom_ribbon(data = refl_tidy[refl_tidy$taxonID == species[1], ],
                #aes(ymin = min_reflectance, ymax = max_reflectance), # min max shading
                #aes(ymin = mean_minus_sd[mean_minus_sd$taxonID == species[1], ], 
                #    ymax = mean_plus_sd[mean_plus_sd$taxonID == species[1], ]), # std dev shading
                aes(ymin = mean_reflectance - sd_reflectance,
                    ymax = mean_reflectance + sd_reflectance), # std dev shading
                colour=NA,
                alpha = shading_alpha,
                fill = shading_colors[1],
                show.legend = F) + 
    
    # PICOL
    geom_ribbon(data = refl_tidy[refl_tidy$taxonID == species[2], ],
                #aes(ymin = min_reflectance, ymax = max_reflectance), # min max shading
                aes(ymin = mean_reflectance - sd_reflectance,
                    ymax = mean_reflectance + sd_reflectance), # std dev shading
                colour=NA,
                alpha = shading_alpha,
                fill = shading_colors[2],
                show.legend = F) + 
    
    # PIEN
    geom_ribbon(data = refl_tidy[refl_tidy$taxonID == species[3], ],
                # aes(ymin = min_reflectance, ymax = max_reflectance), # min max shading
                aes(ymin = mean_reflectance - sd_reflectance,
                    ymax = mean_reflectance + sd_reflectance), # std dev shading
                colour=NA,
                alpha = shading_alpha,
                fill = shading_colors[3],
                show.legend = F) + 
    
    # PIFL2
    geom_ribbon(data = refl_tidy[refl_tidy$taxonID == species[4], ],
                # aes(ymin = min_reflectance, ymax = max_reflectance), # min max shading
                aes(ymin = mean_reflectance - sd_reflectance,
                    ymax = mean_reflectance + sd_reflectance), # std dev shading
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
    ggtitle(paste0("Mean Hyperspectral reflectance per species: ", out_description, " \n",
                   #"(shading shows minimum and maximum refl range per wavelength)")) # min max shading
                   "(shading shows one standard deviation from mean refl range per wavelength)")) # std dev shading
                   
                   
  ggsave(paste0(out_dir,"/figures/","ribbon_plot_", out_description, ".png"), width = 10, height = 6)
  
}

species_table <- data.table("TaxonID" = c("ABLAL",
                                          "PICOL",
                                          "PIEN",
                                          "PIFL2"),
                            "Scientific name" = c("Abies lasiocarpa",
                                                  "Pinus contorta",
                                                  "Picea engelmannii",
                                                  "Pinus flexilis"),
                            "Common Name" = c("Subalpine fir",
                                              "Lodgepole pine",
                                              "Engelmann spruce",
                                              "Limber pine"),
                            "Number of individual trees" = as.data.frame(table(spectra_all$taxonID))$Freq)


knitr::kable(species_table) %>% kable_styling(bootstrap_options = c("striped", "hover"))

# Spectral seperability  --------------------------------------------------

# The [separability](https://www.rdocumentation.org/packages/spatialEco/versions/1.1-0/topics/separability) 
# function within the SpatialEco R package calculates a variety of two-class 
# sample separability metrics:
  
#  * *B* Bhattacharryya distance statistic (Bhattacharyya 1943; Harold 2003) - 
#        Measures the similarity of two discrete or continuous probability distributions.
#   * *JM* Jeffries-Matusita distance statistic (Bruzzone et al., 2005; Swain et al., 1971) - 
#          The J-M distance is a function of separability that directly relates to the 
#          probability of how good a resultant classification will be. 
#   * *M* M-Statistic (Kaufman & Remer 1994) - This is a measure of the difference of 
#          the distributional peaks. A large M-statistic indicates good separation 
#          between the two classes as within-class variance is minimized and 
#          between-class variance maximized (M <1 poor, M >1 good).
#   * *D* Divergence index and *TD* Transformed Divergence index (Du et al., 2004) - 
#         Maximum likelihood approach. Transformed divergence gives an exponentially 
#         decreasing weight to increasing distances between the classes.

# install.packages("spatialEco")
library(spatialEco)
x <- as.numeric(mean_reflectance[1,2:length(wavelengths)]) # first row, ABLAB
y <- as.numeric(mean_reflectance[2,2:length(wavelengths)]) # second row, PICOL
sep <- spatialEco::separability(x, y, plot = TRUE)
print(sep)

# calculate seperability between each pair of species, for each polygon size 


# create a list of each pair of species to compare seperability 
taxon_list <- as.character(mean_reflectance$taxonID)
# each column in taxon_pairs contains the taxon ID's to compare
taxon_pairs <- combn(taxon_list, m = 2) 

# create any empty data frame to populate with seperability metrics
# each row corresponds to a taxon pair in taxon_pairs,
# while each col corresponds to a seperability metric in this order: 
# B, JM, M, mdif, D, TD 
# (since this is the order generated by spatialEco::separability)
seperability_metrics <- data.frame(matrix(NA, 
                                          nrow = ncol(taxon_pairs), 
                                          ncol = ncol(sep)))
# name the columns accordingly to indicate the seperability metric
colnames(seperability_metrics) <- colnames(sep)


print("Taxon pairs to calculate seperability for: ")
print(taxon_pairs)

for(c in 1:ncol(taxon_pairs)){
  print(c)
  
  # print the taxon ID's of the species currently being compared 
  print(taxon_pairs[1,c])
  print(taxon_pairs[2,c])
  
  # get the mean spectral reflectance for species currently being compared.
  # remove the first column with taxonID, since it's not numeric it will 
  # cause an error in the metric calculations if not removed. 
  refl1 <- mean_reflectance %>% filter(taxonID == taxon_pairs[1,c]) %>%
               select(-taxonID) %>% as.numeric()
  refl2 <- mean_reflectance %>% filter(taxonID == taxon_pairs[2,c]) %>%
               select(-taxonID) %>% as.numeric()
  
  # calculate seperability metrics 
  metrics <- spatialEco::separability(refl1, refl2, plot = TRUE)
  print(metrics)
  
  # store seperability metrics in the data frame
  seperability_metrics[c,] <- metrics
}

# visualize using heat map https://www.r-graph-gallery.com/heatmap/ 





# Read high-res RGB data --------------------------------------------------
tile_easting_northing <- "453000_4433000"

rgb_filename <- paste0("../data/NIWO/rgb/2017_NIWO_1_",
                       tile_easting_northing,
                       "_image.tif")
rgb_img <- raster(rgb_filename)
print(paste("ncols of high-res rgb data: ", as.character(rgb_img@ncols)))
print(paste("nrows of high-res rgb data: ", as.character(rgb_img@nrows)))


# Read CHM data 
chm_filename <- paste0("../data/NIWO/chm/NEON_D13_NIWO_DP3_",
                       tile_easting_northing,
                       "_CHM.tif")
chm <- raster(chm_filename)
print(paste("ncols of chm data: ", as.character(chm@ncols)))
print(paste("nrows of chm data: ", as.character(chm@nrows)))

hist(chm, breaks=round(chm@data@max))


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
            paste0(out_dir,"rgb_composite_452000_4432000.tif"), 
            format="GTiff",
            overwrite=TRUE)



# Random Forest Classification --------------------------------------------
# This part of the script might eventually replace the part of the script above,
# so that's why there is a lot of repeated code. 

# set the seed for consistent random generation 
set.seed(14)

# set working directory
setwd("~/github/jubilant-waffle/code/")

# load any local functions in external files 
source("supporting_functions.R")

# code for NEON site 
site_code <- 'NIWO'


# input shapefile scenarios to test 

# directory with shapefiles (tree stem locations and crown polygons)
shapefile_dir <- paste0('../data/', site_code, '/shapefiles/')

# define the layer name for each of the shapefile scenarios. 

# all stem points at NIWO
# "/shapefiles_maxDiameter/mapped_stems_with_crown_diameter.shp"
allStems_layer <- "all_mapped_stems_with_crown_diameter"

# polygons with max diameter, one generated for each stem point at NIWO
# "shapefiles_maxDiameter/polygons_all.shp"
allPolygons_maxDiameter_layer <- "all_polygons_max_diameter"

# polygons with half max diameter, one generated for each stem point at NIWO
# "shapefiles_halfDiameter/polygons_all.shp" 
allPolygons_halfDiameter_layer <- "all_polygons_half_diameter"

# neon_veg workflow polygons generated with max diameter 
# "shapefiles_maxDiameter/polygons_clipped_overlap.shp"
neonvegPolygons_maxDiameter_layer <- "neon-veg_polygons_max_diameter"

# neon_veg workflow polygons generated with half the max diameter diameter 
# "shapefiles_50percentDiameter/polygons_clipped_overlap.shp"
neonvegPolygons_halfDiameter_layer <- "neon-veg_polygons_half_diameter"



# clip all remote sensing data layers using each of the input shapefile scenarios to test.

# hyperspectral imagery directory with .h5 files
# and also .rds rasterstacks, generated earlier to speed up the testing 
h5_dir <- paste0('../data/', site_code, '/hyperspectral/')
# CHM geotiff directory 
chm_dir <- paste0('../data/', site_code, '/chm/')
# slope geotiff  directory 
slope_dir <- paste0('../data/', site_code, '/slope/')
# aspect geotiff directory 
aspect_dir <- paste0('../data/', site_code, '/aspect/')

# define the output directory. If it doesn't exist already, create it.
check_create_dir('../output/') # create top level "output" directory
out_dir <- paste0('../output/', site_code, '/')
check_create_dir(out_dir) # create output folder for site


# read in each of the shapefile layers 
allStems <- rgdal::readOGR(dsn = shapefile_dir,
                           layer = allStems_layer)

allPolygons_maxDiameter <- rgdal::readOGR(dsn = shapefile_dir,
                                          layer = allPolygons_maxDiameter_layer)

allPolygons_halfDiameter <- rgdal::readOGR(dsn = shapefile_dir,
                                           layer = allPolygons_halfDiameter_layer)

neonvegPolygons_maxDiameter <- rgdal::readOGR(dsn = shapefile_dir,
                                              layer = neonvegPolygons_maxDiameter_layer)

neonvegPolygons_halfDiameter <- rgdal::readOGR(dsn = shapefile_dir,
                                               layer = neonvegPolygons_halfDiameter_layer)





# I could instead write a loop that loops through the layer names.
# reads each layer, clips the giant stack of data for each tile. 
shapefileLayerNames <- c("all_mapped_stems_with_crown_diameter",
                         "all_polygons_max_diameter",
                         "all_polygons_half_diameter",
                         "neon-veg_polygons_max_diameter",
                         "neon-veg_polygons_half_diameter")



# make the following into a function to extract training features for each of the shapefiles 

shapefileLayer <- allStems_layer
shapefileLayer <- allPolygons_maxDiameter_layer

# read the shapefile layer 
shp <- rgdal::readOGR(dsn = shapefile_dir,
                      layer = shapefileLayer)

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


# hyperspectral data - list the .h5 files 
h5_list <- list.files(path = h5_dir, full.names = TRUE)
h5_list <- h5_list[grepl("*.h5", h5_list)]

# loop through h5 files 
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
  rasterstack_filename <- paste0(h5_dir, "rasterstack_",
                                 str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][5],"_",
                                 str_split(tail(str_split(h5, "/")[[1]],n=1),"_")[[1]][6],".rds")
  
  print(paste("rasterstack filename: ", rasterstack_filename))
  
  # check to see if a .rds file already exists for the current tile.
  if (file.exists(rasterstack_filename)){
    
    # if it exists, read that instead of re-generating the same rasterstack.
    message("rasterstack was already created for current tile")
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
  
  # figure out which trees are within the current tile by 
  polygons_in <- shp %>% 
    dplyr::filter(xmin >= extent(s)[1] & 
                    xmax < extent(s)[2] & 
                    ymin >= extent(s)[3] & 
                    ymax < extent(s)[4])
  
  print(paste0(as.character(nrow(polygons_in))," polygons in current tile"))
  
  
  
}



# create a list of all the directories within "../data/" to use
# for extracting features 
dataDirectories <- c("hyperspectral",
                            "chm",
                            "slope",
                            "aspect")

# read all remote sensing data products and combine them into a giant stack 

# add a layer with pixel ID's 

# clip using the current shapefile, write these features to file 





