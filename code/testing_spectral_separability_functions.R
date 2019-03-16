

# calculate JM distance - spectral separability ---------------------------

for(i in 1:nrow(shapefileLayerNames)){
  extracted_features_filename <- paste0(out_dir, site_code, "_spectral_reflectance_ALL_",
                                        shapefileLayerNames$description[i],".csv")
  
  distances <- jmDist(extracted_features_filename)
  print(distances)
  
  # Error in solve.default(sigma, m) : 
  # system is computationally singular: reciprocal condition number = 4.4099e-22
  
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


# R MARKDOWN CODE CHUNKS BELOW 

## Spectral seperability analysis

The [separability](https://www.rdocumentation.org/packages/spatialEco/versions/1.1-0/topics/separability) function within the SpatialEco R package calculates a variety of two-class sample separability metrics:
  
  * *B* Bhattacharryya distance statistic (Bhattacharyya 1943; Harold 2003) - Measures the similarity of two discrete or continuous probability distributions.

* *JM* Jeffries-Matusita distance statistic  (Bruzzone et al., 2005; Swain et al., 1971) - The J-M distance is a function of separability that directly relates to the probability of how good a resultant classification will be. 

* *M* M-Statistic (Kaufman & Remer 1994) - This is a measure of the difference of the distributional peaks. A large M-statistic indicates good separation between the two classes as within-class variance is minimized and between-class variance maximized (M <1 poor, M >1 good).

* *D* Divergence index and *TD* Transformed Divergence index (Du et al., 2004) - Maximum likelihood approach. Transformed divergence gives an exponentially decreasing weight to increasing distances between the classes.




```{r results = TRUE}
# install.packages("spatialEco")
# library(spatialEco)
x <- as.numeric(mean_reflectance[1,2:length(wavelengths)]) # first row, ABLAB
y <- as.numeric(mean_reflectance[2,2:length(wavelengths)]) # second row, PICOL
sep <- spatialEco::separability(x, y, plot = TRUE)

# calculate seperability between each pair of species, for each polygon size 


# visualize using heat map https://www.r-graph-gallery.com/heatmap/ 

```

```{r seperability}
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
```
