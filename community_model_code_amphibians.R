##############################################################/############
########   Amphibians of Brazilian Atlantic Rainforest streams #########
##############   Hierarchical multi-species model   ####################
###########################################################################

rm(list=ls(all=TRUE)) # clear memory 
choose.dir(default = "", caption = "Select folder") # choose a folder
getwd() # Check the filepath for the working directory

# Read amphibian occurence data of Brazil's Atlantic Forest
data <- read.csv("amphibian_occ_data.csv", header=T, sep=",", na.strings=TRUE)
data$Occ <- rep(1, dim(data)[1])  # Adding a column with 1s on occurrence data

# See the first ten lines of data
data[1:10,]

# How many times each species was observed
(total.count = tapply(data$Occ, data$Species, sum))

# The number of unique species
(uspecies = as.character(unique(data$Species)))

# nspec is the number of observed species
(nspec = length(uspecies))

# The number of unique sampling locations
(upoints = as.character(unique(data$Stream)))

# nsites is the number of sampled streams
(nsites = length(upoints))

# Reshape the data using the R package "reshape"
library(reshape)

# The detection/non-detection data is reshaped into a three dimensional 
# array X where the first dimension = n.site (streams); the second 
# dimension = nrep (replicates); and the last dimension = nspec (species). 
junk.melt=melt(data,id.var=c("Species", "Stream", "Rep"), measure.var="Occ")
X=cast(junk.melt, Stream ~ Rep ~ Species) # aggregatin the data in a 3D array.

# Add in the missing lines with NAs
for (i in 1: dim(X)[3]) {
  b = which(X[,,i] > 0)
  X[,,i][b] = 1                          # put 1 (presence) if the species were surveyed at least 1 time in the sample
  X[,,i][-b] = 0                         # put 0 if the species was not surveyed
  X[,,i][29,5] = NA;  X[,,i][39,5] = NA; # adding NA's for missing sampling occasions
} # i

# nrep is the number of replicates
(nrep <- dim (X)[2])

# Create all zero encounter histories to add to the detection array X as part of the
# data augmentation to account for additional species (beyond the n observed species). 

# nzeroes is the number of all zero encounter histories to be added
nzeroes = 120

# X.zero is a matrix of zeroes, including the NAs for when a point has not been sampled  
X.zero = matrix(0, nrow = nsites, ncol = nrep)
X.zero[29,5] = NA;  X.zero[39,5] = NA

# Xaug is the augmented version of X.  The first n species were actually observed and the n+1 
# through nzeroes species are all zero encounter histories  
Xaug <- array(0, dim=c(dim(X)[1], dim(X)[2], dim(X)[3]+nzeroes))
Xaug[,,(dim(X)[3]+1):dim(Xaug)[3]] = rep(X.zero, nzeroes)
dimnames(X)=NULL
Xaug[,,1:dim(X)[3]] <-  X

# K is a vector of length nsites indicating the number of reps at each site j  
KK <- X.zero
a=which(KK==0); KK[a] <- 1
K=apply(KK,1,sum, na.rm=TRUE)
K=as.vector(K)

# Create a vector to indicate the methodology type (active (SAVTS) = 1; passive (AAR) = 0)
(Met <- c(0,0,0,1,1))

# Read in the habitat data      
habitat <- read.csv("occupancy_habitat_covariates_anura.csv", header=TRUE, sep=",", na.strings=c("NA"))
head(habitat)

# Standardize the natural forest cover
(forest <- as.vector(habitat$forest))
mforest <- mean(forest, na.rm=TRUE)
sdforest <- sd(forest, na.rm=TRUE)
forest1 <- as.vector((forest-mforest) / sdforest)

# Standardize the agriculture cover
(agriculture <- as.vector(habitat$agriculture))
magriculture <- mean(agriculture, na.rm=TRUE)
sdagriculture <- sd(agriculture, na.rm=TRUE)
agriculture1 <- as.vector((agriculture-magriculture) / sdagriculture)

# Standardize the catchment area
(catchment <- as.vector(habitat$catchment_area))
mcatchment <- mean(catchment, na.rm=TRUE)
sdcatchment<- sd(catchment, na.rm=TRUE)
catchment1 <- as.vector((catchment-mcatchment) / sdcatchment)

# Standardize the stream density (the area of buffer is 12.56 ha)
(density <- as.vector(habitat$stream_length) / 12.56)
mdensity <- mean(density, na.rm=TRUE)
sddensity <- sd(density, na.rm=TRUE)
density1 <- as.vector((density-mdensity) / sddensity)

# Standardize the slope
slope <- as.vector(habitat$slope_mean)
mslope <- mean(slope, na.rm=TRUE)
sdslope <- sd(slope, na.rm=TRUE)
slope1 <- as.vector((slope-mslope) / sdslope)

# Read in the detection data - The sampling dates were converted to Julien dates 
# We assumed the first day as the beginning of southern hemisphere spring (09/23/2015) 
detec <- read.csv("detection_covariates_anura.csv", header=TRUE, sep=",", na.strings=c("NA"))
head(detec)

# Putting the julian date and daily precipitation in a matrix
dates <- as.matrix(detec[,c("date1","date2","date3","date4", "date5")])
rain  <- as.matrix(detec[,c("rain1","rain2","rain3","rain4", "rain5")])

# Standardize the julian date
mdate <- mean(dates, na.rm=TRUE)
sddate <- sqrt(var(dates[1:length(dates)], na.rm=TRUE))
date1 <- (dates-mdate) /  sddate
date2 <- date1*date1     # Quadratic effect of Julian date 

# Standardize daily precipitation
mrain <- mean(rain, na.rm=TRUE)
sdrain <- sqrt(var(rain[1:length(rain)], na.rm=TRUE))
rain1 <- (rain-mrain) /  sdrain

# Adding 0s (mean) to missing data (NA's) - Jags doesn't work with NA's in covariates data  
date1[29,5] = 0;  date1[39,5] = 0;
date1 <- as.matrix(date1)

date2[29,5] = 0;  date2[39,5] = 0;
date2 <- as.matrix(date2)

rain1[29,5] = 0;  rain1[39,5] = 0;
rain1 <- as.matrix(rain1)

