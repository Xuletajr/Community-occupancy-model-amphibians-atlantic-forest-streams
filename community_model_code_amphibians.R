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
