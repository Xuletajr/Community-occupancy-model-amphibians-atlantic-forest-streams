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