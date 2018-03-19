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

# Write the model code to a text file 
cat("
    model{
    
    #Priors
    omega ~ dunif(0,1)
    
    # Hyperpriors - define prior distributions for community-level parameters 
    # in occupancy and detection models
    a0.mean ~ dunif (0,1)
    mu.a0 <- log(a0.mean) - log(1-a0.mean)
    tau.a0 ~ dgamma(0.1,0.1) 
    
    b0.mean ~ dunif(0,1) ##transect
    mu.b0 <- log(b0.mean) - log(1-b0.mean)
    tau.b0 ~ dgamma(0.1,0.1)    
    
    b00.mean ~ dunif(0,1)
    mu.b00 <- log(b00.mean) - log(1-b00.mean)
    tau.b00 ~ dgamma(0.1,0.1) 
    
    mu.a1 ~ dnorm(0, 0.001)
    mu.a2 ~ dnorm(0, 0.001)
    mu.a3 ~ dnorm(0, 0.001)
    mu.a4 ~ dnorm(0, 0.001)
    mu.a5 ~ dnorm(0, 0.001)
    mu.b1 ~ dnorm(0, 0.001)
    mu.b2 ~ dnorm(0, 0.001)
    mu.b3 ~ dnorm(0, 0.001)
    
    tau.a1 ~ dgamma(0.1,0.1)
    tau.a2 ~ dgamma(0.1,0.1)
    tau.a3 ~ dgamma(0.1,0.1)
    tau.a4 ~ dgamma(0.1,0.1) 
    tau.a5 ~ dgamma(0.1,0.1)
    tau.b1 ~ dgamma(0.1,0.1) 
    tau.b2 ~ dgamma(0.1,0.1)
    tau.b3 ~ dgamma(0.1,0.1)
    
    # Create priors for species i from the community-level prior distributions
    for (i in 1:(nspec+nzeroes)) {
    w[i]  ~ dbern(omega)
    a0[i] ~ dnorm(mu.a0,  tau.a0)
    b00[i]~ dnorm(mu.b00, tau.b00) 
    b0[i] ~ dnorm(mu.b0,  tau.b0)
    a1[i] ~ dnorm(mu.a1,  tau.a1)
    a2[i] ~ dnorm(mu.a2,  tau.a2)
    a3[i] ~ dnorm(mu.a3,  tau.a3)
    a4[i] ~ dnorm(mu.a4,  tau.a4)
    a5[i] ~ dnorm(mu.a5,  tau.a5)
    b1[i] ~ dnorm(mu.b1,  tau.b1)    
    b2[i] ~ dnorm(mu.b2,  tau.b2)
    b3[i] ~ dnorm(mu.b3,  tau.b3)
    
    # Estimate the Z matrix (true occurrence for species i at stream j)      
    for (j in 1:nsites) {
    logit(psi[j,i]) <- a0[i] + a1[i] * forest1[j] + a2[i] * agriculture1[j] +
    a3[i] * catchment1[j] + a4[i] * density1[j] + a5[i] * slope1[j]
    
    mu.psi[j,i] <- psi[j,i] * w[i]
    Z[j,i] ~ dbern(mu.psi[j,i])
    
    # Estimate detection for species i at stream j during sampling period k.      
    for (k in 1:K[j]) {  
    logit(p[j,k,i]) <-   b0[i] * Met[k] + b00[i] * (1-Met[k]) + 
    b1[i] * date1[j,k] +  b2[i] * date2[j,k] + b3[i] * rain1[j,k]
    
    mu.p[j,k,i] <- p[j,k,i]*Z[j,i]
    X[j,k,i] ~ dbern(mu.p[j,k,i])
    } # k
    } # j
    } # i
    
    # Derived quantities
    n0 <- sum(w[(nspec+1):(nspec+nzeroes)]) # number of unseen species
    N <- nspec + n0                         # Overall estimated richness
    
    # Determine site level richness estimates for the whole community
    for(j in 1:nsites){
    Nsite[j] <- inprod(Z[j,1:(nspec+nzeroes)],w[1:(nspec+nzeroes)])
    } # j
    
    # Finish writing the text file into a document we call community.model.anura.txt
    }
    ",file="community.model.anura.txt")

# Load all the data
jags.data = list (nspec = nspec, nzeroes = nzeroes, nsites = nsites, K = K, 
                  X = Xaug, forest1 = forest1, agriculture1 = agriculture1, 
                  catchment1 = catchment1, density1 = density1, slope1=slope1, 
                  Met = Met, date1 = date1, date2 = date2, rain1 = rain1)

# Specify the parameters to be monitored
jags.params = c("omega", "mu.a0", "mu.a1", "mu.a2", "mu.a3", "mu.a4", "mu.a5",
                "mu.b0", "mu.b00", "mu.b1", "mu.b2", "mu.b3", "tau.a0", "tau.a1",  
                "tau.a2","tau.a3","tau.a4","tau.a5", "tau.b0", "tau.b00", 
                "tau.b1", "tau.b2", "tau.b3", "a0", "a1", "a2","a3", "a4", "a5",
                "b0", "b00", "b1", "b2", "b3", "N", "Nsite")

# Specify the initial values
zinits <- apply(Xaug,c(1,3),max,na.rm=TRUE)
jags.inits = function (){
  omegaGuess = runif(1, nspec/(nspec+nzeroes), 1)
  psi.meanGuess = runif(1, 0.25,1)
  list(omega=omegaGuess, w=c(rep(1, nspec), rbinom(nzeroes, size=1, prob=omegaGuess)),
       a0=rnorm(nspec+nzeroes), b0=rnorm(nspec+nzeroes),b00=rnorm(nspec+nzeroes),
       Z = zinits)  
}

# MCMC settings
ni <- 50000 # number of total iterations per chain
nt <-    20 # thinning rate 
nb <- 30000 # number of iterations to discard at the beginning
nc <-     3 # number of Markov chains
na <- 10000 # Number of iterations to run in the JAGS adaptive phase.

# Load the jagsUI library
library(jagsUI)

# Run JAGS model
fit <- jagsUI(data = jags.data, inits = jags.inits, jags.params,
              "covar.model.anura.BAF.txt", n.chains = nc, n.thin = nt,
              n.iter = ni, n.burnin = nb, n.adapt = na, parallel=T, store.data=T)