## Community-occupancy-model-amphibians-atlantic-forest-streams
__________________________________________________________________________________________________________________________________________
# Effects of agriculture and topography on tropical amphibian species and communities

### José Wagner Ribeiro Jr, Tadeu Siqueira, Gabriel Lourenço Brejão, and Elise F. Zipkin

Code DOI:  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1218018.svg)](https://doi.org/10.5281/zenodo.1218018)

Please contact the first author for questions about the code or data: José W. Ribeiro Jr (jwribeirojunior@gmail.com)
__________________________________________________________________________________________________________________________________________
# Abstract:
Habitat loss is the greatest threat to the persistence of forest-dependent amphibians, but it is not the only factor influencing species occurrences. The composition of the surrounding matrix, structure of stream networks, and topography are also important landscape characteristics influencing amphibian distributions. Tropical forests have high diversity and endemism of amphibians, but little is known about the specific responses of many of these species to landscape features. In this paper, we quantify the response of amphibian species and communities to landscape-scale characteristics in streams within the fragmented Brazilian Atlantic Forest. We surveyed amphibian communities during a rainy season in 50 independent stream segments using Standardized Acoustic and Visual Transect Sampling (active) and Automated Acoustic Recorders (passive) methods. We developed a hierarchical multi-species occupancy model to quantify the influence of landscape-scale characteristics (forest cover, agriculture, catchment area, stream density and slope) on amphibian occurrence probabilities while accounting for imperfect detection of species using the two survey methods. At the community level, we estimated an overall mean positive relationship between amphibian occurrence probabilities and forest cover, and a negative relationship with agriculture. Catchment area and slope were negatively related with amphibian community structure (95% credible interval [CI] did not overlap zero). The species-level relationships with landscape covariates were highly variable but showed similar patterns to those at the community-level. Species detection probabilities varied widely and were influenced by the sampling method. For most species, the active method resulted in higher detection probabilities than the passive approach. Our findings suggest that small streams and flat topography lead to higher amphibian occurrence probabilities for many species in Brazil’s Atlantic Forest. Our results combined with land use and topographic maps can be used to make predictions of amphibian occurrences and distributions beyond our study area. Such projections can be useful to determine where to conduct future research and prioritize conservation efforts in human-modified landscapes.

Key Words: Anuran, Atlantic Forest, automated acoustic recorders, community model, deforestation, detection error, habitat loss, hierarchical model, landscape, matrix habitat, tropical forest

# Data
## Occupancy data

amphibian_occ_data.csv - contains the data of amphibian occurrence in Brazilian Atlantic Forest streams. The rows are 50 sites in the study area. There are four columns: 1) "stream" - contains sampling site id; 2) "date" - contains the survey date (YYYYMMDD); 3) "species" - is species name code. The names of species with their respectives codes are available in the "amphibian_species.csv"; 4) "rep" - representes a replicate survey event. Rep 1-3 represents the passive method using Automated Acoustic Recorders (AAR) and rep 4-5 represents the active approach of Standardized Acoustic and Visual Transect Sampling (SAVTS).

## Habitat covariates
occupancy_habitat_covariates_anura.csv - contains habitat covariates information for each site. 1) "stream" - contains sampling site id;	2) "long" - is the geographic longitude coordinate as decimal degree; 3)	"lat" - is the geographic latitude coordinate as decimal degree; 4)	"altitude" of each site;	5) "agriculture" -  is the proportion of agriculture area within a buffer; 6) "buildings" - is the proportion of rural buildings area within a buffer ; 7) "forest" - is the proportion of natural forest cover area within a buffer ; 8) "silviculture" - is the proportion of Pinus/Eucalyptus area within a buffer; 9) "lentic_water" - is the proportion of lentic bodywater area within a buffer; 10 "catchment_area" -  is the complete surface area (ha) that contributes to the stream channel in the downstream point from each sampling site ; 11)	"stream_length" - is the stream length network (m) within a buffer; 12)	"slope_mean - " is the mean slope within a buffer, it was derived from the Digital Elevation Model raster image (30-m resolution) from Shuttle Radar Topography Mission (SRTM).

## Detection covariates
detection_covariates_anura.csv - contains information of survey date and daily precipitation for each replicate survey event: 1) stream" - contains sampling site id; 2-6) Each columns contains Julian date for a different replicate survey event (we assumed the first day as the beginning of southern hemisphere spring). ; 7-11) Each columns contains daily precipitation (mm) for a replicate survey event. 

# Code
community_model_code_amphibians.R - R code to run the hierarchical community occupancy model for amphibian in Brazil’s Atlantic Forest streams. Contains code to import and reshape the data and run the model file in JAGS.
