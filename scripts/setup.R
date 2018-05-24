### SETUP ####
# Description: this script loads and installs packages, sources the functions script, and reads input data. This script is not meant to modify or pre-process the input data in any way; this is done by the inputs.R script

# Setup: requires project in root folder with:
# all input data in /inputs directory, 
# all scripts (e.g., functions.R) in /scripts directory, 
# and /outputs directory for optional plots/tables

source('scripts/functions.R')

# Install and load packages
# Tidyverse packages and data.table for data processing, mgcv (GAMs), rgdal (spatial data), doParallel (parallel computing), rstan and coda for Bayesian inference, crayon for colored console outputs (helpful but unnecessary), scales/e1071/MASS/quantreg for miscallaneous functions
ipak(c("tidyverse", "data.table", "car", "psych", "doParallel", "coda", "rstan", "scales", "e1071", "MASS", "crayon", "rgdal", "networkD3", "data.tree", "plotly", "ggridges")) # "mgcv", "quantreg"

options(stringsAsFactors = FALSE)
options(scipen=999)
rstan_options(auto_write = TRUE)

alpha <- scales::alpha
extract <- rstan::extract
select <- dplyr::select # dplyr conflicts with MASS
rescale <- scales::rescale

# options(mc.cores = parallel::detectCores())

### > READ DATA ####
# Most datasets are read into list 'inputs'
inputs = list()
# Invertebrates ####
inv <- read.csv('inputs/invertebrates_2017-04-20.dat', sep='\t', header=TRUE, na.strings="", stringsAsFactors=FALSE)
inputs$taxonomy <- read.csv('inputs/invertebrates_taxonomy_2018-02-23.dat', sep='\t', header=TRUE, stringsAsFactors=FALSE, na.strings=c(""," ","NA"))

# Coordinates for sites and borders of Switzerland
inputs$ch <- readOGR("./inputs", "switzerland", stringsAsFactors = F)
# inputs$ch <- read.csv('inputs/Swiss_border_coordinates.dat', sep='\t', header=TRUE, stringsAsFactors = FALSE)
inputs$xy <- select(inv, SiteId, X, Y)
inputs$xy <- distinct(inputs$xy)

# Env. conditions ####
env <- read.csv('inputs/environmental_data_2017-08-18.dat', sep='\t', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)

# Land use ####
inputs$mp <- read.csv('inputs/site_MP_Strahm_TU.dat', sep='\t', header=TRUE, stringsAsFactors=FALSE)
inputs$land_metadata <- read.csv('inputs/landuse_metadata.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)

# Additional data ####
# idw = list()
# idw$edo_farm <- read.csv('inputs/site_IDW_EDO_1.0_farmland.csv', sep=',', header=TRUE, na.strings= '--', stringsAsFactors=FALSE)
# idw$flo_farm <- read.csv('inputs/site_IDW_FLO_1.0_farmland.csv', sep=',', header=TRUE, na.strings= '--', stringsAsFactors=FALSE)
# idw$edo_forest <- read.csv('inputs/site_IDW_EDO_1.0_forest.csv', sep=',', header=TRUE, na.strings= '--', stringsAsFactors=FALSE)
# idw$flo_forest <- read.csv('inputs/site_IDW_FLO_1.0_forest.csv', sep=',', header=TRUE, na.strings= '--', stringsAsFactors=FALSE)
# idw$ha_flo_forest <- read.csv('inputs/site_IDW_HA-FLO_1.0_forest.csv', sep=',', header=TRUE, na.strings= '--', stringsAsFactors=FALSE)
# idw$ha_flo_farm <- read.csv('inputs/site_IDW_HA-FLO_1.0_farmland.csv', sep=',', header=TRUE, na.strings= '--', stringsAsFactors=FALSE)
as09.metadata <- read.csv('inputs/ArealStatistics_metadata.csv', sep = ',', header = TRUE, stringsAsFactors = F, na.strings = '-')

# inputs$fri <- read.csv('inputs/site_inv_FRI_20171505.dat', sep='\t', header=TRUE, stringsAsFactors = FALSE)
# inputs$bfri <- read.csv('inputs/site_inv_bFRI_20171505.dat', sep='\t', header=TRUE, stringsAsFactors = FALSE)

# srb <- read.csv('inputs/site_RiverBasins.csv', sep = ',', header = TRUE, na.strings = "", stringsAsFactors = F)
# srb$OBJECTID <- NULL
# srb$RiverRhone <- ifelse(srb$RiverBasin == 'Rhone', 1, 0)
# srb$RiverRhein <- ifelse(srb$RiverBasin == 'Rhein', 1, 0)
# srb$RiverInn <- ifelse(srb$RiverBasin == 'Inn', 1, 0)
# srb$RiverPo <- ifelse(srb$RiverBasin == 'Po', 1, 0)

# chem = list()
# chem$chem.m <- read.csv('inputs/chemistry_alldata_2009_2015_2016-06-09.dat', sep='\t', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)
# chem$chem.sd <- read.csv('inputs/Site_chemistry_20161117.csv', sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)
# chem$chem.sites <- readOGR("/inputs", "cc_join", stringsAsFactors = FALSE)
# chem$chem.discharge <- read.csv('inputs/site_chem_Discharge.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
# inputs$fs <- read.csv('inputs/site_FF.csv', sep=',', header=TRUE, stringsAsFactors = FALSE, na.strings = c("N/A"))
# inputs$nutrients <- fread('inputs/site_NP.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)
# sc <- fread('site_Connectivity.txt', sep='\t', header=TRUE, stringsAsFactors = FALSE)
# inputs$inv.pref <- fread('inputs/Invertebrates_functionalGroups.csv', header = TRUE, sep = ';', stringsAsFactors = FALSE)
# inputs$golf <- read.csv('inputs/site_Golf.csv', sep=',', header=TRUE, stringsAsFactors = FALSE)
# inputs$dc <- fread('inputs/site_Discharge_20171505.dat', sep='\t', header=TRUE, stringsAsFactors = FALSE)
# Spatial data (points/polygons)
# sites <- readOGR("/inputs", "SiteId_EZGNR", stringsAsFactors = FALSE) # distinct sites located in sub-catchments
# bg <- readOGR("/inputs", "biogeography", stringsAsFactors = FALSE)

# Land use by buffer zones ####
b <- list()
b$areas <- fread("inputs/catchment_buffer/Area_All_Buffer.csv", sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)

b$forest.sub <- fread("inputs/catchment_buffer/Forest_All_Buffer.csv", sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)
b$urban.sub <- fread("inputs/catchment_buffer/Urban_All_Buffer.csv", sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)

b$forest <- fread("inputs/catchment_buffer/Forest_All_Buffer_Catch.csv", sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)
b$ftype <- fread("inputs/catchment_buffer/FTyp_All_Buffer_Catch.csv", sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)
b$grass <- fread("inputs/catchment_buffer/Grass_All_Buffer_Catch.csv", sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)
b$urban <- fread("inputs/catchment_buffer/Urban_All_Buffer_Catch.csv", sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)
b$arable <- fread("inputs/catchment_buffer/Arable_All_Buffer_Catch.csv", sep=';', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)

