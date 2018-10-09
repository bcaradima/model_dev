### INPUTS ####
# Description: this script prepares the observations and input data. Pre-processing and sampling is performed for observations. Potential inputs for the models are also calculated.
# Taxonomy ####
inputs$taxonomy <- data.table(inputs$taxonomy)

# Add missing taxonomic data for Drusus_manticola, Baetis_niger, Baetis_nexus
species <- rbind(c("Arthropoda", "Insecta", "Trichoptera", "Limnephilidae", "Drusus", "monticola"), 
                 c("Arthropoda", "Insecta", "Ephemeroptera", "Baetidae", "Baetis", "niger"),
                 c("Arthropoda", "Insecta", "Ephemeroptera", "Baetidae", "Baetis", "nexus"))
species <- data.table(species)
colnames(species) <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
inputs$taxonomy <- bind_rows(inputs$taxonomy, species)
rm(species)

# Modify $Species so it's identical to BDM column names
inputs$taxonomy$Species <- ifelse(!is.na(inputs$taxonomy$Species), paste(inputs$taxonomy$Genus, inputs$taxonomy$Species, sep="_"), NA)

# Invertebrate data ####
# Select invertebrate columns of interest
inv <- select(inv, SiteId, SampId, X, Y, Longitude, Latitude, MonitoringProgram, Altitude, Abundance, AbundanceClass, Phylum, Class, Order, Family, Genus, Species, Taxon, Day, Month, Year)

# Replace not-NA species entries with the taxon
inv$Species[!is.na(inv$Species)] <- inv$Taxon[!is.na(inv$Species)]
inv <- select(inv, -Taxon)

# Create presence-absence column: non-NA obs. set to 1, set NA obs. to 0
# inv$PA <- ifelse(inv$Abundance > 0, 1, 0)
# inv$PA[is.na(inv$PA)] <- 1

# SAMPLE BDMf ####
# Read and format column names
bdmf <- read.csv('inputs/BDM.fam_occ_data_2017-07-25.dat', header=T, sep='\t', stringsAsFactors = F)

# Select samples based on sampling window
bdmf <- bdmf[bdmf$samp.wind=='pref' | bdmf$samp.wind=='buf', ]
bdmf <- bdmf[!is.na(bdmf$SiteId), ]

bdmf <- select(bdmf, SiteId, SampId, contains("Occurrence.group."), contains("Occurrence."))
ind <- !(colnames(bdmf) %in% c("SiteId", "SampId"))
cnames <- colnames(bdmf[, ind])
cnames <- gsub("Occurrence.", "", cnames)
cnames <- gsub("group.", "", cnames)
colnames(bdmf) <- c("SiteId","SampId", cnames)

# Set seed and sample
set.seed(2017)
sample.bdmf <- bdmf[sample(nrow(bdmf)),]

n.bdmf <- occur.freq(sample.bdmf)

# # Choose number of folds
# k <- 3
# # Cut the randomly shuffled data into k-folds
# k.folds <- cut(seq(1, nrow(sample.bdmf)), breaks=k, labels=FALSE)
# 
# fold1 <- sample.bdmf[which(k.folds == 1),]
# fold2 <- sample.bdmf[which(k.folds == 2),]
# fold3 <- sample.bdmf[which(k.folds == 3),]
# 
# # Combine the folds manually for the joint model
# train1 <- bind_rows(fold1, fold2)
# test1 <- fold3
# 
# train2 <- bind_rows(fold2, fold3)
# test2 <- fold1
# 
# train3 <- bind_rows(fold1, fold3)
# test3 <- fold2
# 
# # Write data to /outputs 
# write.csv(train1, 'outputs/bdm.family.train1.csv', row.names = F)
# write.csv(test1, 'outputs/bdm.family.test1.csv', row.names = F)
# 
# write.csv(train2, 'outputs/bdm.family.train2.csv', row.names = F)
# write.csv(test2, 'outputs/bdm.family.test2.csv', row.names = F)
# 
# write.csv(train3, 'outputs/bdm.family.train3.csv', row.names = F)
# write.csv(test3, 'outputs/bdm.family.test3.csv', row.names = F)
# 
# write.csv(sample.bdmf, 'outputs/bdm.family.sample.csv', row.names = F)

### SAMPLE BDMs ####
# Read and format column names
bdms <- read.csv('inputs/BDM_occ_data_2017-07-25.dat', sep='\t', header=TRUE, na.strings=c("<Null>", "NA", ""), stringsAsFactors=FALSE)

bdms <- bdms[bdms$samp.wind=='pref' | bdms$samp.wind=='buf', ]

# # Are there repeated samples at a site per year?
# sample.check <- sites %>%
#   select(SiteId, SampId, Day, Month, Year, X, Y) %>%
#   group_by(SiteId) %>%
#   summarise(n.samples=n(), n.years=uniqueN(Year)) 

# Clean BDM dataset column names
bdms <- select(bdms, SiteId, SampId, contains("Occurrence.group."), contains("Occurrence."))
ind <- !(colnames(bdms) %in% c("SiteId", "SampId"))
cnames <- colnames(bdms[, ind])
cnames <- gsub("Occurrence.", "", cnames)
cnames <- gsub("group.", "", cnames)
colnames(bdms) <- c("SiteId","SampId", cnames)
rm(ind, cnames)

# Drop NA SiteIds
bdms <- bdms[!is.na(bdms$SiteId), ]

# Drop two genera because the family occurs as well
bdms <- select(bdms, -Silo, -Stactobia)

# Set seed and sample (i.e., shuffle the data)
set.seed(2017)
sample.bdms <- bdms[sample(nrow(bdms)),]

n.bdms <- occur.freq(sample.bdms)

# Choose number of folds
k <- 3
# Cut the randomly shuffled data into k-folds
k.folds <- cut(seq(1, nrow(sample.bdms)), breaks=k, labels=FALSE)

fold1 <- sample.bdms[which(k.folds == 1),]
fold2 <- sample.bdms[which(k.folds == 2),]
fold3 <- sample.bdms[which(k.folds == 3),]

# Combine the folds manually for the joint model
train1 <- bind_rows(fold1, fold2)
test1 <- fold3

train2 <- bind_rows(fold2, fold3)
test2 <- fold1

train3 <- bind_rows(fold1, fold3)
test3 <- fold2

# # Write data to /outputs 
# write.csv(train1, 'outputs/bdm.species.train1.csv', row.names = F)
# write.csv(test1, 'outputs/bdm.species.test1.csv', row.names = F)
# 
# write.csv(train2, 'outputs/bdm.species.train2.csv', row.names = F)
# write.csv(test2, 'outputs/bdm.species.test2.csv', row.names = F)
# 
# write.csv(train3, 'outputs/bdm.species.train3.csv', row.names = F)
# write.csv(test3, 'outputs/bdm.species.test3.csv', row.names = F)
# 
# write.csv(sample.bdms, 'outputs/bdm.species.sample.csv', row.names = F)

# SAMPLE CFCH ####
invf <- read.csv('inputs/ALL_occ_data_2017-07-25.dat', header = TRUE, sep = '\t', stringsAsFactors=FALSE)
invf <- invf[invf$samp.wind=='pref' | invf$samp.wind=='buf', ]
invf <- invf[!is.na(invf$SiteId), ]

invf <- select(invf, SiteId, SampId, contains("Occurrence.group."), contains("Occurrence."))
ind <- !(colnames(invf) %in% c("SiteId", "SampId"))
cnames <- colnames(invf[, ind])
cnames <- gsub("Occurrence.", "", cnames)
cnames <- gsub("group.", "", cnames)
colnames(invf) <- c("SiteId","SampId", cnames)

# Duplicate sample IDs identified, differ only by Altitude_m value;
# get the unique rows *without* Altitude_m column
invf$SampId[duplicated(invf$SampId)]
invf <- unique(invf)

# Set seed and sample
set.seed(2017)
sample.invf <- invf[sample(nrow(invf)),]

# Count occurrences
n.invf <- occur.freq(sample.invf)

sample.invf <- sample.invf[, c("SiteId", "SampId", names(n.invf[n.invf > 0]))]

# # Choose number of folds
# k <- 3
# # Cut the randomly shuffled data into k-folds
# k.folds <- cut(seq(1, nrow(sample.invf)), breaks=k, labels=FALSE)
# 
# fold1 <- sample.invf[which(k.folds == 1),]
# fold2 <- sample.invf[which(k.folds == 2),]
# fold3 <- sample.invf[which(k.folds == 3),]
# 
# # Combine the folds manually for the joint model
# train1 <- bind_rows(fold1, fold2)
# test1 <- fold3
# 
# train2 <- bind_rows(fold2, fold3)
# test2 <- fold1
# 
# train3 <- bind_rows(fold1, fold3)
# test3 <- fold2
# 
# # Write data to /outputs 
# write.csv(train1, 'outputs/inv.family.train1.csv', row.names = F)
# write.csv(test1, 'outputs/inv.family.test1.csv', row.names = F)
# 
# write.csv(train2, 'outputs/inv.family.train2.csv', row.names = F)
# write.csv(test2, 'outputs/inv.family.test2.csv', row.names = F)
# 
# write.csv(train3, 'outputs/inv.family.train3.csv', row.names = F)
# write.csv(test3, 'outputs/inv.family.test3.csv', row.names = F)
# 
# write.csv(sample.invf, 'outputs/inv.family.sample.csv', row.names = F)

# SAMPLE CFp ####
# # Subset samples in the plateau?
sample.invfp <- sample.invf[sample.invf$SiteId %in% env$SiteId[env$BIOGEO=="Mittelland"],]
n.invfp <- occur.freq(sample.invfp)
sample.invfp <- sample.invfp[, c("SiteId", "SampId", names(n.invfp[n.invfp > 0]))]

# # Choose number of folds
# k <- 3
# # Cut the randomly shuffled data into k-folds
# k.folds <- cut(seq(1, nrow(sample.invfp)), breaks=k, labels=FALSE)
# 
# fold1 <- sample.invfp[which(k.folds == 1),]
# fold2 <- sample.invfp[which(k.folds == 2),]
# fold3 <- sample.invfp[which(k.folds == 3),]
# 
# # Combine the folds manually for the joint model
# train1 <- bind_rows(fold1, fold2)
# test1 <- fold3
# 
# train2 <- bind_rows(fold2, fold3)
# test2 <- fold1
# 
# train3 <- bind_rows(fold1, fold3)
# test3 <- fold2
# 
# # Write data to /outputs 
# write.csv(train1, 'outputs/inv.family.plateau.train1.csv', row.names = F)
# write.csv(test1, 'outputs/inv.family.plateau.test1.csv', row.names = F)
# 
# write.csv(train2, 'outputs/inv.family.plateau.train2.csv', row.names = F)
# write.csv(test2, 'outputs/inv.family.plateau.test2.csv', row.names = F)
# 
# write.csv(train3, 'outputs/inv.family.plateau.train3.csv', row.names = F)
# write.csv(test3, 'outputs/inv.family.plateau.test3.csv', row.names = F)
# 
# write.csv(sample.invfp, 'outputs/inv.family.plateau.sample.csv', row.names = F)

# Rename land use from DE to EN codes
codes <- inputs$land_metadata$Code
codes <- codes[-match("EZG_NR", codes)]
land.use <- env[, c("SiteId", "SampId", codes)]

english_fields <- inputs$land_metadata$DesignatedCode
english_fields <- english_fields[-match("UID", inputs$land_metadata$DesignatedCode)]
english_fields <- c("SiteId", "SampId", english_fields)

colnames(land.use) <- english_fields
inputs$land.use <- as.data.table(land.use)
rm(english_fields)

### > PREPARE PREDICTORS ####
# Create 'p', to contain pre-processed predictors
# Site/sample ID, coordinates, catchment ID, catchment area
p <- select(env, SiteId, SampId, X, Y, EZG_NR, A_EZG)
p$EZG_km <- p$A_EZG/1000000

# Joint land use and select fields of interest
p <- left_join(p, land.use, by = c("SiteId", "SampId"))
p <- select(p, SiteId, SampId, X, Y, EZG_NR, EZG_km, starts_with("p_"), starts_with("v_"), livestock_unit)

# Insecticides ####
# Mean annual spray treatments times the fraction of each cropland in catchment
p <- mutate(p, IAR = 0.44*p_potatoes + 0.03*p_total_grain + 0.01*p_corn + 1.82*p_rapeseed + 0.07*p_beet + 2.66*p_vegetables + 3.10*p_orchard + 0.37*p_vineyard + 0.38*p_legume)

# Urban / Transport ####
p <- mutate(p, Urban = p_settlement + p_railtrack + p_road_highway + p_road_cantonal + p_road_munic)
# p <- mutate(p, HM = p_bldg_env + p_railtrack + p_road_highway + p_road_cantonal + p_road_munic + p_landfill)
# p <- mutate(p, Urban = p_settlement + p_bldg_env + p_urbangreen + p_landfill)
# p <- mutate(p, Transport = l_road_settlement + l_railtrack + l_road_highway + l_road_cantonal + l_road_munic)
# p$Transport <- (p$Transport/1000)/p$CatchArea_km

# Land use: AS 2009 (UI) ####
asm <- as09.metadata[as09.metadata$coefficient == TRUE, ]
asm <- as09.metadata[, c("AS04_72", "AS04_72_Description", "value")]
colnames(asm) <- c("Code_72", "Description", "Coefficient")

files <- list.files(path = 'inputs/ArealStatistics09', pattern="*.csv")

data.list <- lapply(files, function(f){
  x <- read.csv(paste('inputs/ArealStatistics09/', f, sep = ''), sep = ',', header = TRUE, stringsAsFactors = F)
})

data.list <- lapply(data.list, function(table){
  table <- gather(table, Variable, Value, -SITEID)
})

as09.data <- do.call(rbind, data.list)
rm(files, data.list)

# Clean areal statistics data
catchment.area <- as09.data[as09.data$Variable == 'Shape_Area',]
as09.data <- as09.data[!as09.data$Variable %in% c("OID", "Shape_Area"),]
as09.data <- filter(as09.data, !Variable %in% c("OID", "Shape_Area"))
as09.data$Variable <- as.character(as09.data$Variable)
colnames(as09.data) <- c("SiteId", "Code_72", "LandArea")
catchment.area <- catchment.area[, c("SITEID", "Value")]
colnames(catchment.area) <- c("SiteId", "CatchmentArea")
as09.data <- left_join(as09.data, catchment.area, by = "SiteId")

codes <- strsplit(as09.data$Code_72, "_")
codes <- sapply(codes, function(x){ as.integer(x[2])})
as09.data$Code_72 <- codes
as09.data$Percent <- (as09.data$LandArea/as09.data$CatchmentArea)*100
rm(codes, catchment.area)

# Filter land use by AS09_72 urban codes used in Mutzner et al. 2016
as09.data <- as09.data[as09.data$Code_72 %in% asm$Code_72, ]
as09.data <- left_join(as09.data, asm, by = "Code_72")
as09.data$Coefficient[is.na(as09.data$Coefficient)] <- 0
as09.data$UrbanIndex <- as09.data$Percent*as09.data$Coefficient

urban.index <- aggregate(UrbanIndex ~ SiteId, as09.data, sum)
urban.index$UI <- (urban.index$UrbanIndex/max(urban.index$UrbanIndex))
urban.index <- urban.index[,c("SiteId", "UI")]
urban.index[is.na(urban.index$UI)] <- 0

p <- left_join(p, urban.index, by = "SiteId")
rm(urban.index) # remove separate urban index dataset

# Land use: AS 1997 ####
files <- list.files(path = 'inputs/ArealStatistics97', pattern="*.csv")

data.list <- lapply(files, function(f){
  x <- read.csv(paste('inputs/ArealStatistics97/', f, sep = ''), sep = ',', header = TRUE, stringsAsFactors = F)
})

data.list <- lapply(data.list, function(table){
  table <- gather(table, Variable, Value, -SITEID)
})

as97.data <- do.call(rbind, data.list)
rm(files, data.list)

# Clean areal statistics data
catchment.area <- as97.data[as97.data$Variable == 'Shape_Area',]
as97.data <- as97.data[!as97.data$Variable %in% c("OID", "Shape_Area"), ]
as97.data$Variable <- as.character(as97.data$Variable)
colnames(as97.data) <- c("SiteId", "Code_72", "LandArea")
catchment.area <- catchment.area[, c("SITEID", "Value")]
colnames(catchment.area) <- c("SiteId", "CatchmentArea")
as97.data <- left_join(as97.data, catchment.area, by = "SiteId")

codes <- strsplit(as97.data$Code_72, "_")
codes <- sapply(codes, function(x){ as.integer(x[2])})
as97.data$Code_72 <- codes
as97.data$Percent <- (as97.data$LandArea/as97.data$CatchmentArea)*100
rm(codes, catchment.area)

# Land use ~ Buffer ####
# Get the buffer areas for both sub-catchments (AREA) catchments (ASUM)
b$areas <- select(b$areas, -OBJECTID, -h1, -h2, -A_SUBEZG, -A_EZG, -Sub_Area)

# Select land use columns of interest inside catchment b
b$jgrass <- select(b$grass, EZG_NR, Grass_B10, Grass_B20, Grass_B50, Grass_B100, Grass_B500, Grass_B1000)
b$jforest <- select(b$forest, EZG_NR, Forest_B10, Forest_B20, Forest_B50, Forest_B100, Forest_B500, Forest_B1000)
b$jurban <- select(b$urban, EZG_NR, Urban_B10, Urban_B20, Urban_B50, Urban_B100, Urban_B500, Urban_B1000)
b$jarable <- select(b$arable, EZG_NR, Perc_B10, Perc_B20, Perc_B50, Perc_B100, Perc_B500, Perc_B1000)
colnames(b$jarable) <- c("EZG_NR", "A10m", "A20m", "A50m", "A100m", "A500m", "A1km")
b$jftype <- select(b$ftype, EZG_NR, FTyp_B10_VALUE_1, FTyp_B20_VALUE_1, FTyp_B50_VALUE_1, FTyp_B100_VALUE_1, FTyp_B500_VALUE_1, FTyp_B1000_VALUE_1, FTyp_B10_VALUE_4, FTyp_B20_VALUE_4, FTyp_B50_VALUE_4, FTyp_B100_VALUE_4, FTyp_B500_VALUE_4, FTyp_B1000_VALUE_4)
colnames(b$jftype) <- c("EZG_NR", "F1.10", "F1.20", "F1.50", "F1.100", "F1.500", "F1.1km", "F4.10", "F4.20", "F4.50", "F4.100", "F4.500", "F4.1km")

# Calculate proportions of grassland by catchment
b$jgrass <- left_join(b$jgrass, b$areas, by = "EZG_NR")
b$jgrass$G10m <- (b$jgrass$Grass_B10/b$jgrass$ASUM_B10)*100
b$jgrass$G20m <- (b$jgrass$Grass_B20/b$jgrass$ASUM_B20)*100
b$jgrass$G50m <- (b$jgrass$Grass_B50/b$jgrass$ASUM_B50)*100
b$jgrass$G100m <- (b$jgrass$Grass_B100/b$jgrass$ASUM_B100)*100
b$jgrass$G500m <- (b$jgrass$Grass_B500/b$jgrass$ASUM_B500)*100
b$jgrass$G1km <- (b$jgrass$Grass_B1000/b$jgrass$ASUM_B1000)*100
b$jgrass <- select(b$jgrass, EZG_NR, G10m, G20m, G50m, G100m, G500m, G1km)

# Calculate proportions of forest by catchment
b$jforest <- left_join(b$jforest, b$areas, by = "EZG_NR")
b$jforest$F10m <- (b$jforest$Forest_B10/b$jforest$ASUM_B10)*100
b$jforest$F20m <- (b$jforest$Forest_B20/b$jforest$ASUM_B20)*100
b$jforest$F50m <- (b$jforest$Forest_B50/b$jforest$ASUM_B50)*100
b$jforest$F100m <- (b$jforest$Forest_B100/b$jforest$ASUM_B100)*100
b$jforest$F500m <- (b$jforest$Forest_B500/b$jforest$ASUM_B500)*100
b$jforest$F1km <- (b$jforest$Forest_B1000/b$jforest$ASUM_B1000)*100
b$jforest <- select(b$jforest, EZG_NR, F10m, F20m, F50m, F100m, F500m, F1km)

# Calculate proportions of forest (type 1 and type 4) by catchment
b$jftype <- left_join(b$jftype, b$areas, by = "EZG_NR")
b$jftype$F1.10 <- (b$jftype$F1.10/b$jftype$ASUM_B10)*100 
b$jftype$F1.20 <- (b$jftype$F1.20/b$jftype$ASUM_B20)*100
b$jftype$F1.50 <- (b$jftype$F1.50/b$jftype$ASUM_B50)*100
b$jftype$F1.100 <- (b$jftype$F1.100/b$jftype$ASUM_B100)*100
b$jftype$F1.500 <- (b$jftype$F1.500/b$jftype$ASUM_B500)*100
b$jftype$F1.1km <- (b$jftype$F1.1km/b$jftype$ASUM_B1000)*100

b$jftype$F4.10 <- (b$jftype$F4.10/b$jftype$ASUM_B10)*100 
b$jftype$F4.20 <- (b$jftype$F4.20/b$jftype$ASUM_B20)*100
b$jftype$F4.50 <- (b$jftype$F4.50/b$jftype$ASUM_B50)*100
b$jftype$F4.100 <- (b$jftype$F4.100/b$jftype$ASUM_B100)*100
b$jftype$F4.500 <- (b$jftype$F4.500/b$jftype$ASUM_B500)*100
b$jftype$F4.1km <- (b$jftype$F4.1km/b$jftype$ASUM_B1000)*100
b$jftype <- select(b$jftype, EZG_NR, F1.10, F1.20, F1.50, F1.100, F1.500, F1.1km, F4.10, F4.20, F4.50, F4.100, F4.500, F4.1km)

# Calculate proportions of urban by catchment
b$jurban <- left_join(b$jurban, b$areas, by = "EZG_NR")
b$jurban$U10m <- (b$jurban$Urban_B10/b$jurban$ASUM_B10)*100
b$jurban$U20m <- (b$jurban$Urban_B20/b$jurban$ASUM_B20)*100
b$jurban$U50m <- (b$jurban$Urban_B50/b$jurban$ASUM_B50)*100
b$jurban$U100m <- (b$jurban$Urban_B100/b$jurban$ASUM_B100)*100
b$jurban$U500m <- (b$jurban$Urban_B500/b$jurban$ASUM_B500)*100
b$jurban$U1km <- (b$jurban$Urban_B1000/b$jurban$ASUM_B1000)*100
b$jurban <- select(b$jurban, EZG_NR, U10m, U20m, U50m, U100m, U500m, U1km)

# Proportions of urban by sub-catchment
b$urban.sub[is.na(b$urban.sub)] <- 0
b$urban.sub <- left_join(b$urban.sub, b$areas, by = "EZG_NR")
b$urban.sub$Usub10m <- (b$urban.sub$Urban_B10/b$urban.sub$AREA_B10)*100
b$urban.sub$Usub20m <- (b$urban.sub$Urban_B20/b$urban.sub$AREA_B20)*100
b$urban.sub$Usub50m <- (b$urban.sub$Urban_B50/b$urban.sub$AREA_B50)*100
b$urban.sub$Usub100m <- (b$urban.sub$Urban_B100/b$urban.sub$AREA_B100)*100
b$urban.sub$Usub500m <- (b$urban.sub$Urban_B500/b$urban.sub$AREA_B500)*100
b$urban.sub$Usub1km <- (b$urban.sub$Urban_B1000/b$urban.sub$AREA_B1000)*100
b$urban.sub <- select(b$urban.sub, EZG_NR, Usub10m, Usub20m, Usub50m, Usub100m, Usub500m, Usub1km)

# Proportions of agriculture by sub-catchment
b$arable.sub <- select(b$arable, EZG_NR, Arable_B10, Arable_B20, Arable_B50, Arable_B100, Arable_B500, Arable_B1000, AREA_B10, AREA_B20, AREA_B50, AREA_B100, AREA_B500, AREA_B1000)
b$arable.sub$Asub10m <- (b$arable.sub$Arable_B10/b$arable.sub$AREA_B10)*100
b$arable.sub$Asub20m <- (b$arable.sub$Arable_B20/b$arable.sub$AREA_B20)*100
b$arable.sub$Asub50m <- (b$arable.sub$Arable_B50/b$arable.sub$AREA_B50)*100
b$arable.sub$Asub100m <- (b$arable.sub$Arable_B100/b$arable.sub$AREA_B100)*100
b$arable.sub$Asub500m <- (b$arable.sub$Arable_B500/b$arable.sub$AREA_B500)*100
b$arable.sub$Asub1km <- (b$arable.sub$Arable_B1000/b$arable.sub$AREA_B1000)*100
b$arable.sub <- select(b$arable.sub, EZG_NR, Asub10m, Asub20m, Asub50m, Asub100m, Asub500m, Asub1km)

# Calculate proportions of forest by sub-catchment
b$forest.sub[is.na(b$forest.sub)] <- 0
b$forest.sub <- left_join(b$forest.sub, b$areas, by = "EZG_NR")
b$forest.sub$Fsub10m <- (b$forest.sub$Forest_B10/b$forest.sub$AREA_B10)*100
b$forest.sub$Fsub20m <- (b$forest.sub$Forest_B20/b$forest.sub$AREA_B20)*100
b$forest.sub$Fsub50m <- (b$forest.sub$Forest_B50/b$forest.sub$AREA_B50)*100
b$forest.sub$Fsub100m <- (b$forest.sub$Forest_B100/b$forest.sub$AREA_B100)*100
b$forest.sub$Fsub500m <- (b$forest.sub$Forest_B500/b$forest.sub$AREA_B500)*100
b$forest.sub$Fsub1km <- (b$forest.sub$Forest_B1000/b$forest.sub$AREA_B1000)*100
b$forest.sub <- select(b$forest.sub, EZG_NR, Fsub10m, Fsub20m, Fsub50m, Fsub100m, Fsub500m, Fsub1km)

# Join all the buffer proportions to the initial predictor dataset
p <- left_join(p, b$jarable, by = "EZG_NR")
p <- left_join(p, b$jgrass, by = "EZG_NR")
p <- left_join(p, b$jforest, by = "EZG_NR")
p <- left_join(p, b$jurban, by = "EZG_NR")
p <- left_join(p, b$jftype, by = "EZG_NR")

p <- left_join(p, b$forest.sub, by = "EZG_NR")
p <- left_join(p, b$urban.sub, by = "EZG_NR")
p <- left_join(p, b$arable.sub, by = "EZG_NR")
p <- left_join(p, b$areas, by = "EZG_NR")

# Temperature ####
p <- left_join(p, env[, c("SiteId", "SampId", "temperature")], by = c("SiteId", "SampId"))
colnames(p)[match("temperature", colnames(p))] <- "Temp"
# mean maximum summer temperature
# p$Temp <- 17.836 + 1.858*log10(p$CatchArea_km+1) - 0.005*(p$CatchElev)

# Flow velocity ####
flow <- env[ , c("SampId", "v.old", "v.new")]
colnames(flow) <- c("SampId", "FV.old", "FV")
flow$FV[is.na(flow$FV)] <- flow$FV.old[is.na(flow$FV)]
flow$FV.old <- NULL
p <- left_join(p, flow, by = "SampId")
rm(flow)

# Site removed: site on disconnected river in large catchement.
# Small CW and high Q yields flow velocity of ~60 m/s
p <- filter(p, SiteId != "CSCF_69_850_GE")

# Livestock units ####
p$LUD <- p$livestock_unit/p$EZG_km

# FRI / bFRI ####
p <- left_join(p, env[, c("SiteId", "SampId", "bFRI", "FRI")], by = c("SiteId", "SampId"))

p$Forest <- p$p_forest
# Proportion of river length intersected by forests, upstream and within buffer distance
# inputs$bfri$bFRI <- (inputs$bfri$ForestLength/inputs$bfri$RiverLength)*100
# inputs$fri$FRI <- (inputs$fri$ForestLength/inputs$fri$RiverLength)*100
# inputs$bfri <- select(inputs$bfri, SiteId, bFRI)
# inputs$fri <- select(inputs$fri, SiteId, FRI)
# 
# p <- left_join(p, inputs$bfri, by = "SiteId")
# p <- left_join(p, inputs$fri, by = "SiteId")

# IDW land use ####
p <- left_join(p, env[, c("SiteId", "SampId", "A.EDO","F.EDO","A.FLO","F.FLO","A.HAFLO","F.HAFLO")], by = c("SiteId", "SampId"))

# NAs contain zero target land use
# idw$edo_farm$idw.num[is.na(idw$edo_farm$idw.num)] <- 0
# idw$edo_forest$idw.num[is.na(idw$edo_forest$idw.num)] <- 0
# 
# idw$flo_farm$idw.num[is.na(idw$flo_farm$idw.num)] <- 0
# idw$flo_forest$idw.num[is.na(idw$flo_forest$idw.num)] <- 0
# 
# idw$ha_flo_forest$idw.num[is.na(idw$ha_flo_forest$idw.num)] <- 0
# idw$ha_flo_farm$idw.num[is.na(idw$ha_flo_farm$idw.num)] <- 0
# 
# # Calculate and join spatially weighted land use
# idw$edo_farm$A.EDO <- (idw$edo_farm$idw.num/idw$edo_farm$idw.denom)*100
# idw$edo_forest$F.EDO <- (idw$edo_forest$idw.num/idw$edo_forest$idw.denom)*100
# 
# idw$flo_farm$A.FLO <- (idw$flo_farm$idw.num/idw$flo_farm$idw.denom)*100
# idw$flo_forest$F.FLO <- (idw$flo_forest$idw.num/idw$flo_forest$idw.denom)*100
# 
# idw$ha_flo_farm$A.HAFLO <- (idw$ha_flo_farm$idw.num/idw$ha_flo_farm$idw.denom)*100
# idw$ha_flo_forest$F.HAFLO <- (idw$ha_flo_forest$idw.num/idw$ha_flo_forest$idw.denom)*100
# 
# idw$edo_farm <- select(idw$edo_farm, SiteId, A.EDO)
# idw$edo_forest <- select(idw$edo_forest, SiteId, F.EDO)
# 
# idw$flo_farm <- select(idw$flo_farm, SiteId, A.FLO)
# idw$flo_forest <- select(idw$flo_forest, SiteId, F.FLO)
# 
# idw$ha_flo_farm <- select(idw$ha_flo_farm, SiteId, A.HAFLO)
# idw$ha_flo_forest <- select(idw$ha_flo_forest, SiteId, F.HAFLO)
# 
# p <- left_join(p, idw$edo_farm, by = "SiteId")
# p <- left_join(p, idw$edo_forest, by = "SiteId")
# p <- left_join(p, idw$flo_farm, by = "SiteId")
# p <- left_join(p, idw$flo_forest, by = "SiteId")
# p <- left_join(p, idw$ha_flo_farm, by = "SiteId")
# p <- left_join(p, idw$ha_flo_forest, by = "SiteId")

# Toxic units ####
inputs$mp$EZG_EZG_NR <- NULL
colnames(inputs$mp) <- c("SiteId", "TU.Dm", "TU.Cr")
p <- left_join(p, inputs$mp, by = "SiteId")

# Morphology ####
morph <- env[, c("SiteId", "SampId", "width.variability", "bed.modification", "ecomorphology")]
colnames(morph) <- c("SiteId", "SampId", "WV", "BM", "Morph")
p <- left_join(p, morph, by = c("SiteId", "SampId"))
rm(morph) # remove separate morphology dataset

# Wastewater ####
p <- left_join(p, env[, c("SiteId", "SampId", "Discharge")], by = c("SiteId", "SampId"))
p$WW <- p$v_ww # number of connected inhabitants * 400L/person/day * 365 days/year / 1000L/m3
p$Annual_Discharge <- p$Discharge*60*60*24*365 # Annual discharge (m3/yr) = 60 s/min * 60 min/hour * 24 hours/day * 365 days/year
p$WW <- (p$WW/p$Annual_Discharge)*100
p <- p[p$WW <= 100,]# Sites dropped due to WW exceeding Q

# Hydropeaking ####
hp <- select(env, SiteId, SampId, Residual_flow, Residual_flow_Value)
hp$Residual_flow_Value[hp$Residual_flow == 0] <- 1
hp[hp == 9999] <- NA
hp <- select(hp, SiteId, SampId, Residual_flow_Value)
colnames(hp) <- c("SiteId", "SampId", "HP")
p <- left_join(p, hp, by=c("SiteId", "SampId"))
rm(hp)

# Substrates ####
# Substrates are only available for BDM sites
# So perhaps it would be better to just try 1. "mobile blocs>250mm +  coarse sediments (25-250mm)", 2. "mosses+hydrophytes+coarse organic matter" as two classes that are favorable, and 3. "sand and silt <2.5 mm + mud <0.1 mm" as probably rather unfavorable classes. 
substrate <- select(env, SiteId, SampId, starts_with("sfract_"))
substrate <- substrate[substrate$SiteId %in% bdms$SiteId, ]
substrate$SS1 <- substrate$sfract_mobile_blocks + substrate$sfract_coarse_inorganic_sediments
substrate$SS2 <- substrate$sfract_moss + substrate$sfract_hydrophytes + substrate$sfract_coarse_organic_matter
substrate$SS3 <- substrate$sfract_sand_silt + substrate$sfract_fine_sediments

substrate <- substrate[, c("SiteId", "SampId", "SS1", "SS2", "SS3")]
p <- left_join(p, substrate, by = c("SiteId", "SampId"))

inputs$substrate <- substrate; rm(substrate)

# Golf courses ###
# inputs$golf <- select(inputs$golf, SiteId, pgolf)
# colnames(inputs$golf) <- c("SiteId", "p_golf")
# p <- left_join(p, inputs$golf, by = "SiteId")
# p$p_golf[is.na(p$p_golf)] <- 0

# Drop sites ####
# p <- p[!(p$SiteId == "CSCF_50_VD"),] # WW outlier (13.45 WW:Q)

# Random noise ####
p$Noise <- runif(nrow(p), 0, 1)

# Stream order
p <- left_join(p, select(env, SiteId, SampId, FLOZ), by=c("SiteId", "SampId"))

# Delete objects ####
rm(as09.data, as09.metadata, as97.data, asm, b, fold1, fold2, fold3, train1, train2, train3, test1, test2, test3, cnames, ind, k, k.folds)
