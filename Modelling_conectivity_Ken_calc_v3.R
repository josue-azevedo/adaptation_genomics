## -----------------------------------------------------------------------------
#load(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_RDA/RDA_Ken_calcarata.RData") # 


## -----------------------------------------------------------------------------
#save.image(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_RDA/Modelling_Ken_calcarata.RData") # 


## -----------------------------------------------------------------------------
#load(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_RDA/Modelling_Ken_calcarata.RData") # 

## -----------------------------------------------------------------------------
## Loading libraries
library(cowplot)
library(dplyr)
library(ggplot2)
library(LEA)
library(magrittr)
library(plyr)
library(tidyr)
library(ggpubr)
library(vioplot)
library(pbapply)
# Load main SIG libraries
library(dismo)
library(rgdal)
library(maptools)
library(SDMTools)
library(spdep)
library(sdm)
library(raster)
library(maptools)
library(rgeos)
library(igraph)
library(scatterpie)
library(sf)
library(rgdal)
library(rgeos)
library(ggnewscale)
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA
library(usdm)     # Used to run VIF
library(sfheaders)
library(viridis)
library(gdm)
library(StAMPP) # calculate pairwise fst distances
library(RColorBrewer)
#install.packages("gradientForest", repos="http://R-Forge.R-project.org") # gradient forest is not available in CRAN
#library(gradientForest)
#install.packages("sdm")
#library(beepr)
#beep()

# Adding a notin argument to facilitate the analyses
`%notin%` <- Negate(`%in%`)

# Add path
path = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses"

# Data projections #####
oldproj <- paste0(" +init=epsg:4326") #this is WGS84 most commonly used for google earth etc. in decimal degrees
behrmannCRS <- CRS('+proj=cea +lat_ts=30') # This is an equal area projection, i.e., cells have the same area
#.####

## Detailed SAmerican  map
SAshp <- readOGR("/Users/josue/Dropbox/1Doutorado/Chapter_Wer/South_America/South_America.shp")
SAshp2 = spTransform(SAshp,CRSobj = behrmannCRS)

## Shapefile amazonia
bio_SPDF<- readOGR("/Users/josue/Dropbox/1Doutorado/Chapter_2/biodiverse_pipeline-master/wwf_biomes_simp.shp")
bio_SPDF <-spTransform(bio_SPDF, CRSobj=behrmannCRS)

bio_SPDF3 = crop(bio_SPDF,selected_clipped_biomes_3)
# savanna biome=7
# floresta biome= 1
# floresta sub biome 2
# caatinga biome = 13

#refine_ext = drawExtent()
#class      : Extent 
#xmin       : -6732169 
#xmax       : -3112643 
#ymin       : -2019442 
#ymax       : 1242340 

refine_ext = c(-6732169,-3112643,-2019442,1242340) # here I draw a polygon to delimit my study area 

# Cropt to extent including only eastern Amazonia
bio_SPDF2 = crop(bio_SPDF, refine_ext)
plot(bio_SPDF)
plot(bio_SPDF2)

## -----------------------------------------------------------------------------
#install.packages("MigClim",dependencies = TRUE)
library(MigClim)


## -----------------------------------------------------------------------------
#MigClim.userGuide()

## -----------------------------------------------------------------------------

# From the last step of Var Selection GEA
# Molecular data only
classified_ind_molecular = read.csv("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/classified_ind_molecular.csv",header = TRUE)

head(classified_ind_molecular)
classified_ind_molecular[classified_ind_molecular$kmeans==2,]
classified_ind_molecular_pts = classified_ind_molecular
coordinates(classified_ind_molecular_pts) = classified_ind_molecular_pts[,c("long","lat")]
proj4string(classified_ind_molecular_pts) = behrmannCRS
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==2,], add=TRUE,pch=16,cex=.5,col="white")

classified_ind = read.csv("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/class_ind_all_pts.csv",header = TRUE)

classified_ind_6_pts = classified_ind
coordinates(classified_ind_6_pts) = classified_ind_6_pts[,c("long","lat")]
proj4string(classified_ind_6_pts) = behrmannCRS

head(classified_ind)
classified_ind$kind = classified_ind$kmeans
table(classified_ind$kind)
  
classified_ind$kind[classified_ind$kind == 1] <- "warm_wet_stable"
classified_ind$kind[classified_ind$kind == 2] <- "hot_dry_seazonal"

table(classified_ind$kind)

## -----------------------------------------------------------------------------
#### importando shapefile amazonia ####
#bio_SPDF<- readOGR("/Users/josue/Dropbox/1Doutorado/Chapter_2/biodiverse_pipeline-master/wwf_biomes_simp.shp")
#bio_SPDF <-spTransform(bio_SPDF, CRSobj=behrmannCRS)
#1 = Tropical & Subtropical Moist Broadleaf Forests
#2 = Tropical & Subtropical Dry Broadleaf Forests
#3 = Tropical & Subtropical Coniferous Forests
#4 = Temperate Broadleaf & Mixed Forests
#5 = Temperate Conifer Forests
#6 = Boreal Forests/Taiga
#7 = Tropical & Subtropical Grasslands, Savannas & Shrublands
#8 = Temperate Grasslands, Savannas & Shrublands
#9 = Flooded Grasslands & Savannas
#10 = Montane Grasslands & Shrublands
#11 = Tundra
#12 = Mediterranean Forests, Woodlands & Scrub
#13 = Deserts & Xeric Shrublands
#14 = Mangroves

#selected_biomes = bio_SPDF[bio_SPDF@data$BIOME %in% c(1,7,2,9,13),]

#convert sp to sf for use with this function 
#ncsf_sp<-as(selected_biomes, "sf")

#ncmp = st_cast(ncsf_sp,"POLYGON") # Transforming to Multipolygon


## -----------------------------------------------------------------------------
# I will refine the extent with Q-Gis
#convert sf to sp for use with other functions
#ncmp<-as(ncmp, "Spatial")
#writeOGR(ncmp,"/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/layers/selected_biomes.shp" ,"BIOME", driver="ESRI Shapefile",overwrite = TRUE)

## -----------------------------------------------------------------------------
# Taken from the Variable Selection for GEA script
selected_clipped_biomes_2 = readRDS("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/selected_clipped_biomes_2.rds")

plot(selected_clipped_biomes_2)


## -----------------------------------------------------------------------------
selected_clipped_biomes_3 = selected_clipped_biomes_2

selected_clipped_biomes_3 = as(selected_clipped_biomes_3, 'Spatial')
selected_clipped_biomes_3 # Check if needed later

## -----------------------------------------------------------------------------

# Getting Current bioclimatic rasters produced in Variable Selection for GEA
my_bioclim_rasters = readRDS("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/my_bioclim_rasters.rds")

#my_bioclim_rasters[[1]] = CLIMOND
#my_bioclim_rasters[[2]] = CHELSA

preds.climond_stck = my_bioclim_rasters[[2]] # Now with CHELSA

## -----------------------------------------------------------------------------
names(preds.climond_stck) 
preds_modelling_mask = subset(preds.climond_stck,c('bio_4','bio_5','bio_15','bio_18'))#not including bio_5 as it has the same variation as bio_6

# As less variables included, more extreme are the models
preds_modelling_mask = subset(preds.climond_stck,c('bio_2','bio_4','bio_5',"bio_17",'bio_15','bio_18'))#not including bio_5 as it has the same variation as bio_6

# Bio2 is more stable across time
# Bio11-12 also stable
# All other temperature variables are horrible!

names(preds_modelling_mask) 

## -----------------------------------------------------------------------------
plot(preds_modelling_mask)
layerStats(preds_modelling_mask, 'pearson', na.rm=T)

vifstep(raster::extract(preds_modelling_mask,classified_ind_6_pts))

## -----------------------------------------------------------------------------
# Log transforming precipitation variables
#preds_modelling_mask$bio_18 = log(preds_modelling_mask$bio_18)
#preds_modelling_mask$bio_17 = log(preds_modelling_mask$bio_17)

## -----------------------------------------------------------------------------
# Primeiro liste todos os rasters ambientais por nome. path= indica o diretório em que as variáveis estão

# For climond
#future_layers_path <- list.files(path="/Users/josue/Dropbox/4Environmental_layers/future_climond",pattern='zip',full.names = T) 
#future_layers_path

#get_future_bioclim = function(future_layers_path){

list_zip = unzip(zipfile = future_layers_path, list = TRUE)
list_zip2 = grep("Bio", list_zip$Name, value=TRUE)
list_zip2 = list_zip2[1:19]
#list_zip2
link <- paste0("/vsizip/",future_layers_path,"/",list_zip2)
#link
test <- stack(link)
#test
proj4string(test) = oldproj
test = projectRaster(test,crs = behrmannCRS)
test = mask(crop(test,selected_clipped_biomes_2),selected_clipped_biomes_2)
test = raster::resample(test,preds_modelling_mask) # just to have the same resolution as the chelsa rasters
names(test) = c('bio_1','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19')
plot(test[[1]])    
return(test)
}
#list_stack_future_bioclim = lapply(future_layers_path,get_future_bioclim)

## -----------------------------------------------------------------------------

selected_clipped_biomes_oldproj = spTransform(selected_clipped_biomes_3,CRSobj = "+proj=longlat +datum=WGS84 +no_defs ")
plot(selected_clipped_biomes_oldproj)

# For Chelsa (processing 2050 separated from 2100 because of disk space)
#future_layers_path <- list.files(path="/Users/josue/Dropbox/4Environmental_layers/future_chelsa",pattern='2041-2070',full.names = T) 
#future_layers_path = future_layers_path

# For Chelsa (processing 2050 separated from 2100 because of disk space)
future_layers_path <- list.files(path="/Users/josue/Dropbox/4Environmental_layers/future_chelsa",pattern='2011-2040',full.names = T) 
future_layers_path

# Obs - Perhaps a loop may be better at using the computer memory. I needed to have more than 30GB available for this function this time
get_future_bioclim = function(future_layers_path){
  #future_layers_path = future_layers_path[3]
print(gsub('/Users/josue/Dropbox/4Environmental_layers/future_chelsa/',"",future_layers_path))
  link = list.files(path=future_layers_path,pattern="bio",full.names = T)
  test <- stack(link)
  names(test)
  #plot(test[[1:6]])
  #plot(test[[5]])
  # plot(test[[14:19]])
  #proj4string(test) = oldproj
  test = mask(crop(test,selected_clipped_biomes_oldproj),selected_clipped_biomes_oldproj)
  test = raster::aggregate(test,15)
  test = projectRaster(test,crs = behrmannCRS)
  test = raster::resample(test,preds_modelling_mask[[1]]) # just to have the same resolution as the chelsa rasters
  #names(test) = c('bio_1','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9')
  names(test) = c('bio_15','bio_17','bio_18','bio_2','bio_4','bio_5')
  #plot(test)    
  return(test)
}
list_stack_future_bioclim = lapply(future_layers_path,get_future_bioclim)
gc()
#remove(list_stack_future_bioclim)
library(beepr)
beep()
## -----------------------------------------------------------------------------
summary(list_stack_future_bioclim)
#remove(list_stack_future_bioclim)

## -----------------------------------------------------------------------------
list_of_years = as.data.frame(gsub('/Users/josue/Dropbox/4Environmental_layers/future_chelsa/',"",future_layers_path))
colnames(list_of_years) = "model"

list_of_years

## -----------------------------------------------------------------------------
# Calculating the emsemble by median 
# Chelsa 2011-2040 ssp370 ## Parei aqui
list_of_years
model_stacks = list_stack_future_bioclim[c(1,3,5,7,9)]
list_of_layers = 1:6
# First, test the result
mean_future_fun = function(list_of_layers){ # year 1
  layerss = function(model_stacks){
    which_layer = model_stacks[[list_of_layers]] # will return the first layer of each variable
    return(which_layer)
  }
  which_variable1 = lapply(model_stacks,layerss)
  which_variable = calc(stack(which_variable1),median)
  return(which_variable)
}
climond_2040_A1B_list = lapply(list_of_layers,mean_future_fun)
climond_2040_A1B = stack(climond_2040_A1B_list)

names(climond_2040_A1B) = c('bio_15','bio_17','bio_18','bio_2','bio_4','bio_5') # chelsa
plot(climond_2040_A1B)


# Chelsa 2011-2040 ssp585
list_of_years
model_stacks = list_stack_future_bioclim[c(2,4,6,8,10)]
list_of_layers = 1:6

climond_2040_A2_list = lapply(list_of_layers,mean_future_fun)
climond_2040_A2 = stack(climond_2040_A2_list)

names(climond_2040_A2) = c('bio_15','bio_17','bio_18','bio_2','bio_4','bio_5') # chelsa

plot(climond_2040_A2)

current_test = raster::extract(preds_modelling_mask[["bio_4"]],classified_ind_6_pts)

test_AB1 = raster::extract(climond_2040_A1B[["bio_4"]],classified_ind_6_pts)

test_A2 = raster::extract(climond_2040_A2[["bio_4"]],classified_ind_6_pts)

# Temperature increase from 2011 - 2040
boxplot(current_test,test_AB1,test_A2)
range(test_A2)

test_AB1==test_A2

test_within_models1 = raster::extract(model_stacks[[1]][["bio_5"]],classified_ind_6_pts)
test_within_models2 = raster::extract(model_stacks[[2]][["bio_5"]],classified_ind_6_pts)
test_within_models3 = raster::extract(model_stacks[[3]][["bio_5"]],classified_ind_6_pts)
test_within_models4 = raster::extract(model_stacks[[4]][["bio_5"]],classified_ind_6_pts)
test_within_models5 = raster::extract(model_stacks[[5]][["bio_5"]],classified_ind_6_pts)

current_test = raster::extract(preds_modelling_mask[["bio_5"]],classified_ind_6_pts)

boxplot(current_test,test_within_models1,test_within_models2,test_within_models3,test_within_models4,test_within_models5)

## -----------------------------------------------------------------------------
# For 2050 A1B == Chelsa 2041-2070 ssp370
list_of_years
model_stacks = list_stack_future_bioclim[c(1,2,4)]
list_of_layers = 1:19

# First, test the result
mean_future_fun = function(list_of_layers){ # year 1
layerss = function(model_stacks){
which_layer = model_stacks[[list_of_layers]] # will return the first layer of each variable
return(which_layer)
}
which_variable1 = lapply(model_stacks,layerss)
which_variable = calc(stack(which_variable1),median)
return(which_variable)
}
climond_2050_A1B_list = lapply(list_of_layers,mean_future_fun)
climond_2050_A1B = stack(climond_2050_A1B_list)
#names(climond_2050_A1B) = c('bio_1','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19') # climond
names(climond_2050_A1B) = c('bio_1','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9') # chelsa

plot(climond_2050_A1B)

#For climond_2050_A2 = chelsa 2041-2070 ssp585
list_of_years
model_stacks = list_stack_future_bioclim#[c(3,5)]
list_of_layers = 1:19
climond_2050_A2_list = lapply(list_of_layers,mean_future_fun)
climond_2050_A2 = stack(climond_2050_A2_list)
#names(climond_2050_A2) = c('bio_1','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19')
names(climond_2050_A2) = c('bio_1','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9') 
climond_2050_A2

# Skipping chelsa 2100 for now

#For climond_2100_A1B
list_of_years
model_stacks = list_stack_future_bioclim[c(1,3,5)]
list_of_layers = 1:19
climond_2100_A1B_list = lapply(list_of_layers,mean_future_fun)
climond_2100_A1B = stack(climond_2100_A1B_list)
names(climond_2100_A1B) = c('bio_1','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9')
climond_2100_A1B

#For climond_2100_A2
list_of_years
model_stacks = list_stack_future_bioclim[c(2,4,6)]
list_of_layers = 1:19
climond_2100_A2_list = lapply(list_of_layers,mean_future_fun)
climond_2100_A2 = stack(climond_2100_A2_list)
names(climond_2100_A2) = c('bio_1','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9')
climond_2100_A2

par(mfrow=c(1,1))
# Testing selected variables
plot(which_preds[["bio_15"]])
plot(climond_2050_A1B[["bio_15"]]);
plot(climond_2050_A2[["bio_15"]]);
plot(climond_2050_A2[["bio_15"]]-climond_2050_A1B[["bio_15"]])

plot(climond_2100_A1B[["bio_15"]]);
plot(climond_2100_A2[["bio_15"]])
plot(climond_2100_A2[["bio_15"]]-climond_2100_A1B[["bio_15"]])

# Write layers for saving time of processing
#processed_climond_rasters = list(climond_2050_A1B,climond_2050_A2,climond_2100_A1B,climond_2100_A2)
#saveRDS(processed_climond_rasters, "/Users/josue/Dropbox/4Environmental_layers/processed_climond_rasters.rds")

# Read in Future projection processed climond rasters ####
#processed_climond_rasters = readRDS("/Users/josue/Dropbox/4Environmental_layers/processed_climond_rasters.rds")
#processed_climond_rasters[[1]]

# Write layers for saving time of processing
processed_chelsa_rasters = list(climond_2050_A1B,
                                climond_2050_A2,
                                climond_2100_A1B,
                                climond_2100_A2)

#saveRDS(processed_chelsa_rasters, "/Users/josue/Dropbox/4Environmental_layers/processed_chelsa_rasters.rds")

# Read in Future projection processed climond rasters ####
processed_chelsa_rasters = readRDS("/Users/josue/Dropbox/4Environmental_layers/processed_chelsa_rasters.rds")
# 2050 A1B = Chelsa 2041-2070 ssp370
climond_2050_A1B = processed_chelsa_rasters[[1]]
# 2050 A2 = chelsa 2041-2070 ssp585
climond_2050_A2 = processed_chelsa_rasters[[2]] 
# 2100 A1B = Chelsa 2070-2100 ssp370
climond_2100_A1B = processed_chelsa_rasters[[3]]
# 2100 A2 = chelsa 2070-2100 ssp585
climond_2100_A2 = processed_chelsa_rasters[[4]] 


#plot(processed_chelsa_2050_rasters[[1]])
## -----------------------------------------------------------------------------
# For the present
# Primeiro liste todos os rasters ambientais por nome. path= indica o diretório em que as variáveis estão

land_cover_current = raster("/Users/josue/Dropbox/4Environmental_layers/land_cover/MODISLandcover_2010_reclass_1km.tif")
#land_cover_current_b = projectRaster(land_cover_current, crs =behrmannCRS)
#land_cover_current_c = crop(land_cover_current, spTransform(selected_clipped_biomes_3,CRSobj = behrmannCRS))
#land_cover_current_c = aggregate(land_cover_current_c,10)
#land_cover_current_c = raster::resample(land_cover_current_c,preds_modelling_mask) # just to have the same resolution as the chelsa rasters
#remove(land_cover_current)
#remove(land_cover_current_b)
#land_cover_current_c
#plot(preds.climond_stck)


## -----------------------------------------------------------------------------
plot(land_cover_current,col=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33'))


## -----------------------------------------------------------------------------
# Increasing resolution of preds_modelling_maskto match land cover maps (which I will use in low res)
#preds_modelling_mask_d = disaggregate(preds_modelling_mask[[1]],24) 

#current_forest_cover_cl = raster::clump(current_forest_cover)
#raster_to_polygon
#plot(current_forest_cover_cl)
#rc = current_forest_cover_cl
#extract IDs of clumps according to some criteria
#clump9 = data.frame(freq(rc))
#clump9 = clump9[!clump9$count <= 12, ] #Remove clumps smaller than 12 cells, that is, half of the aggregate value
#clump9 = as.vector(clump9$value) # record IDs from clumps which met the criteria in previous step
#head(clump9)
#rc[rc %notin% clump9] = NA #replace cells with IDs which do not belong to the group of interest 
#plot(rc,col="black",legend=FALSE)

#rc_agr = raster::aggregate(rc,24, expand=FALSE)
#plot(rc_agr)
#rc_agr_proj = projectRaster(rc,crs = crs(preds_modelling_mask))

preds_modelling_mask_d = disaggregate(preds_modelling_mask[[1]],17)

preds_modelling_mask_d = raster(resolution=c(1000,1000),ext= extent( preds_modelling_mask[[1]]),crs=behrmannCRS)
values(preds_modelling_mask_d) = 1:ncell(preds_modelling_mask_d)
preds_modelling_mask_d = raster::mask(preds_modelling_mask_d, selected_clipped_biomes_2)

plot(preds_modelling_mask_d)

#current_forest_cover_agr_proj = raster::resample(rc_agr_proj,preds_modelling_mask_d)
#plot(rc_agr_proj)

#current_forest_cover_agr_proj[current_forest_cover_agr_proj > 0] <- 0
#current_forest_cover_agr_proj[is.na(current_forest_cover_agr_proj)] <- 1

#current_forest_cover_agr_proj = raster::mask(current_forest_cover_agr_proj,preds_modelling_mask_d[[1]])

#plot(current_forest_cover_agr_proj)


## -----------------------------------------------------------------------------
# List downloaded files (ex, world_A1B_2050, world_A1B_2100)
land_cover_future_list <- list.files(path="/Users/josue/Dropbox/4Environmental_layers/land_cover/",pattern='tif',full.names = T) 
land_cover_future_list

## -----------------------------------------------------------------------------
get_future_land_cover = function(land_cover_future_list){

# Read raster
land_cover = raster(land_cover_future_list)
# Crop to study area (you do not need sptranform if not changing it)
land_cover_c = mask(crop(land_cover,selected_clipped_biomes_3),selected_clipped_biomes_3)
extent(land_cover_c) <- raster::extent(selected_clipped_biomes_3)
#land_cover_c = projectRaster(land_cover_c, crs =behrmannCRS)
print("cropping")
land_cover_c[land_cover_c < 2] <- NA
land_cover_c[land_cover_c >= 3] <- NA
land_cover_c[land_cover_c == 3] <- NA
land_cover_c[land_cover_c >= 2] <- 1
print("aggregating")
#land_cover_c = aggregate(land_cover_c,10)
#land_cover_c = raster::resample(land_cover_c,preds_modelling_mask) # just to have the same resolution as the chelsa rasters
plot(land_cover_c)
return(land_cover_c)
}
list_rasters_future_land_cover = pblapply(land_cover_future_list,get_future_land_cover,cl = 1)


## -----------------------------------------------------------------------------
rasters_future_land_cover= stack(list_rasters_future_land_cover)

plot(list_rasters_future_land_cover[[1]],col= terrain.colors(255))
rasters_future_land_cover


# Saving processed files to avoid using them again
saveRDS(rasters_future_land_cover,"/Users/josue/Dropbox/4Environmental_layers/rasters_future_land_cover.rds")

#rasters_future_land_cover = readRDS("/Users/josue/Dropbox/4Environmental_layers/rasters_future_land_cover.rds")

## -----------------------------------------------------------------------------
par(bg=NA)
plot(rasters_future_land_cover[[1]],col= "#005a32",legend=FALSE,axes = FALSE,box = FALSE)
plot(SAshp2,col="white",border="black", add=TRUE,lwd=1)
plot(rasters_future_land_cover[[2]],col= "#005a32",legend=FALSE,axes = FALSE,box = FALSE)
plot(SAshp2,col=NA,border="#969696", add=TRUE,lwd=1)
plot(rasters_future_land_cover[[3]],col= "#005a32",legend=FALSE,axes = FALSE,box = FALSE)
plot(SAshp2,col=NA,border="#969696", add=TRUE,lwd=1)
plot(rasters_future_land_cover[[4]],col= "#005a32",legend=FALSE,axes = FALSE,box = FALSE)
plot(SAshp2,col=NA,border="#969696", add=TRUE,lwd=1)
plot(rasters_future_land_cover[[5]],col= "#005a32",legend=FALSE,axes = FALSE,box = FALSE)
plot(SAshp2,col=NA,border="#969696", add=TRUE,lwd=1)


## -----------------------------------------------------------------------------
# Increasing resolution to match modelling layers
list_rasters_future_land_cover = list(rasters_future_land_cover[[1]],
                                      rasters_future_land_cover[[2]],
                                      rasters_future_land_cover[[3]],
                                      rasters_future_land_cover[[4]],
                                      rasters_future_land_cover[[5]])

# Function to change raster resolution
increase_future_lc_resol = function(list_rasters_future_land_cover){
  print("clumping")
  
rasters_future_land_cover_cl = raster::clump(list_rasters_future_land_cover)
#raster_to_polygon
#plot(current_forest_cover_cl)
rc = rasters_future_land_cover_cl
#extract IDs of clumps according to some criteria
#clump9 = data.frame(freq(rc))
#clump9 = clump9[!clump9$count <= 12, ] #Remove clumps smaller than 12 cells, that is, half of the aggregate value
#clump9 = as.vector(clump9$value) # record IDs from clumps which met the criteria in previous step
#head(clump9)
#rc[rc %notin% clump9] = NA #replace cells with IDs which do not belong to the group of interest 
#plot(rc,col="black",legend=FALSE)
#rc_agr = raster::aggregate(rc,24, expand=FALSE)
#plot(rc_agr)
 print("projecting")
rc_agr_proj = projectRaster(rc,crs = crs(preds_modelling_mask))
 print("resampling")
rasters_future_land_cover_agr_proj = raster::resample(rc_agr_proj,preds_modelling_mask_d)

  print("transforming in barrier")
rasters_future_land_cover_agr_proj[rasters_future_land_cover_agr_proj > 0] <- 0
rasters_future_land_cover_agr_proj[is.na(rasters_future_land_cover_agr_proj)] <- 1

return(rasters_future_land_cover_agr_proj)
}

list_rasters_future_land_cover_ag_proj = pblapply(list_rasters_future_land_cover,increase_future_lc_resol,cl = 1)

rasters_future_land_cover_agr_proj = stack(list_rasters_future_land_cover_ag_proj)

#plot(rasters_future_land_cover_agr_proj)

#beepr::beep()
## -----------------------------------------------------------------------------


# Masking and saving the raster to the directory to avoid loosing temp files
rasters_future_land_cover_agr_proj = raster::mask(rasters_future_land_cover_agr_proj,preds_modelling_mask_d[[1]],filename="/Users/josue/Dropbox/4Environmental_layers/rasters_future_land_cover_agr_proj.grd",overwrite=TRUE)

# Get barrier rasters
#rasters_future_land_cover_agr_proj = raster::stack("/Users/josue/Dropbox/4Environmental_layers/rasters_future_land_cover_agr_proj.grd")

names(rasters_future_land_cover_agr_proj) = c("Current","A1B_2050","A1B_2100","A2_2050","A2_2100")

plot(rasters_future_land_cover_agr_proj[[1]])



## ----------------------------------------------------------------------------- #################################################

# Now the modeling really begins

## ----------------------------------------------------------------------------- #################################################

head(classified_ind)
population_pts = classified_ind[,c("kind","long","lat","long3_j","lat3_j")]

# Selecting population to model
population_pts= population_pts[population_pts$kind %in% c('hot_dry_seazonal'),]

# Check selection
unique(population_pts$kind)

nrow(population_pts)

# One point per grid ~ 30, then I added 10km for the thinning
# Thin data
spp_coords_thinned <- spThin::thin( loc.data = population_pts,
                            lat.col = "lat3_j", long.col = "long3_j",
                            spec.col = "kind",
                            thin.par = 40, reps = 1,
                            locs.thinned.list.return = TRUE,
                            write.files = FALSE,
                            write.log.file = FALSE)

spp_coords_thinned = as.data.frame(spp_coords_thinned)
head(spp_coords_thinned)

spp_coords_thinned$kind = 'hot_dry_seazonal'
spp_coords_thinned = spp_coords_thinned[,c(3,1,2)]
colnames(spp_coords_thinned) = c("kind", "long", "lat")
head(spp_coords_thinned)

# Changing to Behrman
coordinates(spp_coords_thinned) = spp_coords_thinned[,c("long","lat")]
proj4string(spp_coords_thinned) = oldproj
population_pts = spTransform(spp_coords_thinned,CRSobj = behrmannCRS)
population_pts@data$long = population_pts@coords[,1]
population_pts@data$lat = population_pts@coords[,2]

# Clean the data to one record per grid cell
bck.na<-preds_modelling_mask[[1]]
bck.na[]<-NA
r<-rasterize(population_pts@data[,c('long', 'lat')],bck.na,fun='count')
occ_data_clean_ok.pa<-rasterToPoints(r,fun=function(x){x>0}, spatial=T) 
population_pts = as.data.frame(occ_data_clean_ok.pa@coords)
population_pts$kind = 'hot_dry_seazonal'
population_pts = population_pts[,c(3,1,2)]
colnames(population_pts) = c("kind", "long", "lat")
coordinates(population_pts) = population_pts[,c("long","lat")]
proj4string(population_pts) = behrmannCRS

population_pts@data <- population_pts@data[,"kind",drop=FALSE] ## Necessary to not consider latitude or longitude as predictors

population_pts_dry_ses = population_pts

## -----------------------------------------------------------------------------
#dir.create("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_pop_modelling")
path_modelling = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_pop_modelling/"

## -----------------------------------------------------------------------------
unique(population_pts@data$kind)
#head(population_pts_all@data$kind)
population_pts@data$Occurrence = 1
#head(population_pts_all@data)
population_pts@data <- population_pts@data[,"Occurrence",drop=FALSE]
population_name= "hot_dry_seazonal"
plot(selected_clipped_biomes_2)
plot(population_pts,add=TRUE)

## -----------------------------------------------------------------------------
# How many unique records intersect with the raster?
print("total number of records for the population")
nrow(population_pts)
n_of_locs= extract(preds_modelling_mask[[1]],population_pts,cellnumbers=TRUE)
n_of_locs=length(unique(n_of_locs[,1]))
print("number of unique points given cell resolution")
n_of_locs
print("proportion of points per variable for the small models (min 5, best>10)")
(n_of_locs/2)-1 


## -----------------------------------------------------------------------------
# Checking variable importance for each modell
var_comb=combn(x=names(preds_modelling_mask), 2, simplify = FALSE)
length(var_comb)


## -----------------------------------------------------------------------------
# function for two by two models
#two_model_fun = function(var_comb){
#  species_sdm <- sdmData(Occurrence~.,train=population_pts,predictors=subset(preds_modelling_mask, var_comb),bg=list(n=n_of_locs*100,method='gRandom',remove=TRUE))
  #getmethodNames() # Check the names of all models
  # First round to test variable importance
  #model_var_sel <- sdm(Occurrence~.,species_sdm,methods = c("rf"), replication=c('cv'), cv.folds=5,n=10) ## Here I can only choose one species each time
#  model_var_sel <- sdm(Occurrence~.,species_sdm,methods = c("rf"),replication=c('boot'),n=10) ## Here I can only choose one species each time
  ## Variable importance calcs ####
#  model_info = getModelInfo(model_var_sel) 
#  modelIDs = model_info$modelID
#  df1=getVarImp(model_var_sel,id=modelIDs,wtest='training')
#  df2=df1@varImportanceMean$AUCtest
#  rownames(df2)= NULL
#return(df2)
#}
#var_imp_list= pbapply::pblapply(var_comb,two_model_fun,cl=6) # pbapply doesn't work on Jupyter Notebooks


## -----------------------------------------------------------------------------
#var_imp_list_bind = bind_rows(var_imp_list)
#var_imp_rf = aggregate(var_imp_list_bind$AUCtest, by=list(Category=var_imp_list_bind$variables), FUN=sum)
#colnames(var_imp_rf)=c("variables", "AUCtest")
#rownames(var_imp_rf)=var_imp_rf$variables
#var_importance_rf=var_imp_rf[order(var_imp_rf$AUCtest),]
#head(var_imp_rf)
#write.csv(var_importance_rf, paste0(path_modelling,population_name,"_var_importance.csv"))


## -----------------------------------------------------------------------------
## If you were to select the variables from the var importance list. Skip it here
#var_importance_rf= var_importance_rf[!var_importance_rf$variables=="elevat",]
#for (i in max_env_var_number:length(names(preds.mask))){
#  sel_best_var_rf= tail(rownames(var_importance_rf),i)
#  preds.no.colinear.rf <- subset(preds.mask,sel_best_var_rf)
#  v_rf = vifstep(preds.no.colinear.rf) 
#  preds.no.colinear.rf.final <- dropLayer(preds.no.colinear.rf,v_rf@excluded)
#  k= length(names(preds.no.colinear.rf.final))
#  if (k ==max_env_var_number) break
#}
##########


## -----------------------------------------------------------------------------
## Now the ensembling of small models###############
#length(var_comb)
# function for two by two models
###################################################
population_name= "hot_dry_seazonal"
var_comb
which_preds=preds_modelling_mask

# If the idea was to model around the records only
#library(adehabitatHR)
#bgExt <- mcp(classified_ind_6_pts, percent = 100) 
#plot(which_preds[[1]])
#plot(classified_ind_6_pts,add=TRUE)
#plot(bgExt, add=T)
#bgExt_buff <- rgeos::gBuffer(bgExt, width = 350000)
#plot(bgExt_buff,border="red", add=T)
#which_preds = raster::mask(which_preds,bgExt_buff)


## gerando pontos de pseudo-ausencias for evalutation
  #psudo <- sampleRandom(which_preds[[1]],n_of_locs*1,xy=TRUE,sp=TRUE)
set.seed(999)  
psudo <- sampleRandom(which_preds[[1]],1000,xy=TRUE,sp=TRUE)
  psudo@data$Occurrence <- 0 # retemos apenas a coluna Occurrence
  psudo@data <- psudo@data[,"Occurrence",drop=FALSE]
  buffer_pts <- gBuffer(population_pts,width=33000) 
  proj4string(buffer_pts) <- proj4string(psudo) # igualando projecao de buffer_pts e de psudo
  pts_outside <- psudo[is.na(over(psudo,buffer_pts)),] 
  pts_outside@data$Occurrence <- 0 
  sp11=population_pts
  sp11@data$Occurrence <- 1 
  sp11@data <- sp11@data[,"Occurrence",drop=FALSE]
  #pts_outside@data <- pts_outside@data[,"Occurrence",drop=F]
  proj4string(sp11) <- proj4string(pts_outside) # igualando projecao de buffer_pts e de psudo
  speciesN <- bind(sp11,pts_outside) # reunindo presencas e pseudoause


## -----------------------------------------------------------------------------
#getwd()
#dir.create(paste0(path,"/sdm_models/"))
#dir.exists(paste0(path,"/sdm_models/"))
population_name= "hot_dry_seazonal"
paste0(path,"/sdm_models/",population_name,"_",paste0(unlist(var_comb[1]),collapse = '_'),"_model",".sdm") # test if the path is correct


## -----------------------------------------------------------------------------

# For ensembling of small models:
var_comb=combn(x=names(preds_modelling_mask), 2, simplify = FALSE)
length(var_comb)

# For traditional modelling
var_comb = list(names(preds_modelling_mask))
length(var_comb)

# I will add a condition for ESM (ensembling of small models x normal modelling)

setwd(paste0(path,"/sdm_models/"))
two_model_fun_final = function(var_comb){
  
#var_comb = var_comb[[1]]
  
  if (length(var_comb)>1) {
    print(paste("running", var_comb))
    species_sdm <- sdmData(Occurrence~.,train=sp11,predictors=subset(which_preds, var_comb),bg=list(n=1000,method='gRandom',remove=TRUE))
    
    replication=c('cv')
    cv.folds=sum(speciesN$Occurrence==1) # equivalent to jackniffe
    n=1
    print("Running Ensemble of Small Models")
    
  } else {
    lista=c(1:2)
    print("Running model with all variables:")
    many_pa_model = function(lista){
    print(paste("pa",lista))  
    species_sdm <- sdmData(Occurrence~.,train=sp11,predictors=subset(which_preds, unique(unlist(var_comb))), bg=list(n=nrow(sp11)*2,method='gRandom',remove=TRUE))
    
    replication=c('cv')
    cv.folds=5 ######## Traditional CV
    n=5
    print(paste(unique(unlist(var_comb))))
    
    invisible(capture.output(model_final <- sdm(Occurrence~.,species_sdm,methods = c("domain.dismo","glm","maxent"), replication=replication, cv.folds=cv.folds,n=n))) ## 
    return(model_final)
  }
    model_final_list = lapply(lista, many_pa_model) 
    model_final=Reduce('+',model_final_list)
  }

#getModelId(model_final)


  evaluation = getEvaluation(model_final)
  good_models = evaluation[evaluation$AUC>0.75 & evaluation$TSS>0.5,]$modelID
  #getModelInfo(model_final)
  model_final = subset(model_final, good_models, drop=TRUE)
  #getModelInfo(model_final)
  
  sdm::write.sdm(model_final, paste0(path,"/sdm_models/",population_name,"_",paste0(unlist(var_comb),collapse = '_'),"_model",".sdm"),overwrite = TRUE)
  
  print("predicting")
  invisible(capture.output(presente <- predict(model_final,subset(which_preds, unlist(var_comb)),filename=paste0(path_modelling, population_name,"_current_pred.grd"),mean=FALSE, overwrite=TRUE)))
  
  # ensemble based on a Weighted averaging that is weighted using TSS statistic with threshold criterion 
  
  #presente= raster::mask(presente,which_preds[[1]]) # This is necessary to remove predictions outside the range
  
  print("rescaling predictions")
  #plot(presente)
  # Funtion to rescale sdm predictions
  x=1:nlayers(presente)
  norm_fun <- function(x) {
    y=presente[[x]]/abs(max(na.exclude(presente[[x]][1:ncell(presente[[x]])])))
    y[1:ncell(y)] = BBmisc::normalize(y[1:ncell(y)], method ="range",range = c(0, 1))
    return(y)
  }
  present_norm_list= lapply(x,norm_fun)
  present_norm=stack(present_norm_list)
  #plot(present_norm)

  # In case I do not change the scale of suitability values
  #present_norm = presente
  
  # Creating my own esenbling function
  # reunindo presencas e pseudoause
  #obs <- speciesN$Occurrence
  #pred <- raster::extract(present_norm,speciesN) # que extrai todos os valores de adequabilidade dos pontos de speciesN
  ## e enfim o ev, que mostra limiares (thresholds) de varias regras diferentes
  # Get evaluations by AUC from the mean models
#  x=1:ncol(pred)
#  ev_fun <- function(x) {
#    ev <- evaluates(obs,pred[,x])
#    ev2 = ev@statistics$AUC
#    return(ev2)
#  }
#  ev_all= lapply(x,ev_fun)
#  AUC_weights = simplify2array(ev_all)

  # Or get it from SDM
  AUC_weights = getEvaluation(model_final)[,2]

  # Calculating weighted mean from AUC
  present_ensemble = weighted.mean(x=present_norm, w=AUC_weights)
  #quartz(12,12)
  #plot(present_ensemble)
  print(paste("ensembling", paste(unlist(var_comb),collapse = "_")))
  
  # Or
  #invisible(capture.output(present_ensemble <- ensemble(model_final,newdata=subset(which_preds, var_comb),filename=paste0(path_modelling, population_name,"_",gsub(":", "-", Sys.time()),"_current_pred.grd"),setting=list(method='weighted',stat='TSS',opt=2),overwrite=TRUE)))
  
  present_ensemble= raster::mask(present_ensemble,which_preds[[1]]) # This is necessary to remove predictions outside the range
  
  plot(present_ensemble)
  #plot(population_pts,add=TRUE)
  return(tryCatch(present_ensemble, error=function(e) NULL))
}

# For small models, use pblapply
#var_imp_list= pblapply(var_comb,two_model_fun_final,cl=1)
# Full model
var_imp_list = two_model_fun_final(var_comb)
plot(var_imp_list)
#beepr::beep()

## -----------------------------------------------------------------------------
# Check if you have generated all the models for the combination of variables (ESM only)
#summary(var_imp_list)

## -----------------------------------------------------------------------------
# From this list, produce the final ENMs models by emsembling them again!
if (length(var_comb)>1) {
present_all_var =stack(var_imp_list)
  plot(present_all_var)
  # Creating my own esenbling function
  obs <- speciesN$Occurrence
  pred <- raster::extract(present_all_var,speciesN) # que extrai todos os valores de adequabilidade dos pontos de speciesN
  ## e enfim o ev, que mostra limiares (thresholds) de varias regras diferentes
  # Get evaluations by AUC from the mean models
  x=1:ncol(pred)
  ev_fun <- function(x) {
    ev <- evaluates(obs,pred[,x])
    ev1 = ev@statistics$AUC
    ev2 = median(ev@threshold_based$TSS)
    ev3 <- ev@threshold_based$threshold[2]
    ev_all = c(ev1,ev2,ev3)  
    return(ev_all)
  }
  ev_all= lapply(x,ev_fun)
  AUC_weights_all_var = simplify2array(ev_all)[1,]
  TSS_weights_all_var = simplify2array(ev_all)[2,]
  sensi_specificity = simplify2array(ev_all)[3,]
  
  # Calculating weighted mean from AUC
present_ensemble_all_var_wer = weighted.mean(x=present_all_var, w=TSS_weights_all_var) 
plot(present_ensemble_all_var_wer)  
plot(population_pts, add=TRUE)
} else {
  present_ensemble_all_var_wer = var_imp_list
}

# Alternative colour
#viridis::turbo(100,direction = 1,alpha = .8)

# Testing colors
# Classic palette BuPu, with 4 colors
coul <- rev(brewer.pal(11, "RdYlBu"))
# Add more colors to this palette :
coul <- colorRampPalette(coul,alpha=.8)(10)


plot(present_ensemble_all_var_wer,col = coul, main="Current",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)

writeRaster(present_ensemble_all_var_wer,
            paste0(path_modelling,"1output_rasters/",
                   "present_ensemble_all_var_wer",
                   ".tif"))

## -----------------------------------------------------------------------------
median(sensi_specificity)

## -----------------------------------------------------------------------------
model.stack.test=present_ensemble_all_var_wer
sp1=population_pts

library(gdistance)
model_cost.ras  <- -1 * log(model.stack.test)    
model_trans.ras <- 1 / model_cost.ras  
trans     <- transition(model_trans.ras, transitionFunction=mean, directions = 8)
trans     <- geoCorrection(trans, type="c", multpl=F)
lin_dist.ras = accCost(trans, sp1)
thisLineagePoints.ras <- rasterize(sp1, model.stack.test, 1)
lin_dist.ras <- raster::distance(thisLineagePoints.ras) / 1000   
plot(lin_dist.ras)
lin_dist.ras[lin_dist.ras < 100 ] <- 100
plot(lin_dist.ras)
lin_weight.ras <- 1 / (lin_dist.ras) # Using inverse of distance after 200km
lin_weight.ras = lin_weight.ras/max(na.omit(lin_weight.ras[1:ncell(lin_weight.ras)]))
plot(lin_weight.ras)
model.stack <- lin_weight.ras * model.stack.test
plot(model.stack)
model.stack.test = model.stack
#
obs <- speciesN$Occurrence
pred <- raster::extract(model.stack,speciesN) # created from pseudoabsence generator
e_sp1 <- evaluates(obs,pred)

#tr_sp1 = median(sensi_specificity)
tr_sp1 <- e_sp1@threshold_based$threshold[2]
tr_sp1 = min(na.omit(raster::extract(model.stack,classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==2,])))
tr_sp1=0.3
#binary_code_sp1 <- c(NA, NA, 0, -10, tr_sp1, 0,tr_sp1, 10, 1)
#matrix_bin_sp1 <- matrix(binary_code_sp1, ncol = 3, byrow = TRUE)
pr.pa2 <- raster(model.stack.test) # creating an empty raster
pr.pa2[] <- ifelse(model.stack.test[] >= tr_sp1,1,0) ## Esse aqui È seu model bin·rio agora

#pr.pa2 = raster(paste0(path_modelling,"1output_rasters/","pr.pa2",".tif"))

# Plotting
plot(pr.pa2, col=c("#f0f0f0","#969696"), axes=FALSE, box=FALSE, legend=FALSE, main="Current distribution")
plot(SAshp2,add=TRUE)
plot(population_pts_dry_ses,pch=21,col="black",bg="white", add=TRUE,cex=1)
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==2,],pch = 21, cex=1.5, col="black", bg="#377eb8",lwd=1.5, add=TRUE)
#plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==1,], pch = 21, cex=1.5, col="black", bg="#e41a1c",lwd=1.5, add=TRUE)


writeRaster(pr.pa2,
            paste0(path_modelling,"1output_rasters/",
                   "pr.pa2",
                   ".tif"))

# Add location of molecular data
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==2,], add=TRUE,pch=16,cex=.5,col="white")

# Testing the match > 50%
present_ensemble_all_var_wer_test = present_ensemble_all_var_wer
present_ensemble_all_var_wer_test[present_ensemble_all_var_wer_test>pr.pa2]=0
plot(present_ensemble_all_var_wer_test)
plot(population_pts, add=TRUE)

## -----------------------------------------------------------------------------
# Getting the path to the variables
fitted_models_path <- list.files(path=paste0(path,"/sdm_models"),pattern=population_name,full.names = T) 
# Selecting which raster to predict to
which_preds = climond_2050_A1B

# If you log present varibles and don't do that below, you'll not not predict anything
#which_preds$bio_18 = log(which_preds$bio_18)
#which_preds$bio_17 = log(which_preds$bio_17)

## -----------------------------------------------------------------------------
pred_two_model_fun = function(var_comb){

# Reading in the models
if (length(var_comb)>1) {
g1 = grep(var_comb[1], fitted_models_path,value = TRUE) 
g2 = grep(var_comb[2], g1,value = TRUE)
print(paste0(g2))
} else {
g2 = grep(paste(unlist(var_comb), collapse = '_'), fitted_models_path,value = TRUE)
print(paste0(g2))
}

model_final = sdm::read.sdm(g2)
    print("predicting")
invisible(capture.output(presente <- predict(model_final,subset(which_preds, unlist(var_comb)),mean=FALSE)))
print("rescaling")
# Funtion to rescale sdm predictions
  x=1:nlayers(presente)
  norm_fun <- function(x) {
    y=presente[[x]]/abs(max(na.exclude(presente[[x]][1:ncell(presente[[x]])])))
    y[1:ncell(y)] = BBmisc::normalize(y[1:ncell(y)], method ="range",range = c(0, 1))
    return(y)
  }
  present_norm_list= lapply(x,norm_fun)
  present_norm=stack(present_norm_list)

to_remove = which(1:nlayers(presente) %notin% getEvaluation(model_final)[,1]) # some models were not evaluated. Considering the number of models, I will subset those

if (length(to_remove)>0) {
present_norm2 = raster::dropLayer(present_norm,to_remove)
} else {
present_norm2 = present_norm
}
    
AUC_weights = getEvaluation(model_final)[,2]

  # Calculating weighted mean from AUC
  present_ensemble = weighted.mean(x=present_norm2, w=AUC_weights)
  
  present_ensemble= raster::mask(present_ensemble,which_preds[[1]]) # This is necessary to remove predictions outside the range
  
  return(tryCatch(present_ensemble, error=function(e) NULL))    
}    

# For 2050 A1B seasonally dry
# If ensemble of small models
which_preds = climond_2050_A1B
if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2050_A1B = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2050_A1B =stack(var_imp_list_2050_A1B) # rem A1b 
ensemble_all_var_2050_A1B = weighted.mean(x=all_var_2050_A1B, w=TSS_weights_all_var)
# If Normal modelling
} else {
ensemble_all_var_2050_A1B = pred_two_model_fun(var_comb)
}


plot(ensemble_all_var_2050_A1B,col = coul, main="2070 ssp370 (high CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)

writeRaster(ensemble_all_var_2050_A1B,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2050_A1B",
                   ".tif"))

# Testing the match > 50%
present_ensemble_all_var_wer_test = ensemble_all_var_2050_A1B
present_ensemble_all_var_wer_test[present_ensemble_all_var_wer_test<0.3]=0
plot(present_ensemble_all_var_wer_test,col = coul)
plot(population_pts, add=TRUE)


## ----------------------------------------------------------------------------- For future 2050 A2
which_preds = climond_2050_A2
if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2050_A2 = pblapply(var_comb,pred_two_model_fun,cl=1)# ESM
all_var_2050_A2 =stack(var_imp_list_2050_A2) #
ensemble_all_var_2050_A2 = weighted.mean(x=all_var_2050_A2, w=TSS_weights_all_var)
} else {
ensemble_all_var_2050_A2 = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2050_A2,col = coul, main="2070 ssp585 (Extreme CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)

writeRaster(ensemble_all_var_2050_A2,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2050_A2",
                   ".tif"))

# Testing the match > 50%
present_ensemble_all_var_wer_test = ensemble_all_var_2050_A2
present_ensemble_all_var_wer_test[present_ensemble_all_var_wer_test<0.3]=0
plot(present_ensemble_all_var_wer_test)
plot(population_pts, add=TRUE)

## ----------------------------------------------------------------------------- For future 2100 A1B
which_preds = climond_2100_A1B

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2100_A1B = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2100_A1B =stack(var_imp_list_2100_A1B) # rem A1b
ensemble_all_var_2100_A1B = weighted.mean(x=all_var_2100_A1B, w=TSS_weights_all_var)
} else {  
ensemble_all_var_2100_A1B = pred_two_model_fun(var_comb)
}

#plot(ensemble_all_var_2100_A1B,col = viridis::turbo(100,direction = 1,alpha = .8), main="2100 ssp370 (high CO2)",axes=FALSE, box=FALSE)

plot(ensemble_all_var_2100_A1B,col = coul, main="2100 ssp370 (high CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)

writeRaster(ensemble_all_var_2100_A1B,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2100_A1B",
                   ".tif"))

# Testing the match > 50%
present_ensemble_all_var_wer_test = ensemble_all_var_2100_A1B
present_ensemble_all_var_wer_test[present_ensemble_all_var_wer_test<0.1]=0
plot(present_ensemble_all_var_wer_test)
plot(population_pts, add=TRUE)


## ----------------------------------------------------------------------------- For future 2100 A2
which_preds = climond_2100_A2

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2100_A2 = pblapply(var_comb,pred_two_model_fun,cl=1)
ensemble_all_var_2100_A2 = weighted.mean(x=all_var_2100_A2, w=TSS_weights_all_var)
} else {
ensemble_all_var_2100_A2 = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2100_A2)
plot(population_pts, add=TRUE)

plot(ensemble_all_var_2100_A2,col = coul, main="2100 ssp585 (extreme CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2100_A2,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2100_A2",
                   ".tif"),overwrite=TRUE)



# Testing the match > 50%
present_ensemble_all_var_wer_test = ensemble_all_var_2100_A2
present_ensemble_all_var_wer_test[present_ensemble_all_var_wer_test<0.1]=0
plot(present_ensemble_all_var_wer_test)
plot(population_pts, add=TRUE)


## ----------------------------------------------------------------------------- For future 2040 A1B
which_preds = climond_2040_A1B
if (length(var_comb)>1) { # If ensemble of small models
  var_imp_list_2040_A1B = pblapply(var_comb,pred_two_model_fun,cl=1)
  ensemble_all_var_2040_A1B = weighted.mean(x=all_var_2040_A1B, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2040_A1B = pred_two_model_fun(var_comb) # SDM
}
#plot(ensemble_all_var_2040_A1B)
#plot(population_pts, add=TRUE)
plot(ensemble_all_var_2040_A1B,col = coul, main="2040 ssp370 (high CO2)",axes=FALSE, box=FALSE)
plot(SAshp2,add=TRUE)

writeRaster(ensemble_all_var_2040_A1B,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2040_A1B",
                   ".tif"),overwrite=TRUE)

## ----------------------------------------------------------------------------- For future 2040 A2
which_preds = climond_2040_A2
if (length(var_comb)>1) { # If ensemble of small models
  var_imp_list_2040_A2 = pblapply(var_comb,pred_two_model_fun,cl=1)
  ensemble_all_var_2040_A2 = weighted.mean(x=all_var_2040_A2, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2040_A2 = pred_two_model_fun(var_comb) # SDM
}
#plot(ensemble_all_var_2040_A2)
#plot(population_pts, add=TRUE)
plot(ensemble_all_var_2040_A2,col = coul, main="2040 ssp585 (extreme CO2)",axes=FALSE, box=FALSE)
# Testing with rcolour brewer
plot(SAshp2, add=TRUE)

writeRaster(ensemble_all_var_2040_A2,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2040_A2",
                   ".tif"),overwrite=TRUE)

## -----------------------------------------------------------------------------warm_wet_stable 
#################################################

head(classified_ind)
population_pts = classified_ind[,c("kind","long","lat","long3_j","lat3_j")]

# Selecting population to model
population_pts= population_pts[population_pts$kind %in% c('warm_wet_stable'),]

# Check selection
unique(population_pts$kind)

nrow(population_pts)

# One point per grid ~ 30, then I added 10km for the thinning
# Thin data
spp_coords_thinned <- spThin::thin( loc.data = population_pts,
                                    lat.col = "lat3_j", long.col = "long3_j",
                                    spec.col = "kind",
                                    thin.par = 40, reps = 1,
                                    locs.thinned.list.return = TRUE,
                                    write.files = FALSE,
                                    write.log.file = FALSE)

spp_coords_thinned = as.data.frame(spp_coords_thinned)
head(spp_coords_thinned)

spp_coords_thinned$kind = 'warm_wet_stable'
spp_coords_thinned = spp_coords_thinned[,c(3,1,2)]
colnames(spp_coords_thinned) = c("kind", "long", "lat")
head(spp_coords_thinned)

# Changing to Behrman
coordinates(spp_coords_thinned) = spp_coords_thinned[,c("long","lat")]
proj4string(spp_coords_thinned) = oldproj
population_pts = spTransform(spp_coords_thinned,CRSobj = behrmannCRS)
population_pts@data$long = population_pts@coords[,1]
population_pts@data$lat = population_pts@coords[,2]

# Clean the data to one record per grid cell
bck.na<-preds_modelling_mask[[1]]
bck.na[]<-NA
r<-rasterize(population_pts@data[,c('long', 'lat')],bck.na,fun='count')
occ_data_clean_ok.pa<-rasterToPoints(r,fun=function(x){x>0}, spatial=T) 
population_pts = as.data.frame(occ_data_clean_ok.pa@coords)
population_pts$kind = 'warm_wet_stable'
population_pts = population_pts[,c(3,1,2)]
colnames(population_pts) = c("kind", "long", "lat")
coordinates(population_pts) = population_pts[,c("long","lat")]
proj4string(population_pts) = behrmannCRS

population_pts@data <- population_pts@data[,"kind",drop=FALSE] ## Necessary to not consider latitude or longitude as predictors

## -----------------------------------------------------------------------------
unique(population_pts@data$kind)

## warm_wet_stable -----------------------------------------------------------------------------
#dir.create("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_pop_modelling")
path_modelling = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_pop_modelling/"


## -----------------------------------------------------------------------------
unique(population_pts@data$kind)
#head(population_pts_all@data$kind)
population_pts@data$Occurrence = 1
#head(population_pts_all@data)
population_pts@data <- population_pts@data[,"Occurrence",drop=FALSE]
population_name= "warm_wet_stable"

population_pts_warm_wet = population_pts

## -----------------------------------------------------------------------------
# How many unique records intersect with the raster?
print("total number of records for the population")
nrow(population_pts)
n_of_locs= extract(preds_modelling_mask[[1]],population_pts,cellnumbers=TRUE)
n_of_locs=length(unique(n_of_locs[,1]))
print("number of unique points given cell resolution")
n_of_locs
print("proportion of points per variable for the small models (min 5, best>10)")
(n_of_locs/2)-1 

## -----------------------------------------------------------------------------
## Now the ensembling of small models###############
#length(var_comb)
# function for two by two models
###################################################
population_name= "warm_wet_stable"
# For ensembling of small models:
which_preds = preds_modelling_mask

var_comb=combn(x=names(preds_modelling_mask), 2, simplify = FALSE)
length(var_comb)
# For traditional modelling
var_comb = list(names(preds_modelling_mask))
length(var_comb)


## gerando pontos de pseudo-ausencias for evalutation
#set.seed(999)
  psudo <- sampleRandom(preds_modelling_mask[[1]],1000,xy=TRUE,sp=TRUE) 
  #plot(psudo)
  psudo@data$Occurrence <- 0 # retemos apenas a coluna Occurrence
  psudo@data <- psudo@data[,"Occurrence",drop=FALSE]
  buffer_pts <- gBuffer(population_pts,width=33000) 
  proj4string(buffer_pts) <- proj4string(psudo) # igualando projecao de buffer_pts e de psudo
  pts_outside <- psudo[is.na(over(psudo,buffer_pts)),] 
  pts_outside@data$Occurrence <- 0 
  sp11=population_pts
  sp11@data$Occurrence <- 1 
  sp11@data <- sp11@data[,"Occurrence",drop=FALSE]
  #pts_outside@data <- pts_outside@data[,"Occurrence",drop=F]
  proj4string(sp11) <- proj4string(pts_outside) # igualando projecao de buffer_pts e de psudo
  speciesN <- bind(sp11,pts_outside) # reunindo presencas e pseudoause
  plot(speciesN[speciesN$Occurrence==0,], col="red",pch=16)
  plot(sp11, col="blue",pch=16,add=TRUE)
plot(selected_clipped_biomes_2,add=TRUE)

## -----------------------------------------------------------------------------
#getwd()
dir.exists(paste0(path,"/sdm_models/"))
population_name= "warm_wet_stable"
paste0(path,"/sdm_models/",population_name,"_",paste0(unlist(var_comb[1]),collapse = '_'),"_model",".sdm") # test if the path is correct


## -----------------------------------------------------------------------------
# Modelling
if (length(var_comb)>1) { # ESM
var_imp_list_warm_wet= pblapply(var_comb,two_model_fun_final,cl=1)
} else { # SDM
  var_imp_list_warm_wet= two_model_fun_final(var_comb)  
}
plot(var_imp_list_warm_wet)

if (length(var_comb)>1) { # ESM
warm_wet_all_var_current =stack(var_imp_list_warm_wet)
plot(warm_wet_all_var_current)
# From this list, produce the final models by emsembling them again!
#var_imp_list_2050
  # Creating my own esenbling function
  obs <- speciesN$Occurrence
  pred <- extract(warm_wet_all_var_current,speciesN) # que extrai todos os valores de adequabilidade dos pontos de speciesN
  ## e enfim o ev, que mostra limiares (thresholds) de varias regras diferentes
  # Get evaluations by AUC from the mean models
  x=1:ncol(pred)
  ev_fun <- function(x) {
    ev <- evaluates(obs,pred[,x])
    ev1 = ev@statistics$AUC
    ev2 = median(ev@threshold_based$TSS)
    ev3 <- ev@threshold_based$threshold[2]
    ev_all = c(ev1,ev2,ev3)  
    return(ev_all)
  }
  ev_all= lapply(x,ev_fun)
  AUC_weights_all_var_wet = simplify2array(ev_all)[1,]
  TSS_weights_all_var_wet = simplify2array(ev_all)[2,]
  sensi_specificity_wet = simplify2array(ev_all)[3,]
  
# Calculating weighted mean from AUC
warm_wet_ensemble_all_var_current = weighted.mean(x=warm_wet_all_var_current, w=TSS_weights_all_var_wet)

} else { SDM
  warm_wet_ensemble_all_var_current = var_imp_list_warm_wet
}

plot(warm_wet_ensemble_all_var_current)
plot(population_pts, add=TRUE)

plot(warm_wet_ensemble_all_var_current,col = coul, main="Current",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(warm_wet_ensemble_all_var_current,
            paste0(path_modelling,"1output_rasters/",
                   "warm_wet_ensemble_all_var_current",
                   ".tif"),overwrite=TRUE)

## -----------------------------------------------------------------------------
model.stack.test=warm_wet_ensemble_all_var_current
sp1=population_pts

library(gdistance)
model_cost.ras  <- -1 * log(model.stack.test)    
model_trans.ras <- 1 / model_cost.ras  
trans     <- transition(model_trans.ras, transitionFunction=mean, directions = 8)
trans     <- geoCorrection(trans, type="c", multpl=F)
lin_dist.ras = accCost(trans, sp1)
thisLineagePoints.ras <- rasterize(sp1, model.stack.test, 1)
lin_dist.ras <- raster::distance(thisLineagePoints.ras) / 1000   
plot(lin_dist.ras)
lin_dist.ras[lin_dist.ras < 300 ] <- 300
plot(lin_dist.ras)
lin_weight.ras <- 1 / (lin_dist.ras) # Using inverse of distance after 200km
lin_weight.ras = lin_weight.ras/max(na.omit(lin_weight.ras[1:ncell(lin_weight.ras)]))
plot(lin_weight.ras)
model.stack <- lin_weight.ras * model.stack.test
plot(model.stack)
#plot(warm_wet_ensemble_all_var_current)
model.stack.test = model.stack
#
obs <- speciesN$Occurrence
pred <- extract(model.stack.test,speciesN) # created from pseudoabsence generator
e_sp1 <- evaluates(obs,pred)
tr_sp1 <- e_sp1@threshold_based$threshold[2]
tr_sp1 = min(na.omit(raster::extract(model.stack,classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==1,])))

# 0.3 is the mean value for threshold for both populations
#(thres_seasonal+thres_warm_wet)
# (0.1395671+0.4601434)/2

#tr_sp1 = median(sensi_specificity)
#binary_code_sp1 <- c(NA, NA, 0, -10, tr_sp1, 0,tr_sp1, 10, 1)
#matrix_bin_sp1 <- matrix(binary_code_sp1, ncol = 3, byrow = TRUE)
tr_sp1=0.3
pr.pa2_warm_wet <- raster(model.stack.test) # creating an empty raster
pr.pa2_warm_wet[] <- ifelse(model.stack.test[] >= tr_sp1,1,0) ## Esse aqui È seu model bin·rio agora


# Plotting
plot(pr.pa2_warm_wet, col=c("#f0f0f0","#969696"), axes=FALSE, box=FALSE, legend=FALSE, main="Current distribution")
plot(SAshp2,add=TRUE)
plot(population_pts_warm_wet,pch=21,col="black",bg="white", add=TRUE,cex=1)
#plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==2,],pch = 21, cex=1.5, col="black", bg="#377eb8",lwd=1.5, add=TRUE)
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==1,], pch = 21, cex=1.5, col="black", bg="#e41a1c",lwd=1.5, add=TRUE)

writeRaster(pr.pa2_warm_wet,
            paste0(path_modelling,"1output_rasters/",
                   "pr.pa2_warm_wet",
                   ".tif"),overwrite=TRUE)

## -----------------------------------------------------------------------------
plot(pr.pa2_warm_wet, main="Warm_wet_stable")
plot(population_pts, add=TRUE)

## -----------------------------------------------------------------------------
# Getting the path to the variables
fitted_models_path <- list.files(path=paste0(path,"/sdm_models/"),pattern=population_name,full.names = T) # Here I had to manually change the population names before the name of the models. Check what's wrong later.


## ----------------------------------------------------------------------------- For future 2050 A1B warm wet
which_preds = climond_2050_A1B

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2050_A1B_w_t = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2050_A1B_w_t =stack(var_imp_list_2050_A1B_w_t) # rem A1b
# Calculating weighted mean from AUC
ensemble_all_var_2050_A1B_w_t = weighted.mean(x=all_var_2050_A1B_w_t, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2050_A1B_w_t = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2050_A1B_w_t)
plot(population_pts, add=TRUE)

# Testing the match > 50%
ensemble_all_var_2050_A2_w_t_test = ensemble_all_var_2050_A1B_w_t
ensemble_all_var_2050_A2_w_t_test[ensemble_all_var_2050_A2_w_t_test<0.3]=0
plot(ensemble_all_var_2050_A2_w_t_test)
plot(population_pts, add=TRUE)

plot(ensemble_all_var_2050_A1B_w_t,col = coul, main="2070 ssp370 (high CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2050_A1B_w_t,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2050_A1B_w_t",
                   ".tif"),overwrite=TRUE)

## ----------------------------------------------------------------------------- For future 2050 A2 warm wet
preds.climond_stck
plot(preds.climond_stck[["bio_12"]])
plot(climond_2050_A1B[["bio_12"]])
plot(climond_2100_A1B[["bio_12"]])

plot(climond_2050_A1B[["bio_4"]]-preds.climond_stck[["bio_4"]])


pre_data = raster::extract(preds.climond_stck[["bio_17"]],speciesN)
future_data = raster::extract(climond_2050_A1B[["bio_17"]],speciesN)
future_data_a2 = raster::extract(climond_2100_A1B[["bio_17"]],speciesN)

# Testing with undowloaded data
#bio_5_ukesm_2050 = raster("/Users/josue/Dropbox/4Environmental_layers/future_chelsa/CHELSA_bio5_2041-2070_ukesm1-0-ll_ssp370_V.2.1.tif")

#bio_5_ukesm_2050 = raster("/Users/josue/Dropbox/4Environmental_layers/future_chelsa/CHELSA_bio5_2041-2070_mpi-esm1-2-hr_ssp370_V.2.1.tif")
#future_data = raster::extract(bio_5_ukesm_2050,sp11)

boxplot(pre_data,future_data,future_data_a2)

which_preds = climond_2050_A2

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2050_A2_w_t = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2050_A2_w_t =stack(var_imp_list_2050_A2_w_t) # rem A1b
# Calculating weighted mean from AUC
ensemble_all_var_2050_A2_w_t = weighted.mean(x=all_var_2050_A2_w_t, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2050_A2_w_t = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2050_A2_w_t,col = coul, main="2070 ssp585 (extreme CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2050_A2_w_t,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2050_A2_w_t",
                   ".tif"),overwrite=TRUE)

# Testing the match > 50%
ensemble_all_var_2050_A2_w_t_test = ensemble_all_var_2050_A2_w_t
ensemble_all_var_2050_A2_w_t_test[ensemble_all_var_2050_A2_w_t_test<0.3]=0
plot(ensemble_all_var_2050_A2_w_t_test)
plot(population_pts, add=TRUE)

## ----------------------------------------------------------------------------- For future 2100 A1B warm wet 
which_preds = climond_2100_A1B

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2100_A1B_w_t = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2100_A1B_w_t =stack(var_imp_list_2100_A1B_w_t) # rem A1b
# Calculating weighted mean from AUC
ensemble_all_var_2100_A1B_w_t = weighted.mean(x=all_var_2100_A1B_w_t, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2100_A1B_w_t = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2100_A1B_w_t,col = coul, main="2100 ssp370 (high CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2100_A1B_w_t,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2100_A1B_w_t",
                   ".tif"),overwrite=TRUE)

# Testing the match > 50%
ensemble_all_var_2050_A2_w_t_test = ensemble_all_var_2100_A1B_w_t
ensemble_all_var_2050_A2_w_t_test[ensemble_all_var_2050_A2_w_t_test<0.3]=0
plot(ensemble_all_var_2050_A2_w_t_test)
plot(population_pts, add=TRUE)


## ----------------------------------------------------------------------------- For future 2100 A1B warm wet
which_preds = climond_2100_A2
if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2100_A2_w_t = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2100_A2_w_t =stack(var_imp_list_2100_A2_w_t) # rem A1b
# Calculating weighted mean from AUC
ensemble_all_var_2100_A2_w_t = weighted.mean(x=all_var_2100_A2_w_t, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2100_A2_w_t = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2100_A2_w_t,col = coul, main="2100 ssp585 (extreme CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2100_A2_w_t,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2100_A2_w_t",
                   ".tif"),overwrite=TRUE)

# Testing the match > 50%
ensemble_all_var_2050_A2_w_t_test = ensemble_all_var_2100_A2_w_t
ensemble_all_var_2050_A2_w_t_test[ensemble_all_var_2050_A2_w_t_test<0.1]=0
plot(ensemble_all_var_2050_A2_w_t_test)
plot(population_pts, add=TRUE)

## ----------------------------------------------------------------------------- For future 2040 A1B warm wet
# Getting the path to the variables
fitted_models_path <- list.files(path=paste0(path,"/sdm_models"),pattern=population_name,full.names = T) 

which_preds = climond_2040_A1B
if (length(var_comb)>1) { # If ensemble of small models
  var_imp_list_2040_A1B_w_t = pblapply(var_comb,pred_two_model_fun,cl=1)
  all_var_2040_A1B_w_t =stack(var_imp_list_2040_A1B_w_t) # rem A1b
  # Calculating weighted mean from AUC
  ensemble_all_var_2040_A1B_w_t = weighted.mean(x=all_var_2040_A1B_w_t, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2040_A1B_w_t = pred_two_model_fun(var_comb) # SDM
}
plot(ensemble_all_var_2040_A1B_w_t,col = coul, main="2040 ssp370 (high CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2040_A1B_w_t,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2040_A1B_w_t",
                   ".tif"),overwrite=TRUE)

## ----------------------------------------------------------------------------- For future 2040 A2 warm wet
which_preds = climond_2040_A2
if (length(var_comb)>1) { # If ensemble of small models
  var_imp_list_2040_A2_w_t = pblapply(var_comb,pred_two_model_fun,cl=1)
  all_var_2040_A2_w_t =stack(var_imp_list_2040_A2_w_t) # rem A1b
  # Calculating weighted mean from AUC
  ensemble_all_var_2040_A2_w_t = weighted.mean(x=all_var_2040_A2_w_t, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2040_A2_w_t = pred_two_model_fun(var_comb) # SDM
}
plot(ensemble_all_var_2040_A2_w_t,col = coul, main="2040 ssp585 (extreme CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2040_A2_w_t,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2040_A2_w_t",
                   ".tif"),overwrite=TRUE)

## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------
## ----------------------------------------------------------------------------- ALLLLLLLL POPPPPPPPPPPPPPPPssssssss
################################################### -----------------------------------------------------------------------------## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------
nrow(classified_ind)
head(classified_ind)
population_pts = classified_ind[,c("kind","long","lat","long3_j","lat3_j")]

# Check selection
unique(population_pts$kind)

nrow(population_pts)

# One point per grid ~ 30, then I added 10km for the thinning
# Thin data
spp_coords_thinned <- spThin::thin( loc.data = population_pts,
                                    lat.col = "lat3_j", long.col = "long3_j",
                                    spec.col = "kind",
                                    thin.par = 40, reps = 1,
                                    locs.thinned.list.return = TRUE,
                                    write.files = FALSE,
                                    write.log.file = FALSE)

spp_coords_thinned = as.data.frame(spp_coords_thinned)
head(spp_coords_thinned)

spp_coords_thinned$kind = 'All_pops'
spp_coords_thinned = spp_coords_thinned[,c(3,1,2)]
colnames(spp_coords_thinned) = c("kind", "long", "lat")
head(spp_coords_thinned)

# Changing to Behrman
coordinates(spp_coords_thinned) = spp_coords_thinned[,c("long","lat")]
proj4string(spp_coords_thinned) = oldproj
population_pts = spTransform(spp_coords_thinned,CRSobj = behrmannCRS)
population_pts@data$long = population_pts@coords[,1]
population_pts@data$lat = population_pts@coords[,2]

writeOGR()

writeOGR(population_pts,"/Users/josue/Downloads/" ,"population_pts", driver="ESRI Shapefile",overwrite = TRUE)

plot(population_pts, add=T)
# Clean the data to one record per grid cell
bck.na<-preds_modelling_mask[[1]]
bck.na[]<-NA
r<-rasterize(population_pts@data[,c('long', 'lat')],bck.na,fun='count')
occ_data_clean_ok.pa<-rasterToPoints(r,fun=function(x){x>0}, spatial=T) 
population_pts = as.data.frame(occ_data_clean_ok.pa@coords)
population_pts$kind = 'All_pops'
population_pts = population_pts[,c(3,1,2)]
colnames(population_pts) = c("kind", "long", "lat")
coordinates(population_pts) = population_pts[,c("long","lat")]
proj4string(population_pts) = behrmannCRS

population_pts@data <- population_pts@data[,"kind",drop=FALSE] ## Necessary to not consider latitude or longitude as predictors

## -----------------------------------------------------------------------------
unique(population_pts@data$kind)


## -----------------------------------------------------------------------------
unique(population_pts@data$kind)
#head(population_pts_all@data$kind)
population_pts@data$Occurrence = 1
#head(population_pts_all@data)
population_pts@data <- population_pts@data[,"Occurrence",drop=FALSE]
population_name= "All_pops"


## -----------------------------------------------------------------------------
# How many unique records intersect with the raster?
print("total number of records for the population")
nrow(population_pts)
n_of_locs= extract(preds_modelling_mask[[1]],population_pts,cellnumbers=TRUE)
n_of_locs=length(unique(n_of_locs[,1]))
print("number of unique points given cell resolution")
n_of_locs
print("proportion of points per variable for the small models (min 5, best>10)")
(n_of_locs/2)-1 


## -----------------------------------------------------------------------------
## Now the ensembling of small models###############
#length(var_comb)
# function for two by two models
###################################################
population_name= "All_pops"
var_comb
which_preds=preds_modelling_mask

## gerando pontos de pseudo-ausencias for evalutation
set.seed(999)
  psudo <- sampleRandom(preds_modelling_mask[[1]],nrow(population_pts),xy=TRUE,sp=TRUE) 
  psudo@data$Occurrence <- 0 # retemos apenas a coluna Occurrence
  psudo@data <- psudo@data[,"Occurrence",drop=FALSE]
  buffer_pts <- gBuffer(population_pts,width=33000) 
  proj4string(buffer_pts) <- proj4string(psudo) # igualando projecao de buffer_pts e de psudo
  pts_outside <- psudo[is.na(over(psudo,buffer_pts)),] 
  pts_outside@data$Occurrence <- 0 
  sp11=population_pts
  sp11@data$Occurrence <- 1 
  sp11@data <- sp11@data[,"Occurrence",drop=FALSE]
  #pts_outside@data <- pts_outside@data[,"Occurrence",drop=F]
  proj4string(sp11) <- proj4string(pts_outside) # igualando projecao de buffer_pts e de psudo
  speciesN <- bind(sp11,pts_outside) # reunindo presencas e pseudoause


## -----------------------------------------------------------------------------
#getwd()
dir.create(paste0(path,"/sdm_models/"))
dir.exists(paste0(path,"/sdm_models/"))
setwd(paste0(path,"/sdm_models/"))

# Important to not skip it
population_name= "All_pops"
paste0(path,"/sdm_models/",population_name,"_",paste0(unlist(var_comb[1]),collapse = '_'),"_model",".sdm") # test if the path is correct

## -----------------------------------------------------------------------------
# Modelling
# Modelling
if (length(var_comb)>1) { # ESM
  var_imp_list_all_pop= pblapply(var_comb,two_model_fun_final,cl=1)
} else { # SDM
  var_imp_list_all_pop= two_model_fun_final(var_comb)  
}

if (length(var_comb)>1) { # ESM
all_pop_all_var_current =stack(var_imp_list_all_pop)
plot(all_pop_all_var_current)
# From this list, produce the final models by emsembling them again!
#var_imp_list_2050
  # Creating my own esenbling function
  obs <- speciesN$Occurrence
  pred <- extract(all_pop_all_var_current,speciesN) # que extrai todos os valores de adequabilidade dos pontos de speciesN
  ## e enfim o ev, que mostra limiares (thresholds) de varias regras diferentes
  # Get evaluations by AUC from the mean models
  x=1:ncol(pred)
  ev_fun <- function(x) {
    ev <- evaluates(obs,pred[,x])
    ev1 = ev@statistics$AUC
    ev2 = median(ev@threshold_based$TSS)
    ev3 <- ev@threshold_based$threshold[2]
    ev_all = c(ev1,ev2,ev3)  
    return(ev_all)
  }
  ev_all= lapply(x,ev_fun)
  AUC_weights_all_var_all_pop = simplify2array(ev_all)[1,]
  TSS_weights_all_var_all_pop = simplify2array(ev_all)[2,]
  sensi_specificity_all_pop = simplify2array(ev_all)[3,]
  
# Calculating weighted mean from AUC
all_pop_ensemble_all_var_current = weighted.mean(x=all_pop_all_var_current, w=TSS_weights_all_var_all_pop)
} else {
  all_pop_ensemble_all_var_current = var_imp_list
}

options(repr.plot.width=7, repr.plot.height=7)
par(mar = c(1, 1, 1, 1))

plot(all_pop_ensemble_all_var_current,col = coul, main="Current",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(all_pop_ensemble_all_var_current,
            paste0(path_modelling,"1output_rasters/",
                   "all_pop_ensemble_all_var_current",
                   ".tif"),overwrite=TRUE)

# To check where the suitability values are smaller than SE+SP
median_SE_SP_all_pop = median(sensi_specificity_all_pop)
median_SE_SP_all_pop
all_pop_ensemble_all_var_current_2 = all_pop_ensemble_all_var_current
all_pop_ensemble_all_var_current_2[all_pop_ensemble_all_var_current_2<median_SE_SP_all_pop]=0
plot(all_pop_ensemble_all_var_current_2)  
plot(population_pts, add=TRUE)



## -----------------------------------------------------------------------------
model.stack.test=all_pop_ensemble_all_var_current
sp1=population_pts

library(gdistance)
model_cost.ras  <- -1 * log(model.stack.test)    
model_trans.ras <- 1 / model_cost.ras  
trans     <- transition(model_trans.ras, transitionFunction=mean, directions = 8)
trans     <- geoCorrection(trans, type="c", multpl=F)
lin_dist.ras = accCost(trans, sp1)
thisLineagePoints.ras <- rasterize(sp1, model.stack.test, 1)
lin_dist.ras <- raster::distance(thisLineagePoints.ras) / 1000   
plot(lin_dist.ras)
lin_dist.ras[lin_dist.ras < 500 ] <- 500
plot(lin_dist.ras)
lin_weight.ras <- 1 / (lin_dist.ras) # Using inverse of distance after 200km
lin_weight.ras = lin_weight.ras/max(na.omit(lin_weight.ras[1:ncell(lin_weight.ras)]))
plot(lin_weight.ras)
model.stack <- lin_weight.ras * model.stack.test
plot(model.stack)
model.stack.test = model.stack
#
obs <- speciesN$Occurrence
pred <- extract(model.stack.test,speciesN) # created from pseudoabsence generator
e_sp1 <- evaluates(obs,pred)
tr_sp1 <- e_sp1@threshold_based$threshold[2]
#tr_sp1 <- 0.5 # mannually changing values to include all points
#binary_code_sp1 <- c(NA, NA, 0, -10, tr_sp1, 0,tr_sp1, 10, 1)
#matrix_bin_sp1 <- matrix(binary_code_sp1, ncol = 3, byrow = TRUE)
min(raster::extract(all_pop_ensemble_all_var_current,classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==1,]))

tr_sp1=0.3
pr.pa2_all <- raster(model.stack.test) # creating an empty raster
pr.pa2_all[] <- ifelse(model.stack.test[] >= tr_sp1,1,0) ## Esse aqui È seu model bin·rio agora


plot(pr.pa2_all, col=c("#f0f0f0","#969696"), axes=FALSE, box=FALSE, legend=FALSE, main="Current distribution")
plot(SAshp2,add=TRUE)
plot(population_pts_dry_ses,pch=21,col="black",bg="white", add=TRUE,cex=1)
plot(population_pts_warm_wet,pch=21,col="black",bg="white", add=TRUE,cex=1)
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==2,],pch = 21, cex=1.5, col="black", bg="#377eb8",lwd=1.5, add=TRUE)
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==1,], pch = 21, cex=1.5, col="black", bg="#e41a1c",lwd=1.5, add=TRUE)

## pr.pa2 esse aqui È o modelo binário ponderado pela distancia e pelo seu threshold

writeRaster(pr.pa2_all,
            paste0(path_modelling,"1output_rasters/",
                   "pr.pa2_all",
                   ".tif"),overwrite=TRUE)

## -----------------------------------------------------------------------------
plot(pr.pa2_all, main="All")
plot(population_pts, add=TRUE)
e_sp1@threshold_based


## -----------------------------------------------------------------------------
median(sensi_specificity)

## -----------------------------------------------------------------------------
fitted_models_path <- list.files(path=paste0(path,"/sdm_models"),pattern=population_name,full.names = T)
fitted_models_path # manually add the population name in front of the layers in the directory (SDM not writing the name correctly)
population_name

## ----------------------------------------------------------------------------- For future 2050 A1B All
which_preds = climond_2050_A1B

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2050_A1B_All = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2050_A1B_All =stack(var_imp_list_2050_A1B_All) # rem A1b
# Calculating weighted mean from AUC
ensemble_all_var_2050_A1B_All = weighted.mean(x=all_var_2050_A1B_All, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2050_A1B_All = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2050_A1B_All,col = coul, main="2070 ssp370 (high CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2050_A1B_All,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2050_A1B_All",
                   ".tif"),overwrite=TRUE)

# Testing the match > 50%
ensemble_all_var_2050_A2_All_test = ensemble_all_var_2050_A1B_All
ensemble_all_var_2050_A2_All_test[ensemble_all_var_2050_A2_All_test<0.5]=0
plot(ensemble_all_var_2050_A2_All_test)
plot(population_pts, add=TRUE)

## ----------------------------------------------------------------------------- For future 2050 A2 All
which_preds = climond_2050_A2

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2050_A2_All = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2050_A2_All =stack(var_imp_list_2050_A2_All) # rem A1b
# Calculating weighted mean from AUC
ensemble_all_var_2050_A2_All = weighted.mean(x=all_var_2050_A2_All, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2050_A2_All = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2050_A2_All,col = coul, main="2070 ssp585 (extreme CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2050_A2_All,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2050_A2_All",
                   ".tif"),overwrite=TRUE)

# Testing the match > 50%
ensemble_all_var_2050_A2_All_test = ensemble_all_var_2050_A2_All
ensemble_all_var_2050_A2_All_test[ensemble_all_var_2050_A2_All_test<0.5]=0
plot(ensemble_all_var_2050_A2_All_test)
plot(population_pts, add=TRUE)

## ----------------------------------------------------------------------------- For future 2100 A1B All 
which_preds = climond_2100_A1B

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2100_A1B_All = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2100_A1B_All =stack(var_imp_list_2100_A1B_All) # rem A1b
# Calculating weighted mean from AUC
ensemble_all_var_2100_A1B_All = weighted.mean(x=all_var_2100_A1B_All, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2100_A1B_All = pred_two_model_fun(var_comb) # SDM
}
plot(ensemble_all_var_2100_A1B_All,col = coul, main="2100 ssp370 (high CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2100_A1B_All,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2100_A1B_All",
                   ".tif"),overwrite=TRUE)

# Testing the match > 50%
ensemble_all_var_2050_A2_All_test = ensemble_all_var_2100_A1B_All
ensemble_all_var_2050_A2_All_test[ensemble_all_var_2050_A2_All_test<0.5]=0
plot(ensemble_all_var_2050_A2_All_test)
plot(population_pts, add=TRUE)

## ----------------------------------------------------------------------------- For future 2100 A12 All
which_preds = climond_2100_A2

if (length(var_comb)>1) { # If ensemble of small models
var_imp_list_2100_A2_All = pblapply(var_comb,pred_two_model_fun,cl=1)
all_var_2100_A2_All =stack(var_imp_list_2100_A2_All) # rem A1b
# Calculating weighted mean from AUC
ensemble_all_var_2100_A2_All = weighted.mean(x=all_var_2100_A2_All, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2100_A2_All = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2100_A2_All,col = coul, main="2100 ssp585 (extreme CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2100_A2_All,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2100_A2_All",
                   ".tif"),overwrite=TRUE)

# Testing the match > 50%
ensemble_all_var_2050_A2_All_test = ensemble_all_var_2100_A2_All
ensemble_all_var_2050_A2_All_test[ensemble_all_var_2050_A2_All_test<0.5]=0
plot(ensemble_all_var_2050_A2_All_test)
plot(population_pts, add=TRUE)

## ----------------------------------------------------------------------------- For future 2040 AB1 All
which_preds = climond_2040_A1B

if (length(var_comb)>1) { # If ensemble of small models
  var_imp_list_2040_A1B_All = pblapply(var_comb,pred_two_model_fun,cl=1)
  all_var_2040_A1B_All =stack(var_imp_list_2040_A1B_All) # rem A1b
  # Calculating weighted mean from AUC
  ensemble_all_var_2040_A1B_All = weighted.mean(x=all_var_2040_A1B_All, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2040_A1B_All = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2040_A1B_All,col = coul, main="2040 ssp370 (high CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)

writeRaster(ensemble_all_var_2040_A1B_All,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2040_A1B_All",
                   ".tif"),overwrite=TRUE)

## ----------------------------------------------------------------------------- For future 2040 AB2 All
which_preds = climond_2040_A2

if (length(var_comb)>1) { # If ensemble of small models
  var_imp_list_2040_A2_All = pblapply(var_comb,pred_two_model_fun,cl=1)
  all_var_2040_A2_All =stack(var_imp_list_2040_A2_All) # rem A1b
  # Calculating weighted mean from AUC
  ensemble_all_var_2040_A2_All = weighted.mean(x=all_var_2040_A2_All, w=TSS_weights_all_var)
} else {
  ensemble_all_var_2040_A2_All = pred_two_model_fun(var_comb) # SDM
}

plot(ensemble_all_var_2040_A2_All,col = coul, main="2040 ssp585 (extreme CO2)",axes=FALSE, box=FALSE)
#plot(population_pts,pch=1,col="black",border="gray", add=TRUE,cex=.6)
plot(SAshp2,add=TRUE)
#ssp370 (high CO2)
#ssp585 (extreme CO2)
ensemble_all_var_2040_A2_All = ensemble_all_var_2040_A2
writeRaster(ensemble_all_var_2040_A2_All,
            paste0(path_modelling,"1output_rasters/",
                   "ensemble_all_var_2040_A2_All",
                   ".tif"),overwrite=TRUE)

## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
# Now, preparing layers and functions for MgClim
# Finally!!!!!!!!!
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## ----------------------------------------------------------------------------- Function to prepare data for dispersal analyses ## -----------------------------------------------------------------------------
# Recover rasters from SDMs above:
# Reload the rasters from their saved paths
pr.pa2 = raster(paste0(path_modelling, "1output_rasters/pr.pa2.tif"))
ensemble_all_var_2040_A1B = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2040_A1B.tif"))
ensemble_all_var_2050_A1B = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2050_A1B.tif"))
ensemble_all_var_2100_A1B = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2100_A1B.tif"))
ensemble_all_var_2040_A2 = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2040_A2.tif"))
ensemble_all_var_2050_A2 = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2050_A2.tif"))
ensemble_all_var_2100_A2 = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2100_A2.tif"))
pr.pa2_warm_wet = raster(paste0(path_modelling, "1output_rasters/pr.pa2_warm_wet.tif"))
ensemble_all_var_2040_A1B_w_t = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2040_A1B_w_t.tif"))
ensemble_all_var_2050_A1B_w_t = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2040_A1B_w_t.tif")) 

# The one above is wrong in the year, review later

ensemble_all_var_2100_A1B_w_t = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2100_A1B_w_t.tif"))
ensemble_all_var_2040_A2_w_t = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2040_A2_w_t.tif"))
ensemble_all_var_2050_A2_w_t = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2050_A2_w_t.tif"))
ensemble_all_var_2100_A2_w_t = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2100_A2_w_t.tif"))
# The rasters with "_All" in their names are not directly found in the provided script. Assuming they follow the same pattern:
ensemble_all_var_2040_A1B_All = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2040_A1B_All.tif"))
ensemble_all_var_2050_A1B_All = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2050_A1B_All.tif"))
ensemble_all_var_2100_A1B_All = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2100_A1B_All.tif"))
ensemble_all_var_2040_A2_All = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2040_A2_All.tif"))
ensemble_all_var_2050_A2_All = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2050_A2_All.tif"))
ensemble_all_var_2100_A2_All = raster(paste0(path_modelling, "1output_rasters/ensemble_all_var_2100_A2_All.tif"))
pr.pa2_all = raster(paste0(path_modelling, "1output_rasters/pr.pa2_all.tif"))



prep_disp_an = function(range_current,
                        model_future_1,model_future_2,
                        model_future_3,
                        stack_barriers,
                        ras_mask,
                        thres_pres){
  # First disaggragate
  presence_forested = disaggregate(range_current,17) %>%
    # Cropping presence layer by forest layer
    raster::crop(.,preds_modelling_mask_d) %>%
    # Resample by model raster
    raster::resample(., preds_modelling_mask_d)
  
  # Create presence forested
  presence_forested = raster::mask(presence_forested,stack_barriers[[1]],maskvalue=1)
  plot(presence_forested)
  presence_forested[is.na(presence_forested)] <- 0
  presence_forested = raster::mask(presence_forested,preds_modelling_mask_d)
  # After resampling, some empty cells received a value
  presence_forested[presence_forested>0] = 1
  
  plot(presence_forested)
  plot(population_pts,add=TRUE, pch=1)
  
  print("preparing rasters for future")
  ensemble_all_var_future = stack(model_future_1,model_future_2,model_future_3)
  ensemble_all_var_future = ensemble_all_var_future*1000
  ensemble_all_var_future[ensemble_all_var_future<thres_pres] = NA
  ensemble_all_var_future=disaggregate(ensemble_all_var_future,17) %>%
    raster::crop(., preds_modelling_mask_d) %>%
    raster::resample(., preds_modelling_mask_d)
  names(ensemble_all_var_future) = c("y2040","y2070","y2100")
  return(list(presence_forested,ensemble_all_var_future))
}
## ----------------------------------------------------------------------------- 

## ----------------------------------------------------------------------------- Function to run sequential dispersal analyses ## -----------------------------------------------------------------------------
dispersal_funtion = function(presence_forested,ensemble_all_var_future,land_cover, deforestation, test_run, path_to_migClim, simulName,rcThreshold,dispSteps, extent_test_mg){
  
  setwd(path_to_migClim)
  
  print(paste("begin of DispAn" ,date()))
  print("Preparing dataframe step 1")
  
  print("Subseting land cover")
  if (deforestation == "strong") {
    land_cover = raster::subset(land_cover,c("Current","A2_2050","A2_2100"))
  } else {
    land_cover = raster::subset(land_cover,c("Current","A1B_2050","A1B_2100"))
  }
  print("land_cover ok")
  
  print("Stacking MigClim_var")
  MigClim_var = raster::stack(presence_forested, ensemble_all_var_future[[1]] ,land_cover[[1]])
  #plot(MigClim_var)
  print("Stacking ok")
  
  if (test_run==TRUE) {
    print("Testing with a smaller dataset")
    MigClim_var = raster::crop(MigClim_var,extent_test_mg)
    plot(MigClim_var)
    print("Cropping worked")
    #plot(MigClim_var,col=c("red","blue"))  
  } else {
    print("That's awsome, you're really trying")
  }
  
  mig_clim_df_praper = function(MigClim_var){
    names(MigClim_var)=c("InitialDist","hsmap1","Barrier")
    #MigClim_var
    coord_final_rasters = coordinates(MigClim_var[[1]])
    #head(coord_final_rasters)
    MigClim_df = as.data.frame(raster::extract(MigClim_var,coord_final_rasters))
    MigClim_df$X = coord_final_rasters[,1]
    MigClim_df$Y = coord_final_rasters[,2]
    MigClim_df = na.exclude(MigClim_df)
    dim(MigClim_df)
    MigClim_df = MigClim_df[,c("X","Y","InitialDist","hsmap1","Barrier")] # It has to be in this order :/
    MigClim_df$hsmap1 = round(MigClim_df$hsmap1,0) # has to be integer
    head(MigClim_df)
    return(MigClim_df)
  }
  
  print("preparing MigClim step1 dataframe")
  MigClim_df = mig_clim_df_praper(MigClim_var)
  print("MigClim DF worked")
  
  print("Almost running")
  N <- MigClim.migrate(iniDist=MigClim_df[,1:3],
                       hsMap=MigClim_df[,4], 
                       rcThreshold=rcThreshold,
                       envChgSteps=1, 
                       dispSteps=30, 
                       dispKernel=c(1.0,0.4,0.16,0.06,0.03),
                       barrier=MigClim_df[,5],
                       barrierType="strong",
                       iniMatAge=1, 
                       propaguleProd=c(0.01,0.08,0.5,0.92),
                       lddFreq=0, #lddMinDist=20, lddMaxDist=40,
                       simulName= paste0(simulName,"_2040"), 
                       replicateNb=1,
                       overWrite=TRUE, 
                       testMode=FALSE, 
                       fullOutput=FALSE, keepTempFiles=FALSE)
  
  print("plotting 2040")
  MigClim.plot(asciiFile=paste0(path_to_migClim,
                                paste0(simulName,"_2040"),"/",
                                paste0(simulName,"_2040"),"_raster.asc"), outDir="", fileFormat="jpeg", fullOutput=FALSE)
  
  print("reading in presence in 2040")
  dist_2050 = raster(paste0(path_to_migClim,paste0(simulName,"_2040"),"/",paste0(simulName,"_2040"),"_raster.asc"))
  
  print("Adjust it to be only a presence absence map")
  presence_maker = function(dist_2050){
    dist_2050_d = dist_2050
    dist_2050_d[dist_2050_d<1] = NA
    dist_2050_d[dist_2050_d>=30000] = NA
    dist_2050_d[dist_2050_d>0] = 1
    proj4string(dist_2050_d) = behrmannCRS
    dist_2050_d = raster::resample(dist_2050_d,preds_modelling_mask_d[[1]])
    dist_2050_d[dist_2050_d>0] = 1
    dist_2050_d[is.na(dist_2050_d)] <- 0
    dist_2050_d = raster::mask(dist_2050_d,preds_modelling_mask_d[[1]])
    return(dist_2050_d)
  }
  
  print("preparing presence in 2040")
  dist_2050_d = presence_maker(raster(paste0(path_to_migClim,paste0(simulName,"_2040"),"/",paste0(simulName,"_2040"),"_raster.asc")))
  #plot(dist_2050_d)
  
  print("Taking a look at the 2040 outputs to see if they make sense")
  plot(stack(crop(MigClim_var,extent_test_mg),
             crop(dist_2050_d,extent_test_mg)))
  
  print("Preparing stack and step2 dataframe")
  MigClim_var = raster::stack(dist_2050_d,
                              ensemble_all_var_future[[2]],
                              land_cover[[2]])
  
  if (test_run==TRUE) {
    MigClim_var = raster::crop(MigClim_var,extent_test_mg)
    plot(MigClim_var)
    #plot(MigClim_var,col=c("red","blue"))  
  } else {
    print("That's awsome, you're going to step 2")
  }
  
  print("Preparing DF to step 2")
  MigClim_df = mig_clim_df_praper(MigClim_var)
  
  print("running MigClim 2070")
  N <- MigClim.migrate(iniDist=MigClim_df[,1:3],
                       hsMap=MigClim_df[,4], 
                       rcThreshold=rcThreshold,
                       envChgSteps=1, 
                       dispSteps=30, 
                       dispKernel=c(1.0,0.4,0.16,0.06,0.03),
                       barrier=MigClim_df[,5],
                       barrierType="strong",
                       iniMatAge=1, 
                       propaguleProd=c(0.01,0.08,0.5,0.92),
                       lddFreq=0, #lddMinDist=20, lddMaxDist=40,
                       simulName= paste0(simulName,"_2070"), 
                       replicateNb=1,
                       overWrite=TRUE, 
                       testMode=FALSE, 
                       fullOutput=FALSE, keepTempFiles=FALSE)
  
  print("preparing folder plot 2070")
  MigClim.plot(asciiFile=paste0(path_to_migClim,
                                paste0(simulName,"_2070"),"/",
                                paste0(simulName,"_2070"),"_raster.asc"),
               outDir="", fileFormat="jpeg", fullOutput=FALSE)
  
  print("reading in presence in 2070")
  dist_2050 = raster(paste0(path_to_migClim,paste0(simulName,"_2070"),"/",paste0(simulName,"_2070"),"_raster.asc"))
  
  # Adjust it to be only a presence absence map
  print("preparing presence in 2070")
  
  dist_2050_d = presence_maker(raster(paste0(path_to_migClim,paste0(simulName,"_2070"),"/",paste0(simulName,"_2070"),"_raster.asc")))
  
  print("Taking a look at 2070 outputs to see if they make sense")
  plot(stack(crop(MigClim_var,extent_test_mg),
             crop(dist_2050_d,extent_test_mg)))
  
  print("Preparing stack and step3 dataframe")
  MigClim_var = raster::stack(dist_2050_d,
                              ensemble_all_var_future[[3]],
                              land_cover[[3]])
  
  if (test_run==TRUE) {
    MigClim_var = raster::crop(MigClim_var,extent_test_mg)
    plot(MigClim_var)
    #plot(MigClim_var,col=c("red","blue"))  
  } else {
    print("That's awsome, you're going to step 3")
  }
  
  print("preparing DF for MigClim 2100")
  MigClim_df = mig_clim_df_praper(MigClim_var)
  
  print("running MigClim 2100")
  N <- MigClim.migrate(iniDist=MigClim_df[,1:3],
                       hsMap=MigClim_df[,4], 
                       rcThreshold=rcThreshold,
                       envChgSteps=1, 
                       dispSteps=30, 
                       dispKernel=c(1.0,0.4,0.16,0.06,0.03),
                       barrier=MigClim_df[,5],
                       barrierType="strong",
                       iniMatAge=1, 
                       propaguleProd=c(0.01,0.08,0.5,0.92),
                       lddFreq=0, #lddMinDist=20, lddMaxDist=40,
                       simulName= paste0(simulName,"_2100"), 
                       replicateNb=1,
                       overWrite=TRUE, 
                       testMode=FALSE, 
                       fullOutput=FALSE, keepTempFiles=FALSE)
  
  print("plotting in local folder")
  MigClim.plot(asciiFile=paste0(path_to_migClim,
                                paste0(simulName,"_2100"),"/",
                                paste0(simulName,"_2100"),"_raster.asc"),
               outDir="", fileFormat="jpeg", fullOutput=FALSE)
  ## -----------------------------------------------------------------------------
  print("reading in presence in 2100")
  dist_2050 = raster(paste0(path_to_migClim,paste0(simulName,"_2100"),"/",paste0(simulName,"_2100"),"_raster.asc"))
  ## -----------------------------------------------------------------------------
  # Adjust it to be only a presence absence map
  print("preparing presence in 2100")
  
  dist_2050_d = presence_maker(raster(paste0(path_to_migClim,paste0(simulName,"_2100"),"/",paste0(simulName,"_2100"),"_raster.asc")))
  # Taking a look at the outputs to see if they make sense
  plot(stack(crop(MigClim_var,extent_test_mg),
             crop(dist_2050_d,extent_test_mg)))
  return(list(dist_2050,dist_2050_d))
  
}

# A function only to get 2040 or 2070 results

get_disp_fun_results = function(path_to_migClim,simulName,year){
  get_disp_fun = function(year){
# Define presence maker function here as well  
  print("Adjust it to be only a presence absence map")
  presence_maker = function(dist_2050){
    dist_2050_d = dist_2050
    dist_2050_d[dist_2050_d<1] = NA
    dist_2050_d[dist_2050_d>=30000] = NA
    dist_2050_d[dist_2050_d>0] = 1
    proj4string(dist_2050_d) = behrmannCRS
    dist_2050_d = raster::resample(dist_2050_d,preds_modelling_mask_d[[1]])
    dist_2050_d[dist_2050_d>0] = 1
    dist_2050_d[is.na(dist_2050_d)] <- 0
    dist_2050_d = raster::mask(dist_2050_d,preds_modelling_mask_d[[1]])
    return(dist_2050_d)
  }
  
  print(paste("reading in presence", gsub("_","",year)))
  
  dist_2050 = raster(paste0(path_to_migClim,paste0(simulName,year),"/",paste0(simulName,year),"_raster.asc"))
  
  # Adjust it to be only a presence absence map
  print(paste("preparing presence in", gsub("_","",year)))
  
  dist_2050_d = presence_maker(raster(paste0(path_to_migClim,paste0(simulName,year),"/",paste0(simulName,year),"_raster.asc")))
  }
  dist_2050_e = lapply(year,get_disp_fun)
return(dist_2050_e)
}

## ----------------------------------------------------------------------------- 
# For testing and displaying figures in a subset of data

extent_test_mg = c(-5256894,
                   -4663615,
                   -535403.8,
                   587589)

#plot(ensemble_all_var_future[[3]])
#extent_test_mg = drawExtent()

## ----------------------------------------------------------------------------- ## ----------------------------------------------------------------------------- 
# MigClim for Dry_ses_AB1
## ----------------------------------------------------------------------------- ## ----------------------------------------------------------------------------- ## -----------------------------------------------------------------------------

# Preparing layers for Dry_ses_AB1
var_Dry_ses_AB1= prep_disp_an(
  range_current=pr.pa2,
  model_future_1=ensemble_all_var_2040_A1B,
  model_future_2=ensemble_all_var_2050_A1B,
  model_future_3=ensemble_all_var_2100_A1B,
  stack_barriers=rasters_future_land_cover_agr_proj,
  ras_mask=preds_modelling_mask_d,
  thres_pres=300)

# Running for Dry_ses AB1
Dry_ses_AB1 = dispersal_funtion(
  presence_forested=var_Dry_ses_AB1[[1]],
  ensemble_all_var_future=var_Dry_ses_AB1[[2]],
  land_cover= rasters_future_land_cover_agr_proj,
  deforestation = "weak",
  test_run=FALSE,
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="Dry_ses_ssp370",
  rcThreshold=0,
  dispSteps=30,
  extent_test_mg=extent_test_mg);date()

Dry_ses_AB1[[1]] # contains several values
Dry_ses_AB1[[2]] # processed presence file
plot(Dry_ses_AB1[[2]], main="Small dataset test")
#raster::zoom(Dry_ses_AB1[[2]],extent_test_mg,new=FALSE)



# Running for Dry_ses AB1 - No threshold
#Dry_ses_AB1_0thres = dispersal_funtion(
#  presence_forested=var_Dry_ses_AB1[[1]],
#  ensemble_all_var_future=var_Dry_ses_AB1[[2]],
#  land_cover= rasters_future_land_cover_agr_proj,
#  deforestation = "weak",
#  test_run=FALSE,
#  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
#  simulName="Dry_ses_ssp370_0thres",
#  rcThreshold=0,
#  dispSteps=30,
#  extent_test_mg=extent_test_mg);date()

#Dry_ses_AB1_0thres[[1]] # contains several values
#Dry_ses_AB1_0thres[[2]] # processed presence file
#plot(Dry_ses_AB1_0thres[[2]], main="Small dataset test")
#raster::zoom(Dry_ses_AB1[[2]],extent_test_mg,new=FALSE)

# MigClim for A2 dry-seas 
# Preparing layers for Dry_ses_A2
var_Dry_ses_A2=prep_disp_an(
  range_current=pr.pa2,
  model_future_1=ensemble_all_var_2040_A2,
  model_future_2=ensemble_all_var_2050_A2,
  model_future_3=ensemble_all_var_2100_A2,
  stack_barriers=rasters_future_land_cover_agr_proj,               ras_mask=preds_modelling_mask_d,
  thres_pres=300)


# Running for Dry_ses_A2
Dry_ses_A2 = dispersal_funtion(
  presence_forested=var_Dry_ses_A2[[1]],
  ensemble_all_var_future=var_Dry_ses_A2[[2]],
  land_cover= rasters_future_land_cover_agr_proj,
  deforestation = "strong",
  test_run=FALSE,
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="Dry_ses_ssp585",
  rcThreshold=0,
  dispSteps=30,
  extent_test_mg=extent_test_mg);date()

Dry_ses_A2[[1]] # contains several values
Dry_ses_A2[[2]] # processed presence file
plot(Dry_ses_A2[[2]], main="Small dataset test Dry_ses_A2")

## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------
## ----------------------------------------------------------------------------- MigClim Warm Wet
## ----------------------------------------------------------------------------- ## -----------------------------------------------------------------------------

# MigClim for AB1 WarmWet 
# Preparing layers for WarmWet_AB1
var_warm_wet_AB1=prep_disp_an(
  range_current=pr.pa2_warm_wet,
  model_future_1=ensemble_all_var_2040_A1B_w_t,
  model_future_2=ensemble_all_var_2050_A1B_w_t,
  model_future_3=ensemble_all_var_2100_A1B_w_t,
  stack_barriers=rasters_future_land_cover_agr_proj,               ras_mask=preds_modelling_mask_d,
  thres_pres=300)


# Dispersal analyses AB1 WarmWet 
Warm_wet_A1B = dispersal_funtion(
  presence_forested=var_warm_wet_AB1[[1]],
  ensemble_all_var_future=var_warm_wet_AB1[[2]],
  land_cover= rasters_future_land_cover_agr_proj,
  deforestation = "weak",
  test_run=FALSE,
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="WarmWet_ssp370", # CHANGE THIS!
  rcThreshold=0,
  dispSteps=30,
  extent_test_mg=extent_test_mg);date()

Warm_wet_A1B[[1]] # contains several values
Warm_wet_A1B[[2]] # processed presence file
plot(Warm_wet_A1B[[2]], main="Small dataset test Warm_wet_A1B")


# MigClim for A2 WarmWet 
# Preparing layers for WarmWet_A2
var_warm_wet_A2=prep_disp_an(
  range_current=pr.pa2_warm_wet,
  model_future_1=ensemble_all_var_2040_A2_w_t,
  model_future_2=ensemble_all_var_2050_A2_w_t,
  model_future_3=ensemble_all_var_2100_A2_w_t,
  stack_barriers=rasters_future_land_cover_agr_proj,               ras_mask=preds_modelling_mask_d,
  thres_pres=300)

Warm_wet_A2 = dispersal_funtion( # CHANGE THIS 
  presence_forested=var_warm_wet_A2[[1]],
  ensemble_all_var_future=var_warm_wet_A2[[2]],
  land_cover= rasters_future_land_cover_agr_proj,
  deforestation = "strong",
  test_run=TRUE,
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="WarmWet_ssp585", # CHANGE THIS!
  rcThreshold=0,
  dispSteps=30,
  extent_test_mg=extent_test_mg);date()

Warm_wet_A2[[1]] # contains several values
Warm_wet_A2[[2]] # processed presence file
plot(Dry_ses_A2[[2]], main="Small dataset test Dry_ses_A2")

## -----------------------------------------------------------------------------## ----------------------------------------------------------------------------- 
## -----------------------------------------------------------------------------## ----------------------------------------------------------------------------- 
## Test to run SSP585 simulations in 2100 under weaker deforestation
## -----------------------------------------------------------------------------## ----------------------------------------------------------------------------- 
## -----------------------------------------------------------------------------## ----------------------------------------------------------------------------- 

Warm_wet_A2_weak = dispersal_funtion( # CHANGE THIS 
  presence_forested=var_warm_wet_A2[[1]],
  ensemble_all_var_future=var_warm_wet_A2[[2]],
  land_cover= rasters_future_land_cover_agr_proj,
  deforestation = "weak",
  test_run=FALSE,
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="WarmWet_ssp585_weak", # CHANGE THIS!
  rcThreshold=0,
  dispSteps=30,
  extent_test_mg=extent_test_mg);date()

Warm_wet_A2_weak[[1]] # contains several values
Warm_wet_A2_weak[[2]] # processed presence file
plot(Warm_wet_A2_weak[[2]], main="Test warm_wet_A2_ssp585")

# Running for Dry_ses_A2
Dry_ses_A2_weak = dispersal_funtion(
  presence_forested=var_Dry_ses_A2[[1]],
  ensemble_all_var_future=var_Dry_ses_A2[[2]],
  land_cover= rasters_future_land_cover_agr_proj,
  deforestation = "weak",
  test_run=FALSE,
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="Dry_ses_ssp585_weak",
  rcThreshold=0,
  dispSteps=30,
  extent_test_mg=extent_test_mg);date()

Dry_ses_A2_weak[[1]] # contains several values
Dry_ses_A2_weak[[2]] # processed presence file
plot(Dry_ses_A2_weak[[2]], main="Small dataset test Dry_ses_A2")


## -----------------------------------------------------------------------------## ----------------------------------------------------------------------------- 
## -----------------------------------------------------------------------------## ----------------------------------------------------------------------------- 

## -----------------------------------------------------------------------------## ----------------------------------------------------------------------------- Dispersal analyses for All_pops
## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------## -----------------------------------------------------------------------------

# MigClim for All_pops_AB1
# Preparing layers for All_pops_AB1
var_All_pops_AB1=prep_disp_an(
  range_current=pr.pa2_all,
  model_future_1=ensemble_all_var_2040_A1B_All, ## Missing
  model_future_2=ensemble_all_var_2050_A1B_All,
  model_future_3=ensemble_all_var_2100_A1B_All,
  stack_barriers=rasters_future_land_cover_agr_proj,               ras_mask=preds_modelling_mask_d,
  thres_pres=300)


All_pops_AB1 = dispersal_funtion( # CHANGE THIS 
  presence_forested=var_All_pops_AB1[[1]],
  ensemble_all_var_future=var_All_pops_AB1[[2]],
  land_cover= rasters_future_land_cover_agr_proj,
  deforestation = "weak",
  test_run=FALSE,
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="All_pops_ssp370", # CHANGE THIS!
  rcThreshold=0,
  dispSteps=30,
  extent_test_mg=extent_test_mg);date()

All_pops_AB1[[1]] # contains several values
All_pops_AB1[[2]] # processed presence file
plot(All_pops_AB1[[2]], main="Small dataset test All_pops_AB1")

# MigClim for All_pops_A2
# Preparing layers for All_pops_A2
var_All_pops_A2=prep_disp_an(
  range_current=pr.pa2_all,
  model_future_1=ensemble_all_var_2040_A2_All, ## Missing
  model_future_2=ensemble_all_var_2050_A2_All,
  model_future_3=ensemble_all_var_2100_A2_All,
  stack_barriers=rasters_future_land_cover_agr_proj,               ras_mask=preds_modelling_mask_d,
  thres_pres=300)

All_pops_A2 = dispersal_funtion( # CHANGE THIS 
  presence_forested=var_All_pops_A2[[1]],
  ensemble_all_var_future=var_All_pops_A2[[2]],
  land_cover= rasters_future_land_cover_agr_proj,
  deforestation = "strong",
  test_run=FALSE,
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="All_pops_ssp585", # CHANGE THIS!
  rcThreshold=0,
  dispSteps=30,
  extent_test_mg=extent_test_mg);date()

All_pops_A2[[1]] # contains several values
All_pops_A2[[2]] # processed presence file
plot(All_pops_A2[[2]], main="Small dataset test All_pops_A2")


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
# Final maps of pop expansion retraction and evo rescue
## --------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
# Raster preparations for plotting range expantion/retract

exp_retract_plot = function(Mig_clim_future,Present_dist,pop_num,title_plot,deforestation){

Mig_clim_future[Mig_clim_future<1]=NA
Present_dist[Present_dist<1]=NA

if (deforestation == "strong_in_2070") {
  Mig_clim_future = raster::mask(Mig_clim_future,rasters_future_land_cover_agr_proj[[5]],maskvalue=1)
}

range_expansion_Dry_ses = raster::mask(Mig_clim_future,Present_dist,inverse=TRUE)

# Plotting expansion / retraction maps
options(repr.plot.width=7, repr.plot.height=7)
par(mar = c(1, 1, 1, 1))
## 

range_expansion_Dry_ses = raster::mask(Mig_clim_future,Present_dist,inverse=TRUE)

plot(bio_SPDF3,col=NA,border=NA,main=title_plot,cex.main=1)#, adj=0)
plot(Present_dist,col=c("white","#b2182b"),legend=FALSE, add=TRUE)
plot(Mig_clim_future,col=c(NA,"#4575b4"),legend=FALSE,add=TRUE)
plot(range_expansion_Dry_ses,col=c(NA,"#fee090"),legend=FALSE, add=TRUE)
plot(bio_SPDF3,col=NA,border="#969696", add=TRUE,lwd=1)
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans %in% pop_num,],pch = 21, cex=1.5, col="black", bg="#35978f",lwd=1.5, add=TRUE)
#plot(bio_SPDF3,col=NA,border="#969696", add=TRUE,lwd=1)
}

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

# Getting rasters for Dry seasonal 2040 & 2070
Dry_ses_AB1_2040_70= get_disp_fun_results(
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="Dry_ses_ssp370",
  year = list("_2040","_2070","_2100"))

plot(Dry_ses_AB1_2040_70[[3]])
Dry_ses_AB1 = stack(Dry_ses_AB1_2040_70[[3]],Dry_ses_AB1_2040_70[[3]])


Dry_ses_A2_2040_70= get_disp_fun_results(
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="Dry_ses_ssp370",
  year = list("_2040","_2070","_2100"))

Dry_ses_A2 = stack(Dry_ses_A2_2040_70[[3]],Dry_ses_A2_2040_70[[3]])

# Getting rasters for WarmWet 2040 & 2070
Warm_Wet_AB1_2040_70= get_disp_fun_results(
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="WarmWet_ssp370",
  year = list("_2040","_2070","_2100"))

Warm_wet_AB1 = stack(Warm_Wet_AB1_2040_70[[3]],Warm_Wet_AB1_2040_70[[3]])

Warm_Wet_A2_2040_70= get_disp_fun_results(
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="WarmWet_ssp585",
  year = list("_2040","_2070","_2100"))

Warm_wet_A2 = stack(Warm_Wet_A2_2040_70[[3]],Warm_Wet_A2_2040_70[[3]])

# Getting rasters for All pops 2040 & 2070
All_pops_AB1_2040_70= get_disp_fun_results(
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="All_pops_ssp370",
  year = list("_2040","_2070","_2100"))

All_pops_AB1 = stack(All_pops_AB1_2040_70[[3]],All_pops_AB1_2040_70[[3]])

All_pops_A2_2040_70= get_disp_fun_results(
  path_to_migClim = paste0(path,"/layers_for_MigClim/"),
  simulName="All_pops_ssp585",
  year = list("_2040","_2070","_2100"))

All_pops_A2 = stack(All_pops_A2_2040_70[[3]],All_pops_A2_2040_70[[3]])



## -----------------------------------------------------------------------------
# # Plotting Range Expansion AB1 2040 a 2100 Dry seasonal
## -----------------------------------------------------------------------------
exp_retract_plot(Mig_clim_future = Dry_ses_AB1[[2]],
                 Present_dist = var_Dry_ses_AB1[[1]],
                 pop_num = 2,
                 title_plot="Dry-seasonal ssp370 2100",
                 deforestation= "default") # or strong_in_2070

exp_retract_plot(Mig_clim_future = Dry_ses_AB1_2040_70[[2]],
                 Present_dist = var_Dry_ses_AB1[[1]],
                 pop_num = 2,
                 title_plot="Dry-seasonal ssp370 2070",
                 deforestation= "default")

exp_retract_plot(Mig_clim_future = Dry_ses_AB1_2040_70[[1]],
                 Present_dist = var_Dry_ses_AB1[[1]],
                 pop_num = 2,
                 title_plot="Dry-seasonal ssp370 2040",
                 deforestation= "default")

# Plotting Range Expansion 2100 Dry seasonal A2

exp_retract_plot(Mig_clim_future = Dry_ses_A2[[2]],
                 Present_dist = var_Dry_ses_A2[[1]],
                 pop_num = 2,
                 title_plot="Dry-seasonal ssp585 2100",
                 deforestation= "default")

exp_retract_plot(Mig_clim_future = Dry_ses_A2_2040_70[[2]],
                 Present_dist = var_Dry_ses_A2[[1]],
                 pop_num = 2,
                 title_plot="Dry-seasonal ssp585 2070",
                 deforestation = "strong_in_2070")

exp_retract_plot(Mig_clim_future = Dry_ses_A2_2040_70[[1]],
                 Present_dist = var_Dry_ses_A2[[1]],
                 pop_num = 2,
                 title_plot="Dry-seasonal ssp585 2040",
                 deforestation= "default")

# Plotting Range Expansion 2100 seasonal no thres
#exp_retract_plot(Mig_clim_future = Dry_ses_AB1_0thres[[2]],
#                 Present_dist = var_Dry_ses_AB1[[1]],
#                 pop_num = 2)

## -----------------------------------------------------------------------------
# # Plotting Range Expansion AB1 2040 a 2100 Warm Wet
## -----------------------------------------------------------------------------

# Plotting Range Expansion 2100 WarmWet AB1
exp_retract_plot(Mig_clim_future = Warm_wet_AB1[[2]],
                 Present_dist = var_warm_wet_AB1[[1]],
                 pop_num = 1,
                 title_plot="Warm Wet ssp370 2100",
                 deforestation="normal")

exp_retract_plot(Mig_clim_future = Warm_Wet_AB1_2040_70[[2]],
                 Present_dist = var_warm_wet_AB1[[1]],
                 pop_num = 1,
                 title_plot="Warm Wet ssp370 2070",
                 deforestation="normal")

exp_retract_plot(Mig_clim_future = Warm_Wet_AB1_2040_70[[1]],
                 Present_dist = var_warm_wet_AB1[[1]],
                 pop_num = 1,
                 title_plot="Warm Wet ssp370 2040",
                 deforestation="normal")

# Plotting Range Expansion 2100 WarmWet A2
exp_retract_plot(Mig_clim_future = Warm_wet_A2[[2]],
                 Present_dist = var_warm_wet_A2[[1]],
                 pop_num = 1,
                 title_plot="Warm Wet ssp585 2100",
                 deforestation="normal")

exp_retract_plot(Mig_clim_future = Warm_Wet_A2_2040_70[[2]],
                 Present_dist = var_warm_wet_AB1[[1]],
                 pop_num = 1,
                 title_plot="Warm Wet ssp585 2070",
                 deforestation = "strong_in_2070")

exp_retract_plot(Mig_clim_future = Warm_Wet_A2_2040_70[[1]],
                 Present_dist = var_warm_wet_AB1[[1]],
                 pop_num = 1,
                 title_plot="Warm Wet ssp585 2040")


## -----------------------------------------------------------------------------
# # Plotting Range Expansion 2040 a 2100 All pops
## -----------------------------------------------------------------------------
exp_retract_plot(Mig_clim_future = All_pops_AB1[[2]],
                 Present_dist = var_All_pops_AB1[[1]],
                 pop_num = c(1,2),
                 title_plot="Warm Wet ssp370 2100")

exp_retract_plot(Mig_clim_future = All_pops_AB1_2040_70[[2]],
                 Present_dist = var_All_pops_AB1[[1]],
                 pop_num = c(1,2),
                 title_plot="Warm Wet ssp370 2070")

exp_retract_plot(Mig_clim_future = All_pops_AB1_2040_70[[1]],
                 Present_dist = var_All_pops_AB1[[1]],
                 pop_num = c(1,2),
                 title_plot="Warm Wet ssp370 2040")
# All pops A2
exp_retract_plot(Mig_clim_future = All_pops_A2[[2]],
                 Present_dist = var_All_pops_A2[[1]],
                 pop_num = c(1,2),
                 title_plot="Warm Wet ssp585 2100")

exp_retract_plot(Mig_clim_future = All_pops_A2_2040_70[[2]],
                 Present_dist = var_All_pops_A2[[1]],
                 pop_num = c(1,2),
                 title_plot="Warm Wet ssp585 2070",
                 deforestation = "strong_in_2070")

exp_retract_plot(Mig_clim_future = All_pops_A2_2040_70[[1]],
                 Present_dist = var_All_pops_A2[[1]],
                 pop_num = c(1,2),
                 title_plot="Warm Wet ssp585 2040")


# Evolutionary rescue function

evol_res_plot = function(Present_dist,Mig_clim_future,Mig_clim_future_B,pop_num,deforestation){

# Distribution in 2100
Mig_clim_future[Mig_clim_future<1]=NA
Present_dist[Present_dist<1]=NA
Mig_clim_future_B[Mig_clim_future_B<1]=NA

if (deforestation == "strong_in_2070") {
  Mig_clim_future = raster::mask(Mig_clim_future,rasters_future_land_cover_agr_proj[[5]],maskvalue=1)
  Mig_clim_future_B = raster::mask(Mig_clim_future_B,rasters_future_land_cover_agr_proj[[5]],maskvalue=1)
} else {
  Mig_clim_future = raster::mask(Mig_clim_future,rasters_future_land_cover_agr_proj[[3]],maskvalue=1)
  Mig_clim_future_B = raster::mask(Mig_clim_future_B,rasters_future_land_cover_agr_proj[[3]],maskvalue=1)
}

dist_2100_d_extinct = mask(Present_dist,Mig_clim_future, inverse=TRUE)

#plot(Mig_clim_future)
#plot(Mig_clim_future_B, col="blue", add=TRUE)
#plot(dist_2100_d_extinct,col="red",add=TRUE)

evol_rescue_warm_to_hot = mask(Mig_clim_future_B,dist_2100_d_extinct)

print(sum(na.exclude(raster::getValues(evol_rescue_warm_to_hot))))

#plot(evol_rescue_warm_to_hot, col="black")

#evol_rescue_hot_to_warm = mask(Mig_clim_future_B,dist_2100_d_extinct)

# Plotting expansion / retraction maps
options(repr.plot.width=7, repr.plot.height=7)
par(mar = c(0, 0, 0, 0))
plot(bio_SPDF3,col=NA,border=NA)
plot(Mig_clim_future,col=c("#053061"),legend=FALSE,add=TRUE)
plot(dist_2100_d_extinct,col=c("#b2182b"),legend=FALSE,add=TRUE)
plot(evol_rescue_warm_to_hot,col="#8073ac",legend=FALSE,add=TRUE)
plot(bio_SPDF3,col=NA,border="#969696", add=TRUE,lwd=1)
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==pop_num,],pch = 21, cex=1.5, col="black", bg="#35978f",lwd=1.5, add=TRUE)
}

# Evolutionary rescue Dry seas to Warm Wet A1B 2100
evol_res_plot(Mig_clim_future = Warm_wet_AB1[[2]],# fut ext
              Present_dist = var_warm_wet_AB1[[1]], #present
              Mig_clim_future_B = Dry_ses_AB1[[2]],# fut sp1 res
              pop_num=1,
              deforestation = "normal")

# Evolutionary rescue Dry seas to Warm Wet A2 2100
evol_res_plot(Mig_clim_future = Warm_wet_A2[[2]], #Warm_wet_A2_weak
              Present_dist = var_warm_wet_A2[[1]],
              Mig_clim_future_B = Dry_ses_A2[[2]], #Dry_ses_A2_weak
              pop_num=1,
              deforestation = "strong_in_2070")

# Evolutionary rescue Dry seas to Warm Wet A1B 2070
evol_res_plot(Mig_clim_future = Warm_Wet_AB1_2040_70[[2]],# fut
              Present_dist = var_warm_wet_AB1[[1]], #present
              Mig_clim_future_B = Dry_ses_AB1_2040_70[[2]],# fut sp1 res
              pop_num=1,
              deforestation = "strong_in_2070")

# Evolutionary rescue Dry seas to Warm Wet A2 2070
evol_res_plot(Mig_clim_future = Warm_Wet_A2_2040_70[[2]],# fut
              Present_dist = var_warm_wet_A2[[1]], #present
              Mig_clim_future_B = Dry_ses_A2_2040_70[[2]],# fut sp1 res
              pop_num=1,
              deforestation = "strong_in_2070")

# --------------------------------
## Calculing some statistics? ----
# --------------------------------

# 1-warm_wet , 2-dry-seas

pop_num = 2
values_ext = raster::extract(rasters_future_land_cover_agr_proj, classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==pop_num,])

sum(values_ext[,1])/length(values_ext[,1]) # 33% of range loss X 0.88%

plot(rasters_future_land_cover_agr_proj[[1]])
plot(classified_ind_molecular_pts[classified_ind_molecular_pts$kmeans==pop_num,],pch=21, add=TRUE)


evo_rescue_calcs = function(Present_dist,Mig_clim_future,Mig_clim_future_B,pop_num,deforestation){
  
  # Distribution in 2100
  Mig_clim_future[Mig_clim_future<1]=NA
  Present_dist[Present_dist<1]=NA
  Mig_clim_future_B[Mig_clim_future_B<1]=NA
  
  if (deforestation == "strong_in_2070") {
    Mig_clim_future = raster::mask(Mig_clim_future,rasters_future_land_cover_agr_proj[[5]],maskvalue=1)
    Mig_clim_future_B = raster::mask(Mig_clim_future_B,rasters_future_land_cover_agr_proj[[5]],maskvalue=1)
  }
  
dist_2100_d_extinct = mask(Present_dist,Mig_clim_future, inverse=TRUE)
  
evol_rescue_warm_to_hot = mask(Mig_clim_future_B,dist_2100_d_extinct)

evol_rescue_sum = sum(na.omit(values(evol_rescue_warm_to_hot)))

extinct_sum = sum(na.omit(values(dist_2100_d_extinct)))-evol_rescue_sum

remaining_sum = sum(na.omit(values(Present_dist))) - abs(extinct_sum-evol_rescue_sum)

return(list(evol_rescue_sum=evol_rescue_sum,remaining_sum=remaining_sum,extinct_sum=extinct_sum))
}

evo_rescue_calcs_A2 = evo_rescue_calcs(
  Mig_clim_future = Warm_Wet_A2_2040_70[[2]],# fut
  Present_dist = var_warm_wet_A2[[1]], #present
  Mig_clim_future_B = Dry_ses_A2_2040_70[[2]],# fut 
  deforestation = "strong_in_2070")




?ecospat::ecospat.CCV.communityEvaluation.prob

