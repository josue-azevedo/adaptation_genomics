# Install packages, version checking ####
#R.Version()

##
#packageVersion("LEA")

## Load RData
#load(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_RDA/RDA_Ken_calcarata.RData")

##
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("LEA")
#.####

##
# Loading libraries ####
library(cowplot)
library(dplyr)
library(ggplot2)
library(magrittr)
library(plyr)
library(tidyr)
library(ggpubr)
library(vioplot)
# Load main SIG libraries
library(dismo)
library(rgdal)
library(maptools)
library(SDMTools)
library(spdep)
library(sdm)
library(raster)
library(rgeos)
library(igraph)
library(scatterpie)
library(sf)
library(ggnewscale)
library(psych)    # Used to investigate correlations among predictors
library(usdm)     # Used to run VIF
library(pegas)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(robust)
library(WMDB)
#library(ggVennDiagram)
#detach("package:ggVennDiagram", unload=TRUE) # errors with sp, rgeos actvie
library(corrplot)
library(maps)
library(corrplot)
library(MASS)
library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(viridis)
library(devtools)
library(NbClust)
#install_github('NESCent/MINOTAUR', build_vignettes=TRUE)
#install.packages('MINOTAUR')
#library(MINOTAUR)

#install.packages(c("robust","WMDB","ggVennDiagram"))
#install.packages("usdm",dependencies = TRUE)
#detach("package:usdm", unload=TRUE)
#.####

##
# Set WD and defining general objects ----
setwd("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata")

#sp = "Ken_striata"
sp = "Ken_calcarata"
#sp_name = "Kentropyx striata"
sp_name = "Kentropyx calcarata"
path = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses"

# Settings for each species (copied from I. Prates script)
if (sp == "Ken_striata") {
  sp.short = "Ken_" # That's how I named species during GBS library sequencing
  n = 11
  bestK = 2 # From sNMF clustering analyses
  #K = 2 # Best K for controlling false-discovery rates based on genomic inflation factor (see below)
} else {
  sp.short = "Ken_" # That's how I named species during GBS library sequencing
  n = 66
  bestK = 6 # From sNMF clustering analyses
  #K = 5 # Best K for controlling false-discovery rates based on genomic inflation factor (see below)
}

#.#####

# Data projections #####
oldproj <- paste0(" +init=epsg:4326") #this is WGS84 most commonly used for google earth etc. in decimal degrees
behrmannCRS <- CRS('+proj=cea +lat_ts=30') # This is an equal area projection, i.e., cells have the same area
#.####

##
# Importing Species points and some bioclimatic data ####
preds.Chelsa <- list.files(path='/Users/josue/Dropbox/1Doutorado/bromeliacea/chelsa_aggregated',pattern='tif',full.names = T) 
#preds.Chelsa # 19 variaveis na lista
#preds <- stack(preds.Chelsa[c(1:3,12:19)]) # selecting temperature only
preds <- stack(preds.Chelsa) # Include all variables

# simplify variable names
names(preds)= gsub("CHELSA_bio10","bio",names(preds))

# Reading species points
kentropyx_RAD_coor = read.csv(paste0("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/kentropyx_analyses/","kentropyx_RAD_coor_",sp, ".csv"), header = TRUE,row.names = 1)
head(kentropyx_RAD_coor)
coordinates(kentropyx_RAD_coor) = kentropyx_RAD_coor[,c("long","lat")]
#plot(kentropyx_RAD_coor)
#head(kentropyx_RAD_coor)

## Addin "Ken_" to the ID list to match the SNPs table # Peculiar to my own dataset
kentropyx_RAD_coor@data$Sample_ID_full = as.factor(paste0("Ken_",kentropyx_RAD_coor@data$Sample_ID_red))

# Creating spatial dataframe with the same projection as the clim raster
proj4string(kentropyx_RAD_coor) = oldproj
kentropyx_RAD_coor <-spTransform(kentropyx_RAD_coor, CRSobj=behrmannCRS)

#head(kentropyx_RAD_coor@coords)
kentropyx_RAD_coor = cbind(kentropyx_RAD_coor$species,as.data.frame(kentropyx_RAD_coor@coords),kentropyx_RAD_coor$Sample_ID_red,kentropyx_RAD_coor@data$Sample_ID_full,kentropyx_RAD_coor$lat.1,kentropyx_RAD_coor$long.1)
coordinates(kentropyx_RAD_coor)=kentropyx_RAD_coor[,2:3]

#head(kentropyx_RAD_coor)
colnames(kentropyx_RAD_coor@data)= c("species","long","lat","Sample_ID_red","Sample_ID_full","lat_wgs","long_wgs")

head(kentropyx_RAD_coor)
#.####
write.csv(kentropyx_RAD_coor,"/Users/josue/Library/CloudStorage/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_1_outfiles_090/VCFtools_SNMF_05ken_calcarata090/coordinates.csv")
##
# Getting shapefiles from distinct directories in my computer ####

## American continent map
worldLowres <- readOGR("/Users/josue/Dropbox/1Doutorado/Chapter_2/biodiverse_pipeline-master/ne_110m_admin_0_countries.shp")
#drawExtent()

## Detailed SAmerican  map
SAshp <- readOGR("/Users/josue/Dropbox/1Doutorado/Chapter_Wer/South_America/South_America.shp")

## Terrestrial ecoregions
ecor <- rgdal::readOGR("~/Dropbox/1Doutorado/Chapter 1/cerrado_biogeo/shapes_cer/wwf_terr_ecos.shp") ## WWF shapefile 

## Cerrado background
cerrado<- ecor[ecor@data$ECO_NAME=="Cerrado",]
#plot(cerrado, add=T, col= "yellow")

## Caatinga background
Caatinga<- ecor[ecor@data$ECO_NAME=="Caatinga",]

## Guiana background
 Guiana<- ecor[ecor@data$ECO_NAME=="Guianan savanna",]

## Llanos background
Llanos<- ecor[ecor@data$ECO_NAME=="Llanos",]

## Major rivers
rivers<- rgdal::readOGR("/Users/josue/Dropbox/1Doutorado/Chapter_Wer/Kentropyx/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp")

## Major lakes
lakes<- rgdal::readOGR("/Users/josue/Dropbox/1Doutorado/Chapter_Wer/Kentropyx/ne_10m_lakes_pluvial/ne_10m_lakes_pluvial.shx")

##
## Shapefile amazonia
bio_SPDF<- readOGR("/Users/josue/Dropbox/1Doutorado/Chapter_2/biodiverse_pipeline-master/wwf_biomes_simp.shp")
bio_SPDF <-spTransform(bio_SPDF, CRSobj=behrmannCRS)
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
plot(bio_SPDF2)

florests = bio_SPDF2[bio_SPDF2@data$BIOME==1,]
plot(florests)

other_biomes = bio_SPDF2[bio_SPDF2@data$BIOME %in% c(7,2,13),]
plot(other_biomes)

##
# Extract individual polygon features (ex: tropical forests)

florests = bio_SPDF2[bio_SPDF2@data$BIOME==1,]
#plot(florests)

other_biomes = bio_SPDF2[bio_SPDF2@data$BIOME %in% c(7,2,13),]
#plot(other_biomes)

# Rasterise
florests_ras= crop(rasterize(florests,preds[[1]]),florests)
#plot(florests_ras)
#plot(kentropyx_RAD_coor, add=TRUE)

other_biomes_ras= crop(rasterize(other_biomes,preds[[1]]),florests)
plot(other_biomes_ras)
plot(kentropyx_RAD_coor, add=TRUE)

##
## Distance to the border of non-forested biomes
rbiome2=crop(mask(preds[[1]],other_biomes),florests)
rborder <- rasterize(st_cast(st_as_sf(other_biomes), "MULTILINESTRING"),rbiome2)
plot(rborder)

##
dborder <- raster::distance(rborder)
plot(dborder)

##
dborder2 = mask(dborder,florests_ras)
dborder3=dborder2
dborder3[is.na(dborder3)] <- 0
plot(dborder3)     ## If you want the log of the distance to the biome border
plot(log(dborder3+100000))     ## If you want the log of the distance to the biome border
plot(kentropyx_RAD_coor, add=TRUE)

##
writeRaster(dborder3,'/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/dborder_non_log.tif')

# Or read it in 
dborder3 = raster('/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/dborder_non_log.tif')

#.####

##
# Read files from VCF tools and filtering SNPs ####

# Data from first iPyrad run
#imputed_genotipes = read.geno("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/LEA_LFMM_Ken_calcarata_66_pc1/genotypes.geno")
#genotipes[1:10,1:10]

##
# Reading genomic data "original version (changing 9 by NA)
#genotipes = read.geno("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/LEA_LFMM_Ken_calcarata_66_pc1/genotypes.geno") # Original dataset

# Original (not processed to geno) - 0.85, MAF = 010
#genotipes = as.matrix(read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/VCFtools_SNMF_ken_calcarata/ken_calcarata_1SNP-locus.012",row.names = 1))

# Stacks 080
#genotipes = as.matrix(read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/stacks_outputs/stacks_3/out_stacks_3_050/VCFtools_SNMF_05populations.snps/populations.snps_1SNP-locus.012",row.names = 1))

# Original (not processed to geno) - 0.85, MAF = 005
#genotipes = as.matrix(read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/VCFtools_SNMF_Nken_calcarata/ken_calcarata_1SNP-locus.012",row.names = 1))

# MAF=0.05, clustering threshold 0.9 in iPyrad - highest value
genotipes = as.matrix(read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_1_outfiles_090/VCFtools_SNMF_05ken_calcarata090/ken_calcarata090_1SNP-locus.012",row.names = 1))

# MAF=0.1, clustering threshold 0.9 in iPyrad - 
#genotipes = as.matrix(read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_1_outfiles_090/VCFtools_SNMF_ken_calcarata090/ken_calcarata090_1SNP-locus.012",row.names = 1))
# Dataset with new choices of clustering threshold and MAF

# Change 9 to NA
genotipes[1:10,1:10]
genotipes[genotipes == 9] = NA
genotipes[1:10,1:10]
##

# Reading list of sample IDs from vcftools (imputed genotipes is the SNP data and individuals is the ID of each column)

# Original data
#individuals = read.table(file = paste0(path, "/VCFtools_SNMF_", sp, "/", sp, "_1SNP-locus.012.indv")) # read list of sample IDs from vcftools
#dim(individuals)

# Stacks 080
#individuals = read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/stacks_outputs/stacks_3/out_stacks_3_050/VCFtools_SNMF_05populations.snps/populations.snps_1SNP-locus.012.indv") # 

# Original data MAF 005
#individuals = read.table(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/VCFtools_SNMF_Nken_calcarata/ken_calcarata_1SNP-locus.012.indv") # read list of sample IDs from vcftools
#dim(individuals)

# New data MAF 0.05, clust 090
individuals = read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_1_outfiles_090/VCFtools_SNMF_05ken_calcarata090/ken_calcarata090_1SNP-locus.012.indv") # From 0.90 threshold

# New data MAF 0.1, clust 090
#individuals = read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_1_outfiles_090/VCFtools_SNMF_ken_calcarata090/ken_calcarata090_1SNP-locus.012.indv") # From 0.90 threshold
nrow(individuals)

##
# change row names to sample IDs
#individuals$V1 %<>% gsub(pattern = sp.short, replacement = "", .) # removes "Ken_" if present. The "." denotes object (individuals) when using "%<>%"
rownames(genotipes) = individuals$V1 # change row names to sample IDs
genotipes[1:10,1:10]
dim(genotipes)

##
# Create list of loci that will serve as a SNP map when retrieving sequences to blast later

# Original MAF 0.1 - 085 
#loci = read.table(file = paste0(path, "/VCFtools_SNMF_", sp, "/", sp, "_1SNP-locus.012.pos")) # read list of locus IDs from VCFtools

# Stacks 080
#loci = read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/stacks_outputs/stacks_3/out_stacks_3_050/VCFtools_SNMF_05populations.snps/populations.snps_1SNP-locus.012.pos") 

# Original MAF 0.05 - 085 
#loci = read.table(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/VCFtools_SNMF_Nken_calcarata/ken_calcarata_1SNP-locus.012.pos") # read list of locus IDs from VCFtools

# MAF 0.05 - 090
loci = read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_1_outfiles_090/VCFtools_SNMF_05ken_calcarata090/ken_calcarata090_1SNP-locus.012.pos") # from 090

# MAF 0.1 - 090
#loci = read.table("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_1_outfiles_090/VCFtools_SNMF_ken_calcarata090/ken_calcarata090_1SNP-locus.012.pos") # from 090

names(loci) = c("locus.name", "SNP.position") # add column names

loci$locus.name %<>% gsub(pattern = "locus_", replacement = "", .) # removes "Ken_" if present. The "." denotes object (individuals) when using "%<>%"

loci$SNP.number = 1:nrow(loci) # adding a column with consecutive numbers to the loci object
#loci # list loci
length(unique(loci$locus.name)) # How many loci?

##
# Adding colnames to object
colnames(genotipes) = loci$locus.name
# colnames(genotipes) = paste0("RAD_", loci$locus.name)
genotipes[1:10,1:10]
dim(genotipes)

# Removing spp with lots of missing data (not needed for K. calcarata)
#genotipes = genotipes[-which(rownames(genotipes)=="Ken_126"),]
# Did not make a difference when removing SNPs with lots of missing data in ipyrad, but it did with Stacks
#individuals = individuals[-which(individuals$V1=="Ken_126"),]
#nrow(individuals)
#kentropyx_RAD_coor = kentropyx_RAD_coor[-which(kentropyx_RAD_coor$Sample_ID_full=="Ken_126"),]
#nrow(kentropyx_RAD_coor)

# % of missing data per individuals
as.data.frame(sort(apply(genotipes, 1, function(x) sum(is.na(x)))/ncol(genotipes)),col.names=c("sp","%"))

summary(apply(genotipes, 1, function(x) sum(is.na(x)))/ncol(genotipes))

# Checking % per SNP
sum(apply(genotipes, 2, function(x) sum(is.na(x))/nrow(genotipes)) <=0.2)
summary(apply(genotipes, 2, function(x) sum(is.na(x)))/nrow(genotipes))

# Removing SNPs with more than 20% of missing data
SNP_to_remove = which(apply(genotipes, 2, function(x) sum(is.na(x))/nrow(genotipes)) <=0.2)

genotipes_t = genotipes[,SNP_to_remove]
dim(genotipes)
dim(genotipes_t)

# Check again missing data per individuals
as.data.frame(sort(apply(genotipes_t, 1, function(x) sum(is.na(x)))/ncol(genotipes_t)),col.names=c("sp","%"))

# I will remove here Ken_126 for testing ##
#genotipes_t[1:10,1:10]
#genotipes_t = genotipes_t[-which(rownames(genotipes_t)=="Ken_126"),]
#dim(genotipes_t)

# If Ken_126 is remove, remove it from coords
#kentropyx_RAD_coor = kentropyx_RAD_coor[-which(kentropyx_RAD_coor$Sample_ID_full=="Ken_126"),]
#kentropyx_RAD_coor$Sample_ID_full == "Ken_87"
###

# Only one individual presented high values of missing data. I will test sNMF with and without it (didnt make any difference)
plot(log(dborder3+100000))     ## If you want the log of the distance to the biome border
plot(kentropyx_RAD_coor, pch=16, add=TRUE)
plot(kentropyx_RAD_coor[kentropyx_RAD_coor$Sample_ID_full=="Ken_126",], pch=16, add=TRUE, col="red")

# Ken_126 is highly removable
# However, didn't change results from the original 0.85 threshold with and without (when removing >0.2 missing data in loci)
#.####

# Run sNMF  #####
# Create a directory to place resulting SNMF bar plots
dir.create(paste0(path, "/SNMF_all_SNPs_new2","_",sp))
dir.exists(paste0(path, "/SNMF_all_SNPs_new2","_",sp))

# Change to the new directory
setwd(paste0(path, "/SNMF_all_SNPs_new2","_",sp))

# Saving and/or loading data
#save.image(file="/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/rdata_beggining_var_sel2.rdata")

#load(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/rdata_beggining_var_sel.rdata")

# Check genetic data, change back NA to 9
genotipes_t[is.na(genotipes_t)] = 9
genotipes_t[1:10,1:10]

# Write down data in LFMM format (LEA saves it all outside R, instead of as objects)
write.lfmm(genotipes_t, paste0(sp, "_genotypes_15.lfmm"))
dim(genotipes_t)

# All project names and trials bellow
#_genotypes_1.lfmm = orginal 085 without filtering
#_genotypes_2.lfmm = orginal 085 with <0.2 missing data
#_genotypes_3.lfmm = orginal 085 with <0.2 missing data without Ken_126 - same k
# genotypes_4.lfmm = 0.90 with <0.2 missing 
# genotypes_5.lfmm = 0.90 with MAF 0.1 <0.2 missing
# genotypes_6.lfmm = 085, filtering, MAF 0.1, < 075 missing data
# genotypes_7.lfmm = 0.90 with MAF 0.1 <0.25 missing
# genotypes_8.lfmm = 085, MAF01, 025 missing
# genotypes_10.lfmm = 085, MAF01, 025 missing # try again
# genotypes_11.lfmm = Stacks 080, 020 missing, no Ken_126 (65 samples)
# genotypes_12.lfmm = Stacks 090, 025 missing, with Ken_126 (65 samples)
# genotypes_14.lfmm = Stacks 085, 010 missing, with Ken_126 (65 samples), MAF 005
# genotypes_15.lfmm = Stacks 085, 020 missing, with Ken_126 (65 samples), MAF 005


##
## Part 2: Running genetic clustering analyses

# Run sNMF (started at 10:30)
#project.snmf_2 = NULL
#project.snmf_4 = NULL
#project.snmf_5 = NULL
#project.snmf_6 = NULL
#project.snmf_7 = NULL
#project.snmf_8 = NULL
#project.snmf_10 = NULL
#project.snmf_11 = NULL
#project.snmf_12 = NULL
#project.snmf_13 = NULL
#project.snmf_14 = NULL
#project.snmf_15 = NULL

project.snmf_15 = snmf("Ken_calcarata_genotypes_15.lfmm", K = 1:10, entropy = TRUE, repetitions = 100, project = "new", CPU = 6)

# If you need to read it in later
#project.snmf_4 = load.snmfProject("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/Ken_calcarata_genotypes_4.snmfProject")

# plot cross-entropy criterion of all runs of the project
#png(filename = paste0(path, "LEA/plots_SNMF_all_SNPs/", sp, "_all_SNPs_entropy_plot_t", t, "_s", s, "_n", n, "_a", a, ".png"))
plot(project.snmf_4, cex = 1.2, col = "blue", pch = 19)
#dev.off() # saving plot as .png image

# Set best K based on the plot
bestK = 6

# get the cross-entropy of the 10 runs for K = best K
ce = cross.entropy(project.snmf_4, K = bestK)
ce

# select the run with the lowest cross-entropy for K = best K
bestrun = which.min(ce)
bestrun
median(ce)
# 5 - 0.502
# 6 - 0.500
# 7 - 0.502

# Entropy value corresponding to best run for best K
e = round(ce[bestrun], digits = 3) # rounding to four digits
e

## SNMF bar plot ####
# Part 3: Plotting ancestry coefficients as a bar plot based on the Q matrix estimated by SNMF
# get the Q matrix, for the best run and K
qmatrix = LEA::Q(project.snmf_4, K = bestK, run = bestrun)
dim(qmatrix)
qmatrix = as.data.frame(qmatrix)
colnames(qmatrix) = c(paste0(rep("cluster", bestK), 1:bestK)) # change column names for cluster name. Will repeat "cluster" K times, then paste0 with second element from 1 to K times
head(qmatrix)
dim(qmatrix)

# Adding ID names for plotting
rownames(qmatrix) = rownames(genotipes_t) # change row names to sample IDs
write.csv(qmatrix, file = paste0("qmatrix_", sp, ".csv"), quote = FALSE, row.names = FALSE)
head(qmatrix)
dim(qmatrix)

#qmatrix = read.csv(paste0("qmatrix_", sp, ".csv"), header = TRUE)

# Caso queira um plot rápido, use a função barchart:
# Colorbrewer colors - Change the order in my colors accordingly to your preference in the maps
#e41a1c - red
#377eb8 - blue
#4daf4a - green
#737373 - grey
#ff7f00 - orange
#984ea3 - purple
##980043- pink

my.colors <- c('#e41a1c','#377eb8','#4daf4a','#737373','#ff7f00','#984ea3')
indi = gsub("Ken_","", rownames(genotipes_t))

#pdf("barchart_norml.pdf",paper = "a4r",width = 15)
LEA::barchart(project.snmf_4, K = bestK, run = bestrun, 
         border = TRUE, space = 0, col = my.colors, 
         xlab = "Individuals", ylab = "Ancestry proportions", 
         main = "Ancestry matrix") -> bp

axis(1, at = 1:length(bp$order), 
     labels = indi[bp$order], las = 3, 
     cex.axis = .4)
#dev.off()

## SNMF Map of genomic clustering with kernel density ####
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

asc.raster="/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/kentropyx_analyses/SA_raster.asc"
grid=createGridFromAsciiRaster(asc.raster)
constraints=getConstraintsFromAsciiRaster(asc.raster, cell_value_min=0)

# Set Qmatrix order as the same in the bp$order
qmatrix_ordered = qmatrix[bp$order,]
dim(qmatrix_ordered)
# Check the order
kentropyx_RAD_coor$Sample_ID_full == rownames(qmatrix_ordered)

# Check name compatibility
setdiff(rownames(qmatrix_ordered),kentropyx_RAD_coor@data$Sample_ID_full)

# Reorder kentropyx_RAD_coor
kentropyx_RAD_coor = kentropyx_RAD_coor[order(match(kentropyx_RAD_coor$Sample_ID_full,rownames(qmatrix_ordered))),]

# Check the order
kentropyx_RAD_coor@data$Sample_ID_full == rownames(qmatrix_ordered)

# Now add the coords
head(kentropyx_RAD_coor@data)
coord.at = kentropyx_RAD_coor@data[,c("long_wgs","lat_wgs")]

# This one has to be in the same order as my.colors from ColourBrewer
lColorGradients = list(
  c("gray95",brewer.pal(9,"Reds")),
  c("gray95",brewer.pal(9,"Blues")),
  c("gray95",brewer.pal(9,"Greens")),
  c("gray95",brewer.pal(9,"Greys")),
  c("gray95",brewer.pal(9,"YlOrBr")),
  c("gray95",brewer.pal(9,"Purples")),
  c("gray95",brewer.pal(9,"RdPu")))

#e41a1c - red
#377eb8 - blue
#4daf4a - green
#737373 - grey
#ff7f00 - orange
#984ea3 - purple

maps(matrix = qmatrix_ordered, coord.at, grid, constraints, method = "max",
     main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5)

map(add = T, interior = F)
points(kentropyx_RAD_coor$long_wgs,kentropyx_RAD_coor$lat_wgs, pch = 21, cex=1.2, col="black", bg="white", lwd=1)

#text(kentropyx_RAD_coor$long_wgs,kentropyx_RAD_coor$lat_wgs,kentropyx_RAD_coor$Sample_ID_red)

# If you want to highlight some of the individuals
#kentropyx_RAD_coor_doido = kentropyx_RAD_coor[kentropyx_RAD_coor$Sample_ID_full=="Ken_67",]
#points(kentropyx_RAD_coor_doido$long_wgs,kentropyx_RAD_coor_doido$lat_wgs, pch = 21, cex=1, col="black", bg="red", lwd=1)

## Pie-charts plot

# Always checking the order of barplot with 
rownames(qmatrix_ordered) == kentropyx_RAD_coor$Sample_ID_full

RAD_seq_selector = cbind.data.frame(qmatrix_ordered,kentropyx_RAD_coor$long_wgs,kentropyx_RAD_coor$lat_wgs)

colnames(RAD_seq_selector) = c(colnames(RAD_seq_selector)[1:ncol(qmatrix_ordered)],"long","lat")

# Creating a column to designate which cluster each individual was classified (based in the maximum value of in each cluster)
RAD_seq_selector$class = colnames(RAD_seq_selector[1:bestK])[max.col(RAD_seq_selector[1:bestK],ties.method="first")]

head(RAD_seq_selector)

#names(RAD_seq_selector[!duplicated(RAD_seq_selector$class),]) # Checked - same order as in barplot

# Then, to match colors in the pie chart, you have to change the order from the qmatrix/RAD_seq_selector to go from cluster1 to cluster 6 in order
RAD_seq_selector = RAD_seq_selector[order(RAD_seq_selector$class),]

head(RAD_seq_selector)

# This one has to be in the same colors as in the barplot (right-to-left)
#my.colors_pie_chart = my.colors[as.numeric(gsub("cluster","" ,)] # reordering my.colors by the cluster order in the reordered qmatrix/rad_seq_selector

#all_objects = list(SAshp,cerrado,Caatinga,Guiana,Llanos,rivers,SAshp)
#saveRDS(all_objects, "/Users/josue/Downloads/all_layers.rds")

ggplot()+
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill="chartreuse4",col=NA)+
  geom_polygon(data=cerrado, aes(long, lat, group=group), fill="gold2")+
  geom_polygon(data=Caatinga, aes(long, lat, group=group), fill="goldenrod3")+
  geom_polygon(data=Guiana, aes(long, lat, group=group), fill="gold2")+
  geom_polygon(data=Llanos, aes(long, lat, group=group), fill="gold2")+
  geom_path(data=rivers, aes(long, lat, group=group), col="dodgerblue3")+
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill=NA,col="grey20",size = .2)+
  #geom_point(data=RAD_seq_selector,aes(x=long ,y=lat), col="black", fill="white", size=4, shape=21,height= 0.4,width=0.4)+
  geom_scatterpie(aes(x=long, y=lat, group = long, r =.8),col = 'grey25', size = .2,
                  data = RAD_seq_selector, cols = colnames(RAD_seq_selector[,c(1:bestK)])) +
  scale_fill_manual(values=my.colors)+
  #geom_text(data=kentropyx_RAD_coor@data,aes(label = Sample_ID_red, x = long, y = lat),size=3)+
  coord_cartesian(xlim=c(-67,-34),ylim = c(-17.5,7))+
  ggsn::scalebar(x.min=-67, x.max=-34, y.min=-17.5, y.max=7,dist = 250,dist_unit = "km", model = 'WGS84',st.size = 2.5,location="topright",transform=TRUE, border.size=.5)+
  #annotation_custom(g, xmin=-50, xmax=-40, ymin=-1.5, ymax=5.57)+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'grey75',size = .1),
        panel.background = element_rect(fill = 'white'))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  #ggtitle(my_title)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(colour = 'black', fill=NA, size=.5))

# To add in future versions
locator()
geom_text(aes(label = rownames(df)),
          size = 3.5)
#.####

# Impute missing data ####
#Now we can impute our missing date based on the number of clusters using LEA instead of just using the most common value for each SNP as in the RDA tutorial #

#genotipes_impute = impute(object=project.snmf_15, input.file="Ken_calcarata_genotypes_15.lfmm", K=bestK, run=bestrun) 

# Now I will convert the newly imputed genotype to geno

#output = lfmm2geno("Ken_calcarata_genotypes_4.lfmm_imputed.lfmm", paste0(sp, "_genotypes_imputed_4.geno"))

# Now I will read back my imputed genotypes list
AllFreq = read.geno(paste0(sp, "_genotypes_imputed_4.geno"))

rownames(AllFreq) = rownames(genotipes_t)
colnames(AllFreq) = colnames(genotipes_t)
AllFreq[1:10,1:10]

###
# Removing NAs of  genotypes assigning maximum value (original way)
#AllFreq <- apply(genotipes_t, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
#sum(is.na(AllFreq)) # No NAs
#AllFreq[1:10,1:10]
dim(AllFreq)
#.####

# Delimiting study area ####
## I will select my study area to clip the polygons
selected_clipped_biomes<- readOGR("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/layers/selected_clipped_biomes3.shp")
# remove invalidities in the county map
gIsValid(selected_clipped_biomes) #returns FALSE: there are invalid self-intersecting geometries in the polygons, which will cause problems
selected_clipped_biomes <- gBuffer(selected_clipped_biomes, byid = TRUE, width = 0)
gIsValid(selected_clipped_biomes) #returns TRUE

# dissolve county map by rating area & fortify to data frame
selected_clipped_biomes@data$union = 1 # create a new field with identical values
selected_clipped_biomes <- unionSpatialPolygons(selected_clipped_biomes, IDs = selected_clipped_biomes@data$union)
#selected_clipped_biomes <- fortify(selected_clipped_biomes)
#plot(selected_clipped_biomes)

# Remove holes in the polygon
selected_clipped_biomes_2 <- sfheaders::sf_remove_holes(st_as_sf(spTransform(selected_clipped_biomes, CRSobj=behrmannCRS)))
plot(selected_clipped_biomes_2)

#rivers_2 = rivers
#rivers_3 = crop(rivers_2,SAshp)
#plot(rivers_3)
#writeOGR(rivers_3,"/Users/josue/Dropbox/1Doutorado/Chapter_Wer/Kentropyx/ne_10m_rivers_lake_centerlines/South_America_rivers.shp" ,driver="ESRI Shapefile",layer=rivers_3@data$name, overwrite_layer = TRUE)
#plot(SAshp, add=TRUE)

#.####

# RDA data prep (swiss knife) ####

# List predictors
## Bioclimatic variables definition ####
#Bio1 = Annual Mean Temperature
#Bio2 = Mean Diurnal Range
#Bio3 = Isothermality
#Bio4 = Temperature Seasonality (sd of monthly temps)
#Bio5 = Max Temperature of Warmest Month
#Bio6 = Min Temperature of Coldest Month
#Bio7 = Temperature Annual Range (min-max mean temps of the year)
#Bio8 = Mean Temperature of Wettest Quarter
#Bio9 = Mean Temperature of Driest Quarter
#Bio10 = Mean Temperature of Warmest Quarter
#Bio11 = Mean Temperature of Coldest Quarter
#Bio12 = Annual Precipitation
#Bio13 = Precipitation of Wettest Month
#Bio14 = Precipitation of Driest Month
#Bio15 = Precipitation Seasonality
#Bio16 = Precipitation of Wettest Quarter
#Bio17 = Precipitation of Driest Quarter
#Bio18 = Precipitation of Warmest Quarter
#Bio19 = Precipitation of Coldest Quarter

## Preparing environmental rasters ####
preds # Preds from Chelsa V.1
preds_chelse_clip = mask(crop(preds,selected_clipped_biomes_2),selected_clipped_biomes_2)
plot(preds_chelse_clip[[1]])
plot(kentropyx_RAD_coor,pch=16,cex=0.5, add=TRUE)
plot(selected_clipped_biomes_2,add=TRUE)

# Preds from Climond
# Primeiro liste todos os rasters ambientais por nome. path= indica o diretório em que as variáveis estão
preds.climond <- list.files(path="/Users/josue/Dropbox/4Environmental_layers/current_climond/CM10_1975H_Bio_V1.2/",pattern='txt',full.names = T) 
#preds.climond
preds.climond = preds.climond[1:19] # selecting only traditional bioclimatic variables
#preds.climond

preds.climond_stck1 = stack(preds.climond)
proj4string(preds.climond_stck1) = oldproj
preds.climond_stck = projectRaster(preds.climond_stck1,res=16700,crs = behrmannCRS)
res(preds.climond_stck) = 16700
preds.climond_stck2 = mask(crop(preds.climond_stck, selected_clipped_biomes_2),selected_clipped_biomes_2)
plot(preds.climond_stck2[[5]])
#preds.climond_stck3 = raster::resample(preds.climond_stck2,preds_chelse_clip[[1]]) # just to have the same resolution as the chelsa rasters
#preds.climond_stck4 = mask(preds.climond_stck3,selected_clipped_biomes_2)

plot(preds.climond_stck2[[1]])
plot(kentropyx_RAD_coor,pch=16,cex=0.5, add=TRUE)
plot(selected_clipped_biomes_2,add=TRUE)

# To make all cells matching across all analyses
preds_modelling_mask = preds.climond_stck2[[1]]
#preds.climond_stck
preds.climond_stck_n = gsub("CM10_1975H_","",names(preds.climond_stck2))
# To check the order of bioclimatic variables before renaming them
gsub("_V1.2","",preds.climond_stck_n)
names(preds.climond_stck2) = c('bio_1','bio_2','bio_3','bio_4','bio_5','bio_6','bio_7','bio_8','bio_9','bio_10','bio_11','bio_12','bio_13','bio_14','bio_15','bio_16','bio_17','bio_18','bio_19')

# Getting preds from Chelsa V.2
preds.chelsa.2 <- list.files(path="/Users/josue/Dropbox/4Environmental_layers/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/",pattern='tif',full.names = T) 
preds.chelsa.2 = grep("CHELSA_bio",preds.chelsa.2, value=TRUE)
preds.chelsa.2_stck1 = stack(preds.chelsa.2)
preds.chelsa.2_crop1 = raster::crop(preds.chelsa.2_stck1, projectRaster(preds_modelling_mask,crs = oldproj))
preds.chelsa.2_project = projectRaster(preds.chelsa.2_crop1,crs = behrmannCRS)
preds.chelsa.2_mask1 = raster::mask(preds.chelsa.2_project,selected_clipped_biomes_2)
preds.chelsa.2_aggregate = raster::aggregate(preds.chelsa.2_mask1,12)
preds.chelsa.2_res = raster::resample(preds.chelsa.2_aggregate,preds_modelling_mask)
preds.chelsa.2_var = raster::mask(preds.chelsa.2_res,selected_clipped_biomes_2)

plot(preds.chelsa.2_var[[1]])

names(preds.chelsa.2_var) = gsub("bio","bio_",gsub("_1981.2010_V.2.1","",gsub("CHELSA_","",names(preds.chelsa.2_var))))
preds.chelsa.2_var

# I now want to scale all rasters to compare them
preds.climond_stck2
preds.chelsa.2_var
preds_chelse_scale = raster::scale(preds.chelsa.2_var)
preds_climond_scale = raster::scale(preds.climond_stck2)

# Difference Chelsa and climond
plot(preds_chelse_scale[["bio_18"]], main="chelsa")
plot(preds_climond_scale[["bio_18"]],main="climond")
plot(preds_chelse_scale[["bio_1"]]-preds_climond_scale[["bio_1"]], main="chelsa - climond")

# With ndvi_3
preds.climond_stck3 = stack(preds.climond_stck2,ndvi)
names(preds.climond_stck3) = gsub("_from_qgis_median","",names(preds.climond_stck3))

preds.chelsa.3_var = stack(preds.chelsa.2_var,ndvi)
names(preds.chelsa.3_var) = gsub("_from_qgis_median","",names(preds.chelsa.3_var))
names(preds.chelsa.3_var)

# Getting processed bioclimatic variables from files already saved in R ####
my_bioclim_rasters = readRDS("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/my_bioclim_rasters.rds")

#my_bioclim_rasters[[1]] = CLIMOND
#my_bioclim_rasters[[2]] = CHELSA

preds.climond_stck3 = my_bioclim_rasters[[2]]
#preds.climond_stck3 = dropLayer(preds.climond_stck3,"bio_19")
#preds.climond_stck3 = dropLayer(preds.climond_stck3,"NDVI")
names(preds.climond_stck3)

## Log transform precipitation variables (not used)
#plot(preds.climond_stck3[["bio_12"]])
#preds.climond_stck3[["bio_12"]] = log(preds.climond_stck3[["bio_12"]]+1)
#plot(preds.climond_stck3[["bio_12"]])

#plot(preds.climond_stck3[["bio_13"]])
#preds.climond_stck3[["bio_13"]] = log(preds.climond_stck3[["bio_13"]]+1)
#plot(preds.climond_stck3[["bio_13"]])

#plot(preds.climond_stck3[["bio_14"]])
#preds.climond_stck3[["bio_14"]] = log(preds.climond_stck3[["bio_14"]]+1)
#plot(preds.climond_stck3[["bio_14"]])

#plot(preds.climond_stck3[["bio_15"]])
#preds.climond_stck3[["bio_15"]] = log(preds.climond_stck3[["bio_15"]]+1)
#plot(preds.climond_stck3[["bio_15"]])

#plot(preds.climond_stck3[["bio_16"]])
#preds.climond_stck3[["bio_16"]] = log(preds.climond_stck3[["bio_16"]]+1)
#plot(preds.climond_stck3[["bio_16"]])

#plot(preds.climond_stck3[["bio_17"]])
#preds.climond_stck3[["bio_17"]] = log(preds.climond_stck3[["bio_17"]]+1)
#plot(preds.climond_stck3[["bio_17"]])

#plot(preds.climond_stck3[["bio_18"]])
#preds.climond_stck3[["bio_18"]] = log(preds.climond_stck3[["bio_18"]]+1)
#plot(preds.climond_stck3[["bio_18"]])

#plot(preds.climond_stck3[["bio_19"]])
#preds.climond_stck3[["bio_19"]] = log(preds.climond_stck3[["bio_19"]]+1)
#plot(preds.climond_stck3[["bio_19"]])

#remove(my_bioclim_rasters)

# Comparativamente, existem muitas diferenças entre Chelsa e Climond. Vou manter o uso do Climond para tudo e ver o que dá na revisão pelo único motivo que os novos predicted land covers são muito otimistas, nem os layers de 2100 parecem representar mais desmatamento do que o que já temos na região.

## Match coordinates x genomic data ####
## I will keep the objects with the same name as in the tutorial
Coordinates <- kentropyx_RAD_coor@data[,c("Sample_ID_full","lat","long")]
colnames(Coordinates) <- c("Population", "Latitude", "Longitude")

# Do never ever forget changing row names order in these kinds of datasets. Coordinates should always match with gendata 
Coordinates = Coordinates[order(match(Coordinates$Population,rownames(AllFreq))),]
Coordinates$Population == rownames(AllFreq)
head(Coordinates)
points(Coordinates$Longitude,Coordinates$Latitude)

# Also changing the order here
kentropyx_RAD_coor = kentropyx_RAD_coor[order(match(kentropyx_RAD_coor$Sample_ID_full,rownames(AllFreq))),]
kentropyx_RAD_coor$Sample_ID_full == rownames(AllFreq)
#.#####

#
# To run a population based analysis #####
# Create a new raster with same res as env to IDing unique localities
r2 = preds.climond_stck3[[1]]
values(r2) <- 1:ncell(r2)
#plot(r2)

# Locality IDs for pops
loc_id = raster::extract(r2,Coordinates[,c("Longitude","Latitude")])
length(loc_id)
unique(loc_id)
loc_id = paste0("pop_", loc_id)
loc_id

# Aggregate pops
AllFreq2 =  aggregate(AllFreq, by = list(loc_id), function(x) mean(x, na.rm = T)/2)
AllFreq2[1:10,1:10]
rownames(AllFreq2) = AllFreq2$Group.1
AllFreq2 =AllFreq2[,-1]
AllFreq2[1:10,1:10]
dim(AllFreq2)

## Filtering on MAF
freq_mean <- colMeans(AllFreq2)
AllFreq2 <- AllFreq2[,-which(freq_mean>=0.95 | freq_mean<=0.05)]
ncol(AllFreq)
ncol(AllFreq2)
AllFreq2[1:10,1:10]

# Then, coordinates should match the new AllFreq
Coordinates$loc_id = loc_id
head(Coordinates)
kentropyx_RAD_coor$loc_id = loc_id # This one is the will keep track of which individual is in which population
head(kentropyx_RAD_coor)

Coordinates = Coordinates[!duplicated(Coordinates$loc_id),]

# Also changing the order here
AllFreq = AllFreq2[order(match(rownames(AllFreq2),Coordinates$loc_id)),]

# Check if population coordinates is matching with 
Coordinates$loc_id == rownames(AllFreq)
AllFreq[1:10,1:10]
#.####
#
## Extracting environmental values from the species points ####
# Climond
Env = raster::extract(preds.climond_stck3, Coordinates[,c("Longitude","Latitude")])

# Plot and check bioclimatic layers variation
plot(preds.climond_stck3[["bio_15"]]);points(Coordinates[,c("Longitude")],Coordinates[,c("Latitude")])

# Check distribution of values
hist(Env[,"bio_10"])

# Check correlation among all variables
cor(Env)

# check for NAs. If the sum is 0, you're good to go!
sum(is.na(Env) == TRUE) 
head(Env)
dim(Env)

## VIF analysis ####
# Diferent from the original RDA tutorial, I will select some variable of interest and run a VIF test

# First I will keep only the max temperature related variables
# Chelsa variables
#Bio3 = Isothermality
#Bio4 = Temperature Seasonality (sd of monthly temps)
#Bio5 = Max Temperature of Warmest Month
#Bio9 = Mean Temperature of Driest Quarter
#Bio10 = Mean Temperature of Warmest Quarter # Não inclui ao final p deixar Bio9
#Bio14 = Precipitation of Driest Month
#Bio15 = Precipitation Seasonality
#Bio17 = Precipitation of Driest Quarter
#Bio18 = Precipitation of Warmest Quarter

# Climond variables 
#Bio01 Annual mean temperature (°C) 
#Bio02 Mean diurnal temperature range (mean(period max-min)) (°C)
#Bio03 Isothermality (Bio02 ÷ Bio07) 
#Bio04 Temperature seasonality (C of V) xxx
#Bio05 Max temperature of warmest week (°C) 
#Bio06 Min temperature of coldest week (°C) 
#Bio07 Temperature annual range (Bio05-Bio06) (°C) ###
#Bio08 Mean temperature of wettest quarter (°C) ###
#Bio09 Mean temperature of driest quarter (°C) ###
#Bio10 Mean temperature of warmest quarter (°C) ###
#Bio11 Mean temperature of coldest quarter (°C)  ###
#Bio12 Annual precipitation (mm) 
#Bio13 Precipitation of wettest week (mm) 
#Bio14 Precipitation of driest week (mm) 
#Bio15 Precipitation seasonality (C of V) ###
#Bio16 Precipitation of wettest quarter (mm) ###
#Bio17 Precipitation of driest quarter (mm) ###
#Bio18 Precipitation of warmest quarter (mm) ###
#Bio19 Precipitation of coldest quarter (mm) ###

##
#Run variation inflaction factor (VIF) for the bioclim variables and prioritise selecting max temperatures

# Vif across raster
general_vif = usdm::vifstep(preds.climond_stck3,th=5)
general_vif

# VIF only for the selected points
general_vif = usdm::vifstep(Env,th=5)
general_vif

# If needed to run variable selection with non correlated variables only
#selected_bio_var = general_vif@results$Variables

# Run a correlogram with all variables
M<-cor(as.data.frame(Env))

# correlogram with hclust reordering
corrplot(M, type="upper", order="hclust", tl.col="black", tl.srt=45)

# Looking at the cor, run vif with pre-selected variables - Climond
usdm::vifstep(as.data.frame(Env[,c("bio_4", "bio_5",  "bio_15", "bio_18")]))

# Finally, selecting only variables that represent extremes 
selected_bio_var = c("bio_4","bio_5","bio_15","bio_18")


## Other non-climatic variables ####
#
# I won't use the following in this version
# ndvi_3 (median from 1999-2012) downloaded from COPERNICUS Global Land Service (Normalized Difference Vegetation Index: Long Term Statistics 1KM: GLOBE 1999-2017 0711)
# The "greenness" of the vegetation.

# Distance to the border of forested biomes (calculated above from the shape of Terrestrial Ecorregions Olson (2001))
# Historical effects of the ecotones and the instability of such areas on the species distribution

# Topographical complexity or rugosity
# More rugosity of topography implies on more microhabitats and possibility for species movement to find the optimal microclimate

# Velocity of climate change
# Measures temperature stability of certain region since the LGM. Generally correlated to Top Complexity on larger geo scales

# Distance to forested biomes border just to delimit the study area right now
#dborder3 = raster('/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/dborder.tif')

##
# Mapping all bioclimatic variation within the datapoits (just to check if there is any variable important not being included)
#masked_preds = crop(mask(preds.climond_stck2,(buffer(kentropyx_RAD_coor,100000))),dborder3)

# For CHELSA
#masked_preds = crop(mask(preds.chelsa.3_var,(buffer(kentropyx_RAD_coor,100000))),dborder3)

##
#names(masked_preds) = c("bio_1 = Annual Mean Temperature",
#"bio_10 = Mean Temp. War Quar",
#"bio_11 = Mean Temp. Cold. Quar",
#"bio_12 = Annual Precipitation",
#"bio_13 = Precip. of Wett Month",
#"bio_14 = Precip. of Dri Month",
#"bio_15 = Precip Seasonality",
#"bio_16 = Precip of Wett Quar",
#"bio_17 = Precip of Dri Quar",
#"bio_18 = Precip of Warm Quar",
#"bio_19 = Precip of Cold Quar",
#"bio_2 = Mean Diurnal Range",
#"bio_3 = Isothermality",
#"bio_4 = Temp. Seasonality",
#"bio_5 = Max Temp. War. Month",
#"bio_6 = Min Temp. Cold. Month",
#"bio_7 = Temp. Annual Range",
#"bio_8 = Mean Temp. Wett. Quar.",
#"bio_9 = Mean Temp. Dri. Quar",
#"NDVI = NDVI")

##
# Plotting precip var
#options(repr.plot.width=17, repr.plot.height=25)
#plot(masked_preds[[5:11]],cex.main = 2)

##
# Plotting Temp var
#options(repr.plot.width=17, repr.plot.height=17)
#plot(masked_preds[[c(1:3,13:19)]],cex.main = 2)

# Plotting Temp var
#options(repr.plot.width=17, repr.plot.height=17)
#plot(masked_preds[[c(20)]],cex.main = 2)

##
#bio_preds_extract_sel = Env[,selected_bio_var]

#bio_preds_extract_sel_h = bio_preds_extract_sel %>% as.data.frame

#ggplot(gather(bio_preds_extract_sel_h), aes(value)) + 
#    geom_histogram(bins = 10) + 
#    facet_wrap(~key, scales = 'free_x')

##
# Distance to forested biomes border (log of the distance, calculated elsewhere)
#options(repr.plot.width=7, repr.plot.height=7)
#dborder3 = raster('/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/dborder_non_log.tif')
#dborder3
#plot(dborder3)

# Maximum EVI
#/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/layers/Maximum_01_05_25km_uint16.tiff
#evi = raster("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/layers/Maximum_01_05_25km_uint16.tiff")
#evi_2 = crop(evi, projectRaster(preds.climond_stck2, crs = oldproj))
#evi_3 = projectRaster(evi_2, crs = crs(preds.climond_stck2))
#plot(log(evi_2))
#remove(evi) # remove these large objects
#remove(evi_2) # remove these large objects
#evi_4 = raster::aggregate(evi_3, 10)
#ndvi = raster::resample(ndvi_4,preds.climond_stck2[[1]])
#plot(ndvi)

##
# Read in and crop NDVI raster
#ndvi = raster("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/kentropyx_analyses/NDVI/NDVI_from_qgis_median.tif")
#ndvi_2 = crop(ndvi, projectRaster(preds.climond_stck2, crs = oldproj))
#ndvi_3 = projectRaster(ndvi_2, crs = crs(preds.climond_stck2))
#plot(ndvi_3)
#remove(ndvi) # remove these large objects
#remove(ndvi_2) # remove these large objects
#ndvi_4 = raster::aggregate(ndvi_3, 10)

#ndvi = raster::resample(ndvi_4,preds.climond_stck2[[1]])
#plot(ndvi)
##
# Also crop and plot selected bioclimatic variables
#options(repr.plot.width=17, repr.plot.height=15)
#bio_preds_sel_crop = crop(subset(preds,c('bio_4','bio_5','bio_6','bio_17','bio_18')),dborder3)
#plot(bio_preds_sel_crop)

##
# Calculating relief roughness
#elevat <- raster("~/Dropbox/1Doutorado/Chapter 1/cerrado_biogeo/rasters_cer/alt_.tif",crs=CRS("+proj=longlat"))
#elevat <- raster::projectRaster(elevat, crs=crs(dborder3))  #spTransform makes the projection
#elevat <- crop(elevat,dborder3)
#plot(elevat)
#elevat
#relief_roughness<- raster::terrain(elevat, opt='roughness', neighbors=8)
#options(repr.plot.width=7, repr.plot.height=7)
#plot(log(relief_roughness+1))

##
# Habitat homogeneity (1/heterogeneity)
# From Tuanmu, M.-N. and W. Jetz. (2015) A global, remote sensing-based characterization of terrestrial habitat heterogeneity for biodiversity and ecosystem modeling. Global Ecology and Biogeography. DOI: 10.1111/geb.12365.
# Homogeneity Similarity of EVI between adjacent pixels
# These metrics are also those which can better capture within-land-cover heterogeneity
# The pixel values of the data layers should be mulitplied by 0.0001 to obtain the actual values of the metrics.
#homogeneity = raster("/Users/josue/Dropbox/4Environmental_layers/Homogeneity_01_05_5km_uint16.tif")
#homogeneity <- raster::projectRaster(homogeneity, crs=crs(dborder3))  #spTransform makes the projection
#homogeneity_c <- crop(homogeneity,dborder3)
#homogeneity_d = homogeneity_c*0.0001
#plot(homogeneity_d)
#remove(homogeneity) # remove these large objects

##
#options(repr.plot.width=18, repr.plot.height=10)
#par(mfrow=c(2,3))
#plot(ndvi_3, main = "NDVI")
#plot(homogeneity_d,main = "Homogeneity")
#plot(relief_roughness,main = "Relief roughness")
#plot(dborder3,main = "Distance to the border")

##
# Lets extract values from all variables again

#bio_preds_extract_sel ## selected bioclim variables
#roughness_extract = raster::extract(relief_roughness, kentropyx_RAD_coor_2)
#ndvi_extract = raster::extract(ndvi_3, kentropyx_RAD_coor_2)
#dborder_extract = raster::extract(dborder3, kentropyx_RAD_coor_2)
#homogeneity_extract =  raster::extract(homogeneity_d, kentropyx_RAD_coor_2)

#all_var_extract = cbind.data.frame(bio_preds_extract_sel,homogeneity_extract)
#all_var_extract = bio_preds_extract_sel

#colnames(all_var_extract) =  gsub("_extract","",colnames(all_var_extract))
#sum(is.na(all_var_extract) == TRUE) # check for NAs. If the sum is 0, you're good to go!
#head(all_var_extract)
#dim(all_var_extract)

##
# Create fake rasters from the extracted values to run VIF
#raster_list = BBmisc::convertColsToList(all_var_extract)
#values_to_raster_fun <- function(raster_list){
#r <- raster(ncol=33, nrow=2)
#values(r) = raster_list
#return(r)
#}
#masked_rasters_df = lapply(raster_list,values_to_raster_fun)
#all_var_extract_fake_ras = raster::stack(masked_rasters_df)
#names(all_var_extract_fake_ras) = colnames(all_var_extract)

##
#usdm::vif(all_var_extract_fake_ras)

##
#all_var_extract = all_var_extract[,-4]
#head(all_var_extract)

##
# Create fake rasters from the extracted values to run VIF
#raster_list = BBmisc::convertColsToList(all_var_extract)
#values_to_raster_fun <- function(raster_list){
#r <- raster(ncol=33, nrow=2)
#values(r) = raster_list
#return(r)
#}
#masked_rasters_df = lapply(raster_list,values_to_raster_fun)
#all_var_extract_fake_ras = raster::stack(masked_rasters_df)
#names(all_var_extract_fake_ras) = colnames(all_var_extract)

##
#usdm::vif(all_var_extract_fake_ras)
#head(all_var_extract,2)

##
#options(repr.plot.width=17, repr.plot.height=17)
#pairs.panels(all_var_extract, scale=TRUE)

##
#all_var_extract2= as.data.frame(all_var_extract)
#all_var_extract2$roughness = log(all_var_extract2$roughness)
#all_var_extract2$bio_17 = log(all_var_extract2$bio_17)
#all_var_extract2$bio_18 = log(all_var_extract2$bio_18)

#bio_preds_extract_sel_h = all_var_extract2 %>% as.data.frame

#ggplot(gather(bio_preds_extract_sel_h), aes(value)) + 
#    geom_histogram(bins = 10) + 
#    facet_wrap(~key, scales = 'free_x')

##
# if you want to use the log values
#all_var_extract = all_var_extract2

##
# Simplifying colnames of variables
#colnames(all_var_extract) = gsub("_extract","", colnames(all_var_extract))
##

## Standardization of the variables ####
#
# List of selected variables
selected_bio_var = c("bio_4","bio_5","bio_15","bio_18")
  
# Extract values from raster
Env = raster::extract(preds.climond_stck3, Coordinates[,c("Longitude","Latitude")])

# To the same scale 
Env <- scale(Env[,selected_bio_var], center=TRUE, scale=TRUE)

head(Env)
colnames(Env) = selected_bio_var
## Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')

## Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- Coordinates$Population
head(Env)
#

## Running a PCA on neutral genetic markers ####
pca <- rda(AllFreq, scale=TRUE) # PCA in vegan uses the rda() call without any predictors
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")
sum_pc = summary(pca)
round(sum_pc$cont$importance[,1:10],2)
## Neutral population structure table
PCs <- scores(pca, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = rownames(AllFreq), PCs)
colnames(PopStruct) <- c("Population", "PC1", "PC2", "PC3")
#

## Table gathering all variables ####
#
Variables <- data.frame(Coordinates, PopStruct[,-1], Env)
head(Variables)

## Administrative boundaries
admin = spTransform(SAshp,CRSobj = behrmannCRS)

## Species range shapefile (download information at beginning of tutorial)
#selected_clipped_biomes_2 = readRDS("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/selected_clipped_biomes_2.rds")
range = selected_clipped_biomes_2

# RDA Variable selection: forward model building procedure ####

## Null model ####
RDA0 <- rda(AllFreq ~ 1 ,  Variables) 

## Full model = Non correlated selected variables  VIF ####
selected_bio_var

# Running with all variables, only Bio18, Bio17 and Bio4 are important
# Then I run with preselected variables and Bio 5 was not important

RDAfull <- rda(formula(paste0("AllFreq","~",paste0(selected_bio_var,collapse=" + "))), Variables) 
RsquareAdj(RDAfull)

## Stepwise procedure with ordiR2step function ####
#mod <- ordistep(RDA0, RDAfull,direction="both",Pin = 0.05, Pout = 0.1, permutations = how(nperm = 199), steps = 50)

## Stepwise procedure with ordiR2step function - Foward 
mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
mod$anova

vif.cca(RDAfull)

## Now considering only pre-selected variables and population structure
## Full model
RDAfull3 <- vegan::rda(formula(paste0("AllFreq","~",paste0(selected_bio_var,collapse=" + ")," + Condition(PC1+PC2+PC3)")), Variables)

mod3 <- ordiR2step(RDA0, RDAfull3, Pin = 0.01, R2permutations = 1000, R2scope = T)

vif.cca(RDAfull3)
##

# RDA Variance partitioning: disentangling the drivers of genetic variation ####
## Full model ####

pRDAfull <- rda(formula(paste0("AllFreq","~",paste0(selected_bio_var,collapse=" + ")," + PC1 + PC2 + Longitude + Latitude ")),  Variables)
RsquareAdj(pRDAfull)
pRDAfull_anova = anova(pRDAfull)

## Pure climate model ####
pRDAclim <- rda(formula(paste0("AllFreq","~",paste0(selected_bio_var,collapse=" + ")," + Condition(Longitude + Latitude + PC1 + PC2)")),  Variables)
RsquareAdj(pRDAclim)
pRDAclim_anova = anova(pRDAclim)

## Pure neutral population structure model ####
pRDAstruct <- rda(formula(paste0("AllFreq","~","PC1 + PC2"," + Condition(Longitude + Latitude +", paste0(selected_bio_var,collapse=" + "),")")),  Variables)
RsquareAdj(pRDAstruct)
pRDAstruct_anova = anova(pRDAstruct)

##Pure geography model ####
pRDAgeog <- rda(formula(paste0("AllFreq","~","Longitude + Latitude"," + Condition(PC1 + PC2 +", paste0(selected_bio_var,collapse=" + "),")")),  Variables)
RsquareAdj(pRDAgeog)
pRDAgeog_anova = anova(pRDAgeog)

# Finally selected = PCs + Variables = Geography doesnt explain too much
# For CHELSA, all 7 selected variables were significant, then I will remove the less important among the temperature var
# For climond
selected_bio_var = c("bio_2","bio_4","bio_5","bio_9","bio_17","bio_18")

# For CHELSA
selected_bio_var = c("bio_4","bio_5","bio_9","bio_14","bio_15","bio_18")

## Correlogram ####

M<-cor(Variables[,c("PC1","PC2","PC3","Latitude", "Longitude",selected_bio_var)])

# correlogram with hclust reordering
corrplot(M, type="upper", order="hclust", tl.col="black", tl.srt=45)

#.####
# Genotype-Environment Associations RDA Part 1: identifying loci under selection ####
RDA_env=NULL
# Run RDA
RDA_env <- rda(formula(paste0("AllFreq","~",paste0(selected_bio_var,collapse=" + ")," + Condition(PC1+PC2+PC3)")),  Variables)
RsquareAdj(RDA_env)
screeplot(RDA_env, main="Eigenvalues of constrained axes")
summary(RDA_env)$concont
# Using >80% of explanation rule

# Run to check axis importance (needs to be in a good computer)
#signif.axis <- anova.cca(wolf.rda, by="axis", parallel=getOption("mc.cores"))

## Function rdadapt
source("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/RDA-landscape-genomics-main/src/rdadapt.R")

## Running the function with K = 3
rdadapt_env<-rdadapt(RDA_env, 3) #ken_rda
sum(rdadapt_env$q.values < 0.05) # FDR of 5%
hist(rdadapt_env$q.values)
hist(rdadapt_env$p.values)

hist(rdadapt_env$p.values, main=NA,family ="Times",xlab="Corrected p-values GIF=1.41") # I got the values below
title("b) RDA", adj =0,family ="Times",cex=1,font=1)


# Below the same rdadapt function to enable modifying parameters
K = 3
zscores<-RDA_env$CCA$v[,1:as.numeric(K)]
resscale <- apply(zscores, 2, scale)
resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
lambda <- median(resmaha)/qchisq(0.5,df=K)
#lambda = 1.2
reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
hist(reschi2test)
#qval <- qvalue(reschi2test)
#q.values_rdadapt<-qval$qvalues
#rdadapt_env = data.frame(p.values=reschi2test, q.values=q.values_rdadapt)
#sum(rdadapt_env$q.values < 0.1)
## Running manually

## P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rdadapt_env$p.values)
#thres_env = 0.0025

# RDA z-scores # Not used anymore
zscores_RDA <-RDA_env$CCA$v[,1:3]

# Running also with no neutral genetic data to check
RsquareAdj(RDAfull)
summary(RDAfull)$concont
rdadapt_env_only <-rdadapt(RDAfull, 3) #ken_rda
sum(rdadapt_env_only$q.values < 0.05)
hist(rdadapt_env_only$q.values)
hist(rdadapt_env_only$p.values)

## Identifying the loci that are below the p-value threshold
outliers_RDA_swiss <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

outliers_RDA_swiss # Final list of SNPs

# Or using a false discovery rate
colnames(AllFreq)[which(rdadapt_env$q.values < 0.05)]

## Top hit outlier per contig
outliers_r <- outliers_RDA_swiss[order(outliers_RDA_swiss$contig, outliers_RDA_swiss$p.value),]

## List of outlier names
outliers_rdadapt_env <- as.character(outliers_r$Loci[!duplicated(outliers_r$contig)])

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers_r$Loci] <- "All outliers"
#TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "Outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
#ggplot() +
  #geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  #geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  #geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  #scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  #geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  #geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
 # xlab("RDA 1") + ylab("RDA 2")+
  #facet_wrap(~"RDA space") +
  #guides(color=guide_legend(title="Locus type")) +
  #theme_minimal(base_size = 11, base_family = "Arial")
  #theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))


ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  #geom_point(data = TAB_sites, aes(x=RDA1*2, y=RDA3*2, colour = cluster_name), size = 4) +
  #scale_color_manual(values = c(c('#e41a1c','#377eb8','#4daf4a','#737373','#ff7f00','#984ea3'))) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.7, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = c("Bio4","Bio5","Bio15","Bio18")), size = 4.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2")+
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Population")) +
  theme_classic(base_size = 12, base_family = "Arial")

##

#.##########################################
# Old way of calculating RDA scores ####
#.##########################################
## ken_rda ####
selected_bio_var
ken_rda = rda(formula(paste0("AllFreq","~",paste0(selected_bio_var,collapse=" + ")," + Condition(PC1+PC2+PC3)")),  Variables)

# Without including neutral genetic info
#ken_rda = RDAfull

##
#ken_rda
class(ken_rda)
RsquareAdj(ken_rda)
summary(ken_rda)$concont

##
#options(repr.plot.width=6, repr.plot.height=5)
screeplot(ken_rda)

##
#test_sig = anova.cca(ken_rda, parallel=getOption("mc.cores"))
#test_sig

##
#signif.axis <- anova.cca(ken_rda, by="axis", parallel=getOption("mc.cores")) ## RUN LATER, takes up to a few hours
#signif.axis

##
#signif.terms <- anova.cca(ken_rda, by="term", parallel=getOption("mc.cores")) ## takes long, go and have a coffee
#signif.terms

##
vif.cca(ken_rda)

## Biplot Individuals/populations ####
# Which samples were classified in which population (SNMF)?
# To check the first individual of each cluster and  against the barplot
#e41a1c - red
#377eb8 - blue
#4daf4a - green
#737373 - grey
#ff7f00 - orange
#984ea3 - purple
# '#e41a1c','#377eb8','#4daf4a','#737373','#ff7f00','#984ea3'

# Here, check if the order of clusters and the colors to be added in the next line is the same as in the barplot
RAD_seq_selector[!duplicated(RAD_seq_selector$class),]

# If colors are in the correct order, add a color name to facilitate recognizing the color scheme from colourbrewer
RAD_seq_selector$color_name = plyr::mapvalues(RAD_seq_selector$class, from=unique(RAD_seq_selector$class), to=c("red","blue","green","gray","orange","purple"))

# Add a name for the clusters
RAD_seq_selector$cluster_name = plyr::mapvalues(RAD_seq_selector$color_name, from=unique(RAD_seq_selector$color_name), to=c("Toc-Ara","Toc","Madeira","North","Xingú","Atl. F."))

# Now assign the colorbrewer pallete to each individual
RAD_seq_selector$map_colour = plyr::mapvalues(RAD_seq_selector$color_name, from=unique(RAD_seq_selector$color_name), to=c('#e41a1c','#377eb8','#4daf4a','#737373','#ff7f00','#984ea3'))

# I will create another object to not mess with column order
ind_pop_assing = RAD_seq_selector

# Match the order of ind_pop_assing individuals with AllFreq (imputed genotypes)
ind_pop_assing = ind_pop_assing[order(match(rownames(ind_pop_assing),rownames(genotipes_t))),]

# Check the order
rownames(ind_pop_assing) == rownames(genotipes_t)
head(ind_pop_assing)
dim(ind_pop_assing)

# If population level analyses
ind_pop_assing$loc_id = loc_id
ind_pop_assing = ind_pop_assing[!duplicated(ind_pop_assing$loc_id),]
ind_pop_assing$loc_id == Variables$loc_id

##
# var_selected2 contains only the selected bioclimatic variables 
var_selected2 = cbind.data.frame(Variables[,selected_bio_var],ind_pop_assing$cluster_name)

# Here I added the colors to the respective populations (pop)
colnames(var_selected2) = c(colnames(Variables[,selected_bio_var]),'pop')
head(var_selected2,10)

##
# This is just for the RDA plotting
levels(var_selected2$pop) <- unique(ind_pop_assing$color_name)
eco <- var_selected2$pop
bg <- unique(ind_pop_assing$map_colour) #  nice colors for our ecotypes
#bg[eco]


##
options(repr.plot.width=7, repr.plot.height=7)
# axes 1 & 2
#plot(RDAfull, type="n", scaling=3,frame=FALSE,xlim=c(-4,5))
#points(RDA_env, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
#points(RDAfull, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, #bg=ind_pop_assing$map_colour) # the wolves
#text(RDAfull, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
#legend("bottomright", legend=unique(var_selected2$pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

##
# axes 1 & 3
#plot(ken_rda, type="n", scaling=3, choices=c(1,3))
#points(ken_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
#points(ken_rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=ind_pop_assing$map_colour, choices=c(1,3))
#text(ken_rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
#legend("bottomright", legend=unique(var_selected2$pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


## Formatting table for ggplot
# For individuals/populations
sites_scores <- scores(RDA_env, choices=c(1:3), display="sites", scaling="none")
TAB_sites <- data.frame(names = row.names(sites_scores), sites_scores)
rownames(TAB_sites) == ind_pop_assing$loc_id
TAB_sites$cluster_name = ind_pop_assing$cluster_name
head(TAB_sites)
TAB_sites$cluster_name <- factor(TAB_sites$cluster_name, levels = c("Toc-Ara","Toc","Madeira","North","Xingú","Atl. F."))
TAB_sites <- TAB_sites[order(TAB_sites$cluster_name),]
# For loci
locus_scores <- scores(RDA_env, choices=c(1:3), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers_r$Loci] <- "All outliers"
#TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "Outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp"))# pull the biplot scores
TAB_var_13 <- as.data.frame(scores(RDA_env, choices=c(1,3), display="bp"))
TAB_var_13_1 = TAB_var_13
TAB_var_13_1[3,2] = -0.13

## Biplot of RDA loci and variables scores RDA 1 e 2
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_sites, aes(x=RDA1*2, y=RDA2*2, colour = cluster_name), size = 4) +
  scale_color_manual(values = c(c('#e41a1c','#377eb8','#4daf4a','#737373','#ff7f00','#984ea3'))) +
  #geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  #scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.7, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = c("Bio4","Bio5","Bio15","Bio18")), size = 4.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2")+
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Population")) +
  theme_classic(base_size = 12, base_family = "Arial")

## Biplot of RDA loci and variables scores RDA 1 e 3
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_sites, aes(x=RDA1*2, y=RDA3*2, colour = cluster_name), size = 4) +
  scale_color_manual(values = c(c('#e41a1c','#377eb8','#4daf4a','#737373','#ff7f00','#984ea3'))) +
  #geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  #scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var_13, aes(xend=RDA1, yend=RDA3, x=0, y=0), colour="black", size=0.7, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_13_1, aes(x=1.1*RDA1, y=1.1*RDA3, label = c("Bio4","Bio5","Bio15","Bio18")), size = 4.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 3")+
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Population")) +
  theme_classic(base_size = 11, base_family = "Arial")


## Formatting table for ggplot Neutral
# For individuals/populations
sites_scores <- scores(RDAfull, choices=c(1:3), display="sites", scaling="none")
TAB_sites <- data.frame(names = row.names(sites_scores), sites_scores)
rownames(TAB_sites) == ind_pop_assing$loc_id
TAB_sites$cluster_name = ind_pop_assing$cluster_name
head(TAB_sites)
TAB_sites$cluster_name <- factor(TAB_sites$cluster_name, levels = c("Toc-Ara","Toc","Madeira","North","Xingú","Atl. F."))
TAB_sites <- TAB_sites[order(TAB_sites$cluster_name),]
# For loci
locus_scores <- scores(RDAfull, choices=c(1:3), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers_r$Loci] <- "All outliers"
#TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "Outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDAfull, choices=c(1,2), display="bp"))# pull the biplot scores
TAB_var_1 = TAB_var
TAB_var_1[2,1] = -0.3
TAB_var_13 <- as.data.frame(scores(RDAfull, choices=c(1,3), display="bp"))
TAB_var_13_1 = TAB_var_13
TAB_var_13_1[3,2] = -0.3

## Biplot of RDA populations and variables scores RDA 1 e 2
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_sites, aes(x=RDA1*2, y=RDA2*2, colour = cluster_name), size = 4) +
  scale_color_manual(values = c(c('#e41a1c','#377eb8','#4daf4a','#737373','#ff7f00','#984ea3'))) +
  #geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  #scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.7, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_1, aes(x=1.1*RDA1, y=1.1*RDA2, label = c("Bio4","Bio5","Bio15","Bio18")), size = 4.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2")+
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Population")) +
  theme_classic(base_size = 12, base_family = "Arial")

## Biplot of RDA populations and variables scores RDA 1 e 3
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_sites, aes(x=RDA1*2, y=RDA3*2, colour = cluster_name), size = 4) +
  scale_color_manual(values = c(c('#e41a1c','#377eb8','#4daf4a','#737373','#ff7f00','#984ea3'))) +
  #geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  #scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var_13, aes(xend=RDA1, yend=RDA3, x=0, y=0), colour="black", size=0.7, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var_13_1, aes(x=1.1*RDA1, y=1.1*RDA3, label = c("Bio4","Bio5","Bio15","Bio18")), size = 4.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 3")+
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Population")) +
  theme_classic(base_size = 12, base_family = "Arial")

##
load.rda <- summary(ken_rda)$species[,1:3]
par(mfrow=c(1, 2))
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")

## Getting outlier loci ####
#
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) # find loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

cand1 <- outliers(load.rda[,1], 3) # 
length(cand1)

##
cand2 <- outliers(load.rda[,2], 3) # 
length(cand2)

##
cand3 <- outliers(load.rda[,3], 3) # 
length(cand3)

#cand4 <- outliers(load.rda[,4], 3) # 
#length(cand4)

##
ken_rda.cand <- c(names(cand1), names(cand2), names(cand3))#,names(cand4)) # just the names of the candidates

##
length(ken_rda.cand[duplicated(ken_rda.cand)]) # 0 duplicate detections (detected on multiple RDA axes)

##
#ken_rda.cand <- ken_rda.cand[!duplicated(ken_rda.cand)] # 137 unique candidates 

##
length(ken_rda.cand)
ncand = length(ken_rda.cand)

##
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
#cand4 <- cbind.data.frame(rep(3,times=length(cand4)), names(cand4), unname(cand4))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
#<- colnames(cand4)

cand <- rbind(cand1, cand2, cand3) #,cand4)
cand$snp <- as.character(cand$snp)
head(cand,3)

##
# Here for some reason I don't remember anymore, I didnt use variables2, instead I created another object with the selected RDA variables
all_var_extract = Variables[,selected_bio_var]

foo <- matrix(nrow=(ncand), ncol= ncol(all_var_extract))  # 8 columns for 8 predictors
colnames(foo) <- colnames(all_var_extract)
head(foo)

##
gen.imp=AllFreq

##
#i=1
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(all_var_extract,2,function(x) cor(x,snp.gen))
}

##
cand <- cbind.data.frame(cand,foo)  
head(cand,3)
dim(cand)
sum(is.na(gen.imp)) # No NAs

##
length(cand$snp[duplicated(cand$snp)])   # 0 
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2])
table(foo[foo[,1]==2,2])
table(foo[foo[,1]==3,2])
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
head(cand,3)

##
#Next, we’ll see to which of the predictors each candidate SNP is most strongly correlated with:
# Check the colnames of "cand" and modify below
cols_cand = ncol(cand)
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,cols_cand+1] <- names(which.max(abs(bar[4:cols_cand]))) # gives the variable
  cand[i,cols_cand+2] <- max(abs(bar[4:cols_cand]))              # gives the correlation
}

##
ncol(cand)
head(cand,5)

##
colnames(cand)[ncol(cand)-1] <- "predictor"
colnames(cand)[ncol(cand)] <- "correlation"
head(cand,3)
table(cand$predictor)

##
sel <- cand$snp
env <- cand$predictor
# '#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02'
env = plyr::mapvalues(env, from=unique(env), to=c('#1b9e77','#7570b3','#e7298a','#e6ab02'))

##
# color by predictor:
col.pred <- rownames(ken_rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("RAD_",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1b9e77','#7570b3','#e7298a','#e6ab02')

##
colnames(all_var_extract)

##
## Plotting candidate SNPs ####
# axes 1 & 2
#dev.off()
plot(ken_rda, type="n", scaling=1, xlim=c(-.2,.2), ylim=c(-.2,.2))
points(ken_rda, display="species", pch=21, cex=.8, col="gray32", bg=col.pred, scaling=3)
points(ken_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(ken_rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=colnames(all_var_extract), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

##
# axes 1 & 3
plot(ken_rda, type="n", scaling=3, xlim=c(-.2,.2), ylim=c(-.2,.2), choices=c(1,3))
points(ken_rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(ken_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(ken_rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=colnames(all_var_extract), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# An additional way of plotting candidate SNPs (new RDA tutorial)
## Biplot of RDA loci and variables scores
## Formatting table for ggplot
locus_scores <- scores(ken_rda, choices=c(1:3), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%cand$snp] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%cand$snp] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(ken_rda, choices=c(1,2), display="bp")) # pull the biplot scores
#
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Arial") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

#################################################
## Running LFMM new way ####
###################################################
# Preparing data (PCA + testing function)
# Carregue os seguintes pacotes

# First, I will run a PCA (function RDA sic) of the environmental variables
pred.pca <- rda(Variables[,selected_bio_var], scale=TRUE)
summary(pred.pca)$cont
screeplot(pred.pca, main = "Screeplot: Eigenvalues of Wolf Predictor Variables")

# correlations between the PC axis and predictors:
round(scores(pred.pca, choices=1:3, display="species", scaling=0), digits=3)

#We’ll store our synthetic PC axis predictor as pred.PC1 for use in LFMM.
pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)
pred.PC2 <- scores(pred.pca, choices=2, display="sites", scaling=0)
pred.PC3 <- scores(pred.pca, choices=3, display="sites", scaling=0)

# Checking geographical distribution of PCs
pca_of_sel_var = RStoolbox::rasterPCA(subset(preds.climond_stck3, selected_bio_var))

plot(pca_of_sel_var$map[[1]], main="PC1")
plot(pca_of_sel_var$map[[2]], main="PC2")
plot(pca_of_sel_var$map[[3]], main="PC3")

# Below is a test version of LFMM2
#LFMM  v.2 computes LFMMs for GEA using a least-squares estimation method (instead of MCMC in v.1). There are two penalty approaches: ridge and lasso. We’ll use ridge today (see ?lfmm_ridge); see Jumentier et al. (2020) for more information on the ridge vs. lasso penalties.

wolf.lfmm <- lfmm_ridge(Y=AllFreq, X=pred.PC1, K=6) # change K as you see fit

#The lfmm package has a nice built-in function to calculate test statistics for the predictor(s), see ?lfmm_test:

wolf.pv <- lfmm_test(Y=AllFreq, X=pred.PC1, lfmm=wolf.lfmm, calibrate="gif")

names(wolf.pv) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values

#Let’s look at the genomic inflation factor (GIF):
wolf.pv$gif

#An appropriately calibrated set of tests will have a GIF of around 1. The elevated GIF for our tests indicates that the results may be overly liberal in identifying candidate SNPs. If the GIF is less than one, the test may be too conservative.

#NOTE: Changing the value of K influences the GIF, so additional tests using the “best” value of K +/- 1 may be needed in some cases.

#Let’s look at how application of the GIF to the p-values impacts the p-value distribution:

hist(wolf.pv$pvalue[,1], main="Unadjusted p-values")        
hist(wolf.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values")

#We want to see a relatively flat histogram (most loci not under selection) with a peak near zero, indicating candidate adaptive markers.

#Note that you can choose a different GIF and readjust the p-values to get a “better” histogram distribution. This process is subjective and can be difficult with empirical data sets, especially those with an IBD signature, such as these wolf data. Remember, you can also change the value of K in your lfmm models and see how this impacts the GIF.

#Below I’ll show you how to manually adjust the GIF correction factor:

# Let's change the GIF and readjust the p-values:
zscore <- wolf.pv$score[,1]   # zscores for first predictor, we only have one in our case...
(gif <- wolf.pv$gif[1])       # default GIF for this predictor

new.gif1 <- 3        # choose your new GIF

# Manual adjustment of the p-values:
adj.pv1 <- pchisq(zscore^2/new.gif1, df=1, lower = FALSE)

#Plot the p-value histograms:
hist(wolf.pv$pvalue[,1], main="Unadjusted p-values")        
hist(wolf.pv$calibrated.pvalue[,1], main="GIF-adjusted p-values (GIF=2.8)")
hist(adj.pv1, main="REadjusted p-values (GIF=2.0)")

#For now, we’ll stick with the default GIF calculated by the lfmm package, though it looks like the application of the GIF may be a bit conservative (e.g. it is compressing the peak at 0).

#Finally, we convert the adjusted p-values to q-values. q-values provide a measure of each SNP’s significance, automatically taking into account the fact that thousands are simultaneously being tested. We can then use an FDR threshold to control the number of false positive detections (given that our p-value distribution is “well-behaved”).

wolf.qv <- qvalue(wolf.pv$calibrated.pvalue)$qvalues
wolf.qv.adj = qvalue(adj.pv1)$qvalues
length(which(wolf.qv < 0.01)) # how many SNPs have an FDR < 1%?
length(which(wolf.qv.adj < 0.05)) # how many SNPs have an FDR < 5%?
(wolf.FDR.1 <- colnames(gen.imp)[which(wolf.qv < 0.05)]) # identify which SNPs these are

## Function to extract candidate SNPs and plotting histogram ####
# I created this function to automatize the search for the best K and extract zscores, pvalues and candidate snps
X = pred.PC2
K = 6
new.gif = NULL

lfmm_fun = function(X,K,FDR,new.gif){
  kazin = c(K-1,K,K+1,K+2)
  kazin_sel = function(kazin){
  wolf.lfmm.PC2 <- lfmm_ridge(Y=AllFreq, X=X, K=kazin) # change K as you see fit
  wolf.pv.PC2 <- lfmm_test(Y=AllFreq, X=X, lfmm=wolf.lfmm.PC2, calibrate="gif")
  return(wolf.pv.PC2)
  }
  wolf.pv.PC1 = pbapply::pblapply(X=kazin,FUN=kazin_sel)
 which_gif =  c(wolf.pv.PC1[[1]]$gif,wolf.pv.PC1[[2]]$gif,wolf.pv.PC1[[3]]$gif,wolf.pv.PC1[[4]]$gif)
smaller_gif =  which.min(which_gif)
wolf.pv.PC2 = wolf.pv.PC1[[smaller_gif]]
zscore.PC2 <- wolf.pv.PC2$score[,1]   # zscores for first predictor
gif.PC2 = ifelse(is.null(new.gif),wolf.pv.PC2$gif[1],new.gif)  
#(gif.PC2 <- wolf.pv.PC2$gif[1])
  adj.pv1.pc2 <- pchisq(zscore.PC2^2/gif.PC2, df=1, lower = FALSE) # ajdusted 
  # plot histograms
  #par(mfrow=c(1,2))
  hist(adj.pv1.pc2, main= paste0("GIF-adjust. p-values GIF=",round(gif.PC2,2), " BestK=", kazin[[smaller_gif]]))
  #hist(adj.pv1.pc2, main=paste0("REadjust. p-values GIF=", paste0(new.gif)))
  # New way getting with BH algorith
  #L = ncol(AllFreq)
  #w = which(sort(adj.pv1.pc2) < q.level*(1:L)/L)
  #length(w)
  #candidates.aov = order(adj.pv1.pc2)[w]
  #FDR.adj.aov = mean(candidates.aov < 901)
  #POW.adj.aov = sum(candidates.aov > 900)/100
  # Old way with qvalues
  wolf.qv.adj = qvalue(adj.pv1.pc2)$qvalues
  length(which(wolf.qv.adj < FDR)) # how many SNPs have an FDR < 10%?
  wolf.FDR.adj.pc2 <- colnames(AllFreq)[which(wolf.qv.adj < FDR)] # 
  return(list(wolf.FDR.adj.pc2,as.data.frame(wolf.qv.adj),as.data.frame(zscore.PC2),as.data.frame(adj.pv1.pc2)))
  #return(list(wolf.FDR.adj.pc2,as.data.frame(wolf.qv.adj),FDR.adj.aov,POW.adj.aov))
  #return(list(w,as.data.frame(adj.pv1.pc2),FDR.adj.aov,POW.adj.aov))
}

## I got basically the same result as in lfmm2() in the LEA package

## Running LFMM across variables ####
# For PC1
cand_spn_lfmm_pc1 = lfmm_fun(X=pred.PC1, K=bestK,new.gif = NULL, FDR = 0.05)
length(cand_spn_lfmm_pc1[[1]]) # lenght of cand SNPs object
head(cand_spn_lfmm_pc1[[2]]) # qvalues
head(cand_spn_lfmm_pc1[[3]]) #  z-scores
head(cand_spn_lfmm_pc1[[4]]) # adj-pvalues
hist(cand_spn_lfmm_pc1[[4]]$adj.pv1.pc2)

# For PC2 ### 
cand_spn_lfmm_pc2 = lfmm_fun(X=pred.PC2, K=bestK, new.gif = NULL, FDR = 0.05)
length(cand_spn_lfmm_pc2[[1]]) # lenght of cand SNPs object
head(cand_spn_lfmm_pc2[[2]]) # qvalues
head(cand_spn_lfmm_pc2[[3]]) #  z-scores
head(cand_spn_lfmm_pc2[[4]]) # adj-pvalues
hist(cand_spn_lfmm_pc2[[4]]$adj.pv1.pc2)

# For PC3
cand_spn_lfmm_pc3 = lfmm_fun(X=pred.PC3, K=bestK,new.gif = NULL, FDR = 0.05)
length(cand_spn_lfmm_pc3[[1]]) # lenght of cand SNPs object
head(cand_spn_lfmm_pc3[[2]]) # qvalues
head(cand_spn_lfmm_pc3[[3]]) #  z-scores
head(cand_spn_lfmm_pc3[[4]]) # adj-pvalues
hist(cand_spn_lfmm_pc3[[4]]$adj.pv1.pc2)

# A way of joining the results of LFMM-PCs and RDA using the qvalues median (not used)
# Joining PCs and RDA
#z2 = cbind(cand_spn_lfmm_pc1[[2]],cand_spn_lfmm_pc2[[2]],cand_spn_lfmm_pc3[[2]])
#z2.median = apply(z2, MARGIN = 1, median)
#z3 = cbind(z2.median,rdadapt_env$q.values)
#z3.median = apply(z3, MARGIN = 1, median)
#hist(cand_spn_lfmm_pc1[[2]]$wolf.qv.adj)
#hist(z3.median, col = "blue")
#lfmm_bios_sel <- colnames(AllFreq)[which(z2.median < 0.1)]
#lfmm_bios_sel

# LFMM for each variable (not used in this version)
#selected_bio_var
#length(selected_bio_var)
#selected_bio_var[[1]]

# LFMM for var 1
# test K = 5:8        # Remember to check the FDR
#cand_spn_lfmm_1 = lfmm_fun(X=Variables[,selected_bio_var[[1]]], K=bestK,new.gif = NULL, FDR = 0.000001)
#length(cand_spn_lfmm_1[[1]])
#head(cand_spn_lfmm_1[[2]])

# LFMM for var 2
#cand_spn_lfmm_2 = lfmm_fun(X=Variables[,selected_bio_var[[2]]], K=bestK, new.gif = NULL,FDR = 0.000001)
#length(cand_spn_lfmm_2[[1]])
#head(cand_spn_lfmm_2[[2]])

# LFMM for var 3 
#cand_spn_lfmm_3 = lfmm_fun(X=Variables[,selected_bio_var[[3]]], K=bestK, new.gif = NULL, FDR = 0.000001)
#length(cand_spn_lfmm_3[[1]])
#head(cand_spn_lfmm_3[[2]])

# LFMM for var 4
#cand_spn_lfmm_4 = lfmm_fun(X=Variables[,selected_bio_var[[4]]], K=bestK, new.gif = NULL, FDR = 0.000001)
#length(cand_spn_lfmm_4[[1]])
#head(cand_spn_lfmm_4[[2]])

# Combining qvalues of multiple tests for LFMM and RDA (not used)
#z2 = cbind(cand_spn_lfmm_1[[3]],cand_spn_lfmm_2[[3]],cand_spn_lfmm_3[[3]],cand_spn_lfmm_4[[3]])#,cand_spn_lfmm_5[[2]],cand_spn_lfmm_6[[2]],rdadapt_env$q.values)
#z2.median = apply(z2, MARGIN = 1, median)
#hist(cand_spn_lfmm_6[[2]]$wolf.qv.adj)
#hist(z2.median, col = "blue")
#lfmm_bios_sel <- colnames(AllFreq)[which(z2.median < 0.1)]
#lfmm_bios_sel

# Intersecting SNPs among methods ####
# Cand SNP PC1 and PC2 tutorials
cand_lfmm = c(cand_spn_lfmm_pc1[[1]],cand_spn_lfmm_pc2[[1]],cand_spn_lfmm_pc3[[1]])
cand_lfmm =cand_lfmm[!duplicated(cand_lfmm)]
cand_lfmm
 
 # RDA Qvalues higher than FDR of 5%
RDA_qvalues <- colnames(AllFreq)[which(rdadapt_env$q.values < 0.05)]
RDA_qvalues

# Union of RDA and LFMM selected with FDR = 5% 
RDA_LFMM_Q_union = c(RDA_qvalues,cand_lfmm)
RDA_LFMM_Q_union = RDA_LFMM_Q_union[!duplicated(RDA_LFMM_Q_union)]
RDA_LFMM_Q_union
 
# Instersect RDAqvalues LFMMqvalues
RDA_LFMM_Q_inter =  intersect(cand_lfmm,RDA_qvalues)
RDA_LFMM_Q_inter

# Union of cand SNPs with top values of each methodology ~10 SNPs for each methodology
top_50_RDA = max(sort(rdadapt_env$q.values)[1:10])
top_50_RDA = colnames(AllFreq)[which(rdadapt_env$q.values <= top_50_RDA)]

top_50_LFMM = max(sort(c(cand_spn_lfmm_pc1[[2]][,1],cand_spn_lfmm_pc2[[2]][,1],cand_spn_lfmm_pc3[[2]][,1]))[1:10])
top_50_LFMM1 = colnames(AllFreq)[which(cand_spn_lfmm_pc1[[2]] <= top_50_LFMM)]
top_50_LFMM2 = colnames(AllFreq)[which(cand_spn_lfmm_pc2[[2]] <= top_50_LFMM)]
top_50_LFMM3 = colnames(AllFreq)[which(cand_spn_lfmm_pc3[[2]] <= top_50_LFMM)]

top_50_LFMM = c(top_50_LFMM1,top_50_LFMM2,top_50_LFMM3)
top_50_LFMM = top_50_LFMM[!duplicated(top_50_LFMM)]

union_RDA_LFMM = c(top_50_LFMM,top_50_RDA)
union_RDA_LFMM

# Cand PC1-2 and RDA-loadings (old tutorial)
cand_snp_lfmm_pc_rda = intersect(cand_lfmm, ken_rda.cand)
cand_snp_lfmm_pc_rda

# Cand individual variables LFMM
#selected_bio_var
#cand_lfmm_bios = c(cand_spn_lfmm_1[[1]],cand_spn_lfmm_2[[1]],cand_spn_lfmm_3[[1]],cand_spn_lfmm_4[[1]])#,cand_spn_lfmm_5[[1]],cand_spn_lfmm_6[[1]])

#cand_lfmm_bios=cand_lfmm_bios[!duplicated(cand_lfmm_bios)]
#cand_lfmm_bios

## LFMM BIOs and RDA old tutorial
#intersect(ken_rda.cand,cand_lfmm_bios)
#cand_snp_lfmm_rda = intersect(cand_lfmm_bios,ken_rda.cand)
#cand_snp_lfmm_rda
#ken_rda.cand = ken_rda.cand[!duplicated(ken_rda.cand)]

## LFMM BIOs and LFMM PCs
#intersect(cand_lfmm_bios,cand_lfmm) # Bom sinal. Mostly intersected
#lfmm_bios_lfmm_pcs = intersect(cand_lfmm_bios,cand_lfmm) 
#lfmm_bios_lfmm_pcs = lfmm_bios_lfmm_pcs[!duplicated(lfmm_bios_lfmm_pcs)]
#lfmm_bios_lfmm_pcs

## LFMM BIOs and RDA new tutorial
#intersect(cand_lfmm_bios,outliers_RDA_swiss$Loci)
#rda_new_lfmm_bios = intersect(cand_lfmm_bios,outliers_RDA_swiss$Loci)
#intersect(cand_lfmm,outliers_RDA_swiss$Loci)

# RDA old tutorial and RDA new tutorial
intersect(outliers_RDA_swiss$Loci,ken_rda.cand)
RDA_RDA_intersect = intersect(outliers_RDA_swiss$Loci,ken_rda.cand)

# Intersecting using DCMS (not used)
# I will combine p-values calculated by LFMM and RDA using DCM
# First create a dataframe with z-statistics of all tests
#statistics_lfmm = cbind.data.frame(z2.median, # Lfmm
#                                   zscores_RDA) # RDA
#head(statistics_lfmm)
# Then a dataframe with the corrected p-values
#pvalues_lfmm = cbind.data.frame(p.values_lfmm_all, # Lfmm
#                                rdadapt_env$p.values) # RDA
#head(pvalues_lfmm)
# (not used) - calculate ranked p-values
#pvalues_lfmm = stat_to_pvalue(dfv = statistics_lfmm,two.tailed=c(FALSE,FALSE,FALSE,FALSE))
#rownames(pvalues_lfmm) = rownames(cand_spn_lfmm_1[[3]])
#res_dcms = Mahalanobis(log10(pvalues_lfmm))
#names(res_dcms) = rownames(cand_spn_lfmm_1[[3]])
#quantile(res_dcms,0.9)
#sum(res_dcms>quantile(res_dcms,0.99))

# Calculate DCMS
#res_dcms = DCMS(dfv = statistics_lfmm,dfp = pvalues_lfmm)
# Bellow I will produce p-values and q-values to filter according to the Benjamini and Hochberg method 
# First fit a rlm model to get values of mean and sd for the pnorm function
#model= MASS::rlm(res_dcms~1)
#dcms_pvalues = pnorm(q=res_dcms, mean=mean(model$residuals), sd=sd(model$residuals), lower.tail=FALSE)
# Check the new p-values distribution
#hist(dcms_pvalues)
# Adjust p-values with BH
#dcms_qvalues = p.adjust(p = dcms_pvalues, method = "BH", n=length(dcms_pvalues))
# List of candidate SNPs
#names(which(dcms_qvalues < 0.00001)) 

#dcms_candidate_selected = names(which(dcms_qvalues < 0.01))[names(which(dcms_qvalues < 0.01)) %in% c(cand_spn_lfmm_1[[1]],colnames(AllFreq)[which(rdadapt_env$q.values < 0.01)])]
#dcms_candidate_selected = names(res_dcms)[res_dcms>quantile(res_dcms,0.999)]

# Using Rdadapt code to calculate one set of zscores and pvalues for lfmm
# Final choice for the manuscript
# Combine zscores
zscores = cbind(cand_spn_lfmm_pc1[[3]],cand_spn_lfmm_pc2[[3]],cand_spn_lfmm_pc3[[3]])
head(zscores)
resscale <- apply(zscores, 2, scale)
# Calculate Mahanolabis distance
resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
# Calculate Genomic Inflaction Factor
lambda <- median(resmaha)/qchisq(0.5,df= ncol(zscores))
#lambda=1.5
# Calculate corrected pvalues
reschi2test <- pchisq(resmaha/lambda,ncol(zscores),lower.tail=FALSE)
hist(reschi2test, main=NA,family ="Times",xlab="Corrected p-values GIF=1.38")
title("a) LFMM", adj =0,family ="Times",cex=1,font=1)

# Calculate Qvalues
qval <- qvalue(reschi2test)
q.values_lfmm<-qval$qvalues
p.values_lfmm_all=reschi2test

# Getting candidate SNPs from FDR < 5% 
SNP.q.values_lfmm = colnames(AllFreq)[which(q.values_lfmm < 0.05)] 
SNP.q.values_lfmm

## As in RDA, identifying the loci that are below the p-value threshold
outliers_RDA_LFMM <- data.frame(Loci = colnames(AllFreq)[which(p.values_lfmm_all<thres_env)], p.value = p.values_lfmm_all[which(p.values_lfmm_all<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(p.values_lfmm_all<thres_env)], split = "_"), function(x) x[1])))

outliers_RDA_swiss # Outliers from RDA p-value corrected
outliers_RDA_LFMM # Outliers from LFMM p-value corrected

# Candidate SNPs intersected from the two methods
intersect(outliers_RDA_swiss$Loci,outliers_RDA_LFMM$Loci)


#.####

# Checking distribution of candidate SNPs across space ####
##

# SNPs in common between LFMM and RDA qvalues
RDA_LFMM_Q_inter = intersect(SNP.q.values_lfmm,RDA_qvalues)
RDA_LFMM_Q_union = c(SNP.q.values_lfmm2,RDA_qvalues)
RDA_LFMM_Q_union = RDA_LFMM_Q_union[!duplicated(RDA_LFMM_Q_union)]
union_RDA_LFMM # Top SNPs by qvalues of each analyses

# SNPs in common LFMM and RDA pvalue-corrected
RDA_LFMM_Q_union_t = c(outliers_RDA_swiss$Loci,outliers_RDA_LFMM$Loci)
RDA_LFMM_Q_union_t = RDA_LFMM_Q_union_t[!duplicated(RDA_LFMM_Q_union_t)]

# Final snp_cand list for comparison
snp_cand = list(outliers_RDA_LFMM$Loci,outliers_RDA_swiss$Loci, intersect(outliers_RDA_swiss$Loci,outliers_RDA_LFMM$Loci),RDA_LFMM_Q_union_t,union_RDA_LFMM)

# Name of the cand_list for plotting
cand_name = c("a) LFMM","b) RDA","c) Intersect", "d) Union","e) Top Union")

#### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# If the analyses is not population based, run the line below
#kentropyx_RAD_coor@data$loc_id = kentropyx_RAD_coor@data$Sample_ID_full
#### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# The function bellow will run the genotipying classification using the candidate SNPs and output a map.

my_SNP_plotter = function(snp_cand,cand_name){

# Just for test
#snp_cand = snp_cand[[1]]
#cand_name = cand_name[[1]]
#old_list_of_snps = c('RAD_12199','RAD_12728','RAD_26653','RAD_35088','RAD_37777','RAD_66046','RAD_67491','RAD_67551','RAD_70231','RAD_71670')

# Create a dataframe with the alele frequency of the CAND-SNPs only
imputed_genotipes_RDA_LFMM = as.data.frame(AllFreq[,snp_cand])
head(imputed_genotipes_RDA_LFMM,3)
colnames(imputed_genotipes_RDA_LFMM) %in% snp_cand
ncol(imputed_genotipes_RDA_LFMM)

# Removing duplicated columns (if any)
imputed_genotipes_RDA_LFMM <- imputed_genotipes_RDA_LFMM[!duplicated(as.list(imputed_genotipes_RDA_LFMM))]
ncol(imputed_genotipes_RDA_LFMM)

rowSums(imputed_genotipes_RDA_LFMM)/ncol(imputed_genotipes_RDA_LFMM)

range(rowSums(imputed_genotipes_RDA_LFMM)/ncol(imputed_genotipes_RDA_LFMM))

##
# Counting the frequency of genotypes per individual (only works for individual based analysis)
levels=c(0,1,2,9) #all unique values in df
gen_freq_bio10 <- sapply(levels,function(x)rowSums(imputed_genotipes_RDA_LFMM==x)) #count occurrences of x in each row
head(gen_freq_bio10)
colnames(gen_freq_bio10) <- c("gen0","gen1","gen2","NA")

# Checking if columns are present between geographic and geno data
kentropyx_RAD_coor@data$loc_id %in% rownames(gen_freq_bio10)

# Then
kentropyx_RAD_coor_bio10 =  kentropyx_RAD_coor[!duplicated(kentropyx_RAD_coor@data$loc_id),]

kentropyx_RAD_coor_bio10 = kentropyx_RAD_coor_bio10[order(match(kentropyx_RAD_coor_bio10$loc_id,rownames(AllFreq))),]

kentropyx_RAD_coor_bio10@data$loc_id == rownames(AllFreq)

# Add the SNP frequency to the geographic dataset
head(kentropyx_RAD_coor_bio10)

# Remove old gen0-2 if any (in case of rerunning)
kentropyx_RAD_coor_bio10 = kentropyx_RAD_coor_bio10[,!names(kentropyx_RAD_coor_bio10) %in% c("0","1","2", "gen0", "gen1", "gen2", "NA")]
head(kentropyx_RAD_coor_bio10)  

# Check if individuals are in the correct order  
kentropyx_RAD_coor_bio10$loc_id == rownames(gen_freq_bio10)

# Bind geographic information to genotype counting
kentropyx_RAD_coor_bio10@data = cbind.data.frame(kentropyx_RAD_coor_bio10@data, gen_freq_bio10)
head(kentropyx_RAD_coor_bio10,3)
dim(kentropyx_RAD_coor_bio10)
                         
# Transforming coordinate system just for plotting
proj4string(kentropyx_RAD_coor_bio10) = behrmannCRS
kentropyx_RAD_coor_bio10_T = spTransform(kentropyx_RAD_coor_bio10,CRSobj=oldproj)
kentropyx_RAD_coor_bio10_T@data$long = kentropyx_RAD_coor_bio10_T@coords[,1]
kentropyx_RAD_coor_bio10_T@data$lat = kentropyx_RAD_coor_bio10_T@coords[,2]
#head(kentropyx_RAD_coor_bio10_T@data,3)
gen_freq_bio10_locs = kentropyx_RAD_coor_bio10_T@data
head(gen_freq_bio10_locs)

# I also added a coordinate column to the genotype table (check if this is needed)
imputed_genotipes_RDA_LFMM$long = kentropyx_RAD_coor_bio10_T@coords[,1]
imputed_genotipes_RDA_LFMM$lat = kentropyx_RAD_coor_bio10_T@coords[,2]
head(imputed_genotipes_RDA_LFMM,3)
colnames(imputed_genotipes_RDA_LFMM)

## Separate individual SNPs for plotting
which_snp = 4
snp_to_plot = imputed_genotipes_RDA_LFMM[,c(which_snp,which(colnames(imputed_genotipes_RDA_LFMM) %in% c("long","lat")))]
snp_to_plot[snp_to_plot == 9] = NA
snp_to_plot=na.exclude(snp_to_plot)
head(snp_to_plot)

# Just to check if certain SNP is exactly the same as other
#imputed_genotipes_RDA_LFMM$RAD_15122 == imputed_genotipes_RDA_LFMM$RAD_52326

## A function to create a list of individual SNPs+Coordinates for plotting
which_snp = 1:(ncol(imputed_genotipes_RDA_LFMM)-2) # just because I have only 10 snps
#which_snp=2
snpslctr = function(which_snp){
snp_to_plot = imputed_genotipes_RDA_LFMM[,c(which_snp,which(colnames(imputed_genotipes_RDA_LFMM) %in% c("long","lat")))]
snp_to_plot[snp_to_plot == 9] = NA
snp_to_plot=na.exclude(snp_to_plot)
snp_to_plot=snp_to_plot[order(snp_to_plot[,1]),]
snp_to_plot$filll= plyr::mapvalues(snp_to_plot[,1], from=c(0,1,2), to=c('#999999','#e41a1c','#377eb8'))
#head(snp_to_plot,3)
return(snp_to_plot)
}
#each_snp_table = lapply(which_snp,snpslctr)
#summary(each_snp_table)

# Preparing PCA raster for the entire area
selected_RDA_var = subset(preds.climond_stck3,selected_bio_var)

# Projecting raster to WGS as well
bio_preds_sel_crop_2 = projectRaster(selected_RDA_var,crs=oldproj)

# Perform raster PC
selected_RDA_masked_PCA= RStoolbox::rasterPCA(selected_RDA_var)

# Para remover ####
#selected_RDA_masked_PCA= preds.climond_stck3[[selected_bio_var]]
#selected_RDA_masked_PCA_t = projectRaster(selected_RDA_masked_PCA$map$PC1,crs=oldproj)
#####

#plot(RStoolbox::rasterPCA(selected_RDA_var))

# Plotting PC axes
plot(selected_RDA_masked_PCA$map[[1]])
plot(selected_RDA_masked_PCA$map[[2]])

# Projecting raster to WGS (just PC1 in this case)
selected_RDA_masked_PCA_t = projectRaster(selected_RDA_masked_PCA$map[[1]],crs=oldproj)

plot(selected_RDA_masked_PCA_t)

# Here, I would manually adjust raster values, but it wasn't efficient
#df_to_plot =data.frame(rasterToPoints(selected_RDA_masked_PCA_t))
#colnames(df_to_plot)=(c("x","y","layer"))
#summary(df_to_plot)
#midvalue= 0
#df_to_plot$layer = scales::rescale(df_to_plot$layer)
#df_to_plot$layer = log(df_to_plot$layer+1 )
#summary(df_to_plot$layer)
#midvalue=mean(na.exclude(df_to_plot$layer))

## Plotting each SNP across PC1 space ####
options(repr.plot.width=7, repr.plot.height=7)

# Original plotting of SNPs #
#ploting_each_snp = function(each_snp_table){
#my_title <- colnames(each_snp_table)[1]
#b=ggplot()+
#  geom_polygon(data=SAshp, aes(long, lat, group=group), fill="white",col=NA)+
  #geom_polygon(data=cerrado, aes(long, lat, group=group), fill="gold2")+
  #geom_polygon(data=Caatinga, aes(long, lat, group=group), fill="goldenrod3")+
  #geom_polygon(data=Guiana, aes(long, lat, group=group), fill="gold2")+
  #geom_polygon(data=Llanos, aes(long, lat, group=group), fill="gold2")+
  #geom_path(data=rivers, aes(long, lat, group=group), col="dodgerblue3")+
  #geom_raster(data =df_to_plot,aes(x=x,y=y,fill= cut(layer, breaks=seq(0, 0.6931, length.out=10))),include.lowest = T)+ 
  #optional pallete
  #scale_fill_viridis(option="viridis",direction=1,name = element_text("PCA values"))+
  #Modified pallete
 # scale_fill_gradient2(low="#08135a" , mid = "#41B6C4", high = "yellow",
#                       name = element_text("PCA values"),midpoint = midvalue,
#                       guide = guide_colourbar(direction = "vertical",
#                                               label.hjust = 0.5,label.vjust = 0.5, raster=FALSE))+  
#  scale_fill_viridis_d(alpha = 0.8,option = "A", limits = c(0, .3),direction = -1)+
                       #, oob = scales::squish) +
#new_scale("fill")+
#geom_raster(data = TAB_PCA, aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
#  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "D", labels = c("Dry Seasonal","","","","Intermediate","","","","Wet aseasonal")) +
#geom_polygon(data=SAshp, aes(long, lat, group=group), fill=NA,col="grey20",size = .2)+
#  geom_point(data=each_snp_table,aes(x=long ,y=lat), col="black", fill=each_snp_table$filll, size=6, shape=21)+
  #geom_scatterpie(aes(x=long, y=lat, group = Sample_ID_full, r =.8),col = NA, size = .2,
  #                data = gen_freq_bio10_locs, cols = colnames(gen_freq_bio10_locs[,c(6:8)])) +
  #scale_fill_manual(values=c('#999999','#e41a1c','#377eb8','#4daf4a'))+
  #geom_text(data=kentropyx_RAD_coor@data,aes(label = Sample_ID_red, x = long, y = lat),size=3)+
#  coord_cartesian(xlim=c(-67,-34),ylim = c(-17.5,7))+
#  ggsn::scalebar(x.min=-67, x.max=-34, y.min=-17.5, y.max=7,dist = 250,dist_unit = "km", model = 'WGS84',st.size = 2.5,location="topright",transform=TRUE, border.size=.5)+
  #annotation_custom(g, xmin=-50, xmax=-40, ymin=-1.5, ymax=5.57)+
#  theme_classic()+
#  theme(panel.grid.major = element_line(colour = 'grey75',size = .1),
#        panel.background = element_rect(fill = 'white'))+
#  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
#  ggtitle(my_title)+
#  theme(plot.title = element_text(hjust = 0.5))+
#  theme(panel.border = element_rect(colour = 'black', fill=NA, size=.5))
#return(b)
#}
#f_plots_of_each_snp = lapply(each_snp_table,ploting_each_snp)
#f_plots_of_each_snp[[2]]
##

## Improved PCA variation to maximize contrast ####
## Vectorization of the climatic rasters for ggplot
PCA_proj <- list(selected_RDA_masked_PCA_t)
#PCA_proj <- list(selected_RDA_masked_PCA_t$PC)
PCA_proj <- lapply(PCA_proj, function(x) rasterToPoints(x))
for(i in 1:length(PCA_proj)){
  PCA_proj[[i]][,3] <- (PCA_proj[[i]][,3]-min(PCA_proj[[i]][,3]))/(max(PCA_proj[[i]][,3])-min(PCA_proj[[i]][,3]))
}

TAB_PCA <- as.data.frame(PCA_proj[1])
colnames(TAB_PCA)[3] <- "value"
# If we were to map 2 PC axis (see RDA knife tutorial)
#TAB_PCA$variable <- factor(c(rep("PC1", nrow(PCA_proj[[1]])), rep("PC2", nrow(PCA_proj[[2]]))), levels = c("PC1","PC2"))

TAB_PCA$variable <- factor(rep("PC1", nrow(TAB_PCA)), levels = c("PC1"))

## Plotting each SNP function ####
ploting_each_snp = function(each_snp_table){
my_title <- colnames(each_snp_table)[1]
b=ggplot()+ 
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill= "white",col=NA)+
  geom_raster(data = TAB_PCA, aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "D", labels = c("Dry Seasonal","","","","Intermediate","","","","Wet aseasonal")) +
  #geom_sf(data = admin, fill=NA, size=0.1) +
  #new_scale("fill")+
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill=NA,col="grey20",size = .2)+
  geom_point(data=each_snp_table,aes(x=long ,y=lat), col="black", fill=each_snp_table$filll, size=6, shape=21)+
  coord_cartesian(xlim=c(-67,-34),ylim = c(-17.5,7))+
  ggsn::scalebar(x.min=-67, x.max=-34, y.min=-17.5, y.max=7,dist = 250,dist_unit = "km", model = 'WGS84',st.size = 2.5,location="topright",transform=TRUE, border.size=.5)+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'grey75',size = .1),
        panel.background = element_rect(fill = 'white'))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  guides(fill=guide_legend(title=my_title)) +
  #ggtitle(my_title)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(colour = 'black', fill=NA, size=.5))
return(b)
}
#f_plots_of_each_snp = lapply(each_snp_table,ploting_each_snp)
#f_plots_of_each_snp[[3]] # Plottig one of the SNPs for checking

## Mapping all SNPs in PDF

#pdf(file=paste0(gsub(":", "-", Sys.time()),"RAD_genotype.pdf"),onefile = TRUE)
#options(repr.plot.width=7, repr.plot.height=10)
#ggarrange(plotlist=f_plots_of_each_snp,
#          ncol = 1, nrow = 2)
#dev.off()

## Now plotting the total Allele Frequency of candidate SNPs ####
# I will group by longitude first
gen_freq_bio10_locs=gen_freq_bio10_locs[order(gen_freq_bio10_locs[,"long_wgs"]),]

options(repr.plot.width=7, repr.plot.height=7)
my_title <- paste("Allele frequency PCA")
ggplot()+
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill="white",col=NA)+
  #geom_polygon(data=cerrado, aes(long, lat, group=group), fill="gold2")+
  #geom_polygon(data=Caatinga, aes(long, lat, group=group), fill="goldenrod3")+
  #geom_polygon(data=Guiana, aes(long, lat, group=group), fill="gold2")+
  #geom_polygon(data=Llanos, aes(long, lat, group=group), fill="gold2")+
  #geom_path(data=rivers, aes(long, lat, group=group), col="dodgerblue3")+
  #geom_raster(data =df_to_plot,aes(x=x,y=y,fill=layer))+ 
  #optional pallete
  # scale_fill_viridis(option="plasma",direction=1)+
  #Modified pallete
  geom_raster(data = TAB_PCA, aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8,direction = 1, option = "D", labels = c("Dry Seasonal","","","","Intermediate","","","","Wet aseasonal"),name ="Climate") +
geom_polygon(data=SAshp, aes(long, lat, group=group), fill=NA,col="grey20",size = .2)+
new_scale("fill")+
  #geom_point(data=kentropyx_RAD_coor@data,aes(x=long ,y=lat), col="black", fill="white", size=4, shape=21,height= 0.4,width=0.4)+
  geom_scatterpie(aes(x=long, y=lat, group = Sample_ID_full, r =.8),col = 'grey75', size = .2,data = gen_freq_bio10_locs, cols = colnames(gen_freq_bio10_locs[,c("gen0", "gen1", "gen2")])) +
  scale_fill_manual(values=c('#999999','#e41a1c','#377eb8'),name = element_text("Genotypes"))+
  geom_text(data=kentropyx_RAD_coor@data,aes(label = Sample_ID_red, x = long, y = lat),size=3)+
  coord_cartesian(xlim=c(-67,-34),ylim = c(-17.5,7))+
  ggsn::scalebar(x.min=-67, x.max=-34, y.min=-17.5, y.max=7,dist = 250,dist_unit = "km", model = 'WGS84',st.size = 2.5,location="topright",transform=TRUE, border.size=.5)+
  #annotation_custom(g, xmin=-50, xmax=-40, ymin=-1.5, ymax=5.57)+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'grey75',size = .1),
        panel.background = element_rect(fill = 'white'))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  ggtitle(my_title)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(colour = 'black', fill=NA, size=.5))

#save.image(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_RDA/RDA_Ken_calcarata.RData") #
#.####

##
# Using RDA to cluster individuals per genotype and climate ####

# Extracting values of PC1 (WGS proj.) for each individual
varyables = as.data.frame(raster::extract(bio_preds_sel_crop_2,imputed_genotipes_RDA_LFMM[,c("long","lat")]))
head(varyables)

# Removing coordinates from genotype table
imputed_genotipes_RDA_LFMM_2 = imputed_genotipes_RDA_LFMM[,1:(ncol(imputed_genotipes_RDA_LFMM)-2)]
head(imputed_genotipes_RDA_LFMM_2)

# Checking the dataframe order
colnames(imputed_genotipes_RDA_LFMM_2) %in% snp_cand
rownames(imputed_genotipes_RDA_LFMM_2) == rownames(AllFreq)

##
## Extracting environmental values for each source population
Variables = raster::extract(preds.climond_stck3, Coordinates[,c("Longitude","Latitude")])
Env <- Variables[,selected_bio_var]
## Standardization of the variables
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()

## Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')

## Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- row.names(imputed_genotipes_RDA_LFMM_2)
colnames(Env) = selected_bio_var
x = rowSums(imputed_genotipes_RDA_LFMM_2)

##
# Run RDA with candidate SNPs only
selected_bio_var # Check selected variables for formula
#new_rda = rda(imputed_genotipes_RDA_LFMM_2~ bio_3  + bio_4 + bio_5 + bio_9  + bio_17 +  bio_18, data=Env)

# Here I only need to change the selected_bio_var in the beggining
RDA_outliers = rda(formula(paste0("imputed_genotipes_RDA_LFMM_2","~",paste0(selected_bio_var,collapse=" + "))),data=Env)

# Check the adaptive RDA
RDA_outliers
RsquareAdj(RDA_outliers)
R_squareAdj = round(RsquareAdj(RDA_outliers)$adj.r.squared,2)

## 
# Check model significance
test_sig_new = anova.cca(RDA_outliers, parallel=getOption("mc.cores"))
test_sig_new

## Check RDA significance
#signif.axis_new <- anova.cca(RDA_outliers, by="axis", parallel=getOption("mc.cores")) ## Takes up to a few hours with a full dataset
#signif.axis_new
#signif.axis_new$Variance

##
#dev.off()
screeplot(RDA_outliers, main="Eigenvalues of constrained axes")

## RDA biplot ####
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))

multi_fact= max(c(max(abs(TAB_var$RDA1))/max(abs(TAB_loci$RDA1)), max(abs(TAB_loci$RDA1))/max(abs(TAB_var$RDA1))))

cca1_varex<-round(summary(RDA_outliers)$cont$importance[2,1]*100,2) #Get percentage of variance explained by first axis
cca2_varex<-round(summary(RDA_outliers)$cont$importance[2,2]*100,2) #Get percentage of variance explained by second axis

RDA1_varex = round(summary(RDA_outliers)$concont[[1]][[2,1]]*100,1)
RDA2_varex = round(summary(RDA_outliers)$concont[[1]][[2,2]]*100,1)

ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*multi_fact, y=RDA2*multi_fact), colour = "#EB8055FF", size = 4, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.7, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.11*RDA1, y=1.11*RDA2, label = c("Bio4","Bio5","Bio15","Bio18")), size = 4.5, family = "Times") +
 annotate("text", -Inf, -Inf, hjust = -0.1, vjust = -1.5, label = paste0("Adj. R-squared = ",R_squareAdj)) +
  xlab(paste0("RDA 1 ","(",RDA1_varex,"%)")) + ylab(paste0("RDA 2 ","(",RDA2_varex,"%)")) +
  ggtitle(cand_name)+
  #facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  #theme_bw(base_size = 11, base_family = "Arial") +
  #theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
  theme_classic(base_size = 14, base_family = "Arial")

##
## Function to predict the adaptive index across the landscape ####
## Using it here to map across space
source("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/RDA-landscape-genomics-main/src/adaptive_index.R")

source("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/RDA-landscape-genomics-main/src/genomic_offset.R")


# Run adaptive genome analyses
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 3, env_pres = subset(preds.climond_stck3,selected_bio_var), range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

plot(res_RDA_proj_current[[1]]) # RDA
plot(selected_RDA_masked_PCA_t) # PCA

# Projecting RDA raster if needed (first Axis only in this case)
res_RDA_1_proj_current_t = projectRaster(res_RDA_proj_current$RDA1,crs=oldproj)

##
## Getting RDA loadings and running Kmeans ####
## Here I will use RDA loadings to cluster individuals with similar genotypes using KMEANS

# Using RDA scores
raw_scores = summary(RDA_outliers)$sites[,1:3] 
x = summary(RDA_outliers)$sites[,1:3] 

# If using predicted values
#x = raster::extract(stack(res_RDA_proj_current),Coordinates[,c("Longitude","Latitude")])
#rownames(x) = Coordinates$loc_id


# weighting RDA axis by axis importance (3 first axis here) 
summary(RDA_outliers)$concont # Axis importance

x = as.matrix(cbind.data.frame(RDA1 = x[,1]*summary(RDA_outliers)$concont[[1]][2,1],
                     RDA2 = x[,2]*summary(RDA_outliers)$concont[[1]][2,2],
              RDA3 = x[,3]*summary(RDA_outliers)$concont[[1]][2,3]))
#              RDA4 = x[,4]*summary(RDA_outliers)$concont[[1]][2,4]))

head(x)

# Or if using PCA of SNPs as scores
#gen_str_cand_snp = rda(imputed_genotipes_RDA_LFMM_2)
#x = scores(gen_str_cand_snp, choices=1:3, display="sites", scaling=0)
#head(x)

# Using NbClust with kmeans and selecting number of clusters with all available indexes (majority rule)
res<-NbClust(x, distance = "euclidean", min.nc=2, max.nc=6, 
             method = c("kmeans"), index = "all")

# Using kmeans classification from nbcluster
kmeans_res = as.data.frame(res$Best.partition)
k_kmeans = max(res$Best.partition)

## RDA biplot ####

sites_scores <- as.data.frame(x)
TAB_sites <- data.frame(names = row.names(sites_scores), sites_scores)
rownames(TAB_sites) == ind_pop_assing$loc_id
TAB_sites$cluster_name = res$Best.partition
head(TAB_sites)
TAB_sites$cluster_name <- factor(TAB_sites$cluster_name, levels = sort(unique(TAB_sites$cluster_name),decreasing = F))
TAB_sites <- TAB_sites[order(TAB_sites$cluster_name),]
head(TAB_sites)

TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))

multi_fact= max(c(max(abs(TAB_var$RDA1))/max(abs(TAB_sites$RDA1)), max(abs(TAB_sites$RDA1))/max(abs(TAB_var$RDA1))))

RDA1_varex = round(summary(RDA_outliers)$concont[[1]][[2,1]]*100,1)
RDA2_varex = round(summary(RDA_outliers)$concont[[1]][[2,2]]*100,1)


c("#4daf4a","#377eb8","#e41a1c")
# volta p mim
ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_sites, aes(x=RDA1, y=RDA2, colour = cluster_name), size = 4) +
  scale_color_manual(values = c('#4daf4a','#e41a1c','#377eb8','#737373','#ff7f00','#984ea3')) +
  geom_segment(data = TAB_var, aes(xend=RDA1*multi_fact, yend=RDA2*multi_fact, x=0, y=0), colour="black", size=0.7, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.11*RDA1, y=1.11*RDA2, label = c("Bio4","Bio5","Bio15","Bio18")), size = 4.5, family = "Times") +
  annotate("text", -Inf, -Inf, hjust = -0.1, vjust = -1.5, label = paste0("Adj. R-squared = ",R_squareAdj)) +
  xlab(paste0("RDA 1 ","(",RDA1_varex,"%)")) + ylab(paste0("RDA 2 ","(",RDA2_varex,"%)")) +
  ggtitle(cand_name)+
  #facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  #theme_bw(base_size = 11, base_family = "Arial") +
  #theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
  theme_classic(base_size = 14, base_family = "Arial")

# Fuzzy cluster (not used)
#library(cluster)
#fannyy <- fanny(x, k=2, memb.exp = 1.5)
#round(fannyy$membership, 2)[1:4,]
#fannyy$clustering 
## Returns multiple cluster memberships for coefficient above a certain 
## value (here >0.1)
#fannyyMA <- round(fannyy$membership, 2) >= 0.3
#apply(fannyyMA, 1, function(x) paste(which(x), collapse="_"))
#table(apply(fannyyMA, 1, function(x) paste(which(x), collapse="_"))
#)
#cl2 = apply(fannyyMA, 1, function(x) paste(which(x), collapse="_"))
#cl2[cl2 == "1_2"] = "3"

# Kmeans classification from fuzzy cluster
#kmeans_res = as.data.frame(cl2)
#kmeans_res$kmeans = as.numeric(cl2)

# kmeans with preselect K value
#set.seed(9999)
#dev.off()
cl <- kmeans(x, k_kmeans)
#plot(x, col = cl$cluster)
#points(cl$centers, col = 1:4, pch = 8, cex = 2)
#abline(h=0)
#abline(v=0)

# Separating in the RDA space

# Using kmeans preselected K value
#kmeans_res = as.data.frame(cl$cluster)

# Add a col-name
colnames(kmeans_res) = "kmeans"

# Adding raw RDA values for selecting a range of predicted values for literature records selection
# Working only for 2 clusters for now
RDA_range_group = cbind.data.frame(raw_scores, kmeans_res)


# Using yumi pts Literature records
yumi = TRUE

literature_records = read.csv("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/yumi_pts.csv", header = TRUE,sep = ";")
literature_records = literature_records[,c("Longitude","Latitude")]
literature_records$Longitude = gsub(",",".",literature_records$Longitude)
literature_records$Latitude = gsub(",",".",literature_records$Latitude)
literature_records$Longitude = as.numeric(literature_records$Longitude)
literature_records$Latitude = as.numeric(literature_records$Latitude)
head(literature_records)
literature_records2 = literature_records
coordinates(literature_records2) = literature_records[,c("Longitude","Latitude")]
proj4string(literature_records2) = oldproj
literature_records2 = spTransform(literature_records2,CRSobj = behrmannCRS)
literature_records2 = as.data.frame(literature_records2@coords)
# Using predicted values
y = as.data.frame(raster::extract(stack(res_RDA_proj_current),literature_records2))
which_na = rowSums(is.na(y)) > 0
y = y[!which_na,]
y$lat = literature_records2$Latitude[!which_na]
y$long = literature_records2$Longitude[!which_na]
head(y)

range(na.omit(y[,1]))
range(raw_scores[,1])
boxplot(y[,1],raw_scores[,1])

# Limiting literature records to be within the range of predicted RDA values from actually sampled data

# For group class = 1
nrow(y)
head(y)
  
RDA_range_group_1 = RDA_range_group[RDA_range_group$kmeans==1,]
RDA_range_group_2 = RDA_range_group[RDA_range_group$kmeans==2,]
y_group_1 = y
y_group_2 = y

# For group one
y_group_1[,1][!between(y_group_1[,1],min(RDA_range_group_1[,1]),max(RDA_range_group_1[,1]))]=NA
#y_group_1[,2][!between(y_group_1[,2],min(RDA_range_group_1[,2]),max(RDA_range_group_1[,2]))]=NA
#y_group_1[,3][!between(y_group_1[,3],min(RDA_range_group_1[,3]),max(RDA_range_group_1[,3]))]=NA

y_group_1

which_na_1 = rowSums(is.na(y_group_1)) > 0
y_group_1 = y_group_1[!which_na_1,]
nrow(y_group_1)
y_group_1$kmeans = 1

# For group 2
y_group_2[,1][!between(y_group_2[,1],min(RDA_range_group_2[,1]),max(RDA_range_group_2[,1]))]=NA
#y_group_2[,2][!between(y_group_2[,2],min(RDA_range_group_2[,2]),max(RDA_range_group_2[,2]))]=NA
#y_group_2[,3][!between(y_group_2[,3],min(RDA_range_group_2[,3]),max(RDA_range_group_2[,3]))]=NA

y_group_2

which_na_2 = rowSums(is.na(y_group_2)) > 0
y_group_2 = y_group_2[!which_na_2,]

nrow(y_group_2)
y_group_2$kmeans = 2

# Check if there are two classified in the same group
if (sum(rownames(y_group_2) %in% rownames(y_group_1)) > 0) {
  print("clusters have shared coordinates, check it")
} else {
  final_y = rbind.data.frame(y_group_2,y_group_1)
  print("clusters are well defined")
}


head(final_y)

# Creating dataframes for plotting
#if (yumi == TRUE) {
  classified_ind_4_yumi = final_y
  print("yumi")
#} else {
  print("not yumi")
  classified_ind_4 = kmeans_res  
  head(classified_ind_4)
  classified_ind_4$Sample_ID_full = rownames(kmeans_res) #
  # Matching Kmeans classification with geographical coordinates
  classified_ind_4$lat = kentropyx_RAD_coor@data$lat[match(classified_ind_4$Sample_ID_full,kentropyx_RAD_coor@data$loc_id)]
  classified_ind_4$long = kentropyx_RAD_coor@data$long[match(classified_ind_4$Sample_ID_full,kentropyx_RAD_coor@data$loc_id)]
#}


head(classified_ind_4)

# For the two yumi pts
classified_ind_4_yumi$kmeans_class_col = classified_ind_4_yumi$kmeans
classified_ind_4_yumi$kmeans_class_col[classified_ind_4_yumi$kmeans_class == 1] <- "#e41a1c"
classified_ind_4_yumi$kmeans_class_col[classified_ind_4_yumi$kmeans_class == 2] <- "#377eb8"

head(classified_ind_4_yumi)

## Using different colors for clusters (here only three, but I have tested up to six)

# 1 Verde = populazinha - verde



classified_ind_4$kmeans_class_col = classified_ind_4$kmeans
classified_ind_4$kmeans_class_col[classified_ind_4$kmeans_class == 1] <- "#4daf4a"
classified_ind_4$kmeans_class_col[classified_ind_4$kmeans_class == 2] <- "#e41a1c"
classified_ind_4$kmeans_class_col[classified_ind_4$kmeans_class == 3] <- "#377eb8"
classified_ind_4$kmeans_class_col[classified_ind_4$kmeans_class == 4] <- "#984ea3"
classified_ind_4$kmeans_class_col[classified_ind_4$kmeans_class == 5] <- "#ff7f00"
classified_ind_4$kmeans_class_col[classified_ind_4$kmeans_class == 6] <- "#ffff33"
head(classified_ind_4)

##
## Jittering individual records for allowing visualization of all
classified_ind_4$long2_j <- classified_ind_4$long # desloca os pontos ligeiramente para o lado
classified_ind_4$lat2_j <- classified_ind_4$lat

## No jittering for yumi
classified_ind_4_yumi$long2_j <- jitter(classified_ind_4_yumi$long,amount = 0) # desloca os pontos ligeiramente para o lado
classified_ind_4_yumi$lat2_j <- jitter(classified_ind_4_yumi$lat,amount = 0)


## Changing species point projection
classified_ind_5 = classified_ind_4
head(classified_ind_5)
coordinates(classified_ind_5) = classified_ind_5[,c("long2_j","lat2_j")]
proj4string(classified_ind_5) = behrmannCRS
classified_ind_5 = spTransform(classified_ind_5,CRSobj = oldproj)
classified_ind_5$long3_j = classified_ind_5@coords[,1]
classified_ind_5$lat3_j = classified_ind_5@coords[,2]
head(classified_ind_5)
classified_ind_6 = classified_ind_5@data
head(classified_ind_5)

# For yumi
classified_ind_5_yumi = classified_ind_4_yumi
head(classified_ind_5_yumi)
coordinates(classified_ind_5_yumi) = classified_ind_5_yumi[,c("long2_j","lat2_j")]
proj4string(classified_ind_5_yumi) = behrmannCRS
classified_ind_5_yumi = spTransform(classified_ind_5_yumi,CRSobj = oldproj)
classified_ind_5_yumi$long3_j = classified_ind_5_yumi@coords[,1]
classified_ind_5_yumi$lat3_j = classified_ind_5_yumi@coords[,2]
head(classified_ind_5_yumi)
classified_ind_6_yumi = classified_ind_5_yumi@data
head(classified_ind_6_yumi)

## 
# Checking the number of individuals per kmeans cluster
#table(cl$cluster)
#table(classified_ind_6$kmeans_class_col)

##
# Add a PCA values per individual column
classified_ind_6_pts = classified_ind_6
coordinates(classified_ind_6_pts) = classified_ind_6_pts[,c("long","lat")]
PCA_value_by_point = raster::extract(projectRaster(selected_RDA_masked_PCA_t,crs = behrmannCRS),classified_ind_6_pts,df=TRUE)

PCA_value_by_point$ID = classified_ind_6$Sample_ID_full
#head(PCA_value_by_point)
classified_ind_6$PCA_value = PCA_value_by_point$PC1[match(classified_ind_6_pts$Sample_ID_full,PCA_value_by_point$ID)]
head(classified_ind_6)

# Save classified ind molecular data
write.csv(classified_ind_6_pts,"/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/classified_ind_molecular.csv")

# Add RDA values column (directly from RDA loadings and not from RDA raster)
RDA_value_by_point = cbind.data.frame(names(summary(RDA_outliers)$sites[,1]),summary(RDA_outliers)$sites[,1])
colnames(RDA_value_by_point) = c("ID","RDA_pres_1")
head(RDA_value_by_point)
#head(RDA_value_by_point)
classified_ind_6$RDAA_value = RDA_value_by_point$RDA_pres_1[match(classified_ind_6$Sample_ID_full,RDA_value_by_point$ID)]
head(classified_ind_6)

# Add RDA2 values column (directly from RDA loadings and not from RDA raster)
RDA2_value_by_point = cbind.data.frame(names(summary(RDA_outliers)$sites[,2]),summary(RDA_outliers)$sites[,2])
colnames(RDA2_value_by_point) = c("ID","RDA_pres_1")
head(RDA2_value_by_point)
#head(RDA_value_by_point)
classified_ind_6$RDAA2_value = RDA2_value_by_point$RDA_pres_1[match(classified_ind_6$Sample_ID_full,RDA2_value_by_point$ID)]
head(classified_ind_6)

##
# Add violin plot column
cl1 = classified_ind_6[classified_ind_6$kmeans_class_col == "#e41a1c",] # 1
cl2 = classified_ind_6[classified_ind_6$kmeans_class_col == "#377eb8",] # 2
cl3 = classified_ind_6[classified_ind_6$kmeans_class_col == "#4daf4a",] # 3
cl4 = classified_ind_6[classified_ind_6$kmeans_class_col == "#984ea3",] # 4
cl5 = classified_ind_6[classified_ind_6$kmeans_class_col == "#ff7f00",] # 5
cl6 = classified_ind_6[classified_ind_6$kmeans_class_col == "#ffff33",] # 6

# Define number of clusters selected with kmeans
k_kmeans = max(classified_ind_6$kmeans)

## Violin plots ####
#(# lines for different numbers of kmean clusters)

# I created a if else condition for plotting according to the number of clusters
if (yumi==FALSE) {
  print("not yumi")
if(k_kmeans == 2){
  # For two clusters
  vioplot(cl1$PCA_value,cl2$PCA_value, col=c('#e41a1c','#377eb8'), main = "Env. distribution PCA",axes=F,names=c("Cl1","Cl2"))
  box(bty="l")
  
  vioplot(cl1$RDAA_value,cl2$RDAA_value, col=c('#e41a1c','#377eb8'), main = "Env. distribution of RDA",axes=F,names=c("Cl1","Cl2"))
  box(bty="l")  
} else  if(k_kmeans == 3){
  vioplot(cl1$PCA_value,cl2$PCA_value,cl3$PCA_value, col=c('#e41a1c','#377eb8','#4daf4a'), main = "Env. distribution PCA",axes=F,names=c("Cl1","Cl2","C3"))
  box(bty="l")
  vioplot(cl1$RDAA_value,cl2$RDAA_value,cl3$RDAA_value, col=c('#e41a1c','#377eb8','#4daf4a'), main = "Env. distribution of RDA",axes=F,names=c("Cl1","Cl2","C3"))
  box(bty="l")
  #boxplot(cl1$RDAA_value,cl2$RDAA_value,cl3$RDAA_value, col=c('#e41a1c','#377eb8','#4daf4a'), main = "Env. distribution of RDA",axes=F,names=c("Cl1","Cl2","C3"))
} else if(k_kmeans == 4){
  vioplot(cl1$PCA_value,cl2$PCA_value,cl3$PCA_value,cl4$PCA_value, col=c('#e41a1c','#377eb8','#4daf4a','#984ea3'), main = "Env. distribution PCA",axes=F,names=c("Cl1","Cl2","C3","Cl4"))
  box(bty="l")
  boxplot(cl1$RDAA_value,cl2$RDAA_value,cl3$RDAA_value,cl4$RDAA_value, col=c('#e41a1c','#377eb8','#4daf4a','#984ea3'), main = "Env. distribution of RDA",axes=F,names=c("Cl1","Cl2","C3","Cl4"))
  box(bty="l")
} else if(k_kmeans == 5){
  vioplot(cl1$RDAA_value,cl2$RDAA_value,cl3$RDAA_value,cl4$RDAA_value,cl5$RDAA_value, col=c('#e41a1c','#377eb8','#4daf4a','#984ea3',"#ff7f00"), main = "Env. distribution of RDA",axes=F,names=c("Cl1","Cl2","C3","Cl4","Cl5"))
}
} else {
  print("yumi has no boxplot")
}
## Plotting genotyping versus PC1 backgeound ####
options(repr.plot.width=7, repr.plot.height=7)

ggplot()+ 
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill= "white",col=NA)+
  geom_raster(data = TAB_PCA, aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = 1, option = "D", labels = c("Positive","","","","Intermediate","","","","Negative")) +
  #geom_sf(data = admin, fill=NA, size=0.1) +
  #new_scale("fill")+
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill=NA,col="grey20",size = .2)+
  geom_point(data=classified_ind_6,aes(x=long3_j ,y=lat3_j), col="black", fill=classified_ind_6$kmeans_class_col, size=6, shape=21)+
  coord_cartesian(xlim=c(-67,-34),ylim = c(-17.5,7))+
  ggsn::scalebar(x.min=-67, x.max=-34, y.min=-17.5, y.max=7,dist = 250,dist_unit = "km", model = 'WGS84',st.size = 2.5,location="topright",transform=TRUE, border.size=.5)+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'grey75',size = .1),
        panel.background = element_rect(fill = 'white'))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  guides(fill=guide_legend(title="PCA Index")) +
  #ggtitle(my_title)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(colour = 'black', fill=NA, size=.5))

# What do I want tomorrow
# 1 - I want to include points of one yumi cluster with values within the range of the raw data. Image to the side represents the total range

##
## Vectorization of the RDA climatic rasters for ggplot ####
RDA_proj <- list(res_RDA_1_proj_current_t)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}
TAB_RDAA <- as.data.frame(RDA_proj[1])
colnames(TAB_RDAA)[3] <- "value"
TAB_RDAA$variable <- factor(rep("RDA1", nrow(TAB_RDAA)), levels = c("RDA1"))

classified_ind_7= classified_ind_6
classified_ind_7$kmeans
classified_ind_7$kmeans = factor(classified_ind_7$kmeans, levels = c(1,3,2))
classified_ind_7 = classified_ind_7[order(classified_ind_7$kmeans),]
classified_ind_7$kmeans

## Plotting kmeans classification per RDA
rda_plot = ggplot()+ 
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill= "white",col=NA)+
  geom_raster(data = TAB_RDAA, aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "D", labels = c("Positive","","","","Intermediate","","","","Negative")) +
  #geom_sf(data = admin, fill=NA, size=0.1) +
  #new_scale("fill")+
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill=NA,col="grey20",size = .2)+
  geom_point(data=classified_ind_7,aes(x=long3_j ,y=lat3_j), col="black", fill=classified_ind_7$kmeans_class_col, size=6, shape=21)+
  coord_cartesian(xlim=c(-67,-34),ylim = c(-17.5,7))+
  ggsn::scalebar(x.min=-67, x.max=-34, y.min=-17.5, y.max=7,dist = 250,dist_unit = "km", model = 'WGS84',st.size = 2.5,location="topright",transform=TRUE, border.size=.5)+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'grey75',size = .1),
        panel.background = element_rect(fill = 'white'))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  guides(fill=guide_legend(title= "Scores RDA 1")) +
  #ggtitle(cand_name)+
  theme(plot.title = element_text(hjust = 0))+
  theme(panel.border = element_rect(colour = 'black', fill=NA, size=.5))
rda_plot

## Plotting kmeans classification per RDA
rda_plot_yumi = ggplot()+ 
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill= "white",col=NA)+
  geom_raster(data = TAB_RDAA, aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "D", labels = c("Positive","","","","Intermediate","","","","Negative")) +
  #geom_sf(data = admin, fill=NA, size=0.1) +
  #new_scale("fill")+
  geom_polygon(data=SAshp, aes(long, lat, group=group), fill=NA,col="grey20",size = .2)+
  geom_point(data=class_ind_all,aes(x=long3_j ,y=lat3_j), col="black", fill=class_ind_all$kmeans_class_col, size=6, shape=21)+
  coord_cartesian(xlim=c(-67,-34),ylim = c(-17.5,7))+
  ggsn::scalebar(x.min=-67, x.max=-34, y.min=-17.5, y.max=7,dist = 250,dist_unit = "km", model = 'WGS84',st.size = 2.5,location="topright",transform=TRUE, border.size=.5)+
  theme_classic()+
  theme(panel.grid.major = element_line(colour = 'grey75',size = .1),
        panel.background = element_rect(fill = 'white'))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  guides(fill=guide_legend(title= "RDA scores")) +
  ggtitle(cand_name)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(colour = 'black', fill=NA, size=.5))
rda_plot_yumi

class_ind_6_pts_g1 = classified_ind_6_pts[classified_ind_6_pts$kmeans==1,]
class_ind_6_pts_g2 = classified_ind_6_pts[classified_ind_6_pts$kmeans==2,]
head(class_ind_6_pts_g2)

dev.off()
class_ind_6_mcp_g1 = adehabitatHR::mcp(class_ind_6_pts_g1,percent=100)
class_ind_6_buffer_g1 <- rgeos::gBuffer(class_ind_6_pts_g1, width = 400000)

class_ind_6_mcp_g2 = adehabitatHR::mcp(class_ind_6_pts_g2,percent=100)
class_ind_6_buffer_g2 <- rgeos::gBuffer(class_ind_6_pts_g2, width = 400000)

# Checking distribution and clipping the records
plot(classified_ind_6_pts,col=classified_ind_6_pts$kmeans_class_col, pch=16,cex=2)
plot(class_ind_6_mcp_g1,col=NA,add=TRUE)
plot(class_ind_6_buffer_g1,col=NA,add=TRUE)
plot(class_ind_6_mcp_g2,col=NA,add=TRUE)
plot(class_ind_6_buffer_g2,col=NA,add=TRUE)

indi_clust_1_pts = classified_ind_6_yumi
indi_clust_1_pts = indi_clust_1_pts[indi_clust_1_pts$kmeans==1,]
coordinates(indi_clust_1_pts) = indi_clust_1_pts[,c("long","lat")]
plot(indi_clust_1_pts,add=TRUE)
indi_clust_2_pts = classified_ind_6_yumi
indi_clust_2_pts = indi_clust_2_pts[indi_clust_2_pts$kmeans==2,]
coordinates(indi_clust_2_pts) = indi_clust_2_pts[,c("long","lat")]
plot(indi_clust_2_pts,pch=1,add=TRUE)

plot(classified_ind_6_pts,col=classified_ind_6_pts$kmeans_class_col, pch=16,cex=2)
indi_clust_1_pts_cr = crop(indi_clust_1_pts,class_ind_6_buffer_g1)
plot(indi_clust_1_pts_cr,add=TRUE)
# Final dataframe will have both molecular as well as the geographical

plot(classified_ind_6_pts,col=classified_ind_6_pts$kmeans_class_col, pch=16,cex=2)
indi_clust_2_pts_cr = crop(indi_clust_2_pts,class_ind_6_buffer_g2)
plot(indi_clust_2_pts_cr,add=TRUE)
# Final dataframe will have both molecular as well as the geographical

indi_clust_1_pts_cr_df = indi_clust_1_pts_cr@data[,c("lat","long","kmeans","kmeans_class_col","long2_j","lat2_j","long3_j","lat3_j")]
class_ind_6_pts_g1 = class_ind_6_pts_g1@data[,c("lat","long","kmeans","kmeans_class_col","long2_j","lat2_j","long3_j","lat3_j")]

indi_clust_2_pts_cr_df = indi_clust_2_pts_cr@data[,c("lat","long","kmeans","kmeans_class_col","long2_j","lat2_j","long3_j","lat3_j")]
class_ind_6_pts_g2 = class_ind_6_pts_g2@data[,c("lat","long","kmeans","kmeans_class_col","long2_j","lat2_j","long3_j","lat3_j")]

class_ind_all = rbind.data.frame(indi_clust_1_pts_cr_df,class_ind_6_pts_g1,indi_clust_2_pts_cr_df,class_ind_6_pts_g2)

class_ind_all_pts = class_ind_all
coordinates(class_ind_all_pts) = class_ind_all_pts[,c("long","lat")]
plot(class_ind_all_pts,col=class_ind_all_pts$kmeans_class_col, pch=16,cex=2)

write.csv(class_ind_all_pts,"/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/class_ind_all_pts.csv",row.names = FALSE)

return(rda_plot)
}

my_maps = pbapply::pbmapply(my_SNP_plotter,snp_cand,cand_name,SIMPLIFY = FALSE)

my_maps[[1]]
my_maps[[2]]
my_maps[[3]]
my_maps[[4]]
my_maps[[5]]

## Ultimo foi valores reais
dev.off()

# Without running the loop, line by line, the object classified_ind_6 contains all information needed for the modeling step. Now it is a matter of running it there

write.csv(classified_ind_6,"/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/SNMF_all_SNPs_new2_Ken_calcarata/classified_ind_6.csv")


##
#plot(res_RDA_proj_current[[1]], col = rev(viridis::inferno(16)))
#res_RDA_proj_current
#points(classified_ind_4[classified_ind_4$Sample_ID_full=="Ken_87",]$long,classified_ind_4[classified_ind_4$Sample_ID_full=="Ken_87",]$lat)
#.####

# Extra ####
## Original Vectorization of the RDA rasters for ggplot (2 RDA axis)
RDA_proj <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}
## Original adaptive genetic turnover projected across  range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
ggplot(data = TAB_RDA) + 
  #geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_tile(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = 1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  #geom_sf(data = admin, fill=NA, size=0.1) +
  #coord_sf(xlim = c(-70, -30), ylim = c(-20, 10), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

##
#save.image(file = "/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_RDA/RDA_Ken_calcarata_2.RData") # 

##
#load("/Users/josue/Dropbox/3Serrapilheira_INPA_postdoc/kentropyx_genomic_analyses/Ken_RDA/RDA_Ken_calcarata_2.RData")

## How to choose PCs
#Another simple approach to decide on the number of principal components is to set a threshold, say 80%, and stop when the first k components account for a percentage of total variation greater than this threshold (Jolliffe 2002). 
#.####