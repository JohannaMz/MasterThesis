
#### Libraries and set up####

library(sf)
library(dplyr)
library(tmap)
library(raster)
library(rgdal)
library(terra)
library(tidyverse)
library(fasterize)
library(stars)
library(ggplot2)
library(maptools)
library(mapview)
library(sp)
library(spatstat)
library(grid)
library(RColorBrewer)
library(lubridate)


setwd("W:/Masterarbeit Daten/Analyse")
#load("WU_datapreparation.RData") #everything for hessen

Sys.setenv(lang = "en_US") #change error messages to english


####preparation for hessen#####
#set up projection to be used thorughout the script
projection <- "+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs" #crs(rehe_WU)
extent <- extent(412000, 587000, 5471000,  5723500) #xmin, xmax, ymin, ymax -> extent of Hessen and a bit outside
resolution <- 50 #meter grids  

#default raster to use for all hessen covariates
r.raster <- raster()
extent(r.raster) <- extent #extent of hessen
res(r.raster) <- 50
crs(r.raster) <- projection




#### Roe deer Wildlife accidents ####

# Load data on wildlife accidents with roe deer in germany 
rehe_WU <- read_sf("SDM/RehWU_shape/WU_rehe.shp") 

names(rehe_WU)
# [1] "rowid_"     "twu_id"     "strassen_i" "land"       "land_lfn"   "tagebuch_n" "funddatum"  "jahr"       "monat"      "tag"       
# [11] "stunde"     "minuten"    "gsb"        
# "kat"
# "abstand_or" "bez"        "wdm"        "strassen_k" "artid"      "geometry"  


# filter for only roe deer data from hessen -> presence points
rehe_WU$land <- as.factor(rehe_WU$land)
rehe_WU_HE <- filter(rehe_WU, land == "HE")
rehe_WU_HE <- st_transform(rehe_WU_HE, projection)



# plot data
germany_sf <- read_sf("SDM/gadm36_DEU_shp/gadm36_DEU_1.shp")
germany_sf <- st_transform(germany_sf, crs = projection)#change CRS to fit WU data

germany <- tm_shape(germany_sf) +
  tm_borders()

germany +
  tm_shape(rehe_WU_HE) +
  tm_symbols(size = 0.5)


#some outliers in RLP and NRW -> crop WU points to Hessen polygon
hessen_sf <- filter(germany_sf, NAME_1 == "Hessen")
hessen_sf$VARNAME_1 <- "HE"
bawu_sf <- filter(germany_sf, NAME_1 == "Baden-Württemberg")
bawu_sf$VARNAME_1 <- "BW"

rehe_WU_HE <- st_crop(rehe_WU_HE, hessen_sf)

#looks better now. 
germany +
  tm_shape(rehe_WU_HE) +
  tm_symbols(size = 0.5)


#### ATKIS Strassen layer (metadatentabelle: ATKIS_Objektartenkatalog_Basis_DLM_7.1) ####
roads_HE <- read_sf("Layers/strassen_HE.shp")

names(roads_HE)
# "bkg_id"     "land"       "objart"     "objid"      "beginn"     "objart_z"   
# "bez" (=bezeichnug)       
# "brf" (= breiteFahrbahn)
summary(roads_HE$brf)
roads_HE$brf[roads_HE$brf == -9998] <- NA #replace the unknowns with NA (only 130 obs without width data)
# "bvb"(=besondere verkehrsbedeutung): 1000 (durchgangsverkehr), 1003 (nahverkehr, zwischen?rtlich), 2000 (ortsverkehr), 2001 (sammelverkehr), 2002 (anliegerverkehr)        
summary(as.factor(roads_HE$bvb)) # all the same: durchgangsverkehr
# "wdm"(=widmung): 1301 (bundesautobahn), 1303 (bundesstrasse = federal), 1305 (landesstrasse, staatsstra?e = state), 1306 (kreisstrasse = district), 1307 (gemeindestrasse = community), 9997 (attribut trifft nicht zu), 9999 (sonstiges)        
roads_HE$wdm <- as.numeric(roads_HE$wdm)
summary(as.factor(roads_HE$wdm))
# 1301 1303 1305 1306 1307 9997 
# 2068 4065 6173 3939 3784   51 
roads_HE <- filter(roads_HE, wdm != 9997) #take out the unknwon roads
roads_HE$wdm_name <- NA
roads_HE$wdm_name <- roads_HE$wdm
roads_HE$wdm_name[roads_HE$wdm == 1301] <- "Highways"
roads_HE$wdm_name[roads_HE$wdm == 1303] <- "Federal roads"
roads_HE$wdm_name[roads_HE$wdm == 1305] <- "State roads"
roads_HE$wdm_name[roads_HE$wdm == 1306] <- "District roads"
roads_HE$wdm_name[roads_HE$wdm == 1307] <- "Community roads"
roads_HE$wdm_name <- factor(roads_HE$wdm_name, 
                            levels = c("Community roads", "District roads", "State roads", "Federal roads", "Highways"))

# "zus"       "bemerkung"  "id"         "strassen_k" "geometry"  
roads_HE <- st_transform(roads_HE, projection)

#plot to check
hessen <- tm_shape(hessen_sf) +
  tm_borders()
hessen +
  tm_shape(roads_HE) +
  tm_lines(size = 0.5)



#### snap coordinates #### 
#change coordinates of the points that are just next to the road lines: snap 


st_snap_points = function(x, y, max_dist = 200) {
  
  if (inherits(x, "sf")) n = nrow(x)
  if (inherits(x, "sfc")) n = length(x)
  
  out = do.call(c,
                lapply(seq(n), function(i) {
                  nrst = st_nearest_points(st_geometry(x)[i], y)
                  nrst_len = st_length(nrst)
                  nrst_mn = which.min(nrst_len)
                  if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
                  return(st_cast(nrst[nrst_mn], "POINT")[2])
                })
  )
  return(out)
}

tst <- st_snap_points(rehe_WU_HE, roads_HE, max_dist = 150) 

rehe_WU_HE$geometry <- tst


#check how many are outside the buffer of 150 meters
roads_HE_buffered150 <- st_buffer(roads_HE, dist = 150)
outside <- sapply(st_intersects(rehe_WU_HE, roads_HE_buffered150),function(x){length(x)==0})

rehe_WU_HE_outside <- rehe_WU_HE[outside, ] #283 points, the same number as indicated by outside of pixels!
rehe_WU_HE <- rehe_WU_HE[!outside, ]

rm(roads_HE_buffered150, outside, rehe_WU_HE_outside)
st_write(rehe_WU_HE, "Layers/rehe_WU_HE.shp")


#### rasterize street map into 50x50m grids: extract wdm and brf####

# after 35 m buffer: points are not included because they are outside the buffer, not because they are missed between pixels. 
roads_HE_buffered <- st_buffer(roads_HE, dist = 35)
#mean width decreases in size 
roads_HE %>% 
  st_drop_geometry() %>% 
  group_by(wdm) %>% 
  summarise(mean_width = mean(brf, na.rm = TRUE))

#rasterize all roads based on their wdm: if multiple take "bigger" road category
roads_HE_raster_wdm <- fasterize(roads_HE_buffered, r.raster, field = "wdm", fun="min")

# rasterize all roads based on their width: if multiple, take largest width
roads_HE_raster_brf <- fasterize(roads_HE_buffered, r.raster, field = "brf", fun="max")

writeRaster(roads_HE_raster_wdm, "Layers/unscaled predictors/HE/roads_HE_raster_wdm.tiff")
writeRaster(roads_HE_raster_brf, "Layers/unscaled predictors/HE/roads_HE_raster_brf.tiff")


#### calculate road density ####
gemeindestrassen_HE <- read_sf("Layers/shp_gep_gemeindestrasse_HE.shp")
names(gemeindestrassen_HE)

gemeindestrassen_HE <- gemeindestrassen_HE %>% 
  dplyr::select(strassen_i, land, objart, objid, beginn, objart_z, bez, brf, bvb, wdm, zus, id, geometry) %>% 
  rename(bkg_id = strassen_i) %>% 
  mutate(bemerkung = NA, 
         strassen_k = NA, 
         wdm_name = "Community roads")
gemeindestrassen_HE <- st_transform(gemeindestrassen_HE, projection)

#adding gemeindestraßen (mostly the small city roads!), as they are not included in the normal roads layer but be wuite important for the density calculation!
#delete the gemeindestrassen in the road layer, so that those that are actually included are not doubled
roads_density <- rbind(filter(roads_HE, wdm != 1307), gemeindestrassen_HE) 

roads_density_sp <- as(roads_density, "Spatial")
roads_density_psp <- as.psp(roads_density_sp)
Sys.time()


# 50 doesnt make sense, because in most cases it is anyways just the one road going through. 
# but this layer can be used to use focal calculations in a  200 meter moving window
roads_density_px <- pixellate(roads_density_psp, eps = 50)  
roads_density_raster <- raster(roads_density_px)
crs(roads_density_raster) <- projection

#adjust to the right resolution, extent and projection!
roads_density_raster <- projectRaster(roads_density_raster, r.raster, method = "ngb")

#for now the unit is meters per 50*50m = 0.0025km^2
# recalculate it into km/km^2: *400/1000 = 0.4

summary(roads_density_raster)
values(roads_density_raster) <- values(roads_density_raster)*0.4


#calculate mean of density in a moving window of 9x9 pixels 
roads_density_raster_focal <- terra::focal(roads_density_raster, w=matrix(1, nrow = 9, ncol = 9), mean)
#values     : 0, 31.85667  (min, max)

roads_density_raster_focal_100m <- terra::focal(roads_density_raster, w=matrix(1, nrow = 5, ncol = 5), mean)
#values     : 0, 42.19112  (min, max)


#then mask on road network:
roads_density_masked <- mask(roads_density_raster_focal, roads_HE_raster_wdm)
roads_density_masked_100m <- mask(roads_density_raster_focal_100m, roads_HE_raster_wdm)
hist(roads_density_masked_100m)


writeRaster(roads_density_masked, "Layers/unscaled predictors/roads_density_masked.tiff", overwrite = TRUE)
writeRaster(roads_density_masked_100m, "Layers/unscaled predictors/roads_density_masked_100m.tiff", overwrite = TRUE)




#### landcover layers ####

landcover <- raster("CLC2018_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif")

#coordinate system: EPSG:3035 (ETRS89, LAEA)

EPSG3035 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
#crs(landcover) <- EPSG3035

#crop only Hessen

germany_transformed <- st_transform(germany_sf, crs = EPSG3035)
#hessen_transformed <- st_transform(hessen_sf, crs = EPSG3035)

landcover_GE <- crop(landcover, germany_transformed)

landcover_HE <- projectRaster(landcover_GE, r.raster, method = "ngb") 
# r.raster: use template raster to get the same extent, projection and resolution. be careful with values, needs method ngb!!

# reclassification of CLC based on Hothorn2012 (S1)
unique(values(landcover_HE))
hist(landcover_HE) 


#clc landcover codes
clc_legend <- read.csv("CLC2018_raster100m/clc_legend.csv")

#add reclassified column to clc_legend
clc_legend$reclass <- NA
# meadows (2.3.1., 2.4.3., 3.2.1)
clc_legend$reclass[clc_legend$GRID_CODE == 18 | clc_legend$GRID_CODE== 21 | clc_legend$GRID_CODE== 26]  <- 1
# swamps (3.2.2., 4.1.1., 4.1.2)
clc_legend$reclass[clc_legend$GRID_CODE == 27 | clc_legend$GRID_CODE== 35 | clc_legend$GRID_CODE== 36]  <- 2
# industry (1.2.1., 1.2.4., 1.3.1., 1.3.2.)
clc_legend$reclass[clc_legend$GRID_CODE == 3 | clc_legend$GRID_CODE== 6 | clc_legend$GRID_CODE== 7 | clc_legend$GRID_CODE== 8]  <- 3
# urban areas (1.1.1., 1.1.2., 1.2.3., 1.3.3., 1.4.1,1.4.2.) 
clc_legend$reclass[clc_legend$GRID_CODE == 1 | clc_legend$GRID_CODE== 2 | clc_legend$GRID_CODE== 5 | clc_legend$GRID_CODE== 9 | clc_legend$GRID_CODE== 10 | clc_legend$GRID_CODE== 11]  <- 4
# complex habitat (2.4.2., 2.2.2., 3.2.4., 2.2.1.) 
clc_legend$reclass[clc_legend$GRID_CODE == 20 | clc_legend$GRID_CODE== 16 | clc_legend$GRID_CODE== 29 | clc_legend$GRID_CODE== 15]  <- 5
# conifer forests (3.1.2.)
clc_legend$reclass[clc_legend$GRID_CODE == 24]  <- 6
# mixed forests (3.1.3.) 
clc_legend$reclass[clc_legend$GRID_CODE == 25]  <- 7
# broadleaf forests (3.1.1.)
clc_legend$reclass[clc_legend$GRID_CODE == 23]  <- 8
# arable areas (2.1.1.)
clc_legend$reclass[clc_legend$GRID_CODE == 12]  <- 9


landcover_HE_reclass <- subs(landcover_HE, dplyr::select(clc_legend, GRID_CODE, reclass), subsWithNA = TRUE) 
hist(landcover_HE_reclass)
unique(values(landcover_HE_reclass))

# calculate landcover percentage 1 to 9 per pixel

for (i in 1:9){
  # write a function that returns percent of classes: values from 1 to 9
  pclass <- function(x, y=i) {
    return(length(which(x %in% y)) / length(x))
  }
  
  # pass pclass function to focal function to return class percent within a 9x9 window (200 meters each direction)
  lc_prop <- terra::focal(landcover_HE_reclass, w=matrix(1, nrow = 9, ncol = 9), pclass)
  
  #assign the names with each landcover value
  assign(paste0("landcover_HE_", i), lc_prop)
  
  #check how long each loop takes
  print(Sys.time())
}


#do the same for a smaller moving window
for (i in 1:9){
  # write a function that returns percent of classes: values from 1 to 9
  pclass <- function(x, y=i) {
    return(length(which(x %in% y)) / length(x))
  }
  
  # pass pclass function to focal function to return class percent within a  5x5 m window (100 meters each direction)
  lc_prop <- terra::focal(landcover_HE_reclass, w=matrix(1, nrow = 5, ncol = 5), pclass)
  
  #assign the names with each landcover value
  assign(paste0("landcover_HE_small_", i), lc_prop)
  
  #check how long each loop takes
  print(Sys.time())
}


#check if all new layers have the same extent, projection and resolution -> yes
landcover_HE_vector <- c(landcover_HE_1, landcover_HE_2, landcover_HE_3, 
                         landcover_HE_4, landcover_HE_5, landcover_HE_6, landcover_HE_7, 
                         landcover_HE_8, landcover_HE_9)


#mask all layers to road network
for (i in 1:9){
  data <- landcover_HE_vector[[i]]
  masked <- mask(data, roads_HE_raster_wdm)
  
  #assign the names with each landcover value
  assign(paste0("landcover_HE_", i, "_masked"), masked)
  
  #check how long each loop takes
  print(Sys.time())
}

#same for small landcover layers
landcover_HE_small_vector <- c(landcover_HE_small_1, landcover_HE_small_2, landcover_HE_small_3, 
                         landcover_HE_small_4, landcover_HE_small_5, landcover_HE_small_6, landcover_HE_small_7, 
                         landcover_HE_small_8, landcover_HE_small_9)

#mask all layers to road network
for (i in 1:9){
  data <- landcover_HE_small_vector[[i]]
  masked <- mask(data, roads_HE_raster_wdm)
  
  #assign the names with each landcover value
  assign(paste0("landcover_HE_small_", i, "_masked"), masked)
  
  #check how long each loop takes
  print(Sys.time())
}

writeRaster(landcover_HE_1_masked, "Layers/unscaled predictors/landcover_HE_1_masked.tiff")
writeRaster(landcover_HE_2_masked, "Layers/unscaled predictors/landcover_HE_2_masked.tiff")
writeRaster(landcover_HE_3_masked, "Layers/unscaled predictors/landcover_HE_3_masked.tiff")
writeRaster(landcover_HE_4_masked, "Layers/unscaled predictors/landcover_HE_4_masked.tiff")
writeRaster(landcover_HE_5_masked, "Layers/unscaled predictors/landcover_HE_5_masked.tiff")
writeRaster(landcover_HE_6_masked, "Layers/unscaled predictors/landcover_HE_6_masked.tiff")
writeRaster(landcover_HE_7_masked, "Layers/unscaled predictors/landcover_HE_7_masked.tiff")
writeRaster(landcover_HE_8_masked, "Layers/unscaled predictors/landcover_HE_8_masked.tiff")
writeRaster(landcover_HE_9_masked, "Layers/unscaled predictors/landcover_HE_9_masked.tiff")

writeRaster(landcover_HE_small_1_masked, "Layers/unscaled predictors/landcover_HE_1_small_masked.tiff")
writeRaster(landcover_HE_small_2_masked, "Layers/unscaled predictors/landcover_HE_2_small_masked.tiff")
writeRaster(landcover_HE_small_3_masked, "Layers/unscaled predictors/landcover_HE_3_small_masked.tiff")
writeRaster(landcover_HE_small_4_masked, "Layers/unscaled predictors/landcover_HE_4_small_masked.tiff")
writeRaster(landcover_HE_small_5_masked, "Layers/unscaled predictors/landcover_HE_5_small_masked.tiff")
writeRaster(landcover_HE_small_6_masked, "Layers/unscaled predictors/landcover_HE_6_small_masked.tiff")
writeRaster(landcover_HE_small_7_masked, "Layers/unscaled predictors/landcover_HE_7_small_masked.tiff")
writeRaster(landcover_HE_small_8_masked, "Layers/unscaled predictors/landcover_HE_8_small_masked.tiff")
writeRaster(landcover_HE_small_9_masked, "Layers/unscaled predictors/landcover_HE_9_small_masked.tiff")


#### stack all raster layers ####
raster_HE_stack <- stack(c(landcover_HE_1_masked, landcover_HE_2_masked, landcover_HE_3_masked, 
                              landcover_HE_4_masked, landcover_HE_5_masked, landcover_HE_6_masked, landcover_HE_7_masked, 
                              landcover_HE_8_masked, landcover_HE_9_masked, 
                           landcover_HE_small_1_masked, landcover_HE_small_2_masked, landcover_HE_small_3_masked, 
                           landcover_HE_small_4_masked, landcover_HE_small_5_masked, landcover_HE_small_6_masked, 
                           landcover_HE_small_7_masked, landcover_HE_small_8_masked, landcover_HE_small_9_masked,
                            roads_HE_raster_brf, roads_density_masked, roads_density_masked_100m, roads_HE_raster_wdm))
names(raster_HE_stack) <- c("Meadows1_200m", "Swamps2_200m", "Industry3_200m", "Urbanareas4_200m", "Complex habitats5_200m", 
                            "ConiferForests6_200m", "MixedForests7_200m", "BroadleafedForests8_200m", "ArableAreas9_200m", 
                            "Meadows1_100m", "Swamps2_100m", "Industry3_100m", "Urbanareas4_100m", "Complex habitats5_100m", 
                            "ConiferForests6_100m", "MixedForests7_100m", "BroadleafedForests8_100m", "ArableAreas9_100m",
                            "RoadWidth", "RoadDensity_200m", "RoadDensity_100m","RoadCategory")

plot(raster_HE_stack)



raster_HE_vector <- c(landcover_HE_1_masked, landcover_HE_2_masked, landcover_HE_3_masked, 
                           landcover_HE_4_masked, landcover_HE_5_masked, landcover_HE_6_masked, landcover_HE_7_masked, 
                           landcover_HE_8_masked, landcover_HE_9_masked, 
                           landcover_HE_small_1_masked, landcover_HE_small_2_masked, landcover_HE_small_3_masked, 
                           landcover_HE_small_4_masked, landcover_HE_small_5_masked, landcover_HE_small_6_masked, 
                           landcover_HE_small_7_masked, landcover_HE_small_8_masked, landcover_HE_small_9_masked,
                           roads_HE_raster_brf, roads_density_masked, roads_density_masked_100m)



#plus standardization of non-categorical covariates (chyn et al), all except for layer 22 road category
#fastest version is looping rather than trying to scale the raster stack...
for (i in 1:21){
  data <- raster_HE_vector[[i]]
  scaled <- scale(data)
  
  #assign the names with each landcover value
  assign(paste0("raster_HE", i, "_scaled"), scaled)
  
  #check how long each loop takes
  print(Sys.time())
}


#final raster stack with road category added: could stack everyhting here, but cannot allocate the vector... so just download all seperatly.
# writeRaster(raster_HE_stack_final_scaled, 
#             filename="Predictor Rasters/raster_stack_scaled.tif")

writeRaster(raster_HE1_scaled, "Predictor Rasters/landcover_HE_1_scaled.tiff")
writeRaster(raster_HE2_scaled, "Predictor Rasters/landcover_HE_2_scaled.tiff")
writeRaster(raster_HE3_scaled, "Predictor Rasters/landcover_HE_3_scaled.tiff")
writeRaster(raster_HE4_scaled, "Predictor Rasters/landcover_HE_4_scaled.tiff")
writeRaster(raster_HE5_scaled, "Predictor Rasters/landcover_HE_5_scaled.tiff")
writeRaster(raster_HE6_scaled, "Predictor Rasters/landcover_HE_6_scaled.tiff")
writeRaster(raster_HE7_scaled, "Predictor Rasters/landcover_HE_7_scaled.tiff")
writeRaster(raster_HE8_scaled, "Predictor Rasters/landcover_HE_8_scaled.tiff")
writeRaster(raster_HE9_scaled, "Predictor Rasters/landcover_HE_9_scaled.tiff")

writeRaster(raster_HE10_scaled, "Predictor Rasters/landcover_HE_1_small_scaled.tiff")
writeRaster(raster_HE11_scaled, "Predictor Rasters/landcover_HE_2_small_scaled.tiff")
writeRaster(raster_HE12_scaled, "Predictor Rasters/landcover_HE_3_small_scaled.tiff")
writeRaster(raster_HE13_scaled, "Predictor Rasters/landcover_HE_4_small_scaled.tiff")
writeRaster(raster_HE14_scaled, "Predictor Rasters/landcover_HE_5_small_scaled.tiff")
writeRaster(raster_HE15_scaled, "Predictor Rasters/landcover_HE_6_small_scaled.tiff")
writeRaster(raster_HE16_scaled, "Predictor Rasters/landcover_HE_7_small_scaled.tiff")
writeRaster(raster_HE17_scaled, "Predictor Rasters/landcover_HE_8_small_scaled.tiff")
writeRaster(raster_HE18_scaled, "Predictor Rasters/landcover_HE_9_small_scaled.tiff")

writeRaster(raster_HE19_scaled, "Predictor Rasters/roads_HE_raster_brf_scaled.tiff")
writeRaster(raster_HE20_scaled, "Predictor Rasters/roads_density_HE_scaled.tiff")
writeRaster(raster_HE21_scaled, "Predictor Rasters/roads_density_HE_100m_scaled.tiff")



#### summary statistics####

#length of road categories
roads_HE$length <- as.numeric(st_length(roads_HE)) #in meters
sum(roads_HE$length)/1000

roadcat_length <- roads_HE %>% 
  st_drop_geometry() %>% 
  group_by(wdm) %>% 
  summarise(length_sum = sum(length)/1000)



#split of road categories per year
rehe_WU_HE$roadcat <- raster::extract(roads_HE_raster_wdm, as.data.frame(st_coordinates(rehe_WU_HE)))
summary(as.factor(rehe_WU_HE$roadcat))

roadcat_count <- rehe_WU_HE %>% 
  st_drop_geometry() %>% 
  count(roadcat, jahr)

roadcat_count <- roadcat_count %>% 
  left_join(roadcat_length, by = c("roadcat" = "wdm")) %>% 
  mutate(perc = n/length_sum)


(WU_per_roadcat <- ggplot(roadcat_count, aes(x = as.factor(roadcat), y = perc, fill = as.factor(jahr))) +
  geom_bar(stat = "identity",
           width = 0.5, position = position_dodge()) +
  theme_minimal() +
  theme(text = element_text(size = 25), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594"), name = "Year") +
  labs(y = "Relative occurence of DVC") +
  scale_x_discrete(name = "Road category", 
                   labels = c("Highways \n(1417 km)", "Federal roads \n(3333 km)", 
                              "State roads \n(7243 km)", "District roads \n(4945 km)", "Community roads \n(2751 km)")))


ggsave(plot = WU_per_roadcat, "Figures/WU_per_roadcat.pdf")


#WVC throughout the year
rehe_WU_HE$date_time <- as.POSIXct(rehe_WU_HE$funddatum, tz = "UTC")


rehe_WU_HE %>% 
  st_drop_geometry() %>% 
  na.omit() %>% 
  group_by(month(date_time), jahr) %>% 
  count() %>% 
  group_by(`month(date_time)`) %>% 
  summarise(mean = mean(n), 
            sd = sd(n)) %>% 
ggplot(aes(x = `month(date_time)`, y = mean)) +
  geom_line(col = "darkblue") +
  geom_ribbon(aes (ymin = mean - sd, ymax = mean +sd), alpha = 0.1) +
  theme_minimal()+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), 
                   labels = c("Jan","Feb","Mar","April","May","Jun","Jul","Aug","Sep","Okt","Nov","Dec")) +
  annotate("text", x = 3, y = 1600, label = "Gestation") +
  geom_vline(xintercept = 5, linetype = "longdash") +
  annotate("text", x = 6.2, y = 1600, label = "Lactation") +
  geom_vline(xintercept = 7.5, linetype = "longdash") +
  annotate("text", x = 8, y = 1600, label = "Rut") +
  geom_vline(xintercept = 8.5, linetype = "longdash") +
  annotate("text", x = 10, y = 1600, label = "Diapause") +
  labs(y = "Mean number of DVC per month", x = "Month of the year")

ggsave("Figures/WU_per_month.pdf")



rehe_WU_HE %>% 
  st_drop_geometry() %>% 
  na.omit() %>% 
  group_by(yday(date_time), jahr) %>% 
  count() %>% 
  group_by(`yday(date_time)`) %>% 
  summarise(mean = mean(n), 
            sd = sd(n)) %>% 
  ggplot(aes(x = `yday(date_time)`, y = mean)) +
  geom_line(col = "darkblue") +
  geom_ribbon(aes (ymin = mean - sd, ymax = mean +sd), alpha = 0.1) +
  theme_minimal()+
  theme(text = element_text(size = 25)) +
  scale_x_continuous(breaks = c(1,32,60,91,121,152,182,213,244,274,305,335), 
                     labels = c("Jan","Feb","Mar","April","May","Jun","Jul","Aug","Sep","Okt","Nov","Dec")) +
  annotate("text", x = 60, y = 70, label = "Gestation", size = 7) +
  geom_vline(xintercept = 121, linetype = "longdash") +
  annotate("text", x = 155, y = 70, label = "Lactation", size = 7) +
  geom_vline(xintercept = 196, linetype = "longdash") +
  annotate("text", x = 210, y = 70, label = "Rut", size = 7) +
  geom_vline(xintercept = 227, linetype = "longdash") +
  annotate("text", x = 300, y = 70, label = "Diapause", size = 7) +
  labs(y = "Mean number of DVC per day", x = "Month of the year")

ggsave("Figures/WU_per_day.pdf")


rehe_WU_HE %>% 
  st_drop_geometry() %>% 
  na.omit() %>% 
  group_by(hour(date_time), jahr) %>% 
  dplyr::count() %>% 
  dplyr::group_by(`hour(date_time)`) %>% 
  dplyr::summarise(mean = mean(n), 
            sd = sd(n)) %>% 
  ggplot(aes(x = `hour(date_time)`, y = mean)) +
  geom_line(col = "darkblue") +
  geom_ribbon(aes (ymin = mean - sd, ymax = mean +sd), alpha = 0.1) +
  theme_minimal()+
  labs(y = "Mean number of DVC per hour", x = "Time of the day")
 
ggsave("Figures/WU_per_hour.pdf")




#time of day/month heatmap 
# seasons comparison plot
rehe_WU_HE$season <- ifelse(rehe_WU_HE$monat >= 1 & rehe_WU_HE$monat <= 4, "Gestation", 
                            ifelse((rehe_WU_HE$monat == 5 | rehe_WU_HE$monat == 6)|
                                     (rehe_WU_HE$monat == 7 & rehe_WU_HE$tag <= 15) ,"Lactation", 
                                   ifelse((rehe_WU_HE$monat == 7 & rehe_WU_HE$tag >= 16)|
                                            (rehe_WU_HE$monat == 8 & rehe_WU_HE$tag <= 15), "Rut", "Diapause")))


library(plyr)
rehe_WU_HE$date <- format(as.Date(rehe_WU_HE$date_time), format = '%d-%m')
ddply_table <- rehe_WU_HE %>% 
  ddply(c("date","stunde"), summarise, N = length(date_time)/5) #because sampling went for 5 years

ddply_table_seasons <- rehe_WU_HE %>% 
  ddply(c("season", "stunde"), summarise, N = length(date_time)/5) #because 5 years

#divided by the length of periods to get per day
ddply_table_seasons$length <- ifelse(ddply_table_seasons$season == "Gestation", time_length(interval(ymd("2015-01-01"), ymd("2015-04-30")), "days"), 
                                     ifelse(ddply_table_seasons$season == "Lactation", time_length(interval(ymd("2015-05-01"), ymd("2015-07-15")), "days"), 
                                            ifelse(ddply_table_seasons$season == "Rut", time_length(interval(ymd("2015-07-16"), ymd("2015-08-15")), "days"), 
                                                   time_length(interval(ymd("2015-08-16"), ymd("2015-12-31")), "days"))))
ddply_table_seasons$N_mean <- ddply_table_seasons$N/ddply_table_seasons$length


detach("package:plyr", unload = TRUE)  #because dplyr and plyr do not like each other... 

dates <- as.factor(rehe_WU_HE$date)
ddply_background <- data.frame(date = rep(levels(dates), each = 24), 
                               stunde = rep(c(0:23), times = 366))

ddply_background %>% 
  dplyr::left_join(ddply_table, by = c("date", "stunde")) %>% 
  dplyr::mutate(N = coalesce(N, 0)) %>% 
 ggplot() + 
  geom_raster(aes(x=as.Date(date, format='%d-%m'), y=stunde, fill=N)) + 
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal()+
  theme(text = element_text(size = 25), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_breaks = "1 month", 
               date_labels =  c("Dec", "Jan","Feb","Mar","April","May","Jun","Jul","Aug","Sep","Okt","Nov","Dec", "Jan")) +
  scale_y_continuous(breaks = seq(0, 23, by = 2)) +
  labs(x = "Time of year (month)", y = "Time of day (hour)", fill = "Mean number \nof DVC")

ggsave("Figures/WU_heatmap.pdf", width = 25, height = 20, units = "cm")


#not a nice plot, but just to check if it is true that during the rut there is also mean more WVC? yes. but only slightly 
ggplot(ddply_table_seasons, aes(x = stunde, y = N_mean, col = season)) +
  geom_line()


seasons_calcN <- rehe_WU_HE %>% 
  st_drop_geometry() %>% 
  na.omit() %>% 
  group_by(season, jahr) %>% 
  dplyr::count()

seasons_calcN$length <- ifelse(seasons_calcN$season == "Gestation", time_length(interval(ymd("2015-01-01"), ymd("2015-04-30")), "days"), 
                                     ifelse(seasons_calcN$season == "Lactation", time_length(interval(ymd("2015-05-01"), ymd("2015-07-15")), "days"), 
                                            ifelse(seasons_calcN$season == "Rut", time_length(interval(ymd("2015-07-16"), ymd("2015-08-15")), "days"), 
                                                   time_length(interval(ymd("2015-08-16"), ymd("2015-12-31")), "days"))))
seasons_calcN$n_mean <- seasons_calcN$n/seasons_calcN$length


#number of WVC per day per season
seasons_calcN %>% 
  dplyr::group_by(season) %>% 
  dplyr::summarise(mean = mean(n_mean), 
                   sd = sd(n_mean))




#check the different accident categories 1 to 5: not very useful in a graph
rehe_WU_HE %>% 
  st_drop_geometry() %>% 
  count(kat, jahr) %>% 
  ggplot(aes(x = as.factor(kat), y = n, fill = as.factor(jahr))) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  theme_minimal() +
  scale_fill_manual(values=c("#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#084594"), name = "Year") +
  labs(y = "Number of DVC") +
  scale_x_discrete(name = "DVC Category")

rehe_WU_HE %>% 
  st_drop_geometry() %>% 
  count(kat)


#width related to road category: most roads are 5 meters with some outliers in all categories. 
# width might be very correlated with wdm... be careful!
ggplot(roads_HE, aes(x = as.factor(wdm), y = brf)) +
  geom_boxplot()
hist(roads_HE$brf)



#Hessen overview map 

insetmap <- germany + 
  tm_shape(hessen_sf) +
  tm_borders(lwd = 3, col = "black") +
  tm_text("VARNAME_1", size = 0.8, bg.color = "white") +
  tm_shape(bawu_sf) +
  tm_borders(lwd = 1) +
  tm_text("VARNAME_1", size = 0.8, bg.color = "white") 
  

(mainmap <- hessen +
  tm_shape(roads_HE) +
  tm_lines(col = "wdm_name", 
           palette=c("#FB8072", "#BEBADA", "#80B1D3", "#FDB462", "#B3DE69"), title.col = "Road category") +
  tm_legend(position=c("right", "top"), frame = FALSE, legend.title.size = 1.5, legend.text.size = 1) +
  tm_layout(inner.margins = c(0, 0, 0, 0.3), frame = FALSE) +
  #extras: nordpfeil & scale bar
  tm_scale_bar(position = c("left", "bottom"), width = 0.15, text.size = 1) +
  tm_compass(position = c("left", "top"), size = 2))


#add insetmap to layout
w <- 0.5
h <- 0.5
vp = viewport(0.8, 0.27, width = w, height = h)
print(insetmap, vp = vp)


tmap_save(mainmap,filename="Figures/overview_map.pdf",
          dpi=100, insets_tm=insetmap, insets_vp=vp)


#overview road density 
roads_density_raster_focal_crop <- crop(roads_density_raster_focal, extent(471900, 472800, 5557000, 5557700))
roads_density_crop <- st_crop(roads_density, extent(471900, 472800, 5557000, 5557700))

(roads_density_map <- tm_shape(roads_density_raster_focal_crop) +
  tm_raster(n = 9.5, title = "Mean road density") +
  tm_legend(position=c(.71, 0.1), frame = FALSE, legend.title.size = 1.5, legend.text.size = 1) +
  tm_layout(inner.margins = c(0, 0, 0, 0.3), frame = FALSE) +
  tm_shape(roads_density_crop) +
  tm_lines(size = 2) +
  #extras: nordpfeil & scale bar
  tm_scale_bar(position = c("left", "bottom"),breaks = c(0, 0.15, 0.3), text.size = 1) +
  tm_compass(position = c("left", "top"), size = 2))


tmap_save(roads_density_map, "Figures/roads_density_map.pdf", asp = 0)



save.image(file='WU_datapreparation.RData')





##### datapreparation for BaWü #####
load("WU_prepration_BW.Rdata")


#set up projection to be used throughout the script!
projection <- "+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs" #crs(rehe_WU)
extent.bw <- extent(388200, 610200, 5261000, 5516000) #xmin, xmax, ymin, ymax -> extent of bawu and a bit outside
resolution <- 50 #meter grids  

#default raster to use for all bawu covariates
r.raster_bw <- raster()
extent(r.raster_bw) <- extent.bw #extent of bawu
res(r.raster_bw) <- 50
crs(r.raster_bw) <- projection


#### ATKIS Strassen layer (metadatentabelle: ATKIS_Objektartenkatalog_Basis_DLM_7.1) ####

roads_BW <- read_sf("Layers/strassen_BW.shp")

names(roads_BW)
summary(roads_BW$brf)
roads_BW$brf[roads_BW$brf == -9998] <- NA #replace the unknowns with NA (only 95 obs without width data)
# "bvb"(=besondere verkehrsbedeutung): 1000 (durchgangsverkehr), 1003 (nahverkehr, zwischen?rtlich), 2000 (ortsverkehr), 2001 (sammelverkehr), 2002 (anliegerverkehr)        
summary(as.factor(roads_BW$bvb)) # all the same: durchgangsverkehr
# "wdm"(=widmung): 1301 (bundesautobahn), 1303 (bundesstrasse = federal), 1305 (landesstrasse, staatsstra?e = state), 1306 (kreisstrasse = district), 1307 (gemeindestrasse = community), 9997 (attribut trifft nicht zu), 9999 (sonstiges)        
roads_BW$wdm <- as.numeric(roads_BW$wdm)
summary(as.factor(roads_BW$wdm))
# 1301 1303 1305 1306 1307 9997 
# 1644 6701 9085 9578 7524   59 

roads_BW <- filter(roads_BW, wdm != 9997) #take out the  unknwon roads
roads_BW$wdm_name <- NA
roads_BW$wdm_name <- roads_BW$wdm
roads_BW$wdm_name[roads_BW$wdm == 1301] <- "Highways"
roads_BW$wdm_name[roads_BW$wdm == 1303] <- "Federal roads"
roads_BW$wdm_name[roads_BW$wdm == 1305] <- "State roads"
roads_BW$wdm_name[roads_BW$wdm == 1306] <- "District roads"
roads_BW$wdm_name[roads_BW$wdm == 1307] <- "Community roads"
roads_BW$wdm_name <- factor(roads_BW$wdm_name, 
                            levels = c("Community roads", "District roads", "State roads", "Federal roads", "Highways"))


roads_BW <- st_transform(roads_BW, projection)


roads_BW_buffered <- st_buffer(roads_BW, dist = 35)


#rasterize all roads based on their wdm: if multiple take "bigger" road category
roads_BW_raster_wdm <- fasterize(roads_BW_buffered, r.raster_bw, field = "wdm", fun="min")

# rasterize all roads based on their width: if multiple, take largest width
roads_BW_raster_brf <- fasterize(roads_BW_buffered, r.raster_bw, field = "brf", fun="max")

writeRaster(roads_BW_raster_wdm, "Layers/unscaled predictors/roads_BW_raster_wdm.tiff")
writeRaster(roads_BW_raster_brf, "Layers/unscaled predictors/roads_BW_raster_brf.tiff")


#### calculate road density for bw####
gemeindestrassen_BW <- read_sf("Layers/shp_gep_gemeindestrasse_BW.shp")
names(gemeindestrassen_BW)

gemeindestrassen_BW <- gemeindestrassen_BW %>% 
  dplyr::select(strassen_i, land, objart, objid, beginn, objart_z, bez, brf, bvb, wdm, zus, id, geometry) %>% 
  rename(bkg_id = strassen_i) %>% 
  mutate(bemerkung = NA, 
         strassen_k = NA, 
         wdm_name = "Community roads")
gemeindestrassen_BW <- st_transform(gemeindestrassen_BW, projection)

#adding gemeindestraßen, as they are not included in the normal roads layer but be wuite important for the density calculation!
#delete the gemeindestrassen in the road layer, so that those that are actually included are not doubled: some zeros now, because gemeindestrassen were deleted but not included in the new layer...
roads_density_BW <- rbind(filter(roads_BW, wdm != 1307), gemeindestrassen_BW) 


roads_density_BW_sp <- as(roads_density_BW, "Spatial")
roads_density_BW_psp <- as.psp(roads_density_BW_sp)
roads_density_BW_px <- pixellate(roads_density_BW_psp, eps = 50)  
roads_density_BW_raster <- raster(roads_density_BW_px)
crs(roads_density_BW_raster) <- projection


#adjust to the right resolution, extent and projection!
roads_density_BW_raster <- projectRaster(roads_density_BW_raster, r.raster_bw, method = "ngb")
#writeRaster(roads_density_BW_raster, "Layers/roads_density_BW_raster_new.tif") #saved, just to be sure. 

#for now the unit is meters per 50*50m = 0.0025km^2
# recalculate it into km/km^2: *400/1000 = 0.4

summary(roads_density_BW_raster)
values(roads_density_BW_raster) <- values(roads_density_BW_raster)*0.4


#calculate mean of density in a moving window of 9x9 pixels 
roads_density_BW_raster_focal <- terra::focal(roads_density_BW_raster, w=matrix(1, nrow = 9, ncol = 9), mean)
#values     : 0, 36.74135  (min, max)

roads_density_BW_raster_focal_100m <- terra::focal(roads_density_BW_raster, w=matrix(1, nrow = 5, ncol = 5), mean)
#values     : 0, 46.66733  (min, max)


#then mask on road network:
roads_density_BW_masked <- mask(roads_density_BW_raster_focal, roads_BW_raster_wdm)
roads_density_BW_masked_100m <- mask(roads_density_BW_raster_focal_100m, roads_BW_raster_wdm)
hist(roads_density_BW_masked_100m)

writeRaster(roads_density_BW_masked, "Layers/unscaled predictors/roads_density_BW_masked.tiff", overwrite = TRUE)
writeRaster(roads_density_BW_masked_100m, "Layers/unscaled predictors/roads_density_BW_masked_100m.tiff", overwrite = TRUE)



#landcover layers BW####

landcover_BW <- projectRaster(landcover_GE, r.raster_bw, method = "ngb") 

# reclassification of CLC based on Hothorn2012 (S1)
unique(values(landcover_BW))
hist(landcover_BW) 

landcover_BW_reclass <- subs(landcover_BW, dplyr::select(clc_legend, GRID_CODE, reclass), subsWithNA = TRUE) 
hist(landcover_BW_reclass)


# calculate landcover percentage 1 to 9 per pixel
# takes about 2 hours
for (i in 1:9){
  # write a function that returns percent of classes: values from 1 to 9
  pclass <- function(x, y=i) {
    return(length(which(x %in% y)) / length(x))
  }
  
  # pass pclass function to focal function to return class percent within a 9x9 window (200 meters each direction)
  lc_prop <- terra::focal(landcover_BW_reclass, w=matrix(1, nrow = 9, ncol = 9), pclass)
  
  #assign the names with each landcover value
  assign(paste0("landcover_BW_", i), lc_prop)
  
  #check how long each loop takes
  print(Sys.time())
}


landcover_BW_vector <- c(landcover_BW_1, landcover_BW_2, landcover_BW_3, 
                         landcover_BW_4, landcover_BW_5, landcover_BW_6, landcover_BW_7, 
                         landcover_BW_8, landcover_BW_9)

#mask all layers to road network
for (i in 1:9){
  data <- landcover_BW_vector[[i]]
  masked <- mask(data, roads_BW_raster_wdm)
  
  #assign the names with each landcover value
  assign(paste0("landcover_BW_", i, "_masked"), masked)
  
  #check how long each loop takes
  print(Sys.time())
}


#do the same for a smaller moving window
for (i in 1:9){
  # write a function that returns percent of classes: values from 1 to 9
  pclass <- function(x, y=i) {
    return(length(which(x %in% y)) / length(x))
  }
  
  # pass pclass function to focal function to return class percent within a  5x5 m window (100 meters each direction)
  lc_prop <- terra::focal(landcover_BW_reclass, w=matrix(1, nrow = 5, ncol = 5), pclass)
  
  #assign the names with each landcover value
  assign(paste0("landcover_BW_small_", i), lc_prop)
  
  #check how long each loop takes
  print(Sys.time())
}



#same for small landcover layers
landcover_BW_small_vector <- c(landcover_BW_small_1, landcover_BW_small_2, landcover_BW_small_3, 
                               landcover_BW_small_4, landcover_BW_small_5, landcover_BW_small_6, landcover_BW_small_7, 
                               landcover_BW_small_8, landcover_BW_small_9)

#mask all layers to road network
for (i in 1:9){
  data <- landcover_BW_small_vector[[i]]
  masked <- mask(data, roads_BW_raster_wdm)
  
  #assign the names with each landcover value
  assign(paste0("landcover_BW_small_", i, "_masked"), masked)
  
  #check how long each loop takes
  print(Sys.time())
}


writeRaster(landcover_BW_1_masked, "Layers/unscaled predictors/landcover_BW_1_masked.tiff")
writeRaster(landcover_BW_2_masked, "Layers/unscaled predictors/landcover_BW_2_masked.tiff")
writeRaster(landcover_BW_3_masked, "Layers/unscaled predictors/landcover_BW_3_masked.tiff")
writeRaster(landcover_BW_4_masked, "Layers/unscaled predictors/landcover_BW_4_masked.tiff")
writeRaster(landcover_BW_5_masked, "Layers/unscaled predictors/landcover_BW_5_masked.tiff")
writeRaster(landcover_BW_6_masked, "Layers/unscaled predictors/landcover_BW_6_masked.tiff")
writeRaster(landcover_BW_7_masked, "Layers/unscaled predictors/landcover_BW_7_masked.tiff")
writeRaster(landcover_BW_8_masked, "Layers/unscaled predictors/landcover_BW_8_masked.tiff")
writeRaster(landcover_BW_9_masked, "Layers/unscaled predictors/landcover_BW_9_masked.tiff")

writeRaster(landcover_BW_small_1_masked, "Layers/unscaled predictors/landcover_BW_1_small_masked.tiff")
writeRaster(landcover_BW_small_2_masked, "Layers/unscaled predictors/landcover_BW_2_small_masked.tiff")
writeRaster(landcover_BW_small_3_masked, "Layers/unscaled predictors/landcover_BW_3_small_masked.tiff")
writeRaster(landcover_BW_small_4_masked, "Layers/unscaled predictors/landcover_BW_4_small_masked.tiff")
writeRaster(landcover_BW_small_5_masked, "Layers/unscaled predictors/landcover_BW_5_small_masked.tiff")
writeRaster(landcover_BW_small_6_masked, "Layers/unscaled predictors/landcover_BW_6_small_masked.tiff")
writeRaster(landcover_BW_small_7_masked, "Layers/unscaled predictors/landcover_BW_7_small_masked.tiff")
writeRaster(landcover_BW_small_8_masked, "Layers/unscaled predictors/landcover_BW_8_small_masked.tiff")
writeRaster(landcover_BW_small_9_masked, "Layers/unscaled predictors/landcover_BW_9_small_masked.tiff")


#### stack all BW raster layers ####
raster_BW_stack <- stack(c(landcover_BW_1_masked, landcover_BW_2_masked, landcover_BW_3_masked, 
                           landcover_BW_4_masked, landcover_BW_5_masked, landcover_BW_6_masked, landcover_BW_7_masked, 
                           landcover_BW_8_masked, landcover_BW_9_masked, 
                           landcover_BW_small_1_masked, landcover_BW_small_2_masked, landcover_BW_small_3_masked, 
                           landcover_BW_small_4_masked, landcover_BW_small_5_masked, landcover_BW_small_6_masked, 
                           landcover_BW_small_7_masked, landcover_BW_small_8_masked, landcover_BW_small_9_masked,
                           roads_BW_raster_brf, roads_density_masked, roads_density_masked_100m, roads_BW_raster_wdm))
names(raster_BW_stack) <- c("Meadows1_200m", "Swamps2_200m", "Industry3_200m", "Urbanareas4_200m", "Complex habitats5_200m", 
                            "ConiferForests6_200m", "MixedForests7_200m", "BroadleafedForests8_200m", "ArableAreas9_200m", 
                            "Meadows1_100m", "Swamps2_100m", "Industry3_100m", "Urbanareas4_100m", "Complex habitats5_100m", 
                            "ConiferForests6_100m", "MixedForests7_100m", "BroadleafedForests8_100m", "ArableAreas9_100m",
                            "RoadWidth", "RoadDensity_200m", "RoadDensity_100m","RoadCategory")

plot(raster_BW_stack)



raster_BW_vector <- c(landcover_BW_1_masked, landcover_BW_2_masked, landcover_BW_3_masked, 
                      landcover_BW_4_masked, landcover_BW_5_masked, landcover_BW_6_masked, landcover_BW_7_masked, 
                      landcover_BW_8_masked, landcover_BW_9_masked, 
                      landcover_BW_small_1_masked, landcover_BW_small_2_masked, landcover_BW_small_3_masked, 
                      landcover_BW_small_4_masked, landcover_BW_small_5_masked, landcover_BW_small_6_masked, 
                      landcover_BW_small_7_masked, landcover_BW_small_8_masked, landcover_BW_small_9_masked,
                      roads_BW_raster_brf, roads_density_BW_masked, roads_density_BW_masked_100m)


#plus standardization of non-categorical covariates (chyn et al), all except for layer 22 road category
#fastest version is looping rather than trying to scale the raster stack...
for (i in 1:21){
  data <- raster_BW_vector[[i]]
  scaled <- scale(data)
  
  #assign the names with each landcover value
  assign(paste0("raster_BW", i, "_scaled"), scaled)
  
  #check how long each loop takes
  print(Sys.time())
}


#final raster stack with road category added: could stack everyhting here, but cannot allocate the vector... so just download all seperatly.
# writeRaster(raster_HE_stack_final_scaled, 
#             filename="Predictor Rasters/raster_stack_scaled.tif")

writeRaster(raster_BW1_scaled, "Predictor Rasters/landcover_BW_1_scaled.tiff")
writeRaster(raster_BW2_scaled, "Predictor Rasters/landcover_BW_2_scaled.tiff")
writeRaster(raster_BW3_scaled, "Predictor Rasters/landcover_BW_3_scaled.tiff")
writeRaster(raster_BW4_scaled, "Predictor Rasters/landcover_BW_4_scaled.tiff")
writeRaster(raster_BW5_scaled, "Predictor Rasters/landcover_BW_5_scaled.tiff")
writeRaster(raster_BW6_scaled, "Predictor Rasters/landcover_BW_6_scaled.tiff")
writeRaster(raster_BW7_scaled, "Predictor Rasters/landcover_BW_7_scaled.tiff")
writeRaster(raster_BW8_scaled, "Predictor Rasters/landcover_BW_8_scaled.tiff")
writeRaster(raster_BW9_scaled, "Predictor Rasters/landcover_BW_9_scaled.tiff")

writeRaster(raster_BW10_scaled, "Predictor Rasters/landcover_BW_1_small_scaled.tiff")
writeRaster(raster_BW11_scaled, "Predictor Rasters/landcover_BW_2_small_scaled.tiff")
writeRaster(raster_BW12_scaled, "Predictor Rasters/landcover_BW_3_small_scaled.tiff")
writeRaster(raster_BW13_scaled, "Predictor Rasters/landcover_BW_4_small_scaled.tiff")
writeRaster(raster_BW14_scaled, "Predictor Rasters/landcover_BW_5_small_scaled.tiff")
writeRaster(raster_BW15_scaled, "Predictor Rasters/landcover_BW_6_small_scaled.tiff")
writeRaster(raster_BW16_scaled, "Predictor Rasters/landcover_BW_7_small_scaled.tiff")
writeRaster(raster_BW17_scaled, "Predictor Rasters/landcover_BW_8_small_scaled.tiff")
writeRaster(raster_BW18_scaled, "Predictor Rasters/landcover_BW_9_small_scaled.tiff")

writeRaster(raster_BW19_scaled, "Predictor Rasters/roads_BW_raster_brf_scaled.tiff")
writeRaster(raster_BW20_scaled, "Predictor Rasters/roads_density_BW_scaled.tiff")
writeRaster(raster_BW21_scaled, "Predictor Rasters/roads_density_BW_100m_scaled.tiff")



save.image(file = "WU_prepration_BW.Rdata")



#compare hessen and baWü landcover

#possible extent with collared roe deer homerange polygons from WU_telemtry Script: roedeer_mcp_agg_sf
compare_HR <- mask(landcover_BW_reclass, roedeer_mcp_agg_sf)

compare_HR <- compare_HR %>% 
  as.data.frame() %>% 
  count(reclass) %>% 
  mutate(county = "HR (within BW)", 
         freq = n/ncell(compare_HR[!is.na(compare_HR$reclass)])) #n of raster cells that have a value


compare_HE <- landcover_HE_reclass %>% 
  as.data.frame() %>% 
  count(reclass) %>% 
  mutate(county = "HE", 
         freq = n/ncell(landcover_HE_reclass)) #n of raster cells

compare_BW <- landcover_BW_reclass %>% 
  as.data.frame() %>% 
  count(reclass) %>% 
  mutate(county = "BW", 
         freq = n/ncell(landcover_BW_reclass)) #n of raster cells

compare <- rbind(compare_HE, compare_BW, compare_HR)
compare <- filter(compare, !is.na(reclass))

(compare_lc_plot <- ggplot(compare, aes(x = as.factor(reclass), y = freq, fill = as.factor(county))) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal() +
    theme(text = element_text(size = 25), axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=c("#E69F00", "#56B4E9",  "lightgreen")) +
    labs(y = "Percentage of total area", fill = "County") +
    scale_x_discrete(name = "Landcover type", 
                     labels = c("Meadows", "Swamps", "Industry", "Urban areas", "Complex habitats", 
                                "Conifer forests", "Mixed forests", "Broad-leafed forests", "Arable lands")))


ggsave(plot = compare_lc_plot, "Figures/compare_lc_plot_withHR.pdf", width = 25, height = 20, units = "cm")
