## Libraries and set up ####
library(sf)
library(dplyr)
library(dismo) #interface with maxEnt
library(rJava)
library(raster) #spatial data manipulation
library(MASS) # 2D kernel density function
library(maptools) ##reading shapefiles
library(SDMtune)
library(ggplot2)
library(ggspatial)
library(mapview)
library(lubridate)
library(tidyr)
library(sftrack)
#library(pscl)
library(glmmTMB)
library(DHARMa)
library(xtable)
library(rstatix)
library(patchwork)
library(MuMIn)

##### Crossings analysis #####

## add roe deer telem data
projection <- "+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"

roedeer_telem <- read.csv2("Roe Deer Telemetry Data/roedeer_telem.csv")

roedeer_telem$ts <- as.POSIXct(roedeer_telem$ts, format = "%Y-%m-%d %H:%M:%S")
summary(roedeer_telem$ts) #some have no timestamps. delete:
roedeer_telem <- filter(roedeer_telem, !is.na(ts))

count(roedeer_telem, gps_validity) #11 is quite bad gps accuracy. but quite low amount. 


#make sf object from dataframe
roedeer_telem_sf <- st_as_sf(roedeer_telem, coords = c("longitude", "latitude"), crs = 4326) #WGS84 system
roedeer_telem_sf <- st_transform(roedeer_telem_sf, projection) #change to right crs



#### find road crossing locations ####

duplicates <- duplicated(roedeer_telem[,c(3,12)]) #check that there are no duplicates in the coluns animal id and ts
roedeer_telem_traj <- as_sftraj(roedeer_telem_sf, group = c(id = "animal_id"), time = "ts") #5 minuten
summary(roedeer_telem_traj)

roedeer_telem_traj_sf <- sf::st_as_sf(roedeer_telem_traj[, c(1:25, 27), drop = TRUE]) #change back to sf object
roedeer_telem_traj_sf$monat <- month(roedeer_telem_traj_sf$ts)
roedeer_telem_traj_sf$tag <- day(roedeer_telem_traj_sf$ts)
roedeer_telem_traj_sf$season <- ifelse(roedeer_telem_traj_sf$monat >= 1 & roedeer_telem_traj_sf$monat <= 4, "Gestation", 
                                       ifelse((roedeer_telem_traj_sf$monat == 5 | roedeer_telem_traj_sf$monat == 6)|
                                                (roedeer_telem_traj_sf$monat == 7 & roedeer_telem_traj_sf$tag <= 15) ,"Lactation", 
                                              ifelse((roedeer_telem_traj_sf$monat == 7 & roedeer_telem_traj_sf$tag >= 16)|
                                                       (roedeer_telem_traj_sf$monat == 8 & roedeer_telem_traj_sf$tag <= 15), "Rut", "Diapause")))



road_intersect <- st_intersection(roedeer_telem_traj_sf, roads_BW) 


#some points are multipoints: deer crossed the road twice within the same trajectory. make 2 points from this.
st_un_multipoint = function(x) {
  g = st_geometry(x)
  i = rep(seq_len(nrow(x)), sapply(g, nrow))
  x = x[i,]
  st_geometry(x) = st_sfc(do.call(c,
                                  lapply(g, function(geom) lapply(1:nrow(geom), function(i) st_point(geom[i,])))))
  x$original_geom_id = i
  x
}

#extract all MULTIPOINTS and split them into multiple POINTs
road_intersect_multipoints <- road_intersect %>% 
  filter(
    st_geometry_type(.)
    %in% c("MULTIPOINT") ) %>% 
  st_un_multipoint()
st_crs(road_intersect_multipoints) <- projection


road_intersect_final <- road_intersect %>% 
  filter(
    st_geometry_type(.)
    %in% c("POINT") ) %>% 
  mutate(original_geom_id = 0) %>% #new colum within the multipoints sf, combining the two points that belonged together
  rbind(., road_intersect_multipoints)

road_intersect_final$stunde <- hour(road_intersect_final$ts)
road_intersect_final$monat <- month(road_intersect_final$ts)

road_intersect_final$season <- ifelse(road_intersect_final$monat >= 1 & road_intersect_final$monat <= 4, "Gestation", 
                                      ifelse((road_intersect_final$monat == 5 | road_intersect_final$monat == 6)|
                                               (road_intersect_final$monat == 7 & road_intersect_final$tag <= 15) ,"Lactation", 
                                             ifelse((road_intersect_final$monat == 7 & road_intersect_final$tag >= 16)|
                                                      (road_intersect_final$monat == 8 & road_intersect_final$tag <= 15), "Rut", "Diapause")))



## define home ranges of roe deer with all telem datapoints
roedeer_telem_all <- read.csv2("Roe Deer Telemetry Data/roedeer_telem_all.csv")
head(roedeer_telem_all)

roedeer_telem_all$ts <- as.POSIXct(roedeer_telem_all$ts, format = "%Y-%m-%d %H:%M:%S")

#make sf object from dataframe
roedeer_telem_all_sf <- st_as_sf(roedeer_telem_all, coords = c("longitude", "latitude"), crs = 4326) #WGS84 system
roedeer_telem_all_sf <- st_transform(roedeer_telem_all_sf, projection) #change to right crs

#create SPdataframe
roedeer_telem_all_sp <- data.frame(ID = roedeer_telem_all_sf$animal_id, 
                                   X = st_coordinates(roedeer_telem_all_sf)[, 1], 
                                   Y = st_coordinates(roedeer_telem_all_sf)[, 2])
library(sp)
coordinates(roedeer_telem_all_sp) <- c("X", "Y")
proj4string(roedeer_telem_all_sp) <- projection

#calculate 99% minimum convex polygon
roedeer_mcp <- mcp(roedeer_telem_all_sp, percent = 99)
roedeer_mcp
#Since the input file used UTM, the area is in hectares by default.

roedeer_mcp_agg <- aggregate(roedeer_mcp)
mapview(roedeer_mcp_agg)
roedeer_mcp_agg_sf <- sf::st_as_sf(roedeer_mcp_agg)

#mask road network and road crossings to the HR of the roe deer
road_intersect_final <- st_intersection(road_intersect_final, roedeer_mcp_agg_sf) #-2 observations outside of HR
roads_BW_HR <- st_intersection(roads_BW, roedeer_mcp_agg_sf)


## identify WVC as crossings ####
                                                                  
######could there be the WVC as a crossing? -> non successful crossing!#
#ID1 death 25.22.2013 -> last collar date in 2011, so no. 
#summary(filter(roedeer_telem_all_sf, animal_id == 1)$ts)

#ID3 death 15.12.2010 -> last collar date on 22-11-2010, so no.
#summary(filter(roedeer_telem_all_sf, animal_id == 3)$ts)

#ID6 death 23.02.2011 -> last collar date on 23-02-2011
summary(filter(roedeer_telem_all_sf, animal_id == 6)$ts)
summary(filter(road_intersect_final, animal_id == 6)$ts) #last road crossing on 18-02-2011
#can you say that tht is the last crossing??? Not sure, but we can be sure that the accident happenend in this area so yes
plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 6 & date(ts) > "2011-01-18")))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 6 & ts ==  "2011-02-18 21:32:10"), add = TRUE, col = "red")
View(filter(road_intersect_final, animal_id == 6 & date(ts) ==  "2011-02-18"))

#new column to indicate the WVC points
road_intersect_final$case <- "crossing"
road_intersect_final$case[road_intersect_final$animal_id == 6 & 
                            road_intersect_final$ts ==  "2011-02-18 21:32:10"] <- "WVC"

#id12, death 21.09.2014 -> last collar date 2012, so no. 
summary(filter(roedeer_telem_all_sf, animal_id == 12)$ts)

#id 18, death 11.05.2011 -> last collar date 2012??? Probably wrong death date. Should be a year later. -> mess up needs to be changed in table<!! 
summary(filter(roedeer_telem_all_sf, animal_id == 18)$ts) #"2012-05-11 17:15:33" 
summary(filter(road_intersect_final, animal_id == 18)$ts) #last road crossing "2012-05-11 08:15:39"
plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 18 & date(ts) < "2012-05-11")))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 18 & ts ==  "2012-05-11 08:15:39"), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 18 & 
                            road_intersect_final$ts ==  "2012-05-11 08:15:39"] <- "WVC"

roedeer_ids$date_of_death[roedeer_ids$animal_id == 18] <- "2012-05-11"


#id27, death 08.02.2012 -> last collar date 08-02-2012
summary(filter(roedeer_telem_all_sf, animal_id == 27)$ts) #"2012-02-08 11:45:13" 
summary(filter(road_intersect_final, animal_id == 27)$ts) #last road crossing "2012-02-08 02:00:12"

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 27 & ts > "2012-02-08 02:00:12")))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 27 & ts ==  "2012-02-08 02:00:12"), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 27 & 
                            road_intersect_final$ts ==  "2012-02-08 02:00:12"] <- "WVC"


#id29, death 16.07.2012 -> last collar date 16-07-2012
summary(filter(roedeer_telem_all_sf, animal_id == 29)$ts) #"2012-07-16 10:15:11"
summary(filter(road_intersect_final, animal_id == 29)$ts) #"2012-07-16 03:15:51"

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 29 & date(ts) == "2012-07-16")))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 29 & ts ==  "2012-07-16 03:15:51"), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 29 & 
                            road_intersect_final$ts ==  "2012-07-16 03:15:51"] <- "WVC"

#id32, death 13.12.2011  -> last collar date 04.04.2012 ---> date do not not fit at all. 
summary(filter(roedeer_telem_all_sf, animal_id == 32)$ts) #"2012-04-04 18:45:53"
summary(filter(road_intersect_final, animal_id == 32)$ts) #"2012-04-02 04:00:49"

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 32 & ts > "2012-04-02 02:00:49")))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 32 & ts ==  "2012-04-02 04:00:49"), add = TRUE, col = "red")

#delete death date as it is completly unclear when she died, but definilty not on that day
roedeer_ids$date_of_death[roedeer_ids$animal_id == 32] <- NA


#id34, death 12.07.2012 -> last collar date 12-07-2012
summary(filter(roedeer_telem_all_sf, animal_id == 34)$ts) #"2012-07-12 22:00:44"
summary(filter(road_intersect_final, animal_id == 34)$ts) #"2012-07-05 23:45:44" 

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 34 & ts > "2012-07-05 23:45:44" )))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 34 & ts == "2012-07-05 23:45:44"), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 34 & 
                            road_intersect_final$ts == "2012-07-05 23:45:44"] <- "WVC"

#id36, death 24.07.2013 -> last collar date 24.07.2013
summary(filter(roedeer_telem_all_sf, animal_id == 36)$ts) #"2013-07-24 02:15:46"
summary(filter(road_intersect_final, animal_id == 36)$ts) #"2013-07-24 00:00:16"

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 36 & date(ts) == "2013-07-24" )))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 36 & ts == "2013-07-24 00:00:16"), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 36 & 
                            road_intersect_final$ts == "2013-07-24 00:00:16"] <- "WVC"


#id39, death 23.07.2013 -> last collar date same
summary(filter(roedeer_telem_all_sf, animal_id == 39)$ts) #"2013-07-23 22:15:20"
summary(filter(road_intersect_final, animal_id == 39)$ts) #"2013-07-22 20:45:15" 

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 39 & ts > "2013-07-22 20:45:15"  )))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 39 & ts == "2013-07-22 20:45:15" ), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 39 & 
                            road_intersect_final$ts == "2013-07-22 20:45:15"] <- "WVC"


#id45, death 26.07.2012, 26-07
summary(filter(roedeer_telem_all_sf, animal_id == 45)$ts) #"2012-07-26 01:45:18"
summary(filter(road_intersect_final, animal_id == 45)$ts) #"2012-07-23 23:30:15" 

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 45 & ts > "2012-07-23 22:30:15")))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 45 & ts == "2012-07-23 23:30:15"  ), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 45 & 
                            road_intersect_final$ts == "2012-07-23 23:30:15"] <- "WVC"

#id47, death 15.10.2012, 15-10-2012
summary(filter(roedeer_telem_all_sf, animal_id == 47)$ts) #""2012-10-15 19:00:13"
summary(filter(road_intersect_final, animal_id == 47)$ts) #"2012-10-15 03:45:10" 

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 47 & date(ts) == "2012-10-15")))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 47 & ts == "2012-10-15 03:45:10" & wdm == 1305), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 47 & 
                            road_intersect_final$ts == "2012-10-15 03:45:10" & road_intersect_final$wdm == 1305] <- "WVC"


#id54, death 16.04.2013, same date
summary(filter(roedeer_telem_all_sf, animal_id == 54)$ts) #"2013-04-16 05:30:10"
summary(filter(road_intersect_final, animal_id == 54)$ts) #"2013-04-15 18:30:50" 

plot(st_geometry(filter(roedeer_telem_traj_sf, animal_id == 54 & ts > "2013-04-15 18:00:50")))
plot(roads_BW_HR, add = TRUE)
plot(filter(road_intersect_final, animal_id == 54 & ts == "2013-04-15 18:30:50"), add = TRUE, col = "red")

road_intersect_final$case[road_intersect_final$animal_id == 54 & 
                            road_intersect_final$ts == "2013-04-15 18:30:50"] <- "WVC"



#### roe deer ID table and summaries ####
#understand roe deer ids and prepare table
roedeer_ids <- read.csv("Roe Deer Telemetry Data/animal_id_roe deer.csv", sep=";", na.strings = "")
roedeer_ids$area <- as.factor(roedeer_ids$area)
roedeer_ids$sex <- as.factor(roedeer_ids$sex)
roedeer_ids$status <- as.factor(roedeer_ids$status)
roedeer_ids$cause_of_death <- as.factor(roedeer_ids$cause_of_death)

roedeer_ids$date.of.end.of.collar...UTC <- as.POSIXct(roedeer_ids$date.of.end.of.collar...UTC, format = "%d.%m.%Y %H:%M")
roedeer_ids$end.date <- as.Date(roedeer_ids$date.of.end.of.collar...UTC, "UTC")

roedeer_ids$date_of_death <- as.POSIXct(roedeer_ids$date_of_death, format =  "%d.%m.%Y")

#delete 21:lieselotte - halsband kaputt. dates for 7 and 37???
roedeer_ids <- roedeer_ids[!(roedeer_ids$animal_id == 21),]

animal_ids <- unique(roedeer_telem_all_sf$animal_id)

telem_min_date <- as.Date(c())
telem_max_date <- as.Date(c())

for (i in animal_ids){
  telem_min_date[i] <- min(filter(roedeer_telem_all_sf, animal_id == i)$ts)
  telem_max_date[i] <- max(filter(roedeer_telem_all_sf, animal_id == i)$ts)
}

telem_min_date <- telem_min_date[!is.na(telem_min_date)]
telem_max_date <- telem_max_date[!is.na(telem_max_date)]

#start and end of collaring peiod based on telemetry data
roedeer_ids$telem_min_date <- telem_min_date
roedeer_ids$telem_max_date <- telem_max_date


#length of sampled days in roe deer seasons

days_count <- roedeer_ids %>%
  rowwise() %>%
  transmute(animal_id,
            date = list(seq(telem_min_date, telem_max_date, by = "day"))) %>%
  unnest(date)

days_count$monat <- month(days_count$date)
days_count$tag <- day(days_count$date)
days_count$season <- ifelse(days_count$monat >= 1 & days_count$monat <= 4, "Gestation", 
                            ifelse((days_count$monat == 5 | days_count$monat == 6)|
                                     (days_count$monat == 7 & days_count$tag <= 15) ,"Lactation", 
                                   ifelse((days_count$monat == 7 & days_count$tag >= 16)|
                                            (days_count$monat == 8 & days_count$tag <= 15), "Rut", "Diapause")))

roedeer_ids <- days_count %>% 
  group_by(animal_id) %>% 
  count() %>% 
  ungroup() %>%
  rename(diff_days = n) %>% 
  left_join(roedeer_ids, .)

roedeer_ids <- days_count %>% 
  filter(season == "Gestation") %>% 
  group_by(animal_id) %>% 
  count() %>% 
  ungroup() %>%
  rename(diff_days_gest = n) %>% 
  left_join(roedeer_ids, .)

roedeer_ids <- days_count %>% 
  filter(season == "Lactation") %>% 
  group_by(animal_id) %>% 
  count() %>% 
  ungroup() %>%
  rename(diff_days_lact = n) %>% 
  left_join(roedeer_ids, .)  

roedeer_ids <- days_count %>% 
  filter(season == "Rut") %>% 
  group_by(animal_id) %>% 
  count() %>% 
  ungroup() %>%
  rename(diff_days_rut = n) %>% 
  left_join(roedeer_ids, .)

roedeer_ids <- days_count %>% 
  filter(season == "Diapause") %>% 
  group_by(animal_id) %>% 
  count() %>% 
  ungroup() %>%
  rename(diff_days_dia = n) %>% 
  left_join(roedeer_ids, .)


roedeer_ids <- roedeer_ids %>% 
  mutate_at(vars(diff_days_gest:diff_days_dia), ~replace_na(., 0))


summary(roedeer_ids$diff_days)
sd(roedeer_ids$diff_days, na.rm = TRUE)


#cause of death
count(roedeer_ids, status)
roedeer_ids$cause_of_death[roedeer_ids$cause_of_death == "vehicle_collision"] <- "vehicle collision"
count(roedeer_ids, cause_of_death)


road_intersect_final <- roedeer_ids %>% 
  dplyr::select(area, animal_id) %>% 
  left_join(road_intersect_final, ., by = "animal_id")


# number of crossings per individual
roedeer_ids <- road_intersect_final %>% 
  st_drop_geometry() %>% 
  group_by(animal_id) %>% 
  count() %>% 
  left_join(roedeer_ids, ., by = "animal_id") %>% 
  rename("no_crossings" = "n")

summary(roedeer_ids$no_crossings) #471.4
sd(roedeer_ids$no_crossings, na.rm= TRUE)


#summary crossings over roads in relation to roads in home ranges 
summary(road_intersect_final$wdm_name)
                                                                  
#in relatino to length of road categories in home ranges
roads_BW_HR$length <- as.numeric(st_length(roads_BW_HR)) #in meters

roads_BW_HR %>% 
  st_drop_geometry() %>% 
  group_by(wdm) %>% 
  summarise(length_sum = sum(length)/1000) #in kilometers
#     wdm length_sum
#   <dbl>      <dbl>
# 1  1301       2.57 -> 134/2.57 -> 52.14 per kilometer
# 2  1305      10.9 -> 7421/10.9 -> 680.03 /km
# 3  1306      18.8 -> 10094/18.8 -> 536.91/km
# 4  1307       1.15 -> 3091/1.15 <- 2687.83/km


#prepare table for latex
roedeer_ids <- dplyr::select(roedeer_ids, area, animal_id, name, sex, 
                              telem_min_date, telem_max_date, diff_days, no_crossings, 
                              date_of_death, cause_of_death)
write.csv(roedeer_ids, "roedeer_ids_latex.csv")


#### decriptive figures ####

## heatmap of crossings

#time of day/month heatmap 
library(plyr)
road_intersect_final$date <- format(as.Date(road_intersect_final$ts), format = '%d-%m')
road_intersect_final$stunde <- hour(road_intersect_final$ts)
ddply_table <- road_intersect_final %>% 
  ddply(c("date","stunde"), summarise, N = length(ts)) #how many crossings in that hour on that day

roedeer_telem$date <- format(as.Date(roedeer_telem$ts), format = '%d-%m')
roedeer_telem$stunde <- hour(roedeer_telem$ts)

ddply_table_telem <- roedeer_telem %>% 
  ddply(c("date","stunde"), summarise, total = length(unique(animal_id))) #how many animals are being followed that day
detach("package:plyr", unload = TRUE)  #because dplyr and plyr do not like each other... 

ddply_table <- ddply_table %>% 
  left_join(ddply_table_telem, by = c("date","stunde")) %>% 
  mutate(n = N/total)

dates <- as.factor(road_intersect_final$date)
ddply_background <- data.frame(date = rep(levels(dates), each = 24), 
                               stunde = rep(c(0:23), times = 366))

ddply_background %>% 
  dplyr::left_join(ddply_table, by = c("date", "stunde")) %>% 
  dplyr::mutate(n = coalesce(n, 0)) %>% 
  ggplot() + 
  geom_raster(aes(x=as.Date(date, format='%d-%m'), y=stunde, fill=n)) + 
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal()+
  theme(text = element_text(size = 25), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_breaks = "1 month", 
               date_labels =  c("Dec", "Jan","Feb","Mar","April","May","Jun","Jul","Aug","Sep","Okt","Nov","Dec", "Jan")) +
  scale_y_continuous(breaks = seq(0, 23, by = 2)) +
  labs(x = "Time of year (month)", y = "Time of day (hour)", fill = "Mean number \nof road crossings")


ggsave("Figures/crossings_heatmap.pdf", width = 25, height = 20, units = "cm")



##decided to not use this graph as it is not corrected for the number of roe deer actually collared and it is not always the same amount so it does not give a good image of the distribution. rather rely on the analysis. 
# # crossings and collisions comparison per day
# crossings_day <- road_intersect_final %>% 
#   st_drop_geometry() %>% 
#   #na.omit() %>% 
#   group_by(yday(ts), year(ts)) %>% 
#   count() %>% 
#   group_by(`yday(ts)`) %>% 
#   summarise(mean = mean(n), 
#             sd = sd(n)) %>% 
#   rename(day = `yday(ts)`)
# 
# collisions_day <- rehe_WU_HE %>% 
#   st_drop_geometry() %>% 
#   na.omit() %>% 
#   group_by(yday(date_time), jahr) %>% 
#   count() %>% 
#   group_by(`yday(date_time)`) %>% 
#   summarise(mean = mean(n), 
#             sd = sd(n)) %>% 
#   rename(day = `yday(date_time)`)
# 
# 
# ggplot(data = crossings_day, aes(x = day, y = mean)) +
#   geom_line(col = "darkblue") +
#   geom_ribbon(data = crossings_day, aes(ymin = mean - sd, ymax = mean +sd), alpha = 0.1, fill = "darkblue") +
#   
#   geom_line(data = collisions_day, mapping = aes(x = day, y = mean), col = "darkred") +
#   geom_ribbon(data = collisions_day, mapping = aes(ymin = mean - sd, ymax = mean +sd), alpha = 0.1, fill = "darkred") +
#   
#   theme_minimal()+
#   theme(text = element_text(size = 25)) +
#   scale_x_continuous(breaks = c(1,32,60,91,121,152,182,213,244,274,305,335), 
#                      labels = c("Jan","Feb","Mar","April","May","Jun","Jul","Aug","Sep","Okt","Nov","Dec")) +
#   annotate("text", x = 60, y = 80, label = "Gestation", size = 7) +
#   geom_vline(xintercept = 121, linetype = "longdash") +
#   annotate("text", x = 155, y = 80, label = "Lactation", size = 7) +
#   geom_vline(xintercept = 196, linetype = "longdash") +
#   annotate("text", x = 210, y = 80, label = "Rut", size = 7) +
#   geom_vline(xintercept = 227, linetype = "longdash") +
#   annotate("text", x = 300, y = 80, label = "Diapause", size = 7) +
#   labs(y = "Mean number per day", x = "Month of the year")
# 
# ggsave("Figures/comp_per_day.pdf")



# ####not output, but gives a comparison between the sexes and their crosing behaviour
#same reason here to not use
# road_intersect_final %>% 
#   st_drop_geometry() %>% 
#   left_join(dplyr::select(roedeer_ids, animal_id, sex), by = c("animal_id" = "animal_id")) %>% 
#   group_by(sex, yday(ts), year(ts)) %>% 
#   count() %>% 
#   group_by(sex, `yday(ts)`) %>% 
#   summarise(mean = mean(n), 
#             sd = sd(n)) %>% 
#   rename(day = `yday(ts)`) %>% 
# ggplot(aes(x = day, y = mean)) +
#   geom_line(aes(col = sex)) +
#   geom_ribbon(aes (ymin = mean - sd, ymax = mean +sd, fill = sex), alpha = 0.1) +
#   theme_minimal()+
#   theme(text = element_text(size = 25)) +
#   scale_x_continuous(breaks = c(1,32,60,91,121,152,182,213,244,274,305,335), 
#                      labels = c("Jan","Feb","Mar","April","May","Jun","Jul","Aug","Sep","Okt","Nov","Dec")) +
#   annotate("text", x = 60, y = 70, label = "Gestation", size = 7) +
#   geom_vline(xintercept = 121, linetype = "longdash") +
#   annotate("text", x = 155, y = 70, label = "Lactation", size = 7) +
#   geom_vline(xintercept = 196, linetype = "longdash") +
#   annotate("text", x = 210, y = 70, label = "Rut", size = 7) +
#   geom_vline(xintercept = 227, linetype = "longdash") +
#   annotate("text", x = 300, y = 70, label = "Diapause", size = 7) +
#   labs(y = "Mean number of WVC per day", x = "Month of the year")





## overview map of bawü with collar regions

germany_sf <- read_sf("SDM/gadm36_DEU_shp/gadm36_DEU_1.shp")
germany_sf <- st_transform(germany_sf, crs = projection)#change CRS to fit WU data

germany <- tm_shape(germany_sf) +
  tm_borders()


hessen_sf <- filter(germany_sf, NAME_1 == "Hessen")
hessen_sf$VARNAME_1 <- "HE"
bawu_sf <- filter(germany_sf, NAME_1 == "Baden-Württemberg")
bawu_sf$VARNAME_1 <- "BW"

#coordinates from one telem id of that area
coords <- data.frame(x_coord = c(427350.28, 432201.76, 421519.86, 477789.35, 418883.02), 
                     y_coord = c(5391062.77, 5398288.23, 5379945.18, 5301569.66, 5383005.76),
                     locations = levels(roedeer_ids$area))

coords_sf <- sf::st_as_sf(coords, coords = c("x_coord", "y_coord"), crs = projection)


#bawü overview map 
insetmap <- germany + 
  tm_shape(hessen_sf) +
  tm_borders(lwd = 1) +
  tm_text("VARNAME_1", size = 0.8, bg.color = "white") +
  tm_shape(bawu_sf) +
  tm_borders(lwd = 3, col = "black") +
  tm_text("VARNAME_1", size = 0.8, bg.color = "white") 

bawu <- tm_shape(bawu_sf) +
  tm_borders()

(mainmap <- bawu +
    tm_shape(roads_BW) +
    tm_lines(col = "wdm_name",
             palette=c("#FB8072", "#BEBADA", "#80B1D3", "#FDB462", "#B3DE69"), title.col = "Road category") +
    tm_legend(position=c("right", "top"), frame = FALSE, legend.title.size = 1.5, legend.text.size = 1) +
    
    tm_shape(coords_sf) +
    tm_symbols(size = 0.5) +
    tm_text("locations", size = 1, col = "black", fontface = "bold", 
            auto.placement=F, xmod = c(2.5, 2.5, 2, 2, 2.5), ymod = c(0, 0.3, -0.2, 0, 0.1)) +#, bg.color = "white"
    
    tm_layout(inner.margins = c(0, 0, 0, 0.3), frame = FALSE) +
    #extras: nordpfeil & scale bar
    tm_scale_bar(position = c("left", "bottom"), width = 0.15, text.size = 1) +
    tm_compass(position = c("left", "top"), size = 2))


#add insetmap to layout
vp = viewport(0.8, 0.27, width = 0.5, height = 0.5)
print(insetmap, vp = vp)


tmap_save(mainmap,filename="Figures/overview_map_BW.pdf",
          dpi=100, insets_tm=insetmap, insets_vp=vp)



####crossing analysis preparation#####

## add risk value to crossings

#road_intersect_final <- st_read("Roe Deer Telemetry Data/road_intersections_final.shp")
#st_crs(road_intersect_final) <- projection

road_intersect_occ <- road_intersect_final %>% 
  st_coordinates() %>% 
  as.data.frame()

#two points fall within now raster: check which ones
extract <- as.data.frame(raster::extract(raster_BW_stack_complete, road_intersect_occ))
index <- stats::complete.cases(extract)
summary(index)
road_intersect_final[!index, ] #15524 and 16216 with no data for any predictor, because just between two pixels. fine for now, as it is only two
# plot(crop(raster_BW_stack_complete[[22]], c(431500, 432500, 5398000, 5398500) ))
# plot(roads_BW, add = TRUE)
# plot(st_geometry(road_intersect_final[16216,]), add = TRUE)
road_intersect_final <- road_intersect_final[-c(15524, 16216), ]

road_intersect_occ <- road_intersect_final %>% 
  st_coordinates() %>% 
  as.data.frame()


# prepare an SWD object for BW with only presence locations. 
crossings_SWD <- prepareSWD(species = "WVC", p = road_intersect_occ, env = raster_BW_stack_complete, categorical = "RoadCategory")

#predict risk with best all seasons model
crossings_risk <- predict(model_allSeasons_final, data = crossings_SWD, type = "cloglog")
hist(crossings_risk, breaks = 20)


#add risk column to data
road_intersect_final$risk <- crossings_risk
road_intersect_final$risk_group <- round(road_intersect_final$risk, 1)

summary(filter(road_intersect_final, case == "WVC")$risk)


summary(filter(road_intersect_final, case == "crossing")$risk)



## distribution map with risk prob

#creating a distribution map for HR of roe deer
predict_allSeasons <- predict(model_allSeasons_final, data = raster_BW_stack_complete,
                              type = "cloglog", extent = extent(418000, 433000, 5379000, 5398500), progress = "text") 

map_allSeasons <- plotPred(predict_allSeasons, lt = "Probability of \nDVC occurrence",colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))


WVC_crossings_sf <- road_intersect_final %>% 
  filter(case == "WVC")%>% 
  rename(Season = season) 

WVC_crossings <- WVC_crossings_sf %>% 
  mutate(X = st_coordinates(.)[,1], 
         Y = st_coordinates(.)[,2]) %>%  
  as.data.frame()

WVC_crossings$Season <- factor(WVC_crossings$Season, levels = c("Gestation", "Lactation", "Rut", "Diapause"))


pdf(file="Figures/PredMap_BWwithWVC.pdf")
map_allSeasons +
  geom_polygon(roedeer_mcp[which(roedeer_ids$area != "Stetten"),], 
               mapping = aes(x = long, y = lat, group = id), alpha = 0.05) +
  annotation_scale(width_hint = 0.1) +
  # geom_point(data = road_intersect_occ[which(road_intersect_final$area != "Stetten"), ], 
  #            aes (x = X, y = Y), shape = 1, alpha = 0.05, size = 0.6)  +
  geom_point(data = WVC_crossings[-c(1,2),], aes (x = X, y = Y), shape = 8, col = "black") #WVC points -ids 29, 47 in Stetten

dev.off()


#plot for each individual the road crossing behaviour: only for those that actually have crossed roads
roedeer_ids_vector <- roedeer_ids %>% 
  filter(no_crossings > 5)

roedeer_ids_vector <- unique(roedeer_ids_vector$animal_id)

for (i in roedeer_ids_vector) {
  
  extent_id <- extent(roedeer_mcp[which(roedeer_ids$animal_id == i),])
  
  
  if(nrow(filter(WVC_crossings, animal_id == i)) > 0){
    WVC <- filter(WVC_crossings, animal_id == i)
  } else {
    WVC <- NULL
  }
  
  plot <- roads_BW %>% 
    st_crop(extent_id) %>% 
    ggplot() +
    geom_sf()
  
  data <- samplepoints_sf_ids_analysis %>%
    filter(id == i)  %>%
    st_centroid() %>%
    mutate(freq = no.crossings/days,
           X = st_coordinates(.)[,1],
           Y = st_coordinates(.)[,2]) %>%
    as.data.frame()
  
  area <- as.character(data$area.y[1])
  sex <- ifelse(data$sex[1] == "f", "female", "male")
  
  plot <- plot + data %>%
    geom_point(mapping = aes(x = X, y = Y, size = freq, col = risk), fill = "black") +
    facet_wrap(~Season) +
    scale_size(range = c(1.5, 6), name = "Crossings/day") + 
    scale_color_gradientn(colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"), 
                          limits = c(0,1), name = "Probability of \nDVC occurrence") +
    geom_polygon(roedeer_mcp[which(roedeer_ids$animal_id == i),],
                 mapping = aes(x = long, y = lat, group = id), alpha = 0.1) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), text = element_text(size=15)) +
    annotation_scale(width_hint = 0.1)  +
    labs(title = paste("ID", i, "- Sex:", sex, "- Area:", area), x = "", y = "")
  
  
  if(!is.null(WVC)){
    plot <- plot +
      geom_point(data = WVC, aes(x = X, y = Y), col = "black", shape = 8, size = 5, 
                 inherit.aes = FALSE) 
  }
  
  
  path <- file.path(paste("Figures/PredMaps_ID", i,".pdf", sep = ""))
  pdf(file = path)
  print(plot)
  dev.off()
  
  
}



## Analysis of risk####

## data preparation samplepoints_sf_ids

#inforation for rut 
samplepoints_sf_rut_list <- list()
for (i in roedeer_ids_vector) {
  data <- samplepoints_sf %>% 
    dplyr::select(pointID, risk.rut) %>% 
    st_join(filter(roedeer_mcp_sf, id == i)) %>% 
    dplyr::select(pointID, risk.rut, id, area, diff_days_rut) %>% 
    filter(!is.na(id)) %>% 
    rename(days = diff_days_rut, 
           risk = risk.rut)
  
  #data$no.crossings <- NA
  data$no.crossings <- 
    lengths(st_intersects(data, 
                          road_intersect_final[which(road_intersect_final$animal_id == i & road_intersect_final$season == "Rut"), ]))
  samplepoints_sf_rut_list[[i]] <- data
  print(i)
}

samplepoints_sf_rut_ids <- do.call(rbind, samplepoints_sf_rut_list)
samplepoints_sf_rut_ids$Season <- "Rut"
samplepoints_sf_rut_ids <- filter(samplepoints_sf_rut_ids, days != 0)#if number of days is zero delte them 


#inforation for lactation
samplepoints_sf_lact_list <- list()
for (i in roedeer_ids_vector) {
  data <- samplepoints_sf %>% 
    dplyr::select(pointID, risk.lact) %>% 
    st_join(filter(roedeer_mcp_sf, id == i)) %>% 
    dplyr::select(pointID, risk.lact, id, area, diff_days_lact) %>% 
    filter(!is.na(id)) %>% 
    rename(days = diff_days_lact, 
           risk = risk.lact)
  
  #data$no.crossings <- NA
  data$no.crossings <- 
    lengths(st_intersects(data, 
                          road_intersect_final[which(road_intersect_final$animal_id == i & road_intersect_final$season == "Lactation"), ]))
  samplepoints_sf_lact_list[[i]] <- data
  print(i)
}

samplepoints_sf_lact_ids <- do.call(rbind, samplepoints_sf_lact_list)
samplepoints_sf_lact_ids$Season <- "Lactation"
samplepoints_sf_lact_ids <- filter(samplepoints_sf_lact_ids, days != 0)#if number of days is zero delte them 


#inforation for gestation
samplepoints_sf_gest_list <- list()
for (i in roedeer_ids_vector) {
  data <- samplepoints_sf %>% 
    dplyr::select(pointID, risk.gest) %>% 
    st_join(filter(roedeer_mcp_sf, id == i)) %>% 
    dplyr::select(pointID, risk.gest, id, area, diff_days_gest) %>% 
    filter(!is.na(id)) %>% 
    rename(days = diff_days_gest, 
           risk = risk.gest)
  
  #data$no.crossings <- NA
  data$no.crossings <- 
    lengths(st_intersects(data, 
                          road_intersect_final[which(road_intersect_final$animal_id == i & road_intersect_final$season == "Gestation"), ]))
  samplepoints_sf_gest_list[[i]] <- data
  print(i)
}

samplepoints_sf_gest_ids <- do.call(rbind, samplepoints_sf_gest_list)
samplepoints_sf_gest_ids$Season <- "Gestation"
samplepoints_sf_gest_ids <- filter(samplepoints_sf_gest_ids, days != 0)#if number of days is zero delte them   


#inforation for diaation
samplepoints_sf_dia_list <- list()
for (i in roedeer_ids_vector) {
  data <- samplepoints_sf %>% 
    dplyr::select(pointID, risk.dia) %>% 
    st_join(filter(roedeer_mcp_sf, id == i)) %>% 
    dplyr::select(pointID, risk.dia, id, area, diff_days_dia) %>% 
    filter(!is.na(id)) %>% 
    rename(days = diff_days_dia, 
           risk = risk.dia)
  
  #data$no.crossings <- NA
  data$no.crossings <- 
    lengths(st_intersects(data, 
                          road_intersect_final[which(road_intersect_final$animal_id == i & road_intersect_final$season == "Diapause"), ]))
  samplepoints_sf_dia_list[[i]] <- data
  print(i)
}

samplepoints_sf_dia_ids <- do.call(rbind, samplepoints_sf_dia_list)
samplepoints_sf_dia_ids$Season <- "Diapause"
samplepoints_sf_dia_ids <- filter(samplepoints_sf_dia_ids, days != 0)#if number of days is zero delete them  

#add them together to dataframe
samplepoints_sf_ids <- rbind(samplepoints_sf_rut_ids, samplepoints_sf_dia_ids, samplepoints_sf_gest_ids, samplepoints_sf_lact_ids)

#join some id information to the table 
samplepoints_sf_ids <- roedeer_ids %>% 
  dplyr::select(area, animal_id, sex) %>% 
  left_join(samplepoints_sf_ids, ., by = c("id" = "animal_id"))
samplepoints_sf_ids$Season <- factor(samplepoints_sf_ids$Season, levels = c("Gestation", "Lactation", "Rut", "Diapause"))



# #visual example for methods 
# samplepoints_ex <- filter(samplepoints_sf_ids, (pointID > 140 & pointID < 150) & Season == "Rut")
# road_intersect_ex <- filter(road_intersect_final, season == "Rut")
# road_intersect_ex$animal_id <- droplevels(road_intersect_ex$animal_id)
# 
# #example image
# tm_shape(samplepoints_ex) +
#  tm_polygons(col = "risk") +
#  tm_shape(road_intersect_ex) +
#  tm_dots(col = "animal_id", size = 2, shape = 21, palette = "Set2") + 
#   tm_shape(roads_BW) +
#   tm_lines()+
#   tm_legend(show = FALSE)


## analysis of the WVC locations

View(WVC_crossings_sf)


WVC_crossings_summary <- data.frame(matrix(data = NA, nrow = nrow(WVC_crossings), ncol = 12))
colnames(WVC_crossings_summary) <- c("animal_id", "area", "wdm_name", "risk", "Season","sex", "days", "no_crossings",  "min_allCrossings", "max_allCrossings",  "mean_allCrossings", "no_samplingbuffers")

for (i in 1:nrow(WVC_crossings)) {
  y <- unique(WVC_crossings$animal_id)[i]
  WVC_crossings_id <- filter(WVC_crossings_sf, animal_id == y)
  
  samplepoints_sf_filter <- filter(samplepoints_sf_ids, id == y & Season == WVC_crossings_id$Season)
  WVC_crossings_join <- WVC_crossings_id %>% 
    st_join(samplepoints_sf_filter, join = st_within)
  
  WVC_crossings_summary[i, 1] <- WVC_crossings_join$animal_id
  WVC_crossings_summary[i, 2] <- as.character(WVC_crossings_join$area)
  WVC_crossings_summary[i, 3] <- as.character(WVC_crossings_join$wdm_name)
  WVC_crossings_summary[i, 4] <- WVC_crossings_join$risk.y
  WVC_crossings_summary[i, 5] <- WVC_crossings_join$Season.x
  WVC_crossings_summary[i, 6] <- as.character(WVC_crossings_join$sex)
  WVC_crossings_summary[i, 7] <- WVC_crossings_join$days
  WVC_crossings_summary[i, 8] <- WVC_crossings_join$freq
  WVC_crossings_summary[i, 9] <- min(samplepoints_sf_filter$freq)
  WVC_crossings_summary[i, 10] <- max(samplepoints_sf_filter$freq)
  WVC_crossings_summary[i, 11] <- mean(samplepoints_sf_filter$freq)
  WVC_crossings_summary[i, 12] <- nrow(samplepoints_sf_filter)
}

WVC_crossings_summary <- WVC_crossings_summary[order(WVC_crossings_summary$animal_id), ]
WVC_crossings_summary <- WVC_crossings_summary %>% 
  mutate_if(is.numeric, round, digits = 2)
write.csv(WVC_crossings_summary, "WVC_crossings_summary.csv")
#info for ID 34, as the points are actually not included because they are outside of her home range
join <- st_join(filter(road_intersect_final, animal_id == 34 & season == "Rut"), samplepoints_sf, join = st_within)

mean(WVC_crossings$risk)
sd(WVC_crossings$risk)


## analysis
#individuals that were only measured for less than 7 days during a prediod should be deleted for analysis -> bias by eg ID 29 during rut
samplepoints_sf_ids_analysis <- filter(samplepoints_sf_ids, days > 7)

#sex could influence more risky behaviour when crossing = interaction
# area not included as area influences the number of roads available and thus biases the number of crossings
# season included as during eg rut, more risky behaviour possible. Probably from males, but avoidig three way interactions, so one more model with all same but without sex*risk. 
# offset of days, as each season is different in length and therefore less possibily for crossing

# zero inflation could be changed by risk: more zeros in low risk zones because if no crossing no risk
# season: in certain seasons, less crossings generally due to ecological processes
# possibly add a mixed effect of individuals nested in area to account for differences: not possible as each individual can only have one sex. so extra analysis loooking at the differences of road crossings in individuals

summary(samplepoints_sf_ids_analysis)

# library(lattice)
# xyplot(no.crossings ~ risk|Season, group = sex, data = samplepoints_sf_ids_analysis )

samplepoints_sf_ids_analysis  %>% 
  st_drop_geometry() %>% 
  mutate(freq = no.crossings/days, #freqeuncy of crossings per day
         risk_group = as.factor(round(risk, 1))) %>% 
  group_by(sex, risk_group, Season) %>% 
  summarize(mean_freq = mean(freq), 
            sd_freq = sd(freq)) %>% 
  ggplot(aes(x = risk_group, y = mean_freq, fill = sex)) + 
  geom_bar(stat = "identity",position = position_dodge()) +
  facet_grid(~Season) + #, scales = "free_y"
  theme_minimal() +
  labs(x = "Modelled DVC risk", y = "Mean number of crossings per day", fill = "Sex") +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_errorbar(aes(x = risk_group, ymin = 0, ymax = mean_freq + sd_freq), position=position_dodge(width=0.9), width = .1) #extramely large sd!
ggsave("Figures/crossings_descriptives_figure.pdf", width = 25, height = 15, units = "cm")


#these are the two extremly high crossings.
#View(filter(samplepoints_sf_ids_analysis , sex == "m" & round(risk, 1) == 0.3 & Season == "Rut")) #males 22, 29, 47, mostly 29 resp
#View(filter(samplepoints_sf_ids_analysis , sex == "m" & round(risk, 1) == 0.4 & Season == "Rut")) #males 22, 29, 47, mostly 29 & 47


View(roedeer_ids) #not every individual has been surveyed in every season, so reason to again not include id as mixed effect 


#we know that males generally cross more often. influenced differently per season than females.

#dorsnt seem like males are signif more risky 
samplepoints_sf_ids_analysis  %>% 
  #filter(no.crossings == 0) %>% 
  st_drop_geometry() %>% 
  group_by(sex) %>% 
  summarize(sum_crossings = sum(no.crossings))

samplepoints_sf_ids_analysis  %>% 
  filter(no.crossings != 0) %>% 
  ggplot(aes(y = risk, x = sex)) + #, fill = area.y
  geom_boxplot() +
  theme_minimal() +
  theme(text = element_text(size = 25))  +
  scale_x_discrete(labels = c("f" = "Females (n = 12796)", "m" = "Males (n = 8158)")) +
  labs(x = "Sex", y = "Modelled DVC risk")
ggsave("Figures/crossings_descriptives_sex.pdf", width = 15, height = 15, units = "cm")

#enough points per sex and season
samplepoints_sf_ids_analysis  %>% 
  st_drop_geometry() %>% 
  count(sex, Season)

#enough poits per area and season
samplepoints_sf_ids_analysis  %>% 
  st_drop_geometry() %>% 
  count(area.y, Season)
#whats with oberbruch? always the same nmber of sampling points within the area in all seasons but should be enough points
View(filter(samplepoints_sf_ids_analysis , area.y == "Oberbruch"))


#males in rut are quite little points available 
samplepoints_sf_ids_analysis  %>% 
  st_drop_geometry() %>% 
  group_by(sex, Season) %>% 
  summarise(no.crossings = sum(no.crossings), 
            sum.days = sum(days))



#analysis with both sexes, but then ID has to be taken out as random effect. do own anova with only ID to show the differences between individuals? change to all seasons, offset per day

#who is the outlier? just a female with the point in the middle of her home range... 
#samplepoints_sf_ids_analysis[2041,]
#road_intersect_ID40 <- filter(road_intersect_final, animal_id == 40 & season == "Gestation")
# intersect <- st_intersects(samplepoints_sf[which(samplepoints_sf$pointID == 109), ], 
#     road_intersect_ID40)
# plot(st_geometry(samplepoints_sf[which(samplepoints_sf$pointID == 109), ]))
# plot(road_intersect_ID40[intersect[[1]], ], add = T)


analysis.ids.uni.1 <- glmmTMB(no.crossings ~ as.factor(id), 
                              zi = ~ 1, 
                              data = samplepoints_sf_ids_analysis, #[-2041,], #outlier of no crossings above 500influencing the dispersion test significantly
                              family = "nbinom1")
summary(analysis.ids.uni.1)

plot(simulateResiduals(analysis.ids.uni.1)) #Hartig: Specifically, if you have a lot of data points, residual diagnostics will nearly inevitably become significant, because having a perfectly fitting model is very unlikely. That, however, doesn’t necessarily mean that you need to change your model. The p-values confirm that there is a deviation from your null hypothesis. It is, however, in your discretion to decide whether this deviation is worth worrying about. If you see a dispersion parameter of 1.01, I would not worry, even if the test is significant. A significant value of 5, however, is clearly a reason to move to a model that accounts for overdispersion.
testDispersion(analysis.ids.uni.1)
plot(residuals(analysis.ids.uni.1)~fitted(analysis.ids.uni.1))

#significant influence of id on no of ccrossings
#car::Anova(analysis.ids.uni.1)
analysis.ids.uni.1.aov <- aov(analysis.ids.uni.1)
summary(analysis.ids.uni.1.aov)
#tukey <- TukeyHSD(analysis.ids.uni.1.aov) #tukey post-hoc test for pairwise comparison
#plot(tukey, las = 1) #too many comparisons to plot it nicely

#intercept only model with random effect ID: how much variation is provided by individual based differences: individual clustering: difference to fixed model that asks is each individual different to the intercept individual. 
samplepoints_sf_ids_analysis$id <- as.factor(samplepoints_sf_ids_analysis$id)
analysis.ids.uni.2 <- glmmTMB(no.crossings ~ 1 + (1|id), 
                              zi = ~ 1, 
                              data = samplepoints_sf_ids_analysis, #[-2041,], #outlier of no crossings above 500influencing the dispersion test significantly
                              family = "nbinom1")
summary(analysis.ids.uni.2)
coef(analysis.ids.uni.2)



samplepoints_sf_ids_analysis$freq <- samplepoints_sf_ids_analysis$no.crossings/samplepoints_sf_ids_analysis$days

ggplot(samplepoints_sf_ids_analysis, aes(x = risk, y = freq, col = Season)) +
  geom_point(size = 0.7) +
  facet_wrap(~id, ncol = 9) +
  scale_color_manual(values = c("#FFED6F", "#FDC086", "#7FC97F", "#BEAED4")) +
  labs(x = "Modelled DVC risk", y = "Number of crossings per day") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 15),
        legend.title= element_text(size = 15),
        legend.text = element_text(size = 10)) 
ggsave("Figures/crossings_ID.pdf", width = 20, height = 26, units = "cm")



# ggplot(samplepoints_sf_ids_analysis, aes(x = id, y = freq)) +
#   geom_boxplot() +
#   theme_minimal()


# analysis looking into the further factors influencing number of crossings, without ID in mixed effects
analysis.ids.glmm <- glmmTMB(no.crossings ~ risk*as.factor(sex)*Season + offset(log(days)) + (1|area.y), 
                             zi = ~ risk + Season + as.factor(sex), 
                             data = samplepoints_sf_ids_analysis, 
                             family = "nbinom1") #complexest

analysis.ids.glmm1 <- glmmTMB(no.crossings ~ risk + as.factor(sex)*Season + offset(log(days)) + (1|area.y), 
                              zi = ~ risk + Season + as.factor(sex), 
                              data = samplepoints_sf_ids_analysis, 
                              family = "nbinom1") #complex model with sex*season
analysis.ids.glmm2 <- glmmTMB(no.crossings ~ risk*as.factor(sex) + Season + offset(log(days)) + (1|area.y), 
                              zi = ~ risk + Season + as.factor(sex), 
                              data = samplepoints_sf_ids_analysis, 
                              family = "nbinom1") #complex model with risk*sex


#summary(analysis.ids.glmm2)
anova(analysis.ids.glmm1, analysis.ids.glmm2) #no model is much better than the other, 2 lower AIC value
# sim <- simulateResiduals(analysis.ids.glmm1)
# plot(sim) #no influence on dispersion by the outlier ID 40. 


#simplifying the count model
analysis.ids.glmm1.1 <- glmmTMB(no.crossings ~ risk + as.factor(sex) + Season + offset(log(days)) + (1|area.y), zi = ~ risk + Season + as.factor(sex), data = samplepoints_sf_ids_analysis, family = "nbinom1") #simplified without interaction
analysis.ids.glmm1.2 <- glmmTMB(no.crossings ~ risk + as.factor(sex) + offset(log(days)) + (1|area.y), 
                                zi = ~ risk + Season + as.factor(sex), data = samplepoints_sf_ids_analysis, family = "nbinom1") #simplified without Season in count model
analysis.ids.glmm1.3 <- glmmTMB(no.crossings ~ as.factor(sex)*Season + offset(log(days)) + (1|area.y), 
                                zi = ~ risk + Season + as.factor(sex), 
                                data = samplepoints_sf_ids_analysis, 
                                family = "nbinom1") #simplified without risk in count model, but interaction
analysis.ids.glmm1.4 <- glmmTMB(no.crossings ~ as.factor(sex) + Season + offset(log(days)) + (1|area.y), 
                                zi = ~ risk + Season + as.factor(sex), 
                                data = samplepoints_sf_ids_analysis, 
                                family = "nbinom1") #simplified without risk in count model, without interaction


anova(analysis.ids.glmm1, analysis.ids.glmm1.1) # 1 is significantly better, so interaction better thn additive
anova(analysis.ids.glmm1, analysis.ids.glmm1.2) # 1 significantly better, than 1.2, so keep season in
anova(analysis.ids.glmm1, analysis.ids.glmm1.3) #1 is significantly better, so keep risk in
anova(analysis.ids.glmm1, analysis.ids.glmm1.4) #1 sig better


analysis.ids.glmm2.2 <- glmmTMB(no.crossings ~ risk*as.factor(sex) + offset(log(days)) + (1|area.y), 
                                zi = ~ risk + Season + as.factor(sex), 
                                data = samplepoints_sf_ids_analysis, 
                                family = "nbinom1") #simplified without risk in count model, but interaction





#going on with 1, now adjusting the zi part
analysis.ids.glmm1.0.1 <- glmmTMB(no.crossings ~ risk + as.factor(sex)*Season + offset(log(days)) + (1|area.y),
                                  zi = ~ Season + as.factor(sex), 
                                  data = samplepoints_sf_ids_analysis, 
                                  family = "nbinom1") #simplified without risk in zi model
analysis.ids.glmm1.0.2 <- glmmTMB(no.crossings ~ risk + as.factor(sex)*Season + offset(log(days)) + (1|area.y),
                                  zi = ~ risk + as.factor(sex), 
                                  data = samplepoints_sf_ids_analysis, 
                                  family = "nbinom1") #simplified without season in zi model
analysis.ids.glmm1.0.3 <- glmmTMB(no.crossings ~ risk + as.factor(sex)*Season + offset(log(days)) + (1|area.y),
                                  zi = ~ risk + Season, 
                                  data = samplepoints_sf_ids_analysis, 
                                  family = "nbinom1") #simplified without risk in zi model


anova(analysis.ids.glmm1, analysis.ids.glmm1.0.1) #1.0.1 better than 1
anova(analysis.ids.glmm1, analysis.ids.glmm1.0.2) #1 better
anova(analysis.ids.glmm1, analysis.ids.glmm1.0.3) #1 better


analysis.ids.glmm1.1.1 <- glmmTMB(no.crossings ~ risk + as.factor(sex) + Season + offset(log(days)) + (1|area.y), zi = ~ risk + Season, data = samplepoints_sf_ids_analysis, family = "nbinom1") #without sex in zi model
analysis.ids.glmm1.1.2 <- glmmTMB(no.crossings ~ risk + as.factor(sex) + Season + offset(log(days)) + (1|area.y), zi = ~ risk + as.factor(sex), data = samplepoints_sf_ids_analysis, family = "nbinom1") #without season in zi model
analysis.ids.glmm1.1.3 <- glmmTMB(no.crossings ~ risk + as.factor(sex) + Season + offset(log(days)) + (1|area.y), zi = ~ Season + as.factor(sex), data = samplepoints_sf_ids_analysis, family = "nbinom1") #without risk in zi model


aic_table <- AIC(analysis.ids.glmm1, analysis.ids.glmm1.1, analysis.ids.glmm1.2, analysis.ids.glmm1.3, analysis.ids.glmm1.4, analysis.ids.glmm1.0.1, analysis.ids.glmm1.0.2, analysis.ids.glmm1.0.3, analysis.ids.glmm2, analysis.ids.glmm2.2, analysis.ids.glmm1.1.1, analysis.ids.glmm1.1.2, analysis.ids.glmm1.1.3)
rownames(aic_table) <- c("count: Sex*Season + Risk, zi: Season + Risk + Sex", 
                         "count: Risk + Season + Sex, zi: Season + Risk + Sex", 
                         "count: Risk + Sex, zi: Season + Risk + Sex", 
                         "count: Sex*Season, zi: Season + Risk + Sex", 
                         "count: Sex + Season, zi: Season + Risk + Sex", 
                         "count: Sex*Season + Risk, zi: Season + Sex",  
                         "count: Sex*Season + Risk, zi: Risk + Sex",  
                         "count: Sex*Season + Risk, zi: Season + Risk", 
                         "count: Risk*Sex + Season, zi: Season + Risk + Sex",
                         "count: Risk*Sex, zi: Season + Risk + Sex", 
                         "count: Risk + Season + Sex, zi: Season + Risk",
                         "count: Risk + Season + Sex, zi: Risk + Sex",
                         "count: Risk + Season + Sex, zi: Season + Sex")


aic_table <- aic_table[order(aic_table$AIC), ]
#save for latex
aic_table_ltx <- xtable(aic_table, caption = "AIC comparison of all models with number of crossings as response variable.", label = "aic_table")
print(aic_table_ltx,file="Figures/aic_table_ltx.tex",table.placement = "h", include.rownames = TRUE,
      caption.placement="bottom")

rm(analysis.ids.glmm1, analysis.ids.glmm2, analysis.ids.glmm1.1, analysis.ids.glmm1.2, analysis.ids.glmm1.3, analysis.ids.glmm1.4, analysis.ids.glmm1.0.1, analysis.ids.glmm1.0.2, analysis.ids.glmm1.0.3, analysis.ids.glmm, analysis.ids.glmm1.1.1, analysis.ids.glmm1.1.2, analysis.ids.glmm1.1.3, analysis.ids.glmm2.2, analysis.ids.glmm2.4, aic_table, aic_table_ltx)
#so analysis.ids.glmm1.0.1 best

analysis.ids.glmm.final <- analysis.ids.glmm1.1.3
#tried witout influencial point [-2041,]: but little difference in risk significance

summary(analysis.ids.glmm.final) #count model log link, zi model logit link
coef(analysis.ids.glmm.final)
car::Anova(analysis.ids.glmm.final)


plot(fitted(analysis.ids.glmm.final), residuals(analysis.ids.glmm.final)) #residulas are not very nice. 
sim <- simulateResiduals(analysis.ids.glmm.final)
pdf("Figures/DHARMa_glmm.final.pdf", width = 7, height = 5) 
plot(sim) #KS test significant: however qq plot nearly linear, overall distribution roughly OK. re~fitted also.
dev.off()


#count model coefficients
exp(-2.87069) #0.06 intercept number of crossings
exp(0.41596) # +1.515825  per unit increased risk (sig)
exp(0.21853) # + 1.244 for males (sig)
exp(0.27760) # + 1.32 during lactation (.)

exp(-2.87 + 0.5*0.42) #0.06994822 calc for no of crossings for females during gestation at 0.5 risk
exp(-2.87 + 0.5*0.42 + 1*0.228) #0.08786093 calc for no of crossings for females during lactation at 0.5 risk

risk_seq = rep(seq(0, 1, len = 100))
effect_plot_df <- data.frame(risk = rep(seq(0, 1, len = 100), times = 8), 
                             sex = rep(c("m", "f"), each = 400),   
                             Season = rep(rep(levels(samplepoints_sf_ids_analysis$Season), each = 100), times = 2), 
                             days = 1, 
                             area.y = "Stetten")
effect_plot_df$Season <- factor(effect_plot_df$Season, levels = c("Gestation", "Lactation", "Rut", "Diapause"))

predict <- predict(analysis.ids.glmm.final, newdata = effect_plot_df, type = "response", se.fit = TRUE)

ggplot(effect_plot_df, aes(x = risk, y = predict$fit, color = Season)) + #linetype = sex, 
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = predict$fit - 1.98*predict$se.fit, ymax = predict$fit + 1.98*predict$se.fit, fill = Season, color = NULL), alpha = 0.1) + #se times 1.98 is confidence interval. 
  facet_grid(~ sex, labeller = as_labeller(c(
    `f` = "Females",
    `m` = "Males"))) +
  theme_minimal() +
  theme(text = element_text(size = 25), axis.text.x = element_text(angle = 45, hjust = 1))  +
  #scale_linetype_manual(values = c("twodash", "solid"), name = "Sex") +
  #ylim(0, 30) +
  scale_color_manual(values=c("#FFED6F", "#FDC086", "#7FC97F", "#BEAED4")) +
  scale_fill_manual(values=c("#FFED6F", "#FDC086", "#7FC97F", "#BEAED4")) +
  labs(x = "Modelled DVC risk", y = "Predicted number of crossings per day")

ggsave("Figures/crossings_analysis_pred.pdf", width = 25, height = 20, units = "cm")


#The zero-inflation model estimates the probability of an extra zero such that a positive contrast indicates a higher chance of absence
#zi model coefficients as odds ratio: one unit increase the odds ration changes by ...
exp( -1.1359) # 0.321133 intercept of prob of zeros occuring
exp(1.3480) # + 3.849718 during lact (sig)
exp(-0.5953) # + 0.5513971 for males (sig) ??? less zeros in males than in females... 
# zeros <- filter(samplepoints_sf_ids_analysis, no.crossings == 0)
# summary(as.factor(zeros$sex))
#    f    m 
# 1318  656 


# library(texreg)
# analysis.ids.glmm.final_ltx <- texreg(analysis.ids.glmm.final, label = "table: analysis.ids.glmm.final",
#         caption = "Model output for the final negative-binomial zero inflated GLMM with number of crossings as response variable. Count model coefficients are from a negative binomial distribution with log link, while Zero-inflation model coefficients are from a binomial distribution with logit link.",
#         single.row = TRUE,
#        scriptsize = TRUE, #smaller font size
#       #sideways = TRUE, #rotates the table by 90 degrees
#         custom.names = c(),  #customize variable names
#        custom.model.names = c(""), include.aic = TRUE, include.bic = FALSE,include.aicc = FALSE, return.string = TRUE)
# print(analysis.ids.glmm.final_ltx, file="figures/analysis.ids.glmm.final_ltx.tex")
