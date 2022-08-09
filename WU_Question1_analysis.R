

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

setwd("W:/Masterarbeit Daten/Analyse")


#### SDM modelling prep ####
# add raster layers if not in environment
setwd("W:/Masterarbeit Daten/Analyse/Predictor Rasters/HE")

raslist <- list.files(getwd(), pattern=".tiff$", full.names=FALSE)

for(r in raslist){
  assign(substr(r, 1, nchar(r)-5), raster(r))
  print(r)
}


setwd("W:/Masterarbeit Daten/Analyse")

#make a raster stack: scaled variables! 
raster_HE_stack_complete <- stack(landcover_HE_1_scaled, landcover_HE_2_scaled, landcover_HE_3_scaled, landcover_HE_4_scaled, 
                                  landcover_HE_5_scaled, landcover_HE_6_scaled, landcover_HE_7_scaled, landcover_HE_8_scaled,
                                  landcover_HE_9_scaled, landcover_HE_1_small_scaled, landcover_HE_2_small_scaled, landcover_HE_3_small_scaled,
                                  landcover_HE_4_small_scaled, landcover_HE_5_small_scaled, landcover_HE_6_small_scaled,
                                  landcover_HE_7_small_scaled, landcover_HE_8_small_scaled, landcover_HE_9_small_scaled,
                                  roads_density_HE_scaled, roads_density_HE_100m_scaled, roads_HE_raster_brf_scaled, roads_HE_raster_wdm)

names(raster_HE_stack_complete) <- c("Meadows1_200m", "Swamps2_200m", "Industry3_200m", "Urbanareas4_200m", "Complex habitats5_200m", 
                                     "ConiferForests6_200m", "MixedForests7_200m", "BroadleafedForests8_200m", "ArableAreas9_200m", 
                                     "Meadows1_100m", "Swamps2_100m", "Industry3_100m", "Urbanareas4_100m", "Complex habitats5_100m", 
                                     "ConiferForests6_100m", "MixedForests7_100m", "BroadleafedForests8_100m", "ArableAreas9_100m", 
                                     "RoadDensity_200m", "RoadDensity_100m","RoadWidth","RoadCategory")


# load unscaled layers
setwd("W:/Masterarbeit Daten/Analyse/Layers/unscaled predictors/HE")
raslist <- list.files(getwd(), pattern=".tiff$", full.names=FALSE)

for(r in raslist){
  assign(substr(r, 1, nchar(r)-5), raster(r))
  print(r)
}


setwd("W:/Masterarbeit Daten/Analyse")


raster_HE_stack_complete_unscaled <- stack(landcover_HE_1_masked, landcover_HE_2_masked, landcover_HE_3_masked, landcover_HE_4_masked, 
                                           landcover_HE_5_masked, landcover_HE_6_masked, landcover_HE_7_masked, landcover_HE_8_masked,
                                           landcover_HE_9_masked, landcover_HE_1_small_masked, landcover_HE_2_small_masked, landcover_HE_3_small_masked,
                                           landcover_HE_4_small_masked, landcover_HE_5_small_masked, landcover_HE_6_small_masked,
                                           landcover_HE_7_small_masked, landcover_HE_8_small_masked, landcover_HE_9_small_masked,
                                           roads_density_HE_masked, roads_density_HE_masked_100m, roads_HE_raster_brf, roads_HE_raster_wdm)
names(raster_HE_stack_complete_unscaled) <- c("Meadows1_200m", "Swamps2_200m", "Industry3_200m", "Urbanareas4_200m", "Complex habitats5_200m", 
                                              "ConiferForests6_200m", "MixedForests7_200m", "BroadleafedForests8_200m", "ArableAreas9_200m", 
                                              "Meadows1_100m", "Swamps2_100m", "Industry3_100m", "Urbanareas4_100m", "Complex habitats5_100m", 
                                              "ConiferForests6_100m", "MixedForests7_100m", "BroadleafedForests8_100m", "ArableAreas9_100m", 
                                              "RoadDensity_200m", "RoadDensity_100m","RoadWidth","RoadCategory")




### prepare presence (WVC) and background locations ###

#select columns that contain the coordinnates of the locations
projection <- "+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"

rehe_WU_HE <- st_read("Layers/rehe_WU_HE.shp")
rehe_WU_HE <- st_transform(rehe_WU_HE, projection)

#extract coordinates in matrix
rehe_WU_HE_occ <- rehe_WU_HE %>% 
  st_coordinates() %>% 
  as.data.frame()


# construct background locations
occur.ras <- rasterize(rehe_WU_HE_occ, raster_HE_stack_complete, 1) #fill cells with data with 1

presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

rehe_WU_HE_dens <- kde2d(x = pres.locs[,1], y = pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras))) 
rehe_WU_HE_dens_raster <- raster(rehe_WU_HE_dens)
plot(rehe_WU_HE_dens_raster)


# default raster to use for all hessen covariates
r.raster <- raster()
extent(r.raster) <- extent(412000, 587000, 5471000,  5723500) #xmin, xmax, ymin, ymax -> extent of Hessen and a bit outside
res(r.raster) <- 50
crs(r.raster) <- projection


crs(rehe_WU_HE_dens_raster) <- projection
bias_file <- projectRaster(rehe_WU_HE_dens_raster, r.raster)

bias_file_masked <- mask(bias_file, raster_HE_stack_complete[[22]]) #mask the bias file with wdm raster as all others too

# sample 10000 background points with prob = bias file
bg_HE <- xyFromCell(bias_file_masked, sample(which(!is.na(values(bias_file_masked))), 10000, 
                                             prob=values(bias_file_masked)[!is.na(values(bias_file_masked))])) 

bg_HE_df <- bg_HE %>% 
  as.data.frame() %>% 
  rename(X = x, 
         Y = y)
#bg_sf <- st_as_sf(bg_HE_df, coords = c("X", "Y"), crs = projection)



#### SDMtune models ####

# All seasons model

# prepare an SWD object
data_all <- prepareSWD(species = "WVC", p = rehe_WU_HE_occ, a = bg_HE_df, env = raster_HE_stack_complete, categorical = "RoadCategory")

#explore SWD object and export (if necessary)
data_all@species #ect


#train first model #

# using cross validation instead of training and testing datasets: generally better 
# makes a SDMmodelCV object that hosts all the models trained during the cross validation
# k = 5 seems agreeable as the dataset is very large. 
folds = randomFolds(data_all, k = 5, only_presence = TRUE)
default_model_crossval <- train(method = "Maxent", data = data_all, fc = "l", folds = folds)

#variable importance is base on the way how the variables are ordered, but permutation is updated between the models 
vi_cross <- varImp(default_model_crossval, permut = 5)

plotVarImp(vi_cross) #only permutation importance (actually the interesting one)



# variable selection
# Prepare background locations to test autocorrelation
bg_extra_coords <- dismo::randomPoints(raster_HE_stack_complete, 10000)
bg_extra <- prepareSWD(species = "WVC", a = bg_extra_coords,
                       env = raster_HE_stack_complete, categorical = "RoadCategory")


pdf(file="Figures/plot_cor.pdf")
plotCor(bg_extra, method = "spearman", cor_th = 0.7)
dev.off()

corVar(bg_extra, method = "spearman", cor_th = 0.7) #print pairs of correlated variables


# Remove variables with correlation higher than 0.7 accounting for the AUC,
# in the following example the variable importance is computed as permutation
# importance

# vs is the model after variableselection
model_crossval_vs <- varSel(default_model_crossval, metric = "auc", bg4cor = bg_extra, cor_th = 0.7,
                            permut = 10) #takes about 45 minutes (doesnt really change after 3 permutations. )

#Removed variables: Meadows1_200m, Industry3_200m, Urbanareas4_200m, ArableAreas9_200m, 
#Urbanareas4_100m, Complex.habitats5_100m, ConiferForests6_100m, MixedForests7_100m, 
#BroadleafedForests8_100m, RoadDensity_100m
vi_cross_vs <- varImp(model_crossval_vs)


# optimize hyperparameters

# Define the hyperparameters to test
getTunableArgs(model_crossval_vs)
#reg: amount of regularization, controls overfitting by shrinking some parameters towards zero which penalizes model compelxity
#feature classes: linear (L), quadratic (Q), product (P), threshold (T), hinge (H), and category indicator (C) 
#-> more presence points >80 lead to all feature classes being able to be used. merow2013
#number of iterations: default 500

args <- list(reg = seq(0.2, 5, 0.2),
             fc = c("l", "lq", "lh", "lp", "lqp", "lqph")) 

crossval_vs_optimized <- optimizeModel(model_crossval_vs, hypers = args, metric = "auc", seed = 789) #takes 1:52 hours

crossval_vs_optimized@results
model_crossval_vs_optimized <- crossval_vs_optimized@models[[1]] #best model collection after optimizing process

#with cross validation: 
# fc: lqph
# reg: 0.4
# iter: 500

plot(crossval_vs_optimized) #plot the optmization process

vs_optim <- varImp(model_crossval_vs_optimized) #the output is the average of the variable importance of each model trained during the cross validation
plotVarImp(vs_optim)


# reduce variation (optimize model parsimony) 
# discards variables with low contribution

# reduce model, with permutation importance lower than 2% only of removing the variables 
# the model performance does not decrease according to the AUC
model_allSeasons_final <- reduceVar(model_crossval_vs_optimized, th = 2, 
                                    metric = "auc", permut = 10, use_jk = TRUE)
#Removed variables: Swamps2_200m, Swamps2_100m
varImp_allSeasons_final <- varImp(model_allSeasons_final)
plotVarImp(varImp_allSeasons_final)


ggplot(varImp_allSeasons_final, aes(x = reorder(Variable,Permutation_importance), y = Permutation_importance)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Permutation_importance -sd, ymax = Permutation_importance +sd), width = .2) +
  theme_minimal()+
  theme(text = element_text(size = 25))  +
  labs(y = "Permutation importance in %", x = "") +
  scale_x_discrete(labels = c("BroadleafedForests8_200m" = "Broadleafed forests (w/200m)", 
                              "RoadDensity_200m" = "Road density (w/200m)", 
                              "RoadCategory" = "Road category", 
                              "ArableAreas9_100m" = "Arable Areas (w/100m)", 
                              "Meadows1_100m" = "Meadows (w/100m)", 
                              "RoadWidth" = "Roadwidth", 
                              "ConiferForests6_200m" = "Conifer forests (w/200m)", 
                              "MixedForests7_200m" = "Mixed forests (w/200m)", 
                              "Complex.habitats5_200m" = "Complex habitats (w/200m)", 
                              "Industry3_100m" = "Industry (w/100m)")) +
  coord_flip()
ggsave("Figures/VarImp_varImp_allSeasons_final.pdf")




### evaluation statistics ###
auc(model_allSeasons_final@models[[1]]) #the cv data, that makes an average AUC out of the 5 models

# calculate mean and sd for all models 
# make tabe to fill with all seasons
model_performance <- data.frame(matrix(NA, nrow = 5, ncol = 6))
colnames(model_performance) <- c("AUC_mean", "AUC_sd", "TSS_mean", "TSS_sd", "COR_mean", "COR_sd")
rownames(model_performance) <- c("All","Gestation", "Lactation", "Rut", "Diapause")


auc_list <- c()
for(i in 1:5){
  auc_list[i] <- auc(model_allSeasons_final@models[[i]])
}
model_performance[1,1] <- mean(auc_list) #0.7416079
model_performance[1,2] <- sd(auc_list) # 0.0003050889

#TSS: true skill statistics (Allouche, tsoar and kadmon 20061)

tss_list <- c()
for(i in 1:5){
  tss_list[i] <- tss(model_allSeasons_final@models[[i]])
}

model_performance[1,3] <- mean(tss_list) #0.3689115
model_performance[1,4] <- sd(tss_list) #0.001044798


# plot ROC curve
library(plotROC)
plotROC_JM <- function(model, test = NULL, text = 25, legend = TRUE) {
  
  if (!requireNamespace("plotROC", quietly = TRUE)) {
    stop("You need the packege \"plotROC\" to run this function,",
         " please install it.",
         call. = FALSE)
  }
  
  if (class(model@model) == "Maxent") {
    type <- "raw"
  } else {
    type <- "link"
  }
  
  df <- data.frame(set = "Train", pa = model@data@pa,
                   pred = predict(model, data = model@data, type = type),
                   stringsAsFactors = FALSE)
  auc <- auc(model)
  labels <- paste("Train", round(auc, 3))
  
  if (!is.null(test)) {
    df_test <- data.frame(set = "Test", pa = test@pa,
                          pred = predict(model, data = test, type = type),
                          stringsAsFactors = FALSE)
    df <- rbind(df, df_test)
    auc <- auc(model, test = test)
    labels <- c(paste("Test", round(auc, 3)), labels)
  }
  
  my_plot <- ggplot(df, aes(m = .data$pred, d = .data$pa, group = .data$set)) +
    plotROC::geom_roc(n.cuts = 0, aes(color = .data$set), size = 0.5) +
    ggplot2::scale_color_discrete(name = "AUC", labels = labels) +
    ggplot2::geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "grey",
                          linetype = 2) +
    ggplot2::labs(x = "False Positive Rate", y = "True Positive Rate") +
    ggplot2::coord_fixed() +
    scale_x_continuous(breaks=seq(0, 1, 0.5)) +
    theme_minimal() +
    theme(text = element_text(size = text))  +
    theme(text = ggplot2::element_text(colour = "#666666"))
  
  if (!is.null(test)) {
    my_plot <- my_plot +
      ggplot2::guides(colour = ggplot2::guide_legend(reverse = TRUE))
  }
  
  if (legend == FALSE) {
    my_plot <- my_plot +
      theme(legend.position = "none")
    
  }
  
  return(my_plot)
}

#both curves also very similar!
pdf(file="Figures/ROC_varImp_allSeasons_final.pdf")
plotROC_JM(model_allSeasons_final@models[[1]]) 
dev.off()

#training line shows the fit of the model data, 
#testing indicates the fit of the model to the testing data: real test of the models predictive power



### response curves for all-seasons model ###

##raster layers need to be loaded for this.
#open source code: https://rdrr.io/cran/SDMtune/src/R/utils.R to be able to change the code
.get_train_args <- function(model) {
  
  args <- list(data = model@data)
  
  if (class(model) == "SDMmodelCV") {
    args$folds <- model@folds
    model <- model@models[[1]]@model
  } else {
    args$folds <- NULL
    model <- model@model
  }
  
  args$method <- class(model)
  
  if (args$method == "Maxent") {
    args$fc <- model@fc
    args$reg <- model@reg
    args$iter <- model@iter
  } else if (args$method == "Maxnet") {
    args$fc <- model@fc
    args$reg <- model@reg
  } else if (args$method == "ANN") {
    args$size <- model@size
    args$decay <- model@decay
    args$rang <- model@rang
    args$maxit <- model@maxit
  } else if (args$method == "RF") {
    args$mtry <- model@mtry
    args$ntree <- model@ntree
    args$nodesize <- model@nodesize
  } else {
    args$distribution <- model@distribution
    args$n.trees <- model@n.trees
    args$interaction.depth <- model@interaction.depth
    args$shrinkage <- model@shrinkage
    args$bag.fraction <- model@bag.fraction
  }
  return(args)
}

.create_model_from_settings <- function(model, settings, verbose = FALSE) {
  
  args <- .get_train_args(model)
  args[names(settings)] <- settings
  args$verbose <- verbose
  output <- suppressMessages(do.call("train", args))
  
  return(output)
}

#changed the functions get_plot_data and plotResponse, so that x axis is backtransformed. 

# Needs now furthermore: 
#- backtransformRaster: the raster layer with the unscaled variables
#- a SWDdata with the set of unscaled predictors that were in the final model 

data_unscaled <- prepareSWD(species = "WVC", p = rehe_WU_HE_occ, a = bg_HE_df, 
                            env = raster_HE_stack_complete_unscaled[[c(5,6,7,8,10,12, 18,19, 21,22)]], #only those layers that are also included in the model
                            categorical = "RoadCategory")

.get_plot_data_JM <- function(model, var, df, cont_vars, cat_vars, n_rows, fun,
                              marginal, type, categ,backtransformRaster) {
  
  data <- data.frame(matrix(NA, nrow = 1, ncol = ncol(df)))
  colnames(data) <- colnames(df)
  data[cont_vars] <- apply(df[cont_vars], 2, fun)
  data[cat_vars] <- as.factor(apply(df[cat_vars], 2, raster::modal))
  data <- do.call("rbind", replicate(n_rows, data, simplify = FALSE))
  
  if (var %in% cont_vars) {
    var_min <- min(model@data@data[var])
    var_max <- max(model@data@data[var])
    data[var] <- seq(var_min, var_max, length.out = n_rows)
  } else {
    data[var] <- factor(categ)
  }
  
  if (!marginal) {
    new_data <- model@data
    new_data@data <- new_data@data[, var, drop = FALSE]
    settings <- list(data = new_data)
    model <- .create_model_from_settings(model, settings)
  }
  
  pred <- predict(model, data, type = type)
  
  if (var %in% cont_vars) {
    backtransform <-  data[, var]*sd(values(backtransformRaster), na.rm = TRUE) + mean(values(backtransformRaster), na.rm = TRUE)
    plot_data <- data.frame(x = data[, var], x_backtransformed = backtransform, y = pred)
  } else {
    plot_data <- data.frame(x = data[, var], y = pred)  
  }
  
  return(plot_data)
}


#https://rdrr.io/cran/SDMtune/src/R/plotResponse.R

plotResponse_JM <- function(model, var, type = NULL, only_presence = FALSE,
                            marginal = TRUE, fun = mean, rug = FALSE,
                            color = "darkblue", backtransformRaster = NULL, backtransformData = NULL, varname) {
  
  if (!var %in% names(model@data@data))
    stop(paste(var, "is not used to train the model!"))
  
  
  
  p_scaled <- model@data@data[model@data@pa == 1, ]
  
  #either df has all poits or only presence points
  if (only_presence) {
    df <- p_scaled
  } else {
    df <- model@data@data
  }
  
  #names of cont and cat variables
  cont_vars <- names(Filter(is.numeric, p_scaled))
  cat_vars <- names(Filter(is.factor, p_scaled))
  
  
  if (var %in% cat_vars) {
    categ <- as.numeric(levels(df[, var]))
    n_rows <- length(categ)
  } else {
    n_rows <- 100
  }
  
  
  
  
  if (class(model) == "SDMmodel") {
    plot_data <- .get_plot_data_JM(model, var, df, cont_vars, cat_vars, n_rows,
                                   fun, marginal, type, categ, backtransformRaster)
    
    if (var %in% cont_vars) {
      my_plot <- ggplot(plot_data, aes(x = .data$x_backtransformed, y = .data$y)) +
        ggplot2::geom_line(colour = color)
      
    } else {
      my_plot <- ggplot(plot_data, aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_bar(stat = "identity", fill = color)
    }
  } else {
    nf <- length(model@models)
    plot_data <- .get_plot_data_JM(model@models[[1]], var, df, cont_vars, cat_vars,
                                   n_rows, fun, marginal, type, categ, backtransformRaster)
    if (var %in% cont_vars) {
      colnames(plot_data) <- c("x", "x_backtransformed", "y_1")
    } else {
      colnames(plot_data) <- c("x", "y_1")
    }
    for (i in 2:nf)
      plot_data[paste0("y_", i)] <- .get_plot_data_JM(model@models[[i]], var, df,
                                                      cont_vars, cat_vars, n_rows,
                                                      fun, marginal, type, categ, backtransformRaster)$y
    if (var %in% cont_vars) {
      plot_data$y <- rowMeans(plot_data[, -c(1,2)])
      plot_data$sd <- apply(plot_data[, 3:(nf + 2)], 1, sd, na.rm = TRUE)
    } else {
      plot_data$y <- rowMeans(plot_data[, -1])
      plot_data$sd <- apply(plot_data[, 2:(nf + 2)], 1, sd, na.rm = TRUE)  
    }
    plot_data$y_min <- plot_data$y - plot_data$sd
    plot_data$y_max <- plot_data$y + plot_data$sd
    
    if (var %in% cont_vars) {
      my_plot <- ggplot(plot_data,
                        aes(x = .data$x_backtransformed, y = .data$y, ymin = .data$y_min,
                            ymax = .data$y_max)) +
        ggplot2::geom_line(colour = color) +
        ggplot2:: geom_ribbon(fill = color, alpha = 0.2)
      
    } else {
      my_plot <- ggplot(plot_data, aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_bar(stat = "identity", fill = color) +
        ggplot2::geom_errorbar(aes(ymin = .data$y_min, ymax = .data$y_max),
                               width = 0.2, size = 0.3)
    }
  }
  
  my_plot <- my_plot +
    ggplot2::labs(x = varname, y = "Probability of occurrence") +
    ggplot2::theme_minimal() +
    theme(text = element_text(size = 25))  +
    ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))
  
  if (rug == TRUE & var %in% cont_vars) {
    #get presence and absence points
    p <- backtransformData@data[backtransformData@pa == 1, ]  #.get_presence(model@data)
    a <- backtransformData@data[backtransformData@pa == 0, ]   #.get_absence(model@data)
    
    #rug for presences and absences
    p_rug <- data.frame(x = p[, var])
    a_rug <- data.frame(x = a[, var])
    
    my_plot <- my_plot +
      ggplot2::geom_rug(data = p_rug, inherit.aes = FALSE, aes(.data$x),
                        sides = "t", color = "#4C4C4C") +
      ggplot2::geom_rug(data = a_rug, inherit.aes = FALSE, aes(.data$x),
                        sides = "b", color = "#4C4C4C")
  }
  
  return(my_plot)
}



pdf(file="Figures/ResponseAllBroadleafedForests.pdf")
plotResponse_JM(model = model_allSeasons_final, var = "BroadleafedForests8_200m", type = "cloglog", 
                only_presence = FALSE,
                marginal = FALSE, fun = mean, rug = TRUE,
                color = "#386CB0", backtransformRaster =  landcover_HE_8_masked, backtransformData = data_unscaled, varname = "Broad-leafed forests (w/200m)")
dev.off()

pdf(file="Figures/ResponseAllRoadDensity.pdf")
plotResponse_JM(model = model_allSeasons_final, var = "RoadDensity_200m", type = "cloglog", 
                only_presence = FALSE,
                marginal = FALSE, fun = mean, rug = TRUE,
                color = "#386CB0",backtransformRaster =  roads_density_masked, 
                backtransformData = data_unscaled, varname = "Road density (w/200m)")
dev.off()



ResponseAllRoadCat <- plotResponse_JM(model_allSeasons_final, var = "RoadCategory", type = "cloglog", 
                                      only_presence =FALSE,
                                      color = "#386CB0",marginal = FALSE, rug = TRUE, varname = "Road category")
pdf(file="Figures/ResponseAllRoadCat.pdf")
ResponseAllRoadCat +
  scale_x_discrete(labels = c("1301" = "Highways", 
                              "1303" = "Federal roads", 
                              "1305" = "State roads", 
                              "1306" = "District roads", 
                              "1307" = "Community roads")) +
  labs(y = "Probability of occurrence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   
dev.off()


pdf(file="Figures/ResponseAllMeadows.pdf")
plotResponse_JM(model = model_allSeasons_final, var = "Meadows1_100m", type = "cloglog", 
                only_presence = FALSE,
                marginal = FALSE, fun = mean, rug = TRUE,
                color = "#386CB0",backtransformRaster =  landcover_HE_1_small_masked, 
                backtransformData = data_unscaled, varname = "Meadows (w/100m)")
dev.off()

pdf(file="Figures/ResponseAllRoadwidth.pdf")
plotResponse_JM(model = model_allSeasons_final, var = "RoadWidth", type = "cloglog", 
                only_presence = FALSE,
                marginal = FALSE, fun = mean, rug = TRUE,
                color = "#386CB0",backtransformRaster =  roads_HE_raster_brf, 
                backtransformData = data_unscaled, varname = "Roadwidth")
dev.off()


pdf(file="Figures/ResponseAllArable.pdf")
plotResponse_JM(model = model_allSeasons_final, var = "ArableAreas9_100m", type = "cloglog", 
                only_presence = FALSE,
                marginal = FALSE, fun = mean, rug = TRUE,
                color = "#386CB0",backtransformRaster =  landcover_HE_9_small_masked, 
                backtransformData = data_unscaled, varname = "Arable areas (w/100m)")
dev.off()

pdf(file="Figures/ResponseAllConifer.pdf")
plotResponse_JM(model = model_allSeasons_final, var = "ConiferForests6_200m", type = "cloglog", 
                only_presence = FALSE,
                marginal = FALSE, fun = mean, rug = TRUE,
                color = "#386CB0",backtransformRaster =  landcover_HE_6_masked, 
                backtransformData = data_unscaled, varname = "Conifer forests (w/200m)")
dev.off()

pdf(file="Figures/ResponseAllMixed.pdf")
plotResponse_JM(model = model_allSeasons_final, var = "MixedForests7_200m", type = "cloglog", 
                only_presence = FALSE,
                marginal = FALSE, fun = mean, rug = TRUE,
                color = "#386CB0",backtransformRaster =  landcover_HE_7_masked, 
                backtransformData = data_unscaled, varname = "Mixed forests (w/200m)")
dev.off()



### Model predictions with all-seasons model ###

#predict: output is a vecotr containing all the predicted values for the training locations
predict_all <- predict(model_allSeasons_final, data = data_all, type = "cloglog")

#cor statistic
cor_list <- c()
for(i in 1:5){
  predict <- predict(model_allSeasons_final@models[[i]], data = data_all, type = "cloglog")
  cor_list[i] <- cor(predict, data_all@pa)
}

model_performance[1,5] <- mean(cor_list) #0.3689115
model_performance[1,6] <- sd(cor_list) #0.001044798

#pocplot(predict_all[data_all@pa == 1], predict_all[data_all@pa == 0], linearize = TRUE)

# predictionts only for the presence locations
p <- data_all@data[data_all@pa == 1, ]
hist(predict(model_allSeasons_final, data = p, type = "cloglog")) #prediction for all presence records i.e. something like residual analysis 


#creating a distribution map
map_cv <- predict(model_allSeasons_final, data = raster_HE_stack_complete[[c(5,6,7,8,10,12, 18,19, 21,22)]],
                  type = "cloglog", extent = extent(470000, 480000, 5555000, 5565000), progress = "text") #add rasterstack as data, can take long: check out arguemnts progress, parallel and extent!

pdf(file="Figures/PredMap_All.pdf")
plotPred(map_cv, lt = "Probability of \nDVC occurrence",colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")) #can be costumized, more on https://consbiol-unibern.github.io/SDMtune/articles/articles/make_predictions.html
dev.off()

map_large <- predict(model_allSeasons_final, data = raster_HE_stack_complete[[c(5,6,7,8,10,12, 18,19, 21,22)]],
                     type = "cloglog", extent = extent(470000, 500000, 5555000, 5595000), progress = "text")

plotPred(map_large, lt = "Probability of \nDVC occurrence",colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")) 



#### Season models ####

rehe_WU_HE$season <- ifelse(rehe_WU_HE$monat >= 1 & rehe_WU_HE$monat <= 4, "Gestation", 
                            ifelse((rehe_WU_HE$monat == 5 | rehe_WU_HE$monat == 6)|
                                     (rehe_WU_HE$monat == 7 & rehe_WU_HE$tag <= 15) ,"Lactation", 
                                   ifelse((rehe_WU_HE$monat == 7 & rehe_WU_HE$tag >= 16)|
                                            (rehe_WU_HE$monat == 8 & rehe_WU_HE$tag <= 15), "Rut", "Diapause")))

rehe_WU_HE %>% 
  st_drop_geometry() %>% 
  group_by(season, jahr) %>% 
  count() %>% 
  group_by(season) %>% 
  summarise(mean = mean(n), 
            sd = sd(n))



### Rut model 

rehe_WU_HE_occ_rut <- rehe_WU_HE %>% 
  filter(season == "Rut") %>% 
  st_coordinates() %>% 
  as.data.frame()

data_rut <- prepareSWD(species = "WVC", p = rehe_WU_HE_occ_rut, a = bg_HE_df, env = raster_HE_stack_complete, categorical = "RoadCategory")

folds = randomFolds(data_rut, k = 5, only_presence = TRUE)
rut_model_crossval <- train(method = "Maxent", data = data_rut, fc = "l", folds = folds)


model_rut_crossval_vs <- varSel(rut_model_crossval, metric = "auc", bg4cor = bg_extra, cor_th = 0.7,
                                permut = 10) #11 minutes
#Removed variables: Meadows1_200m, Urbanareas4_200m, Industry3_100m, Urbanareas4_100m, 
#Complex.habitats5_100m, ConiferForests6_100m, MixedForests7_100m, BroadleafedForests8_100m, 
#ArableAreas9_100m, RoadDensity_100m
vi_rut_cross_vs <- varImp(model_rut_crossval_vs)

rut_crossval_vs_optimized <- optimizeModel(model_rut_crossval_vs, hypers = args, 
                                           metric = "auc", seed = 789) #30 min

rut_crossval_vs_optimized@results
model_rut_crossval_vs_optimized <- rut_crossval_vs_optimized@models[[1]]#best model collection after optimizing process
# fc: lqph
# reg: 0.4
# iter: 500
model_rut_final <- reduceVar(model_rut_crossval_vs_optimized, th = 2, 
                             metric = "auc", permut = 10, use_jk = TRUE)
#Removed variables: Swamps2_200m, Swamps2_100m
VarImp_rut_final <- varImp(model_rut_final)
plotVarImp(vi_rut_vs_optim_reduced)


auc_list <- c()
for(i in 1:5){
  auc_list[i] <- auc(model_rut_final@models[[i]])
}
model_performance[4,1] <- mean(auc_list) #0.7402395
model_performance[4,2] <- sd(auc_list) # 0.001194891


#TSS: true skill statistics (Allouche, tsoar and kadmon 20061)

tss_list <- c()
for(i in 1:5){
  tss_list[i] <- tss(model_rut_final@models[[i]])
}
model_performance[4,3] <- mean(tss_list) #0.3640745
model_performance[4,4] <- sd(tss_list) #0.002544266

#cor statistic
cor_list <- c()
for(i in 1:5){
  predict <- predict(model_rut_final@models[[i]], data = data_rut, type = "cloglog")
  cor_list[i] <- cor(predict, data_rut@pa)
}

model_performance[4,5] <- mean(cor_list) #0.3928272
model_performance[4,6] <- sd(cor_list) #0.001044798



### Gestation Model


rehe_WU_HE_occ_gest <- rehe_WU_HE %>% 
  filter(season == "Gestation") %>% 
  st_coordinates() %>% 
  as.data.frame()

data_gest <- prepareSWD(species = "WVC", p = rehe_WU_HE_occ_gest, a = bg_HE_df, env = raster_HE_stack_complete, categorical = "RoadCategory")

folds = randomFolds(data_gest, k = 5, only_presence = TRUE)
gest_model_crossval <- train(method = "Maxent", data = data_gest, fc = "l", folds = folds)

model_gest_crossval_vs <- varSel(gest_model_crossval, metric = "auc", bg4cor = bg_extra, cor_th = 0.7,
                                 permut = 10) #19 min
#Removed variables: Meadows1_200m, Urbanareas4_200m, ConiferForests6_200m, ArableAreas9_200m, 
#Industry3_100m, Urbanareas4_100m, Complex.habitats5_100m, MixedForests7_100m, 
#BroadleafedForests8_100m, RoadDensity_100m

vi_gest_cross_vs <- varImp(model_gest_crossval_vs)


gest_crossval_vs_optimized <- optimizeModel(model_gest_crossval_vs, hypers = args, 
                                            metric = "auc", seed = 789) 


gest_crossval_vs_optimized@results
model_gest_crossval_vs_optimized <- gest_crossval_vs_optimized@models[[1]]#best model collection after optimizing process
# fc: lqph
# reg: 0.4
# iter: 500

model_gest_final <- reduceVar(model_gest_crossval_vs_optimized, th = 2, 
                              metric = "auc", permut = 10, use_jk = TRUE)

VarImp_gest_final <- varImp(model_gest_final)
plotVarImp(VarImp_gest_final)


auc_list <- c()
for(i in 1:5){
  auc_list[i] <- auc(model_gest_final@models[[i]])
}
model_performance[2,1] <- mean(auc_list) #0.7523003
model_performance[2,2] <- sd(auc_list) # 0.0004652207


#TSS: true skill statistics (Allouche, tsoar and kadmon 20061)

tss_list <- c()
for(i in 1:5){
  tss_list[i] <- tss(model_gest_final@models[[i]])
}
model_performance[2,3] <- mean(tss_list) # 0.3807672
model_performance[2,4] <- sd(tss_list) #0.00137268

#cor statistic
cor_list <- c()
for(i in 1:5){
  predict <- predict(model_gest_final@models[[i]], data = data_gest, type = "cloglog")
  cor_list[i] <- cor(predict, data_gest@pa)
}

model_performance[2,5] <- mean(cor_list) #0.4521806
model_performance[2,6] <- sd(cor_list) #0.000311345



### Lactation Model

rehe_WU_HE_occ_lact <- rehe_WU_HE %>% 
  filter(season == "Lactation") %>% 
  st_coordinates() %>% 
  as.data.frame()

data_lact <- prepareSWD(species = "WVC", p = rehe_WU_HE_occ_lact, a = bg_HE_df, env = raster_HE_stack_complete, categorical = "RoadCategory")

folds = randomFolds(data_lact, k = 5, only_presence = TRUE)
lact_model_crossval <- train(method = "Maxent", data = data_lact, fc = "l", folds = folds)

model_lact_crossval_vs <- varSel(lact_model_crossval, metric = "auc", bg4cor = bg_extra, cor_th = 0.7,
                                 permut = 10) #15 min
#Removed variables: Meadows1_200m, Industry3_200m, Urbanareas4_200m, ArableAreas9_200m, 
#Urbanareas4_100m, Complex.habitats5_100m, ConiferForests6_100m, MixedForests7_100m, 
#BroadleafedForests8_100m, RoadDensity_100m
vi_lact_cross_vs <- varImp(model_lact_crossval_vs)

lact_crossval_vs_optimized <- optimizeModel(model_lact_crossval_vs, hypers = args, 
                                            metric = "auc", seed = 789) #41 min 


lact_crossval_vs_optimized@results
model_lact_crossval_vs_optimized <- lact_crossval_vs_optimized@models[[1]]#best model collection after optimizing process
# fc: lqph
# reg: 0.4
# iter: 500

model_lact_final <- reduceVar(model_lact_crossval_vs_optimized, th = 2, 
                              metric = "auc", permut = 10, use_jk = TRUE)

#Removed variables: Swamps2_200m

VarImp_lact_final <- varImp(model_lact_final)
plotVarImp(VarImp_lact_final)



auc_list <- c()
for(i in 1:5){
  auc_list[i] <- auc(model_lact_final@models[[i]])
}

model_performance[3,1] <- mean(auc_list) #0.7367191
model_performance[3,2] <- sd(auc_list) # 0.0004898628


#TSS: true skill statistics (Allouche, tsoar and kadmon 20061)

tss_list <- c()
for(i in 1:5){
  tss_list[i] <- tss(model_lact_final@models[[i]])
}
model_performance[3,3] <- mean(tss_list) # 0.3578278
model_performance[3,4] <- sd(tss_list) #0.001622327


#cor statistic
cor_list <- c()
for(i in 1:5){
  predict <- predict(model_lact_final@models[[i]], data = data_lact, type = "cloglog")
  cor_list[i] <- cor(predict, data_lact@pa)
}

model_performance[3,5] <- mean(cor_list) #0.4346337
model_performance[3,6] <- sd(cor_list) #0.0002692445


### Diapause Model 


rehe_WU_HE_occ_dia <- rehe_WU_HE %>% 
  filter(season == "Diapause") %>% 
  st_coordinates() %>% 
  as.data.frame()

data_dia <- prepareSWD(species = "WVC", p = rehe_WU_HE_occ_dia, a = bg_HE_df, env = raster_HE_stack_complete, categorical = "RoadCategory")

folds = randomFolds(data_dia, k = 5, only_presence = TRUE)
dia_model_crossval <- train(method = "Maxent", data = data_dia, fc = "l", folds = folds)

model_dia_crossval_vs <- varSel(dia_model_crossval, metric = "auc", bg4cor = bg_extra, cor_th = 0.7,
                                permut = 10) #15 min
#Removed variables: Meadows1_200m, Urbanareas4_200m, ArableAreas9_200m, 
#Industry3_100m, Urbanareas4_100m, Complex.habitats5_100m, ConiferForests6_100m, 
#MixedForests7_100m, BroadleafedForests8_100m, RoadDensity_100m

vi_dia_cross_vs <- varImp(model_dia_crossval_vs)



dia_crossval_vs_optimized <- optimizeModel(model_dia_crossval_vs, hypers = args, 
                                           metric = "auc", seed = 789) #58 min 


dia_crossval_vs_optimized@results
model_dia_crossval_vs_optimized <- dia_crossval_vs_optimized@models[[1]]#best model collection after optimizing process
# fc: lqph
# reg: 0.4
# iter: 500

model_dia_final <- reduceVar(model_dia_crossval_vs_optimized, th = 2, 
                             metric = "auc", permut = 10, use_jk = TRUE)

#Removed variables: Swamps2_200m, Swamps2_100m
VarImp_dia_final <- varImp(model_dia_final)
plotVarImp(VarImp_dia_final)



auc_list <- c()
for(i in 1:5){
  auc_list[i] <- auc(model_dia_final@models[[i]])
}

model_performance[5,1] <- mean(auc_list) #0.7466063
model_performance[5,2] <- sd(auc_list) #0.0002898015


#TSS: true skill statistics (Allouche, tsoar and kadmon 20061)

tss_list <- c()
for(i in 1:5){
  tss_list[i] <- tss(model_dia_final@models[[i]])
}
model_performance[5,3] <- mean(tss_list) #0.3825857
model_performance[5,4] <- sd(tss_list)  # 0.001241072


#cor statistic
cor_list <- c()
for(i in 1:5){
  predict <- predict(model_dia_final@models[[i]], data = data_dia, type = "cloglog")
  cor_list[i] <- cor(predict, data_dia@pa)
}

model_performance[5,5] <- mean(cor_list) #0.4543712
model_performance[5,6] <- sd(cor_list) #0.0006027654



### Combined response curves for season models ###

#function to extract the data to be plotted
plotResponseData_JM <- function(model, var, type = NULL, only_presence = FALSE,
                                marginal = TRUE, fun = mean,
                                backtransformRaster = NULL) {
  
  if (!var %in% names(model@data@data))
    stop(paste(var, "is not used to train the model!"))
  
  
  
  p_scaled <- model@data@data[model@data@pa == 1, ]
  
  #either df has all poits or only presence points
  if (only_presence) {
    df <- p_scaled
  } else {
    df <- model@data@data
  }
  
  #names of cont and cat variables
  cont_vars <- names(Filter(is.numeric, p_scaled))
  cat_vars <- names(Filter(is.factor, p_scaled))
  
  
  if (var %in% cat_vars) {
    categ <- as.numeric(levels(df[, var]))
    n_rows <- length(categ)
  } else {
    n_rows <- 100
  }
  
  
  
  
  if (class(model) == "SDMmodel") {
    plot_data <- .get_plot_data_JM(model, var, df, cont_vars, cat_vars, n_rows,
                                   fun, marginal, type, categ, backtransformRaster)
    
    if (var %in% cont_vars) {
      my_plot <- ggplot(plot_data, aes(x = .data$x_backtransformed, y = .data$y)) +
        ggplot2::geom_line(colour = color)
      
    } else {
      my_plot <- ggplot(plot_data, aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_bar(stat = "identity", fill = color)
    }
  } else {
    nf <- length(model@models)
    plot_data <- .get_plot_data_JM(model@models[[1]], var, df, cont_vars, cat_vars,
                                   n_rows, fun, marginal, type, categ, backtransformRaster)
    if (var %in% cont_vars) {
      colnames(plot_data) <- c("x", "x_backtransformed", "y_1")
    } else {
      colnames(plot_data) <- c("x", "y_1")
    }
    for (i in 2:nf)
      plot_data[paste0("y_", i)] <- .get_plot_data_JM(model@models[[i]], var, df,
                                                      cont_vars, cat_vars, n_rows,
                                                      fun, marginal, type, categ, backtransformRaster)$y
    if (var %in% cont_vars) {
      plot_data$y <- rowMeans(plot_data[, -c(1,2)])
      plot_data$sd <- apply(plot_data[, 3:(nf + 2)], 1, sd, na.rm = TRUE)
    } else {
      plot_data$y <- rowMeans(plot_data[, -1])
      plot_data$sd <- apply(plot_data[, 2:(nf + 2)], 1, sd, na.rm = TRUE)  
    }
    plot_data$y_min <- plot_data$y - plot_data$sd
    plot_data$y_max <- plot_data$y + plot_data$sd
  }
  
  return(plot_data)
}

# function to bind all response datasets together to be able to plot them in one graph
BindResponseData <- function(var, backtransformRaster = NULL){
  
  ResponseRut <- plotResponseData_JM(model = model_rut_final, 
                                     var = var, 
                                     type = "cloglog", only_presence = FALSE,
                                     marginal = FALSE, fun = mean,
                                     backtransformRaster =  backtransformRaster)
  ResponseRut$Season <- "Rut"
  ResponseGest <- plotResponseData_JM(model = model_gest_final, 
                                      var = var, 
                                      type = "cloglog", only_presence = FALSE,
                                      marginal = FALSE, fun = mean,
                                      backtransformRaster =  backtransformRaster)
  ResponseGest$Season <- "Gestation"
  ResponseLact <- plotResponseData_JM(model = model_lact_final, 
                                      var = var, 
                                      type = "cloglog", only_presence = FALSE,
                                      marginal = FALSE, fun = mean,
                                      backtransformRaster =  backtransformRaster)
  ResponseLact$Season <- "Lactation"
  ResponseDia <- plotResponseData_JM(model = model_dia_final, 
                                     var = var, 
                                     type = "cloglog", only_presence = FALSE,
                                     marginal = FALSE, fun = mean,
                                     backtransformRaster =  backtransformRaster)
  ResponseDia$Season <- "Diapause"
  
  Response <- rbind(ResponseDia, ResponseGest, 
                    ResponseLact, ResponseRut)
  Response$Season <- factor(Response$Season, levels = c("Gestation", "Lactation", "Rut", "Diapause"))
  return(Response)
  
}


#Broadleafed Forests
Response_broadleafed <- BindResponseData(var = "BroadleafedForests8_200m", 
                                         backtransformRaster = landcover_HE_8_masked) 

pdf(file = "Figures/ResponseSeasonsBroadleafed.pdf")
ggplot(Response_broadleafed,
       aes(x = x_backtransformed, y = y, ymin = y_min, ymax = y_max, 
           fill = Season, linetype = Season)) +
  ggplot2::geom_line() +
  ggplot2:: geom_ribbon(alpha = 0.5)+
  scale_fill_manual(values=c("#FFFF99", "#FDC086", "#7FC97F", "#BEAED4")) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  ggplot2::labs(x = "Broad-leafed forests (w/200m)", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))

dev.off()

#Road density
Response_RoadDensity <- BindResponseData(var = "RoadDensity_200m", 
                                         backtransformRaster = roads_density_masked) 


pdf(file = "Figures/ResponseSeasonsRoadDensity.pdf")
ggplot(Response_RoadDensity,
       aes(x = x_backtransformed, y = y, ymin = y_min, ymax = y_max, 
           fill = Season, linetype = Season)) +
  ggplot2::geom_line() +
  ggplot2:: geom_ribbon(alpha = 0.5)+
  scale_fill_manual(values=c("#FFFF99", "#FDC086", "#7FC97F", "#BEAED4")) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  ggplot2::labs(x = "Road density (w/200m)", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))

dev.off()



#Road category
Response_RoadCat <- BindResponseData(var = "RoadCategory") 

pdf(file = "Figures/ResponseSeasonsRoadCat.pdf")
ggplot(Response_RoadCat, aes(x = x, y = y, fill = Season)) +
  ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(position=position_dodge(width=0.9), aes(ymin = y_min, ymax = y_max), width = .2, size = 0.3) +
  scale_fill_manual(values=c("#FFFF99", "#FDC086", "#7FC97F", "#BEAED4")) +
  ggplot2::labs(x = "Road category", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  scale_x_discrete(labels = c("1301" = "Highways", 
                              "1303" = "Federal roads", 
                              "1305" = "State roads", 
                              "1306" = "District roads", 
                              "1307" = "Community roads")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "#666666"))
dev.off()


#Meadows 
Response_Meadows <- BindResponseData(var = "Meadows1_100m", 
                                     backtransformRaster = landcover_HE_1_small_masked) 

pdf(file = "Figures/ResponseSeasonsMeadows.pdf")
ggplot(Response_Meadows,
       aes(x = x_backtransformed, y = y, ymin = y_min, ymax = y_max, 
           fill = Season, linetype = Season)) +
  ggplot2::geom_line() +
  ggplot2:: geom_ribbon(alpha = 0.5)+
  scale_fill_manual(values=c("#FFFF99", "#FDC086", "#7FC97F", "#BEAED4")) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  ggplot2::labs(x = "Meadows (w/100m)", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))

dev.off()

#roadwidth
Response_Roadwidth <- BindResponseData(var = "RoadWidth", 
                                       backtransformRaster = roads_HE_raster_brf) 

pdf(file = "Figures/ResponseSeasonRoadwidth.pdf")
ggplot(Response_Roadwidth,
       aes(x = x_backtransformed, y = y, ymin = y_min, ymax = y_max, 
           fill = Season, linetype = Season)) +
  ggplot2::geom_line() +
  ggplot2:: geom_ribbon(alpha = 0.5)+
  scale_fill_manual(values=c("#FFFF99", "#FDC086", "#7FC97F", "#BEAED4")) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  ggplot2::labs(x = "Roadwidth", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))

dev.off()

#Conifer
Response_conifer <- BindResponseData(var = "ConiferForests6_200m", 
                                     backtransformRaster = landcover_HE_6_masked) 

pdf(file = "Figures/ResponseSeasonConifer.pdf")
ggplot(Response_conifer,
       aes(x = x_backtransformed, y = y, ymin = y_min, ymax = y_max, 
           fill = Season, linetype = Season)) +
  ggplot2::geom_line() +
  ggplot2:: geom_ribbon(alpha = 0.5)+
  scale_fill_manual(values=c("#FFFF99", "#FDC086", "#7FC97F", "#BEAED4")) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  ggplot2::labs(x = "Conifer forests (w/200m)", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))

dev.off()


#mixed forests
Response_mixed <- BindResponseData(var = "MixedForests7_200m", 
                                   backtransformRaster = landcover_HE_7_masked) 

pdf(file = "Figures/ResponseSeasonmixed.pdf")
ggplot(Response_mixed,
       aes(x = x_backtransformed, y = y, ymin = y_min, ymax = y_max, 
           fill = Season, linetype = Season)) +
  ggplot2::geom_line() +
  ggplot2:: geom_ribbon(alpha = 0.5)+
  scale_fill_manual(values=c("#FFFF99", "#FDC086", "#7FC97F", "#BEAED4")) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  ggplot2::labs(x = "Mixed forests (w/200m)", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))

dev.off()


#arable areas: Rut model has varibale (w/200m) while all others have (w/100)
ResponseRutArable200 <- plotResponseData_JM(model = model_rut_final, 
                                            var = "ArableAreas9_200m", 
                                            type = "cloglog", only_presence = FALSE,
                                            marginal = FALSE, fun = mean,
                                            backtransformRaster =  landcover_HE_9_masked)
ResponseRutArable200$Season <- "Rut"

pdf(file = "Figures/ResponseRutArable200.pdf")

ggplot(ResponseRutArable200,
       aes(x = x_backtransformed, y = y, ymin = y_min, ymax = y_max, fill = Season, linetype = Season)) +
  ggplot2::geom_line(linetype =  "dotted") +
  ggplot2:: geom_ribbon(alpha = 0.5, fill = "#7FC97F")+
  ggplot2::labs(x = "Arable areas (w/200m)", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))

dev.off()


#arable areas 100m
ResponseGestArable100 <- plotResponseData_JM(model = model_gest_final, 
                                             var = "ArableAreas9_100m", 
                                             type = "cloglog", only_presence = FALSE,
                                             marginal = FALSE, fun = mean,
                                             backtransformRaster =  landcover_HE_9_small_masked)
ResponseGestArable100$Season <- "Gestation"
ResponseLactArable100 <- plotResponseData_JM(model = model_lact_final, 
                                             var = "ArableAreas9_100m", 
                                             type = "cloglog", only_presence = FALSE,
                                             marginal = FALSE, fun = mean,
                                             backtransformRaster = landcover_HE_9_small_masked)
ResponseLactArable100$Season <- "Lactation"
ResponseDiaArable100 <- plotResponseData_JM(model = model_dia_final, 
                                            var = "ArableAreas9_100m", 
                                            type = "cloglog", only_presence = FALSE,
                                            marginal = FALSE, fun = mean,
                                            backtransformRaster =  landcover_HE_9_small_masked)
ResponseDiaArable100$Season <- "Diapause"

ResponseArable100 <- rbind(ResponseDiaArable100, ResponseGestArable100, 
                           ResponseLactArable100)
ResponseArable100$Season <- factor(ResponseArable100$Season, levels = c("Gestation", "Lactation", "Diapause"))

pdf(file = "Figures/ResponseSeasonArable100.pdf")
ggplot(ResponseArable100,
       aes(x = x_backtransformed, y = y, ymin = y_min, ymax = y_max, 
           fill = Season, linetype = Season)) +
  ggplot2::geom_line() +
  ggplot2:: geom_ribbon(alpha = 0.5)+
  scale_fill_manual(values=c("#FFFF99", "#FDC086", "#BEAED4")) +
  scale_linetype_manual(values=c("solid", "twodash",  "longdash")) +
  ggplot2::labs(x = "Arable areas (w/100m)", y = "Probability of occurrence") +
  ggplot2::theme_minimal() +
  theme(text = element_text(size = 25))  +
  ggplot2::theme(text = ggplot2::element_text(colour = "#666666"))

dev.off()



# combined Variable Importance Plot
VarImp_allSeasons_final$Season <- "All"
VarImp_lact_final$Season <- "Lactation"
VarImp_rut_final$Season <- "Rut"
VarImp_gest_final$Season <- "Gestation"
VarImp_dia_final$Season <- "Diapause"


vi_combined2 <- rbind(VarImp_allSeasons_final ,VarImp_lact_final, VarImp_rut_final, VarImp_gest_final, VarImp_dia_final)

newdata <- data.frame(Variable = rep(levels(as.factor(vi_combined2$Variable)), times = 5), 
                      Season = rep(c("All", "Lactation", "Rut", "Gestation", "Diapause"), each = 13))


vi_combined <- newdata %>% 
  left_join(vi_combined2, by = c("Variable", "Season")) %>% 
  dplyr::mutate(Permutation_importance = coalesce(Permutation_importance, 0), 
                sd = coalesce(sd, 0))

vi_combined$Season <- factor(vi_combined$Season, levels = c("All", "Gestation", "Lactation", "Rut", "Diapause"))

ggplot(vi_combined, aes(x = reorder(Variable, Permutation_importance), y = Permutation_importance, fill = Season)) +
  geom_bar(stat = "identity",  position = position_dodge()) +
  geom_errorbar(position=position_dodge(width=0.9), aes(ymin = Permutation_importance -sd, ymax = Permutation_importance +sd), width = .2) +
  theme_minimal()+
  theme(text = element_text(size = 20))  +
  labs(y = "Permutation importance in %", x = "") +
  scale_fill_manual(values=c("#386CB0", "#FFFF99", "#FDC086", "#7FC97F", "#BEAED4")) +
  scale_x_discrete(labels = c("BroadleafedForests8_200m" = "Broad-leafed forests (w/200m)", 
                              "RoadDensity_200m" = "Road density (w/200m)", 
                              "RoadCategory" = "Road category", 
                              "ArableAreas9_100m" = "Arable areas (w/100m)", 
                              "Meadows1_100m" = "Meadows (w/100m)", 
                              "RoadWidth" = "Roadwidth", 
                              "ConiferForests6_200m" = "Conifer forests (w/200m)", 
                              "MixedForests7_200m" = "Mixed forests (w/200m)", 
                              "Complex.habitats5_200m" = "Complex habitats (w/200m)", 
                              "Industry3_100m" = "Industry (w/100m)", 
                              "ArableAreas9_200m" = "Arable areas (w/200m)", 
                              "Industry3_200m" = "Industry (w/200m)", 
                              "Swamps2_100m" = "Swamps (w/100m)")) +
  coord_flip()
ggsave("Figures/VarImp_combined.pdf")


#export model performance table for latex format
library(xtable)
model_performance
model_performance_ltx <- xtable(model_performance, caption = "Model performance")
print(model_performance_ltx,file="Figures/model_performance_ltx.tex",table.placement = "h", 
      include.rownames = TRUE, caption.placement="bottom")





#### Baden-Württemberg Data Preparation for validation ####


## read in environmental variables for BW#

setwd("W:/Masterarbeit Daten/Analyse/Predictor Rasters/BW")

raslist <- list.files(getwd(), pattern=".tiff$", full.names=FALSE)

for(r in raslist){
  assign(substr(r, 1, nchar(r)-5), raster(r))
  print(r)
}


setwd("W:/Masterarbeit Daten/Analyse")


#make a raster stack: scaled variables! 
raster_BW_stack_complete <- stack(landcover_BW_1_scaled, landcover_BW_2_scaled, landcover_BW_3_scaled, landcover_BW_4_scaled, 
                                  landcover_BW_5_scaled, landcover_BW_6_scaled, landcover_BW_7_scaled, landcover_BW_8_scaled,
                                  landcover_BW_9_scaled, landcover_BW_1_small_scaled, landcover_BW_2_small_scaled, landcover_BW_3_small_scaled,
                                  landcover_BW_4_small_scaled, landcover_BW_5_small_scaled, landcover_BW_6_small_scaled,
                                  landcover_BW_7_small_scaled, landcover_BW_8_small_scaled, landcover_BW_9_small_scaled,
                                  roads_density_BW_scaled, roads_density_BW_100m_scaled, roads_BW_raster_brf_scaled, roads_BW_raster_wdm)

names(raster_BW_stack_complete) <- c("Meadows1_200m", "Swamps2_200m", "Industry3_200m", "Urbanareas4_200m", "Complex habitats5_200m", 
                                     "ConiferForests6_200m", "MixedForests7_200m", "BroadleafedForests8_200m", "ArableAreas9_200m", 
                                     "Meadows1_100m", "Swamps2_100m", "Industry3_100m", "Urbanareas4_100m", "Complex habitats5_100m", 
                                     "ConiferForests6_100m", "MixedForests7_100m", "BroadleafedForests8_100m", "ArableAreas9_100m", 
                                     "RoadDensity_200m", "RoadDensity_100m","RoadWidth","RoadCategory")


rm(landcover_BW_1_scaled, landcover_BW_2_scaled, landcover_BW_3_scaled, landcover_BW_4_scaled, 
   landcover_BW_5_scaled, landcover_BW_6_scaled, landcover_BW_7_scaled, landcover_BW_8_scaled,
   landcover_BW_9_scaled, landcover_BW_1_small_scaled, landcover_BW_2_small_scaled, landcover_BW_3_small_scaled,
   landcover_BW_4_small_scaled, landcover_BW_5_small_scaled, landcover_BW_6_small_scaled,
   landcover_BW_7_small_scaled, landcover_BW_8_small_scaled, landcover_BW_9_small_scaled,
   roads_density_BW_scaled, roads_density_BW_100m_scaled, roads_BW_raster_brf_scaled, roads_BW_raster_wdm)

## add BW road data: same script as in WU_preparation


#### ATKIS Strassen layer (metadatentabelle: ATKIS_Objektartenkatalog_Basis_DLM_7.1) 
#copied from W_datapreparation!!
roads_BW <- read_sf("Layers/strassen_BW.shp")
roads_BW$brf[roads_BW$brf == -9998] <- NA 
roads_BW$wdm <- as.numeric(roads_BW$wdm)
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


## set up libraries and rdata
library(sf)
library(dplyr)
library(tmap)
library(raster)
library(ggplot2)
library(mapview)
library(grid)
library(lubridate)
library(adehabitatHR)

setwd("W:/Masterarbeit Daten/Analyse")


# WVC validation with WU_BW data 
WU_BW <- read.csv2("W:/Masterarbeit Daten/Analyse/SDM/BW_Wildunfall_1.5.21-30.4.22.csv")
#summary(WU_BW)

WU_BW <- filter(WU_BW, !is.na(US_GEO_X) & !is.na(US_GEO_Y))
WU_BW <- dplyr::select(WU_BW, UN_KEY, KO_UDATUM, KO_UMONAT, KO_UZEIT, KO_USTDE, KO_WOTAG, US_GDE, US_STRA1, SS_EDATUM, SS_EZEIT, US_GEO_X, US_GEO_Y, UN_KAT, Hergangstext) #unterschied zwischen KO_UDATM und SS_EDATUM verstehen -> mathias jost

# extract rows that have Reh in the "Hergansgtext"
rehe_WU_BW <- WU_BW[(grep("Reh", WU_BW[,14])), ] #12907 rows. 

rehe_WU_BW <- st_as_sf(rehe_WU_BW, coords = c("US_GEO_X", "US_GEO_Y"), crs = "EPSG:4326")
rehe_WU_BW <- st_transform(rehe_WU_BW, projection)

bawu_sf <- filter(germany_sf, NAME_1 == "Baden-Württemberg")
bawu_sf$VARNAME_1 <- "BW"

rehe_WU_BW <- st_crop(rehe_WU_BW, bawu_sf) #all points within BW. 

#summary(rehe_WU_BW$UN_KAT) # unfallkategory mostly 5 as well. 

#### snap coordinates #### same as in WU_datapreparation for the WU in hessen
#change coordinates of the points that are just next to the road lines: snap 
tst <- st_snap_points(rehe_WU_BW, roads_BW, max_dist = 150) #takes about 3.5 hours
rehe_WU_BW$geometry <- tst


#check how many are outside the buffer of 150 meters
roads_BW_buffered150 <- st_buffer(roads_BW, dist = 150)
outside <- sapply(st_intersects(rehe_WU_BW, roads_BW_buffered150),function(x){length(x)==0})

rehe_WU_BW_outside <- rehe_WU_BW[outside, ] #272 points
rehe_WU_BW <- rehe_WU_BW[!outside, ]


#split the points per season
rehe_WU_BW$KO_USTD <- paste(rehe_WU_BW$KO_USTD, ":00", sep = "")
rehe_WU_BW$date_time <- with(rehe_WU_BW, ymd(KO_UDAT) + hm(KO_USTD))
rehe_WU_BW$date_time <- as.POSIXct(rehe_WU_BW$date_time) #perice on the hour
rehe_WU_BW <- rename(rehe_WU_BW, monat = KO_UMON)
rehe_WU_BW$tag <- day(rehe_WU_BW$date_time)

rehe_WU_BW$season <- ifelse(rehe_WU_BW$monat >= 1 & rehe_WU_BW$monat <= 4, "Gestation", 
                            ifelse((rehe_WU_BW$monat == 5 | rehe_WU_BW$monat == 6)|
                                     (rehe_WU_BW$monat == 7 & rehe_WU_BW$tag <= 15) ,"Lactation", 
                                   ifelse((rehe_WU_BW$monat == 7 & rehe_WU_BW$tag >= 16)|
                                            (rehe_WU_BW$monat == 8 & rehe_WU_BW$tag <= 15), "Rut", "Diapause")))




rehe_WU_BW_occ <- rehe_WU_BW %>% 
  st_coordinates() %>% 
  as.data.frame()


#construct background locations just a it was done for the WU in hessen, because most points are again from accidant category 5!
occur.ras <- rasterize(rehe_WU_BW_occ, raster_BW_stack_complete, 1) #fill cells with data with 1

presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

rehe_WU_BW_dens <- kde2d(x = pres.locs[,1], y = pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras))) 
rehe_WU_BW_dens_raster <- raster(rehe_WU_BW_dens)
plot(rehe_WU_BW_dens_raster)


#set up projection to be used thorughout the script!
extent.bw <- extent(388200, 610200, 5261000, 5516000) #xmin, xmax, ymin, ymax -> extent of bawu and a bit outside
#default raster to use for all bawu covariates
r.raster_bw <- raster()
extent(r.raster_bw) <- extent.bw #extent of bawu
res(r.raster_bw) <- 50
crs(r.raster_bw) <- projection


crs(rehe_WU_BW_dens_raster) <- projection
bias_file <- projectRaster(rehe_WU_BW_dens_raster, r.raster_bw)

bias_file_masked <- mask(bias_file, raster_BW_stack_complete[[22]]) #mask the bias file with wdm raster as all others too

#extract 10000 background points with probability of being sampled based on the bias raster 
bg_BW <- xyFromCell(bias_file_masked, sample(which(!is.na(values(bias_file_masked))), 10000, 
                                             prob=values(bias_file_masked)[!is.na(values(bias_file_masked))])) 

bg_BW_df <- bg_BW %>% 
  as.data.frame() %>% 
  rename(X = x, 
         Y = y)
#bg_sf <- st_as_sf(bg_BW_df, coords = c("X", "Y"), crs = projection)


### SDM validation with BW ####

# make SWD dataframe
# check which points do not have predicotr information: they are spread randomly (lie just between two points), so okay to delete them
extract <- as.data.frame(raster::extract(raster_BW_stack_complete, rehe_WU_BW_occ))
index <- stats::complete.cases(extract)
summary(index)
rehe_WU_BW[!index, ]
# plot(crop(raster_BW_stack_complete[[22]], c(389000, 389150, 5283900, 5284150) ))
# plot(roads_BW, add = TRUE)
# plot(st_geometry(rehe_WU_BW[which(rehe_WU_BW$KO_DAT == 20210521 & rehe_WU_BW$KO_UZEI == 310), ], add = TRUE))
rehe_WU_BW <- rehe_WU_BW[index, ]
# 18 locations without data, deleted. 

# extract risk value from all seasons model for all points
rehe_WU_BW_occ <- rehe_WU_BW %>% 
  st_coordinates() %>% 
  as.data.frame()



# prepare an SWD object for BW with only presence locations. 
rehe_WU_BW_SWD <- prepareSWD(species = "WVC", p = rehe_WU_BW_occ, a = bg_BW_df, env = raster_BW_stack_complete, categorical = "RoadCategory")
#Info: 14 background locations are NA for some environmental variables, they are discarded!

#compute auc value for all seasons model
pdf(file="Figures/ROC_allSeasons_compBW.pdf")
plotROC_JM(model_allSeasons_final@models[[1]], test = rehe_WU_BW_SWD, legend = FALSE) #testing data from BW performs equally well as training data!
dev.off()

#predict risk with best all seasons model: split pa == 0 und pa == 1
# rehe_WU_BW$risk_allSeasons <- predict(model_allSeasons_final, data = rehe_WU_BW_SWD@data[rehe_WU_BW_SWD@pa == 1, ], type = "cloglog") #predict for all true presence points
pred_all_df <- data.frame(ID = seq(1:length(rehe_WU_BW_SWD@pa)),
                          p_num = rehe_WU_BW_SWD@pa,
                          p = as.factor(ifelse(rehe_WU_BW_SWD@pa == 1, "DVC Presence", "Background")), 
                          pred_risk = predict(model_allSeasons_final, data = rehe_WU_BW_SWD, type = "cloglog"))


# compuite COR value
cor(pred_all_df$pred_risk, pred_all_df$p_num) #0.4484805


# histogramm of presence and absence risk distribution
ggplot(pred_all_df , aes(x=pred_risk)) +
  geom_histogram(aes(color = p, fill = p), position = 'identity', alpha = 0.4, bins = 30) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  labs(x = "Predicted DVC risk", y = "Frequency of occurence", fill = "") +
  guides(color = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 20))

ggsave("Figures/WVC_BW_presabs.pdf", width = 25, height = 20, units = "cm")


#two sample test between both
mean(filter(pred_all_df, p == "DVC Presence")$pred_risk)
mean(filter(pred_all_df, p == "Background")$pred_risk)

sd(filter(pred_all_df, p == "DVC Presence")$pred_risk)
sd(filter(pred_all_df, p == "Background")$pred_risk) #double as large variance

#variance not equal: use Welch t-test (default)
t.test(pred_risk ~ p, data = pred_all_df)

pred_all_df %>% 
  cohens_d(pred_risk ~ p, var.equal = FALSE) #ignore vorzeichen -> means high effect size.



#load functions from Phillips2010 for POC curves
source("SDM/POC_Phillips2010/plots.R")


pdf(file="Figures/POC_allSeasons.pdf")
pocplot(pred_all_df[which(pred_all_df$p_num == 1), ]$pred_risk, pred_all_df[which(pred_all_df$p_num == 0), ]$pred_risk, linearize = TRUE) #extremely small SD, but that is due to the large amount of datapoints. taking samples of the data, increases the SD and also more black rugs appear. 
#using the pocplot command, because we have biased background data as samples from the entire area which can also be known presence locations.  (violates the assumption made using ecalp)
dev.off()



# look at spatial distribution of residuals 
mapview(rehe_WU_BW, zcol = "res")
#rehe_WU_BW$res <- 1 - rehe_WU_BW$risk_allSeasons



#auc for gestation model
rehe_WU_BW_SWD_gest <- prepareSWD(species = "WVC", p = rehe_WU_BW_occ[which(rehe_WU_BW$season == "Gestation"), ], a = bg_BW_df, env = raster_BW_stack_complete, categorical = "RoadCategory")

ROC_gest_compBW <- plotROC_JM(model_gest_final@models[[1]], test = rehe_WU_BW_SWD_gest, text = 15) + 
  labs(title = "Gestation")#testing data from BW performs equally well as training data!
ROC_gest_compBW 


pred_gest_df <- data.frame(ID = seq(1:length(rehe_WU_BW_SWD_gest@pa)),
                           p_num = rehe_WU_BW_SWD_gest@pa,
                           p = as.factor(ifelse(rehe_WU_BW_SWD_gest@pa == 1, "DVC Presence", "Background")),
                           pred_risk = predict(model_gest_final, data = rehe_WU_BW_SWD_gest, type = "cloglog"))


# compuite COR value
cor(pred_gest_df$pred_risk, pred_gest_df$p_num) # 0.3954328


# histogramm of presence and absence risk distribution
ggplot(pred_gest_df , aes(x=pred_risk)) +
  geom_histogram(aes(color = p, fill = p), position = 'identity', alpha = 0.4, bins = 30) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  labs(x = "Predicted DVC risk", y = "Frequency of occurence", fill = "") +
  guides(color = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 20))

#ggsave("Figures/WVC_BW_presabs.pdf", width = 25, height = 20, units = "cm")

#variance not equal: use Welch t-test (default)
t.test(pred_risk ~ p, data = pred_gest_df)

pdf(file="Figures/poc_gest_compBW.pdf")
pocplot(pred_gest_df[which(pred_gest_df$p_num == 1), ]$pred_risk, pred_gest_df[which(pred_gest_df$p_num == 0), ]$pred_risk, linearize = TRUE, title_name = "Gestation")
dev.off()

#auc for lactation model
rehe_WU_BW_SWD_lact <- prepareSWD(species = "WVC", p = rehe_WU_BW_occ[which(rehe_WU_BW$season == "Lactation"), ], a = bg_BW_df, env = raster_BW_stack_complete, categorical = "RoadCategory")

ROC_lact_compBW <- plotROC_JM(model_lact_final@models[[1]], test = rehe_WU_BW_SWD_lact, text = 15)  + 
  labs(title = "Lactation")#testing data from BW performs equally well as training data!



pred_lact_df <- data.frame(ID = seq(1:length(rehe_WU_BW_SWD_lact@pa)),
                           p_num = rehe_WU_BW_SWD_lact@pa,
                           p = as.factor(ifelse(rehe_WU_BW_SWD_lact@pa == 1, "DVC Presence", "Background")),
                           pred_risk = predict(model_lact_final, data = rehe_WU_BW_SWD_lact, type = "cloglog"))


# compuite COR value
cor(pred_lact_df$pred_risk, pred_lact_df$p_num) # 0.3339332


# histogramm of presence and absence risk distribution
ggplot(pred_lact_df , aes(x=pred_risk)) +
  geom_histogram(aes(color = p, fill = p), position = 'identity', alpha = 0.4, bins = 30) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  labs(x = "Predicted DVC risk", y = "Frequency of occurence", fill = "") +
  guides(color = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 20))

#ggsave("Figures/WVC_BW_presabs.pdf", width = 25, height = 20, units = "cm")

#variance not equal: use Welch t-test (default)
t.test(pred_risk ~ p, data = pred_lact_df)


pdf(file="Figures/poc_lact_compBW.pdf")
pocplot(pred_lact_df[which(pred_lact_df$p_num == 1), ]$pred_risk, pred_lact_df[which(pred_lact_df$p_num == 0), ]$pred_risk, linearize = TRUE, title_name = "Lactation")
dev.off()


#auc for rut model
rehe_WU_BW_SWD_rut <- prepareSWD(species = "WVC", p = rehe_WU_BW_occ[which(rehe_WU_BW$season == "Rut"), ], a = bg_BW_df, env = raster_BW_stack_complete, categorical = "RoadCategory")

ROC_rut_compBW <- plotROC_JM(model_rut_final@models[[1]], test = rehe_WU_BW_SWD_rut, text = 15)  + 
  labs(title = "Rut")#testing data from BW performs equally well as training data!



pred_rut_df <- data.frame(ID = seq(1:length(rehe_WU_BW_SWD_rut@pa)),
                          p_num = rehe_WU_BW_SWD_rut@pa,
                          p = as.factor(ifelse(rehe_WU_BW_SWD_rut@pa == 1, "DVC Presence", "Background")),
                          pred_risk = predict(model_rut_final, data = rehe_WU_BW_SWD_rut, type = "cloglog"))


# compuite COR value
cor(pred_rut_df$pred_risk, pred_rut_df$p_num) # 0.2560288


# histogramm of presence and absence risk distribution
ggplot(pred_rut_df , aes(x=pred_risk)) +
  geom_histogram(aes(color = p, fill = p), position = 'identity', alpha = 0.4, bins = 30) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  labs(x = "Predicted DVC risk", y = "Frequency of occurence", fill = "") +
  guides(color = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 20))

#ggsave("Figures/WVC_BW_presabs.pdf", width = 25, height = 20, units = "cm")

#variance not equal: use Welch t-test (default)
t.test(pred_risk ~ p, data = pred_rut_df)


pdf(file="Figures/poc_rut_compBW.pdf")
pocplot(pred_rut_df[which(pred_rut_df$p_num == 1), ]$pred_risk, pred_rut_df[which(pred_rut_df$p_num == 0), ]$pred_risk, linearize = TRUE, title_name = "Rut")
dev.off()

#auc for diapause model
rehe_WU_BW_SWD_dia <- prepareSWD(species = "WVC", p = rehe_WU_BW_occ[which(rehe_WU_BW$season == "Diapause"), ], a = bg_BW_df, env = raster_BW_stack_complete, categorical = "RoadCategory")

ROC_dia_compBW <- plotROC_JM(model_dia_final@models[[1]], test = rehe_WU_BW_SWD_dia, text = 15)  + 
  labs(title = "Diapause")#testing data from BW performs equally well as training data!


pred_dia_df <- data.frame(ID = seq(1:length(rehe_WU_BW_SWD_dia@pa)),
                          p_num = rehe_WU_BW_SWD_dia@pa,
                          p = as.factor(ifelse(rehe_WU_BW_SWD_dia@pa == 1, "DVC Presence", "Background")),
                          pred_risk = predict(model_dia_final, data = rehe_WU_BW_SWD_dia, type = "cloglog"))


# compuite COR value
cor(pred_dia_df$pred_risk, pred_dia_df$p_num) # 0.4070859


# histogramm of presence and absence risk distribution
ggplot(pred_dia_df , aes(x=pred_risk)) +
  geom_histogram(aes(color = p, fill = p), position = 'identity', alpha = 0.4, bins = 30) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  labs(x = "Predicted DVC risk", y = "Frequency of occurence", fill = "") +
  guides(color = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 20))

#ggsave("Figures/WVC_BW_presabs.pdf", width = 25, height = 20, units = "cm")

#variance not equal: use Welch t-test (default)
t.test(pred_risk ~ p, data = pred_dia_df)

pdf(file="Figures/poc_dia_compBW.pdf")
pocplot(pred_dia_df[which(pred_dia_df$p_num == 1), ]$pred_risk, pred_dia_df[which(pred_dia_df$p_num == 0), ]$pred_risk, linearize = TRUE, title_name = "Diapause")
dev.off()



#combine all ROC  plots
pdf(file="Figures/ROC_SeasonModels_compBW.pdf")
ROC_gest_compBW + ROC_lact_compBW + ROC_rut_compBW + ROC_dia_compBW + plot_layout(ncol = 2)
dev.off()


#compare all risk probs in bg points between seasons

bg_risk_dia <- filter(pred_dia_df, p == "Background") %>% 
  ggplot() +
  geom_boxplot(aes(y = pred_risk)) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size = 25)) +
  ylim(0,1) +
  labs(title = "Diapause", y = "Modelled DVC risk")

bg_risk_gest <- filter(pred_gest_df, p == "Background") %>% 
  ggplot() +
  geom_boxplot(aes(y = pred_risk))+
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        text = element_text(size = 25)) +
  ylim(0,1) +
  labs(title = "Gestation", y = "Modelled DVC risk")

bg_risk_lact <- filter(pred_lact_df, p == "Background") %>% 
  ggplot() +
  geom_boxplot(aes(y = pred_risk))+
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size = 25)) +
  ylim(0,1) +
  labs(title = "Lactation", y = "Modelled DVC risk")

bg_risk_rut <- filter(pred_rut_df, p == "Background") %>% 
  ggplot() +
  geom_boxplot(aes(y = pred_risk))+
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), 
        text = element_text(size = 25)) +
  ylim(0,1) +
  labs(title = "Rut", y = "Modelled DVC risk")

pdf(file="Figures/bgpoints_comp.pdf", width = 15, height = 10)
bg_risk_gest + bg_risk_lact + bg_risk_rut + bg_risk_dia + plot_layout(nrow = 1)
dev.off()




