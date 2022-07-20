

## Libraries and set up
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

load("WU_analysis.RData")
Sys.setenv(lang = "en_US") #change error messages to english



### SDM modelling prep ###

setwd("W:/Masterarbeit Daten/Analyse/Predictor Rasters/HE")
#setwd("W:/Masterarbeit Daten/Analyse/Layers/unscaled predictors")

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




## prepare presence (WVC) and background locations ###

#select columns that contain the coordinnates of the locations
#presence WU records

projection <- "+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"

rehe_WU_HE <- st_read("Layers/rehe_WU_HE.shp")
rehe_WU_HE <- st_transform(rehe_WU_HE, projection)

#extract coordinates in matrix
rehe_WU_HE_occ <- rehe_WU_HE %>% 
  st_coordinates() %>% 
  as.data.frame()


#construct background locations
occur.ras <- rasterize(rehe_WU_HE_occ, raster_HE_stack_complete, 1) #fill cells with data with 1

presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

rehe_WU_HE_dens <- kde2d(x = pres.locs[,1], y = pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras))) 
rehe_WU_HE_dens_raster <- raster(rehe_WU_HE_dens)
plot(rehe_WU_HE_dens_raster)


#default raster to use for all hessen covariates
r.raster <- raster()
extent(r.raster) <- extent(412000, 587000, 5471000,  5723500) #xmin, xmax, ymin, ymax -> extent of Hessen and a bit outside
res(r.raster) <- 50
crs(r.raster) <- projection


crs(rehe_WU_HE_dens_raster) <- projection
bias_file <- projectRaster(rehe_WU_HE_dens_raster, r.raster)

bias_file_masked <- mask(bias_file, raster_HE_stack_complete[[22]]) #mask the bias file with wdm raster as all others too


bg_HE <- xyFromCell(bias_file_masked, sample(which(!is.na(values(bias_file_masked))), 10000, 
                                             prob=values(bias_file_masked)[!is.na(values(bias_file_masked))])) 

bg_HE_df <- bg_HE %>% 
  as.data.frame() %>% 
  rename(X = x, 
         Y = y)
#bg_sf <- st_as_sf(bg_HE_df, coords = c("X", "Y"), crs = projection)





##### SDMtune models #####

## All seasons model

# prepare an SWD object: Info: 14 presence locations are NA for some environmental variables, they are discarded! from density road layer
data_all <- prepareSWD(species = "WVC", p = rehe_WU_HE_occ, a = bg_HE_df, env = raster_HE_stack_complete, categorical = "RoadCategory")

#explore SWD object and export (if necessary)
data_all@species #ect


#train first model ####

# #one option: split into training and testing data 
# datasets <- trainValTest(data_all, test = 0.2, only_presence = TRUE)
# train <- datasets[[1]]
# test <- datasets[[2]]
# #all feature combinations possible, regularizazion multiplier equal to 1 and 500 algorithm iterations 
# default_model <- train(method = "Maxent", data = train, fc = "l") 
# default_model
# #data: object with presence/background locations used to train the model, model: model configurations
# slotNames(default_model) 


#other option: using cross validation instead of training and testing datasets: generally better 
#makes a SDMmodelCV object that hosts all the models trained during the cross validation
#k = 5 seems agreeable as the dataset is very large. 
folds = randomFolds(data_all, k = 5, only_presence = TRUE)
default_model_crossval <- train(method = "Maxent", data = data_all, fc = "l", folds = folds)


# Compute variable importance: similar results for both methods. 
#variable importance is base on the way how the variables are ordered, but permutation is updated between the models 
vi <- varImp(default_model, permut = 5)
vi_cross <- varImp(default_model_crossval, permut = 5)

plotVarImp(vi_cross) #only permutation importance (actually the interesting one)


#but from here on work with cross validated model: more robust modeling technique. 



#variable selection####
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

#vs is the model after variableselection
model_crossval_vs <- varSel(default_model_crossval, metric = "auc", bg4cor = bg_extra, cor_th = 0.7,
                            permut = 10) #takes about 45 minutes (doesnt really change after 3 permutations. )

#Removed variables: Meadows1_200m, Industry3_200m, Urbanareas4_200m, ArableAreas9_200m, 
#Urbanareas4_100m, Complex.habitats5_100m, ConiferForests6_100m, MixedForests7_100m, 
#BroadleafedForests8_100m, RoadDensity_100m
vi_cross_vs <- varImp(model_crossval_vs)


#optimize hyperparameters####

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
model_crossval_vs_optimized <- crossval_vs_optimized@models[[1]]#best model collection after optimizing process

#with cross validation: 
# fc: lqph
# reg: 0.4
# iter: 500

plot(crossval_vs_optimized) #plot the optmization process, little info in the plot

vs_optim <- varImp(model_crossval_vs_optimized) #the output is the average of the variable importance of each model trained during the cross validation
plotVarImp(vs_optim)


#reduce variation (optimize model parsimony) ####
#discards variables with low contribution

#reduce model, with permutation importance lower than 2% only of removing the variables 
#the model performance does not decrease according to the AUC
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


auc(model_allSeasons_final@models[[1]]) #the cv data, that makes an avera AUC out of the 5 models
#calculate mean and sd for all models 

#make tabe to fill with all seasons
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

#plotROC_JM(model_final, test = test_optim) 
#plots also the testing AUC, training is usually higher than testing line
#training line shows the fit of the model data, 
#testing indicates the fit of the model to the testing data: real test of the models predictive power



### response curves ###

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

#Needs now furthermore: 
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



### Model predictions ###

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



#predictionts only for the presence locations
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



## Season models ####

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



### Combined response curves for season models 

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

#function to bind all response datasets together to be able to plot them in one graph
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

#Road denisty

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


library(xtable)
model_performance
model_performance_ltx <- xtable(model_performance, caption = "Model performance")
print(model_performance_ltx,file="Figures/model_performance_ltx.tex",table.placement = "h", 
      include.rownames = TRUE, caption.placement="bottom")





#### Baden-Württemberg Preparation ####


## read in environmental variables for BW#####

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


#### ATKIS Strassen layer (metadatentabelle: ATKIS_Objektartenkatalog_Basis_DLM_7.1) ####
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
#load("WU_telemetry.RData")
Sys.setenv(lang = "en_US") #change error messages to english


# WVC validation with WU_BW


WU_BW <- read.csv2("W:/Masterarbeit Daten/Analyse/SDM/BW_Wildunfall_1.5.21-30.4.22.csv")
#summary(WU_BW)

WU_BW <- filter(WU_BW, !is.na(US_GEO_X) & !is.na(US_GEO_Y))
WU_BW <- dplyr::select(WU_BW, UN_KEY, KO_UDATUM, KO_UMONAT, KO_UZEIT, KO_USTDE, KO_WOTAG, US_GDE, US_STRA1, SS_EDATUM, SS_EZEIT, US_GEO_X, US_GEO_Y, UN_KAT, Hergangstext) #unterschied zwischen KO_UDATM und SS_EDATUM verstehen -> mathias jost

#View(WU_BW[(grep("Hirsch", WU_BW[,14])), ]) #nothing more included. a lot of names with hirsch. 


# extract rows that have Reh in the "Hergansgtext"
rehe_WU_BW <- WU_BW[(grep("Reh", WU_BW[,14])), ] #12907 rows. 
rm(WU_BW)

rehe_WU_BW <- st_as_sf(rehe_WU_BW, coords = c("US_GEO_X", "US_GEO_Y"), crs = "EPSG:4326")
rehe_WU_BW <- st_transform(rehe_WU_BW, projection)

bawu_sf <- filter(germany_sf, NAME_1 == "Baden-Württemberg")
bawu_sf$VARNAME_1 <- "BW"

rehe_WU_BW <- st_crop(rehe_WU_BW, bawu_sf) #all points within BW. 

#summary(rehe_WU_BW$UN_KAT) # unfallkategory mostly 5 as well. 

#### snap coordinates #### same as in WU_datapreparation for the WU in hessen
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

tst <- st_snap_points(rehe_WU_BW, roads_BW, max_dist = 150) #takes about 3.5 hours
rehe_WU_BW$geometry <- tst


#check how many are outside the buffer of 150 meters
roads_BW_buffered150 <- st_buffer(roads_BW, dist = 150)
outside <- sapply(st_intersects(rehe_WU_BW, roads_BW_buffered150),function(x){length(x)==0})

rehe_WU_BW_outside <- rehe_WU_BW[outside, ] #272 points
rehe_WU_BW <- rehe_WU_BW[!outside, ]

rm(roads_BW_buffered150, outside, rehe_WU_BW_outside)

st_write(rehe_WU_BW, "Layers/rehe_WU_BW.shp", append = FALSE)


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


#construct background locations just s it was done for the WU in hessen, because most points are again from accidant category 5!
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


### validation ####

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





# Crossings analysis #####

## add roe deer telem data


projection <- "+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs"

roedeer_telem <- read.csv2("Roe Deer Telemetry Data/roedeer_telem.csv")
head(roedeer_telem)

roedeer_telem$ts <- as.POSIXct(roedeer_telem$ts, format = "%Y-%m-%d %H:%M:%S")
summary(roedeer_telem$ts) #some have no timestamps. delete:
roedeer_telem <- filter(roedeer_telem, !is.na(ts))

count(roedeer_telem, gps_validity) #11 is quite bad gps accuracy. but quite low amount. 


#make sf object from dataframe
roedeer_telem_sf <- st_as_sf(roedeer_telem, coords = c("longitude", "latitude"), crs = 4326) #WGS84 system
roedeer_telem_sf <- st_transform(roedeer_telem_sf, projection) #change to right crs



## find road crossing locations


# #example for one id
# id_22 <- filter(roedeer_telem_sf, (animal_id == 22| animal_id == 19) & month(ts) == 2 & year(ts) == 2012)
# id_22_traj <- as_sftraj(id_22, group = c(id = "animal_id"), time = "ts")
# id_22_sf <- sf::st_as_sf(id_22_traj[, c(1:25, 27), drop = TRUE])
# 
# roads_crop <- st_crop(roads_BW, extent(418200, 420000, 5382000, 5390000))
# intersections22 <- st_intersection(id_22_sf, roads_crop)
# 
# plot(st_geometry(id_22_sf))
# plot(roads_crop, add = TRUE)
# plot(intersections22, add = TRUE, col = "red")
# rm(id_22, id_22_df, id_22_sf, id_22_traj)


#for entire dataset
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


## WVC as crossings 


######could there be the WVC as a crossing? -> non successful crossing!#####
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



## roe deer ID table and summaries

#understand roe deer ids and prepare table #####
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

rm(telem_min_date, telem_max_date)

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


#summary crossings over roads
summary(road_intersect_final$wdm_name)
#in relation to roads in home ranges 


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
# roedeer_ids <- dplyr::select(roedeer_ids, area, animal_id, name, sex, 
#                              telem_min_date, telem_max_date, diff_days, no_crossings, 
#                              date_of_death, cause_of_death)
# write.csv(roedeer_ids, "roedeer_ids_latex.csv")


## save new layers and Rdata

#save layers: some warnings that however do not influence the shapefile
#st_write(road_intersect_final, "Roe Deer Telemetry Data/road_intersections_final.shp")
st_write(roads_BW_HR, "Roe Deer Telemetry Data/roads_BW_HR.shp")

save.image(file='WU_telemetry.RData')



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



#### Actual crossing analysis #####

## add risk value to crossings

#road_intersect_final <- st_read("Roe Deer Telemetry Data/road_intersections_final.shp")
#st_crs(road_intersect_final) <- projection

road_intersect_occ <- road_intersect_final %>% 
  st_coordinates() %>% 
  as.data.frame()

#two points fall within now raster:check which ones
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


#only for those that actually have crossed roads
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


#analysis of risk including animal ID

  
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
colnames(WVC_crossings_summary) <- c("animal_id", "area", "wdm_name", "risk", "Season","sex", "days", "no_crossings",           "min_allCrossings", "max_allCrossings",  "mean_allCrossings", "no_samplingbuffers")

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
#save old data just t be sure 
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


#these are the two extremly high crossings. two males that 
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

#enough poits per sex and season
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

#who is the outlier? just a female with the ppoint in the middle of her home range... 
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


# Try out an SDM for crossings ####

#rasterize mcp polygons. count if multiple are on top of each other as bias file
mcp.raster <- rasterize(roedeer_mcp, r.raster_bw, 1, fun = "count")

raster_HR_wdm <- mask(raster_BW_stack_complete[[22]], roedeer_mcp_agg)
mcp.raster_masked <- mask(mcp.raster, raster_HR_wdm) #mask the bias file with wdm raster for only home range roads

#extract 10000 background points with probability of being sampled based on the bias raster 
bg.HR <- xyFromCell(mcp.raster_masked, sample(which(!is.na(values(mcp.raster_masked))), 900, 
                                              prob=values(mcp.raster_masked)[!is.na(values(mcp.raster_masked))])) 


bg.HR_df <- bg.HR %>% 
  as.data.frame() %>% 
  rename(X = x, 
         Y = y)
#bg_sf <- st_as_sf(bg_df, coords = c("x", "y"), crs = projection)



crossings_SWD_withbg <- prepareSWD(species = "WVC", p = road_intersect_occ, a = bg.HR_df, env = raster_BW_stack_complete, categorical = "RoadCategory")

folds = randomFolds(crossings_SWD_withbg, k = 5, only_presence = TRUE)
crossings_model_crossval <- train(method = "Maxent", data = crossings_SWD_withbg, fc = "l", folds = folds)


# Prepare background locations to test autocorrelation
bg_extra_HR_coords <- dismo::randomPoints(raster_BW_stack_complete, 900)
bg_extra_HR <- prepareSWD(species = "WVC", a = bg_extra_HR_coords,
                          env = raster_BW_stack_complete, categorical = "RoadCategory")

model_crossings_crossval_vs <- varSel(crossings_model_crossval, metric = "auc",
                                      bg4cor = bg_extra_HR, cor_th = 0.7, permut = 10) #7 minutes for 3 permut
# Removed variables: Swamps2_200m, Urbanareas4_200m, Complex.habitats5_200m, MixedForests7_200m, Meadows1_100m, Industry3_100m, ConiferForests6_100m, BroadleafedForests8_100m, ArableAreas9_100m, RoadDensity_200m


vi_crossings_cross_vs <- varImp(model_crossings_crossval_vs)

crossings_crossval_vs_optimized <- optimizeModel(model_crossings_crossval_vs, hypers = args, 
                                                 metric = "auc", seed = 789) #22 min

crossings_crossval_vs_optimized@results
model_crossings_crossval_vs_optimized <- crossings_crossval_vs_optimized@models[[1]]#best model collection after optimizing process
# fc: lqph
# reg: 0.6
# iter: 500
model_crossings_final <- reduceVar(model_crossings_crossval_vs_optimized, th = 2, 
                                   metric = "auc", permut = 10, use_jk = TRUE)
#Removed variables: none. 

VarImp_crossings_final <- varImp(model_crossings_final)
plotVarImp(VarImp_crossings_final)


auc_list <- c()
for(i in 1:5){
  auc_list[i] <- auc(model_crossings_final@models[[i]])
}
mean(auc_list) #0.8547028
sd(auc_list) #0.001143279



#TSS: true skill statistics (Allouche, tsoar and kadmon 20061)
tss_list <- c()
for(i in 1:5){
  tss_list[i] <- tss(model_crossings_final@models[[i]])
}
mean(tss_list) #0.5531266
sd(tss_list) #0.006355024




## response curves


#response curves still on scaled scale. needs to have the unscaled predictors again
plotResponse(model_crossings_final, "Meadows1_200m")
plotResponse(model_crossings_final, "ArableAreas9_200m")
plotResponse(model_crossings_final, "RoadCategory")

library(plotROC)
plotROC_JM(model_crossings_final@models[[1]])



predict_crossingrisk <- predict(model_crossings_final, data = raster_BW_stack_complete,
                                type = "cloglog", extent = extent(418000, 433000, 5379000, 5398500), progress = "text") 

map_crossingrisk <- plotPred(predict_crossingrisk, lt = "Probability of \nDVC occurrence",colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))



#make map of collision risk minus crossing risk
#means: aounrd 0 the roe deer cross at a certian prob and the collision risk is similar
#in negative areas: the crossing prob is higher than collision risk
#positive values: collision risk is higher than crossing prob
dis_colcross <- predict_crossingrisk + predict_allSeasons
plot(crop(dis_colcross, c(418000, 433000, 5379000, 5398500)))



save.image("WU_analysis.RData") 






