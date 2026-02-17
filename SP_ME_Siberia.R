# setwd('C:/Romanchuk_R/worms/')
# 
# library(rgbif)
# library(geodata)
# library(readxl)
# library(readr)
# library(dplyr)
# library(dismo)
# library(sf)
# library(maps)
# library(corrplot)
# library(spThin)
# library(writexl)
# 
# geodata_path('C:/Romanchuk_R/')
# 
# ########### Границы моделирования ###########
# # Определяем границы Сибири
# my_coords <- c(60.0, 48.0, 95.0, 75.0)  # xmin, ymin, xmax, ymax
# WS_vect <- vect( # Шире, чтобы точки не лежали на границе?
#   matrix(c(
#     59.0, 76.0,
#     96.0, 76.0,
#     96.0, 47.0,
#     59.0, 47.0,
#     59.0, 76.0
#   ), ncol = 2, byrow = TRUE),
#   type = "polygon",
#   crs = "EPSG:4326"
# )
# world_map <- map("world", plot = TRUE, fill = TRUE)
# world_sf <- st_as_sf(world_map, crs = 4326)
# world_vect <- vect(world_sf)
# worldC <- crop(world_vect, WS_vect)

# gadm_path <- "C:/Romanchuk_R/"
# 
# if (!file.exists(paste0(gadm_path, "gadm41_RUS_1.rds"))) {
#   RU_gadm <- gadm(country = "RUS", level = 1, path = gadm_path)
# } else {
#   RU_gadm <- readRDS(paste0(gadm_path, "gadm41_RUS_1.rds"))
# }
# 
# if (!file.exists(paste0(gadm_path, "gadm41_RUS_1.rds"))) {
#   KZ_gadm <- gadm(country = "KAZ", level = 1, path = gadm_path)
# } else {
#   KZ_gadm <- readRDS(paste0(gadm_path, "gadm41_RUS_1.rds"))
# }
# 
# if (!file.exists(paste0(gadm_path, "gadm41_RUS_1.rds"))) {
#   CN_gadm <- gadm(country = "CHN", level = 1, path = gadm_path)
# } else {
#   CN_gadm <- readRDS(paste0(gadm_path, "gadm41_RUS_1.rds"))
# }
# 
# if (!file.exists(paste0(gadm_path, "gadm41_RUS_1.rds"))) {
#   MN_gadm <- gadm(country = "MNG", level = 1, path = gadm_path)
# } else {
#   MN_gadm <- readRDS(paste0(gadm_path, "gadm41_RUS_1.rds"))
# }
# 
# RU_c <- crop(RU_gadm, worldC)
# KZ_c <- crop(KZ_gadm, worldC)
# CN_c <- crop(CN_gadm, worldC)
# MN_c <- crop(MN_gadm, worldC)
# 
# ########### Находки ###########
# # Данные из GBIF
# focalSp <- 'Aporrectodea caliginosa (Savigny, 1826)' # МЕНЯТЬ ВИД
# focalSpKey <- name_backbone(focalSp)$usageKey
# Species <- occ_search(
#   taxonKey = focalSpKey,
#   limit = 20000,
#   hasCoordinate = TRUE,
#   geometry = paste0("POLYGON((",
#                     my_coords[1], " ", my_coords[2], ",",
#                     my_coords[3], " ", my_coords[2], ",",
#                     my_coords[3], " ", my_coords[4], ",",
#                     my_coords[1], " ", my_coords[4], ",",
#                     my_coords[1], " ", my_coords[2], "))")
# )
# 
# ifelse(Species$meta$count != 0,
#        SpC <- (data.frame(longy = Species$data$decimalLongitude,
#                           laty = Species$data$decimalLatitude,
#                           Specie = "Yes")),
#        SpC <- 0)
# 
# if(nrow(SpC) != 0) {
#   SpC <- SpC[!duplicated(SpC$longy), ]
#   SpC <- SpC[!duplicated(SpC$laty), ]
#   d_SpC <- vect(SpC, geom = c('longy', 'laty'), crs = 'EPSG:4326')
# } else {
#   d_SpC <- 0
# }
# 
# 
# # наши данные
# data_geo <- read_xlsx('Our_Data.xlsx', sheet = 1)
# write_excel_csv(data_geo, 'data_geo.csv')
# data_Mx <- read.csv('data_geo.csv')
# data_Mx_S <- data.frame(longy = data_Mx$Longitude, laty = data_Mx$Latitude,
#                         Specie = data_Mx$A.c) # МЕНЯТЬ ВИД
# write.csv(data_Mx_S, 'data_Mx_S.csv')
# d_filt <- filter((data_Mx_S), Specie == 'Yes')
# d_filt <- d_filt[!duplicated(d_filt$longy), ]
# d_filt <- d_filt[!duplicated(d_filt$laty), ]
# d_Sp <- vect(d_filt, geom = c('longy', 'laty'), crs = 'EPSG: 4326')
# 
# # Объединяем
# if (nrow(d_SpC) != 0) {
#   Summ_V <- rbind(d_Sp, d_SpC)
# } else {
#   Summ_V <- d_Sp
# }
# 
# ######### Прореживаем находки ###########
# # Преобразовать обратно в data.frame для thin()
# Summ_df <- as.data.frame(Summ_V, geom = "XY")
# names(Summ_df) <- c("Specie", "longy", "laty")
# # Проредить точки с минимальным расстоянием 5 км
# thinned <- thin(
#   loc.data = Summ_df,
#   lat.col = "laty",  long.col = "longy",   spec.col = "Specie",
#   thin.par = 5,  reps = 10,  locs.thinned.list.return = TRUE,
#   write.files = FALSE,  write.log.file = FALSE)
# # Взять первый результат
# Summ_thinned <- thinned[[5]]
# names(Summ_thinned) <- c("longy", "laty")
# # Добавить столбец Вида
# Summ_thinned$Specie <- "Yes"
# # Преобразовать обратно в vect
# Summ_V <- vect(data.frame(longy = Summ_thinned$longy,
#                           laty = Summ_thinned$laty,
#                           Specie = Summ_thinned$Specie),
#                geom = c('longy', 'laty'), crs = 'EPSG:4326')
# 
# ########### Сбор предикторов ###########
# lrs <- list.files("C:/Romanchuk_R/predictors/World_Clim/Bioclim/wc2.1_2.5m/",
#                   pattern = "\\.tif$", full.names = TRUE)
# 
# predicts <- c(rast(lrs))
# predicts <- crop(predicts, WS_vect)
# 
# varsSum <- extract(predicts, Summ_V)
# varsSum$ID = NULL
# write_xlsx(varsSum, 'varsSum.csv', row.names = F)
# 
# svg("histogram_Ac.svg")
# par(mfrow = c(5, 4), mar = c(2, 2, 2, 1))
# for(col_name in names(varsSum)) {
#   hist(varsSum[[col_name]],
#        main = paste("Гистограмма", col_name),
#        xlab = col_name,
#        col = "lightblue")
# }
# dev.off()
# 
# # убираем строки с NA
# rowsNA <- apply(is.na(varsSum), 1, any)
# varsSum1 <- varsSum[!rowsNA, ]
# Summ_V1 <- Summ_V[!rowsNA, ]
# # убираем переменные с нулевыми значениями
# colSums(round(varsSum1,2))
# # убираем переменные с малым числом ненулевых значений
# colSums(varsSum1 != 0)

# ########### Отсеивание предикторов ###########
# varCorS <- round(cor(varsSum1), 3)
# colnames(varCorS) = c('BIO1','BIO2','BIO3','BIO4','BIO5',
#                       'BIO6','BIO7','BIO8','BIO9','BIO10',
#                       'BIO11','BIO12','BIO13','BIO14','BIO15',
#                       'BIO16','BIO17','BIO18','BIO19','SRTM')
# write.csv(varCorS, 'varCorS1.csv')
# # corrplot.mixed(varCorS, order = 'AOE')
# 
# # varsSum1$ID = NULL
# 
# # ищем через excel и убираем не интересующие переменные -
# # оставляем из тех, что коррелируют (-0.7, +0,7)
# # наиболее подходящие нам
# 
# varsSum1$wc_bio01 = NULL
# varsSum1$wc_bio02 = NULL
# varsSum1$wc_bio03 = NULL
# varsSum1$wc_bio04 = NULL
# varsSum1$wc_bio05 = NULL
# varsSum1$wc_bio06 = NULL
# varsSum1$wc_bio07 = NULL
# varsSum1$wc_bio08 = NULL
# # varsSum1$wc_bio09 = NULL
# # varsSum1$wc_bio10 = NULL
# # varsSum1$wc_bio11 = NULL
# varsSum1$wc_bio12 = NULL
# varsSum1$wc_bio13 = NULL
# varsSum1$wc_bio14 = NULL
# # varsSum1$wc_bio15 = NULL
# varsSum1$wc_bio16 = NULL
# # varsSum1$wc_bio17 = NULL
# varsSum1$wc_bio18 = NULL
# varsSum1$wc_bio19 = NULL
# varsSum1$wc_srtm = NULL
# 
# varCorS <- round(cor(varsSum1), 3)
# write.csv(varCorS, 'varCorS2')
# # svg("CorrPlot_Ac.svg")
# # corrplot.mixed(varCorS, order = 'AOE')
# # dev.off()
# 
# ########### MaxEnt ############
# Summ_V1$x <- geom(Summ_V1)[,3]
# Summ_V1$y <- geom(Summ_V1)[,4]
# 
# inputMxS <- Summ_V1[,c('x','y')]
# inputMxS <- as.data.frame(inputMxS)
# 
# # делим на обучающую и тестовую выборку
# inputMxS$num <- sample(1:5, nrow(inputMxS), replace = T)
# inTrainS <- inputMxS[inputMxS$num != 5,]
# inTestS  <- inputMxS[inputMxS$num == 5,]
# inTrainS$num  <- NULL
# inTestS$num  <- NULL
# 
# names(varsSum1)
# names(varsSum1) # менять столбец
# predictsMxS <- predicts[[names(varsSum1)]]
# names(predictsMxS)
# 
# # добавляем фоновые точки, в 5 раз больше
# bgpoints <- spatSample(worldC, (nrow(Summ_V1))*5)
# bgpoints$x <- geom(bgpoints)[,3]
# bgpoints$y <- geom(bgpoints)[,4]
# bgpnt <- as.data.frame(cbind(x = bgpoints$x, y = bgpoints$y))
# 
# # Запуск MaxEnt с тестовыми данными
# mxModelS <- maxent(
#   x = stack(predictsMxS),
#   p = inTrainS,
#   a = bgpnt,
#   args = c(
#     "randomtestpoints=20",         # процент тестовых данных
#     "replicates=15",               # прогонов модели
#     "responsecurves=true",
#     "jackknife=true",
#     "plots=true",
#     "pictures=true"
#   )
# )
# 
# mxModelS
# 
# MX_maps <- predict(mxModelS, stack(predictsMxS))
# MX_maps_mean <- mean(MX_maps)
# writeRaster(MX_maps_mean, 'Ac_ME_map.tif', overwrite = TRUE)
# 
# ########## Растры будущего ###########
# # Создание растров с прогнозными данными и объндиняем с высотой
# elev <- rast('C:/Romanchuk_R/predictors/World_Clim/wc_srtm_2.5m.tif')
# BCC_245_40 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/BCC_CSM2_MR/BCC_245_40.tif"), elev)
# BCC_245_60 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/BCC_CSM2_MR/BCC_245_60.tif"), elev)
# BCC_585_40 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/BCC_CSM2_MR/BCC_585_40.tif"), elev)
# BCC_585_60 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/BCC_CSM2_MR/BCC_585_60.tif"), elev)
# CMCC_245_40 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/CMCC_ESM2/CMCC_245_40.tif"), elev)
# CMCC_245_60 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/CMCC_ESM2/CMCC_245_60.tif"), elev)
# CMCC_585_40 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/CMCC_ESM2/CMCC_585_40.tif"), elev)
# CMCC_585_60 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/CMCC_ESM2/CMCC_585_60.tif"), elev)
# MPI_245_40 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/MPI_ESM1_2_HR/MPI_245_40.tif"), elev)
# MPI_245_60 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/MPI_ESM1_2_HR/MPI_245_60.tif"), elev)
# MPI_585_40 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/MPI_ESM1_2_HR/MPI_585_40.tif"), elev)
# MPI_585_60 <- c(rast("C:/Romanchuk_R/predictors/World_Clim/MPI_ESM1_2_HR/MPI_585_60.tif"), elev)
# 
# # обрезаем под сибирь
# prBCC_245_40 <- crop(BCC_245_40, WS_vect)
# prBCC_245_60 <- crop(BCC_245_60, WS_vect)
# prBCC_585_40 <- crop(BCC_585_40, WS_vect)
# prBCC_585_60 <- crop(BCC_585_60, WS_vect)
# prCMCC_245_40 <- crop(CMCC_245_40, WS_vect)
# prCMCC_245_60 <- crop(CMCC_245_60, WS_vect)
# prCMCC_585_40 <- crop(CMCC_585_40, WS_vect)
# prCMCC_585_60 <- crop(CMCC_585_60, WS_vect)
# prMPI_245_40 <- crop(MPI_245_40, WS_vect)
# prMPI_245_60 <- crop(MPI_245_60, WS_vect)
# prMPI_585_40 <- crop(MPI_585_40, WS_vect)
# prMPI_585_60 <- crop(MPI_585_60, WS_vect)
# 
# # Функция для переименования всех растров
# rename_rast <- function(rast_list) {
#   correct_names <- c(
#     "wc_bio01", "wc_bio02", "wc_bio03", "wc_bio04", "wc_bio05",
#     "wc_bio06", "wc_bio07", "wc_bio08", "wc_bio09", "wc_bio10",
#     "wc_bio11", "wc_bio12", "wc_bio13", "wc_bio14", "wc_bio15",
#     "wc_bio16", "wc_bio17", "wc_bio18", "wc_bio19", "wc_srtm"
#   )
# 
#   renamed_list <- lapply(rast_list, function(raster) {
#     names(raster) <- correct_names
#     return(raster)
#   })
# 
#   return(renamed_list)
# }
# 
# # Список всех растров будущего
# future_rast <- list(
#   prBCC_245_40 = prBCC_245_40,
#   prBCC_245_60 = prBCC_245_60,
#   prBCC_585_40 = prBCC_585_40,
#   prBCC_585_60 = prBCC_585_60,
#   prCMCC_245_40 = prCMCC_245_40,
#   prCMCC_245_60 = prCMCC_245_60,
#   prCMCC_585_40 = prCMCC_585_40,
#   prCMCC_585_60 = prCMCC_585_60,
#   prMPI_245_40 = prMPI_245_40,
#   prMPI_245_60 = prMPI_245_60,
#   prMPI_585_40 = prMPI_585_40,
#   prMPI_585_60 = prMPI_585_60
# )
# # Переименовываем все
# rast_renamed <- rename_rast(future_rast)
# # names(rast_renamed$prBCC_245_40)
# 
# ########## Для статистики ###########
# BCC1_Sum <- extract(rast_renamed$prBCC_245_40, Summ_V)
# BCC2_Sum <- extract(rast_renamed$prBCC_245_60, Summ_V)
# BCC3_Sum <- extract(rast_renamed$prBCC_585_40, Summ_V)
# BCC4_Sum <- extract(rast_renamed$prBCC_585_60, Summ_V)
# CMCC1_Sum <- extract(rast_renamed$prCMCC_245_40, Summ_V)
# CMCC2_Sum <- extract(rast_renamed$prCMCC_245_60, Summ_V)
# CMCC3_Sum <- extract(rast_renamed$prCMCC_585_40, Summ_V)
# CMCC4_Sum <- extract(rast_renamed$prCMCC_585_60, Summ_V)
# MPI1_Sum <- extract(rast_renamed$prMPI_245_40, Summ_V)
# MPI2_Sum <- extract(rast_renamed$prMPI_245_60, Summ_V)
# MPI3_Sum <- extract(rast_renamed$prMPI_585_40, Summ_V)
# MPI4_Sum <- extract(rast_renamed$prMPI_585_60, Summ_V)
# 
# vars_245_40 = rbind(BCC1_Sum, CMCC1_Sum, MPI1_Sum)
# vars_585_40 = rbind(BCC3_Sum, CMCC3_Sum, MPI3_Sum)
# vars_245_60 = rbind(BCC2_Sum, CMCC2_Sum, MPI2_Sum)
# vars_585_60 = rbind(BCC4_Sum, CMCC4_Sum, MPI4_Sum)
# stat_list <- list(
#   "Actual" = varsSum,
#   "vars_245_40" = vars_245_40,
#   "vars_585_40" = vars_585_40,
#   "vars_245_60" = vars_245_60,
#   "vars_585_60" = vars_585_60
# )
# write_xlsx(stat_list, "for_stat.xlsx")
# 
# ########## Прогноз ###########
# BCC1 <- mean(predict(mxModelS, stack(rast_renamed$prBCC_245_40)))
# BCC2 <- mean(predict(mxModelS, stack(rast_renamed$prBCC_245_60)))
# BCC3 <- mean(predict(mxModelS, stack(rast_renamed$prBCC_585_40)))
# BCC4 <- mean(predict(mxModelS, stack(rast_renamed$prBCC_585_60)))
# CMCC1 <- mean(predict(mxModelS, stack(rast_renamed$prCMCC_245_40)))
# CMCC2 <- mean(predict(mxModelS, stack(rast_renamed$prCMCC_245_60)))
# CMCC3 <- mean(predict(mxModelS, stack(rast_renamed$prCMCC_585_40)))
# CMCC4 <- mean(predict(mxModelS, stack(rast_renamed$prCMCC_585_60)))
# MPI1 <- mean(predict(mxModelS, stack(rast_renamed$prMPI_245_40)))
# MPI2 <- mean(predict(mxModelS, stack(rast_renamed$prMPI_245_60)))
# MPI3 <- mean(predict(mxModelS, stack(rast_renamed$prMPI_585_40)))
# MPI4 <- mean(predict(mxModelS, stack(rast_renamed$prMPI_585_60)))
# 
# # Сохраняем растры для дальнейшего исопользования
# # Создаем список всех растров
# rast_list <- list(
#   BCC1 = BCC1, BCC2 = BCC2, BCC3 = BCC3,BCC4 = BCC4,
#   CMCC1 = CMCC1,CMCC2 = CMCC2,CMCC3 = CMCC3,CMCC4 = CMCC4,
#   MPI1 = MPI1, MPI2 = MPI2, MPI3 = MPI3,MPI4 = MPI4
# )
# # Сохранение каждого растра отдельно
# for (name in names(rast_list)) {
#   filename <- paste0("rasters/Сибирь/Ac/", name, ".tif")
#   writeRaster(rast_list[[name]], filename, overwrite = TRUE)
#   cat("Сохранен:", filename, "\n")
# }
# 
# ########### Усредненные карты ##########
# L_prAc <- list.files("C:/Romanchuk_R/worms/rasters/Сибирь/Ac",
#                      pattern = "\\.tif$", full.names = TRUE)
# prAc <- rast(L_prAc)
# 
# # Объединение в raster stack и усредняем
# mean_245_40 <- prAc[[c("BCC1", "CMCC1", "MPI1")]]
# stack_245_40 <- stack(mean_245_40)
# map_245_40 <- mean(stack_245_40, na.rm = TRUE)
# 
# mean_585_40 <- prAc[[c("BCC3", "CMCC3", "MPI3")]]
# stack_585_40 <- stack(mean_585_40)
# map_585_40 <- mean(stack_585_40, na.rm = TRUE)
# 
# mean_245_60 <- prAc[[c("BCC2", "CMCC2", "MPI2")]]
# stack_245_60 <- stack(mean_245_60)
# map_245_60 <- mean(stack_245_60, na.rm = TRUE)
# 
# mean_585_60 <- prAc[[c("BCC4", "CMCC4", "MPI4")]]
# stack_585_60 <- stack(mean_585_60)
# map_585_60 <- mean(stack_585_60, na.rm = TRUE)

# ########### Сохранение карты ##########
# border <- function() {
# plot(RU_c, border = 'black', lwd = 0.5, col = NA, add = T)
# plot(KZ_c, border = 'blue',lwd = 0.5, col = NA,add = T)
# plot(MN_c, border = 'purple',lwd = 0.5, col = NA,add = T)
# plot(CN_c, border = 'brown',lwd = 0.5, col = NA,add = T)
# }
# 
# plot_245_40 <- function() {
#   plot(map_245_40, main = "plot_245_40",
#        col = rev(terrain.colors(5)),
#        xlab = "Долгота",
#        ylab = "Широта")
#   border()
# }
# 
# plot_585_40 <- function() {
#   plot(map_585_40, main = "plot_585_40",
#        col = rev(terrain.colors(5)),
#        xlab = "Долгота",
#        ylab = "Широта")
#   border()
# }
# 
# plot_245_60 <- function() {
#   plot(map_245_60, main = "plot_245_60",
#        col = rev(terrain.colors(5)),
#        xlab = "Долгота",
#        ylab = "Широта")
#   border()
# }
# 
# plot_585_60 <- function() {
#   plot(map_585_60, main = "plot_585_60",
#        col = rev(terrain.colors(5)),
#        xlab = "Долгота",
#        ylab = "Широта")
#   border()
# }
# 
# svg("Predicts_Mean_Ac.svg", width = 8, height = 12) # Менять название
# par(mfrow = c(2, 2))
# plot_245_40()
# plot_585_40()
# plot_245_60()
# plot_585_60()
# dev.off()
# par(mfrow = c(1, 1))
# 
# MX_maps_rast <- rast('C:/Romanchuk_R/worms/Ac_ME_map.tif')
# MX_plot <- function() {
#   plot(MX_maps_rast,
#        main = "Западная Сибирь",
# col = rev(terrain.colors(5)),
# xlab = "Долгота",
# ylab = "Широта")
#   plot(Summ_V, col = 'red', add = T)
#   border()
# }
# svg("Maxent_plot_Ac.svg", width = 6, height = 8) # Менять название
# MX_plot()
# dev.off()
#
# ######### вероятности ########
# avg_0 <- freq(MX_maps_rast > 0.6)
# avg_2_40 <- freq(map_245_40 > 0.6, useNA = "no")
# avg_2_60 <- freq(map_245_60 > 0.6, useNA = "no")
# avg_5_40 <- freq(map_585_40 > 0.6,  useNA = "no")
# avg_5_60 <- freq(map_585_60 > 0.6,  useNA = "no")
# low_0 <- freq(MX_maps_rast > 0.2)
# low_2_40 <- freq(map_245_40 > 0.2, useNA = "no")
# low_2_60 <- freq(map_245_60 > 0.2, useNA = "no")
# low_5_40 <- freq(map_585_40 > 0.2,  useNA = "no")
# low_5_60 <- freq(map_585_60 > 0.2,  useNA = "no")
# L_hal <- list(
#   Actual_a = avg_0[2, 3],
#   a245_40 = avg_2_40[2, 2],
#   a245_60 = avg_2_60[2, 2],
#   a585_40 = avg_5_40[2, 2],
#   a585_60 = avg_5_60[2, 2],
#   Actual_l = low_0[2, 3],
#   l245_40 = low_2_40[2, 2],
#   l245_60 = low_2_60[2, 2],
#   l585_40 = low_5_40[2, 2],
#   l585_60 = low_5_60[2, 2]
# )
# Tabl2 <- cbind(L_hal)
#
# ########### сравнение растров #########
# vals_act <- values(MX_maps_rast)
# vals_f240 <- values(map_245_40)
# vals_f540 <- values(map_585_40)
# vals_f260 <- values(map_245_60)
# vals_f560 <- values(map_585_60)
# 
# # Удаление NA
# vals_act <- na.omit(vals_act)
# vals_f240 <- na.omit(vals_f240)
# vals_f540 <- na.omit(vals_f540)
# vals_f260 <- na.omit(vals_f260)
# vals_f560 <- na.omit(vals_f560)
# 
# df_anova <- data.frame(
#   value = c(vals_act, vals_f240, vals_f540, vals_f260, vals_f560),
#   group = factor(rep(c('act', 'f240', 'f540', 'f260', 'f560'),
#                      times=c(length(vals_act), length(vals_f240), length(vals_f540),
#                              length(vals_f260), length(vals_f560))))
# )
# 
# # Однофакторный ANOVA
# anova_result <- aov(value ~ group, data=df_anova)
# summary(anova_result)
# 
# # Пост-хок тест Тьюки
# tukey <- TukeyHSD(anova_result)
#
# ########## Разница растров ##########
# MX_maps_raster <- raster(MX_maps_rast)
# act_f240 <- map_245_40 - MX_maps_raster
# act_f540 <- map_585_40 - MX_maps_raster
# act_f260 <- map_245_60 - MX_maps_raster
# act_f560 <- map_585_60 - MX_maps_raster
# 
# svg("plot_разница_new.svg", width = 8, height = 12) # Менять название
# par(mfrow = c(2, 2)) #, mar = c(1.5, 1.5, 1.5, 0.5))
# plot(act_f240, main = "act_f240",
#      col = colorRampPalette(c("red", "lightgrey", "blue"))(9),  # красный-серый-синий
#      xlab = "Долгота",
#      ylab = "Широта",
#      zlim = c(-0.5, 0.5)
# )
# plot(act_f540, main = "act_f540",
#      col = colorRampPalette(c("red", "lightgrey", "blue"))(9),
#      xlab = "Долгота",
#      ylab = "Широта",
#      zlim = c(-0.5, 0.5)
# )
# plot(act_f260, main = "act_f260",
#      col = colorRampPalette(c("red", "lightgrey", "blue"))(9),  # красный-серый-синий
#      xlab = "Долгота",
#      ylab = "Широта",
#      zlim = c(-0.5, 0.5)
# )
# plot(act_f560, main = "act_f560",
#      col = colorRampPalette(c("red", "lightgrey", "blue"))(9),  # красный-серый-синий
#      xlab = "Долгота",
#      ylab = "Широта",
#      zlim = c(-0.5, 0.5)
# )
# dev.off()
# par(mfrow = c(1, 1))
# 