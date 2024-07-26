install.packages("tidyverse")


# Paquetes ----------------------------------------------------------------
library(sf)
library(tidyverse)
library(terra)
library(tmap)
library(tmaptools)
library(kuenm)

library(terra)
library(flexsdm)
library(stringr)

# Mapa base para chequeo de capas -----------------------------------------
# Cargo limites políticos de Sudamerica
library(rnaturalearth)
library(rnaturalearthdata)
# Países
Sudam_raw <- ne_countries(scale = 10, 
                          continent = "south america", 
                          returnclass = "sf")

French_guiana <- ne_countries(scale = 10,
                              type = "map_units",
                              geounit = "French Guiana",
                              returnclass = "sf")
Sudam <- st_union(Sudam_raw, French_guiana)


bb <- c(-83, # longitud mínima 
        -30, # latitud mínima
        -33, # longitud máxima
        12) # latitud máxima


# Linea de costa del LGM
LGMcoast <- st_read("C:/Users/ale_f/OneDrive/Documentos/Capas_GIS/Capas_Sudam/lgm_LINEA_DE_COSTA/LGMcoastline_RayAdams_2001.gpkg")


mapa_base <- tm_shape(Sudam, bbox = bb) + tm_borders(lwd = 2) +
  tmap_options(check.and.fix = TRUE)



# 0. Working directory ----------------------------------------------------
setwd ("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm")


#1 PREPARO DATOS DE OCURRENCIA ------
## 1.1. Partición de los datos ---------------------------------------------

# Para poder evaluar los modelos necesito datos independientes que me permitan
# evaluar que tan bien funcionan. Como no tengo fuentes independientes voy
# a particionar los datos de modo que calibre los modelos con una parte y con
# la otra los evalúo. Usaré la función part_sblock() delpaquete flexsdm para
# particionar en bloques espaciales los datos.

# Previamente debo agregar una columna que identifique las presencias -----
Occs_pres <- Occs_t |> 
  mutate(pr_ab = rep("1", 1963))   #mutate=agrega o quita columnas


#cargo capa raster para usar como env_layer

pc <- rast("M_variables/set_1/PC1.asc")


# Corro la función para seleccionar la mejor partición ----
set.seed(10)

occ_part <- part_sblock(
  data = Occs_pres,
  env_layer = pc,
  pr_ab = "pr_ab",
  x = "longitude",
  y = "latitude",
  n_part = 4,
  min_res_mult = 3,
  max_res_mult = 50,
  num_grids = 30,
  min_occ = 50,
  prop = 0.5 #reviser!!!! puncher
)

Occs_p <- occ_part$part

occ_part$grid

# Transform best block partition to a raster layer with same resolution and extent than
# predictor variables
block_layer <- get_block(env_layer = pc, best_grid = occ_part$grid)


# Plot best spatial blocking partitions ----
cl <- c("#64146D", "#9E2962", "#F47C15", "#FCFFA4")
plot(block_layer, col=cl, legend=FALSE, axes=FALSE)
points(Occs_p[,c("x", "y")])

# Number of presences per block
Occs_p %>%
  dplyr::group_by(.part) %>%
  dplyr::count()

#creamos la tabla de datos occ_train.csv ----

occ_train_block <- Occs_p %>% 
  filter(!.part=="1") %>% #selecciono de la columna .part las filas con el valor deseado. (!) hace lo contrario
  mutate(species = rep("nicotiana_glauca","1449"), #mutate crea columnas
         .before = "x") %>% 
  rename(longitude="x", #rename asigno nombre a cada columna
         latitude="y") %>% 
  select(1:3) #me quedo con las 3 primeras columnas

#guardo csv de occ_train_block
write.csv(occ_train_block, 
          file = "Occs_train_block.csv",
          row.names = F)

# 1.2 Preparing data  # otra forma de preparar datos pero esta vez a partir de randon¿m----
## splitting training and testing records 
help(kuenm_occsplit)

Occs_t <- read.csv("palan_kuenm.csv") 
set.seed(1234) # Para reproducibilidad
split <- kuenm_occsplit(occ = Occs_t, train.proportion = 0.6, method = "random",
                        save = TRUE, name = "Occs")





#1.3 calib_area()  generar área M ------
pc <- rast("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/current/PC1.tif")


area_M <- calib_area(data= Occs_t, x= "longitude", y="latitude", 
                     method =  c('buffer', width=200000),
                     crs = crs(pc))

plot(area_M)

#Guardo el vector con el area M
area_M_sf <- st_as_sf(area_M)
st_write(area_M_sf,
         "M_palan.gpkg")


plot(area_M_sf)

# Cortar raster para calibrar -----

#Genero un objeto con el path hacia las capas que quiero cargar
paths_capas <- list.files (path = "E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/current",
  pattern = ".tif$",full.names = TRUE)

#Genero un stack a partir de los archivos del objeto path - raster apilados
pc <- raster::stack(paths_capas)


#!!!!!Cargo las mascara que usaremos
Mask = area_M_sf

#El crop me recorta los rasters a partir de un shape. Es recomendable antes de usar la funcion mask.
StackCropped <- raster::crop(x = pc,y = Mask)

#Ahora uso la funcion mask para cortar el stack cropeado, nuevamente usando el shapefile de mi mascara
StackMasked <- raster::mask(StackCropped, Mask)


plot(StackMasked)


#Creo una carpeta para guardar mis raster dentro de mi directorio
dir.create("M_variables/set_1")

#Creo los rasters en mi carpeta
lapply(names(StackMasked), function(x){
  writeRaster(StackMasked[[x]], paste0("M_variables/set_1/",
                                       x,".asc"),overwrite=TRUE)})



# 2 Model calibration --------------------------------------------------
## 2.1 candidate model creation ----- 
help(kuenm_cal)

oj <- "Occs_joint.csv"
otr <- "Occs_train.csv"
mvars <- "M_variables"
bcal <- "batch_cal1"
candir <- "Candidate_models"
regm <- c(0.5, 2)
mxpath <- "C:/maxent"

kuenm_cal(occ.joint = oj, occ.tra = otr, M.var.dir = mvars, 
          batch = bcal, out.dir = candir, max.memory = 1000, 
          reg.mult = regm, f.clas = "basic", args = NULL, 
          maxent.path = mxpath, wait = FALSE, run = TRUE)


## 2.2 evaluation of candidate models ------
help(kuenm_ceval)

ote <- "Occs_test.csv"
cresdir <- "Calibration_results1"

kuenm_ceval(path = candir, occ.joint = oj, occ.tra = otr, 
            occ.test = ote, batch = bcal, out.eval = cresdir,
            threshold = 5, rand.percent = 50, iterations = 500,
            kept = F, selection = "OR_AICc", parallel.proc = FALSE)


# 3 Model projections --------------------------------------------------
# 3.1 Create G varibles ----------------


paths_capas_ipsl245 <- list.files (path = "E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/IPSL-CM6A-LR_ssp245",
                           pattern = ".tif$",full.names = TRUE)


IPSL_126 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/IPSL-CM6A-LR_ssp126/pcs.tif")

IPSL_245 <- raster::stack(paths_capas_ipsl245)

IPSL_370 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/IPSL-CM6A-LR_ssp370/pcs.tif")

IPSL_585 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/IPSL-CM6A-LR_ssp585/pcs.tif")


MIROC_126 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/MIROC6_ssp126/pcs.tif")

MIROC_245 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/MIROC6_ssp245/pcs.tif")

MIROC_370 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/MIROC6_ssp370/pcs.tif")

MIROC_585 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/MIROC6_ssp585/pcs.tif")


BCC_CSM2_126 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/BCC-CSM2-MR_ssp126/pcs.tif")

BCC_CSM2_245 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/BCC-CSM2-MR_ssp245/pcs.tif")

BCC_CSM2_370 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/BCC-CSM2-MR_ssp370/pcs.tif")

BCC_CSM2_585 <- raster::stack("E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/BCC-CSM2-MR_ssp585/pcs.tif")





paths_CurrentBios <- list.files(path = "E:/DOCTORADO/MNE/MNE/2024/1_MNE PALAN_COMPLETO/TODO_JUNTO/kuenm/current",
                                pattern = ".tif$",full.names = TRUE)

CurrentBios <- raster::stack (paths_CurrentBios)


#Creo una carpeta para guardar mis raster dentro de mi directorio
dir.create("G_variables/set_1")

dir.create("G_variables/set_1/IPSL_126")
dir.create("G_variables/set_1/IPSL_245")
dir.create("G_variables/set_1/IPSL_370")
dir.create("G_variables/set_1/IPSL_585")

dir.create("G_variables/set_1/MIROC_126")
dir.create("G_variables/set_1/MIROC_245")
dir.create("G_variables/set_1/MIROC_370")
dir.create("G_variables/set_1/MIROC_585")

dir.create("G_variables/set_1/BCC_CSM2_126")
dir.create("G_variables/set_1/BCC_CSM2_245")
dir.create("G_variables/set_1/BCC_CSM2_370")
dir.create("G_variables/set_1/BCC_CSM2_585")

dir.create("G_variables/set_1/current")


#Creo los rasters en mi carpeta
lapply(names(IPSL_126), function(x){
  writeRaster(IPSL_126[[x]], paste0("G_variables/set_1/IPSL_126/",
                                    x,".asc"),overwrite=TRUE)})
lapply(names(IPSL_245), function(x){
  writeRaster(IPSL_245[[x]], paste0("G_variables/set_1/IPSL_245/",
                                       x,".asc"),overwrite=TRUE)})
lapply(names(IPSL_370), function(x){
  writeRaster(IPSL_370[[x]], paste0("G_variables/set_1/IPSL_370/",
                                    x,".asc"),overwrite=TRUE)})
lapply(names(IPSL_585), function(x){
  writeRaster(IPSL_585[[x]], paste0("G_variables/set_1/IPSL_585/",
                                    x,".asc"),overwrite=TRUE)})


lapply(names(MIROC_126), function(x){
  writeRaster(MIROC_126[[x]], paste0("G_variables/set_1/MIROC_126/",
                                     x,".asc"),overwrite=TRUE)})
lapply(names(MIROC_245), function(x){
  writeRaster(MIROC_245[[x]], paste0("G_variables/set_1/MIROC_245/",
                                    x,".asc"),overwrite=TRUE)})
lapply(names(MIROC_370), function(x){
  writeRaster(MIROC_370[[x]], paste0("G_variables/set_1/MIROC_370/",
                                     x,".asc"),overwrite=TRUE)})
lapply(names(MIROC_585), function(x){
  writeRaster(MIROC_585[[x]], paste0("G_variables/set_1/MIROC_585/",
                                     x,".asc"),overwrite=TRUE)})


lapply(names(BCC_CSM2_126), function(x){
  writeRaster(BCC_CSM2_126[[x]], paste0("G_variables/set_1/BCC_CSM2_126/",
                                     x,".asc"),overwrite=TRUE)})
lapply(names(BCC_CSM2_245), function(x){
  writeRaster(BCC_CSM2_245[[x]], paste0("G_variables/set_1/BCC_CSM2_245/",
                                        x,".asc"),overwrite=TRUE)})
lapply(names(BCC_CSM2_370), function(x){
  writeRaster(BCC_CSM2_370[[x]], paste0("G_variables/set_1/BCC_CSM2_370/",
                                        x,".asc"),overwrite=TRUE)})
lapply(names(BCC_CSM2_585), function(x){
  writeRaster(BCC_CSM2_585[[x]], paste0("G_variables/set_1/BCC_CSM2_585/",
                                        x,".asc"),overwrite=TRUE)})


lapply(names(CurrentBios), function(x){
  writeRaster(CurrentBios[[x]], paste0("G_variables/set_1/current/",
                                     x,".asc"),overwrite=TRUE)})



#3.2 creo modelos con proyeccion ------------

help(kuenm_mod)

bfmod <- "batch_model"
moddir <- "Final_models1"
gvars <- "G_variables"


kuenm_mod(occ.joint = oj, M.var.dir = mvars, out.eval = cresdir, 
          batch = bfmod, rep.n = 2, rep.type = "Bootstrap",
          jackknife = TRUE, out.dir = moddir, max.memory = 1000, 
          out.format = "cloglog", project = TRUE, G.var.dir = gvars, 
          ext.type = "ext_clam", write.mess = TRUE, write.clamp = FALSE, 
          maxent.path = mxpath, args = NULL, wait = FALSE, run = TRUE)





# 3.3 Descriptive statistics of results --------------------------------------------------------------------------
help(kuenm_modstats)

spname <- "Nicotiana_glauca"
modstats <- "2_Final_Model_Stats"
proj_scenarios <- c("current", 
                    "IPSL_245", "IPSL_370", "IPSL_585",
                    "MIROC_245", "MIROC_370", "MIROC_585")

kuenm_modstats(sp.name = spname, fmod.dir = moddir, format = "asc", 
               project = TRUE, 
               statistics = c("med", "range","mean"), 
               replicated = TRUE,
               proj.scenarios = proj_scenarios, 
               ext.type = "EC", out.dir = modstats)

# 3.4 Umbralizacion----------------------------------------------------------

# Calculo Umbral para definir limites de la distribucion
# Cargo proyeccion del modelo sobre variables de calibracion, y ocurrencias de calibracion
SDMPresent = rast("Final_Model_Stats1/Statistics_EC/Present_med.tif")
Occs_LatLong = Occs_t[,-1]
# Uso funcion que robe de "https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/"
# pero la modifico para primero calcular el valor de thrshold que quiero:
# Esta funcion tiene para calcular tres tipos de thrEsholds: MTP, P5 Y P10.
thresh = function(sdm, occs, type = "mtp"){
  sdm1 <- raster::raster(sdm)
  occPredVals <- raster::extract(sdm1, occs) # Extraigo valores de idoneidad de presencias
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p5"){
    if(length(occPredVals) < 10){
      p5 <- floor(length(occPredVals) * 0.95)
    } else {
      p5 <- ceiling(length(occPredVals) * 0.95)
    }
    thresh <- rev(sort(occPredVals))[p5]
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  return(thresh)
}

# Calculo el umbral que quiero
Umbral_mtp = thresh(SDMPresent, Occs_LatLong, type = "mtp")
Umbral_p5 = thresh(SDMPresent, Occs_LatLong, type = "p5")
Umbral_p10 = thresh(SDMPresent, Occs_LatLong, type = "p10")


# Funcion para binarizar los mapas de acuerdo al umbral calculado previamente
sdm_threshold <- function(sdm, Umbral, binary = FALSE){
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < Umbral] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}

SDMPresent_bin = sdm_threshold(SDMPresent, Umbral = Umbral_p5, binary = F)

SDMMH = rast("Final_Model_Stats1/Statistics_EC/MH_med.tif")
SDMMH_bin = sdm_threshold(SDMMH, Umbral = Umbral_p5, binary = F)

SDMLGM = rast("Final_Model_Stats1/Statistics_EC/LGM_med.tif")
SDMLGM_bin = sdm_threshold(SDMLGM, Umbral = Umbral_p5, binary = F)

SDMLIG = rast("Final_Model_Stats1/Statistics_EC/LIG_med.tif")
SDMLIG_bin = sdm_threshold(SDMLIG, Umbral = Umbral_p5, binary = F)





# 4 Mapas -------------------------------------------------------------
# vuelvo al directorio de trabajo original del proyecto
setwd("..")

Occs_sf = Occs_t |>  
  st_as_sf(coords = c("Longitude", "Latitude")) %>% 
  st_set_crs(4326)
# Cambiarle la leyenda, y que comience desdel el umbral,

## 4.1. Presente ---------------------------------------------------------
SDM_Present_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(SDMPresent_bin) + tm_raster(style = "cont",
                                       palette = viridisLite::viridis(20, direction = -1),
                                       title= "Suitability",
                                       legend.format =list(text.separator="-"),
                                       legend.reverse = T) +
  tm_shape(Occs_t_sf) + tm_symbols(size = 0.1, col = "red", alpha = 0.8) +
  mapa_base +
  tm_layout(title = "Present",
            title.position = c("left", "top"),
            legend.position = c("right", "bottom"))

tmap_save(SDM_Present_map, "plots/ENM_Present.pdf")
## 4.2 MH --------------------------------------------------------------
SDM_MH_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(SDMMH_bin) + tm_raster(style = "cont",
                                  palette = viridisLite::viridis(20, direction = -1),
                                  title= "Suitability",
                                  legend.format =list(text.separator="-"),
                                  legend.reverse = T) +
  mapa_base +
  tm_layout(title = "MH(6 ka)",
            title.position = c("left", "top"),
            legend.position = c("right", "bottom"))

tmap_save(SDM_MH_map, "plots/ENM_MidHolocene.pdf")

## 4.3 LGM --------------------------------------------------------------
SDM_LGM_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(SDMLGM_bin) + tm_raster(style = "cont",
                                   palette = viridisLite::viridis(20, direction = -1),
                                   title= "Suitability",
                                   legend.format =list(text.separator="-"),
                                   legend.reverse = T) +
  mapa_base +
  tm_shape(LGMcoast) + tm_borders(lty = "dashed") +
  tm_layout(title = "LGM(20 ka)",
            title.position = c("left", "top"),
            legend.position = c("right", "bottom"))

tmap_save(SDM_LGM_map, "plots/ENM_LGM.pdf")

## 4.4 LIG --------------------------------------------------------------
SDM_LIG_map <-  tm_shape(Sudam, bbox = bb) + tm_polygons(col = "lightgrey") +
  tm_shape(SDMLIG_bin) + tm_raster(style = "cont",
                                   palette = viridisLite::viridis(20, direction = -1),
                                   title= "Suitability",
                                   legend.format =list(text.separator="-"),
                                   legend.reverse = T) +
  mapa_base +
  tm_layout(title = "LIG(130 ka)",
            title.position = c("left", "top"),
            legend.position = c("right", "bottom"))

tmap_save(SDM_LIG_map, "plots/ENM_LIG.pdf")
