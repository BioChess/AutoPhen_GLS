#### 0.1 Working directory ####
Sys.setenv(TZ = 'GMT') ### !!! IMPORTANT: IF NOT SET LIKE THIS, MAKE PROBLEMS 
#TO CREATE DATE_TIME FROM PASTING DATE plus TIME

place <- 'PCUB'  # pc in the UB
# place <- "home"  # laptop at home


if (place == "PCUB") WD <- 'D:/Dropbox/'
if (place == 'home') WD <- 'F:/Dropbox/'  # pc in the UB
setwd(paste0(WD,'Diego Vicente/TESIS/01_GLS_analysis'))

wd<-getwd()

#### 0.2 Install packages ####
# Installing all the packages and load their libraries we need to run the script

pacman::p_load('devtools',
               'tidyverse', 'plyr', 'lubridate', 'readxl', 'readr', # data carpenter
               'sp', 'sf', 'spData', 'trip', # spatial data
               'adehabitatHR', 'track2KBA', # migrate behaviour and kernels
               'ggplot2', 'ggpubr', 'ggspatial',# plots
               'patchwork', 'cowplot', # composing plots
               'caTools'# runmean and runmin calcs
)

source(paste0(wd,'/RCode/Migratory patterns/06_FUNCTIONS_automated_phenology_v9.5b.R'))
source(paste0(wd,'/RCode/Migratory patterns/06B_FUNCTIONS_Manual_interactive_reclassification_v1.R'))
source(paste0(wd,'/RCode/Migratory patterns/06C_FUNCTIONS_Phenology_dates_v2.R'))
#### 0.3 Spatial environment ####
world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
load(paste0(wd, '/Raw Data/Projections.Rdata')) # projections
projection <- projections$WGS84
sf_use_s2(FALSE)

#### 1. Load files  and information ####
GLS.files <- read_csv('Raw Data/GLS.files3.csv', col_names = TRUE)
colonies <- read_excel('Raw Data/Colony_coordinates.xlsx')
equinox<-read_excel('Raw Data/Equinox_by_year.xlsx')
colony_days <- read_csv('Raw Data/colony_days_75_25_v2.csv')

tracks_phen <- read_csv('Results/Phenology/Carryover_effects/CALspp_phenology_tracks_corrected.csv')
win_areas <- read_sf(dsn = "Raw Data/wintering_areas")
ecoregions <- win_areas %>% dplyr::select(ECOREGION)

ids <-unique(tracks_phen$ID)


track.list <- list()
phen.list <- list()
pre_env<-ls()
for (i in 1:length(ids)) {
  # i = which(ids == 'BH061001001')
  # i = 14
  track.ph <- tracks_phen %>%
    filter(ID == ids[i])
  ####2. Track information #####
  # track.ph <- tracks_phen %>% filter(ID == ids[i])
  info<- GLS.files%>%dplyr::filter(geo_trip == ids[i])
  geoID <- info$geo_alfanum
  geo_trip <- info$geo_trip
  Ring<- info$Ring #Bird ID
  Sp <- info$Sp
  resucitat<-info$Resucitat
  
  cat('Processing geo ', geo_trip, '(', i, ') \n')
  
  start_year <- year(track.ph$date[1])
  end_year <- year(track.ph$date[nrow(track.ph)]) #Trip end
  # Colony info
  col.name <- info$Colony
  
  if(col.name == 'M.Clara') col.name <- 'MClara'
  if(col.name == 'Cala Morell') col.name <- 'CalaMorell'
  if(col.name == 'Curral Velho') col.name <- 'CVelho'
  if(col.name == 'Ilheu de Cima') col.name <- 'ICima'
  
  colony<- colonies%>%
    filter(Colony == col.name)%>%
    dplyr::select(Longitude = Lon, Latitude = Lat)
  col.lon <- colony$Longitude
  col.lat <- colony$Latitude
  colony_sf<-colony%>%
    st_as_sf(coords = c("Longitude", "Latitude"),
             crs = projections$WGS84, agr = "constant")
  
  
  #### 2.1. Correct possible manual mistakes ####
  
  track.ph2 <- check_phen(track.ph, "Breed2")
  
  track.ph2$Breed3 <- phen_periods(track.ph2$Breed2)
  
  #### 2.2. Identify unique B and W stages ####
  BW_periods <- unique(track.ph2$Breed3[grep("^B|^W", track.ph2$Breed3)])
  B_periods <- grep("^B", BW_periods, value = TRUE)
  W_periods <- grep("^W", BW_periods, value = TRUE)
  
  if(B_periods[1] == 'B1' & track.ph2$Breed2[1] != 'B'){
    Bperiods <- 'B2'
    track.ph2 <- track.ph2 %>%
      mutate(Breed3 = case_when(Breed3 == 'B1' ~ 'B2',
                                T ~ Breed3))
  }
  if(length(W_periods) == 0){
    Wperiods <- NULL
  }
  
  
  df_format <- formatFields(dataGroup = track.ph2, 
                            fieldID   = "ID", 
                            fieldDateTime = "date", 
                            fieldLon  = "x", 
                            fieldLat  = "y" )
  
  df_sf <- df_format%>%
    st_as_sf(coords = c("Longitude", "Latitude"),
             crs = projections$WGS84, agr = "constant", remove = FALSE)
  
  ####3. Phenology plot 1  #####
  breed_df <- track.ph2 %>%
    mutate(
      Shape = case_when(Equinox == 'Sep' ~ 24,
                        Equinox == 'Mar' ~ 25,
                        Equinox == 'No' ~ 21
      ),
      Color = case_when(Breed2 == 'B' ~ 'orange2',
                        Breed2 == 'mig1' ~ 'grey40',
                        Breed2 == 'mig2' ~ 'grey70',
                        Breed2 == 'W' ~ 'steelblue',
      ),
      lab = case_when(Breed2 == 'B' ~ 'Breeding',
                      Breed2 == 'mig1' ~ 'Migration post',
                      Breed2 == 'mig2' ~ 'Migration pre',
                      Breed2 == 'W' ~ 'Wintering',
      )
    )
  breaks <- unique(breed_df$Color)
  labels <- unique(breed_df$lab)
  
  # Breeding vs non-breeding plot
  breed_plot <- ggplot()+
    # World
    {if(!is.null(world))
      geom_sf(data = world, fill=grey(0.9), color=grey(0.6), lwd = 0.2)
    } +
    
    coord_sf(xlim = c(min(breed_df$x), max(breed_df$x)), 
             ylim = c(min(breed_df$y), max(breed_df$y)), expand = TRUE) +
    # Trip
    geom_path(data= breed_df, aes(x, y))+
    {if(length(unique(breed_df$Shape)) > 1)
      geom_point(data= breed_df, aes(x, y, fill = Color,
                                     shape = Shape), size = 3)
    } +
    {if(length(unique(breed_df$Shape)) == 1)
      geom_point(data= breed_df, aes(x, y, fill = Color),
                 shape = 21, size = 3)
    }+
    # Colors
    scale_fill_identity(guide = 'legend',
                        name = 'Track classification',
                        breaks = breaks,
                        labels = labels
    )+
    # Shapes
    scale_shape_identity(breaks = c(21, 25, 22),
                         labels = c('No', 'March', 'September'),
                         name = 'Equinox',
                         guide = 'legend')+
    # Labels
    labs(title = 'Resultating segmentation')+
    xlab('Longitude') + ylab('Latitude')+
    # Correct legends colors
    guides(fill = guide_legend(override.aes = list(fill = breaks,
                                                   shape = 21)))+
    #Theme
    my_theme
  
  ggsave(plot = breed_plot,
         filename = paste0(wd,'/Plots/Phenology/Carryover_effects/Segclust_corrected/', 
                           geo_trip, '_', end_year,'_', Sp, '_', col.name, '_',
                           'segclust_corrected.jpg'),
         width = 15, height = 18, units = 'cm')
  
  
  cat('First plot saved \n')
  
  #### 4. Assign wintering areas ####
  
  df_win <- st_join(df_sf, ecoregions, left = FALSE)
  
  WAs <- df_win %>%
    filter(str_detect(Breed3, "^W\\d+$")) %>%  # Seleccionar filas donde Breed3 es W1, W2, W3, etc.
    group_by(Breed3) %>%  # Agrupar por Breed3 (W1, W2, W3, ...)
    dplyr::summarise(
      Ecoregion = names(which.max(table(ECOREGION)))  # Calcular la Ecoregion con más frecuencia
    ) %>%
    right_join(
      tibble(Breed3 = unique(df_win$Breed3[str_detect(df_win$Breed3, "^W\\d+$")])), 
      by = "Breed3"
    )
  
  if(nrow(WAs) == 0){
    WAs <- ecoregions %>%
      filter(ECOREGION == 'Saharan_Upwelling')%>%
      add_column(Breed3 = 'W1')%>%
      dplyr::select(Breed3, Ecoregion = ECOREGION)
   
  }
  
  ####5. Resident vs migrant classification  #####
  
  if(all(df_format$Breed3 == 'B1')){
    Migrant <- 'Resident'
  } else if(!any(df_format$Breed3 %in% c('mig1', 'mig2'))){
    Migrant <- 'Resident2'
  } else if(any(df_format$Breed3 == 'mig1') & all(df_format$Breed2 != 'W')){
    Migrant <- 'Resident_mig'
  } else if(!any(df_format$Breed3 == 'mig1')){
    Migrant <- 'Resident2_mig'
  } else { 
    Migrant <- 'Migrant'
  }
  
  if(geo_trip %in% c('X462001001')) Migrant <- 'Resident2'
  
  ####6. KDE phenology determination  #####
  
  ph_seg.list <- subset_to_kde(df_format, "Breed3")
  
  # Main processing loop
  results_list <- list()
  for (k in names(ph_seg.list)) {
    # k = 'B1'
    print(k)
    ph_seg <- ph_seg.list[[k]]
    
    rle_breed3 <- rle(df_format[['Breed3']])
    
    if(Migrant == 'Resident2'){
      ph_seg <- ph_seg %>%
        filter(Breed3 == k)
    }
    if(length(rle_breed3[2]) > 1){
      if(rle_breed3$values[2] == 'W1' & rle_breed3$values[2] != 'B2'){
        if(k == 'B1'){
          ph_seg <- ph_seg %>%
            filter(Breed3 == k)
        } else if (k == 'W1'){
          ph_seg <- ph_seg %>%
            filter(Breed3 != 'B1')
        } else if (!(k %in% c('B1', 'W1'))) {
          ph_seg <- ph_seg
        }
      }
    }

    
    
    stg <- ph_seg[!(ph_seg$Breed3 %in% c("mig1", "mig2")), ]
    
    if (is.null(ph_seg) || nrow(stg) == 0) {
      # If no data, initialize as NA
      results_list[[k]] <- list(KDE90_sf = NULL,
                                arr_point = data.frame(Longitude = NA, Latitude = NA), arr_date = NA,
                                dep_point = data.frame(Longitude = NA, Latitude = NA), dep_date = NA
                                )
      next
    }
    
    # Calculate KDE for non-migratory periods
    KDE90_sf <- calculate_kde(ph_seg, proj = projections$WGS84)
    
    # Find first and last points inside KDE
    phen_dates <- find_dates(ph_seg, KDE90_sf, proj = projections$WGS84)
    
    # Store results
    results_list[[k]] <- c(list(KDE90_sf = KDE90_sf), phen_dates)
  }

  if(Migrant == 'Resident'){
    results_list[["B1"]][["dep_date"]] <- as.Date(NA)
  } else if(Migrant == 'Resident2'){
    results_list[["B1"]][["dep_date"]] <- results_list[["W1"]][["arr_date"]]
    
    results_list[["W1"]][["dep_date"]] <- results_list[["B1"]][["arr_date"]]
  }
  
  
  ####7. Phenology plot 2  #####
  
  map_trip <- plot_phen_map(df_format, colony, results_list, world, my_theme)
  
  ggsave(plot = map_trip,
         filename = paste0(wd,'/Plots/Phenology/Carryover_effects/', 
                           geo_trip, '_', end_year,'_', Sp, '_', col.name, '_',
                           'phenology.jpg'),
         width = 15, height = 18, units = 'cm')
  
  
  cat('Second plot saved \n')
  
  ####8. Saving phenology dates  #####

  phen_dates <- extract_phen_dates(df_format)

  if (Migrant %in% c("Resident2", "Resident2_mig")) {
    phen_dates <- phen_dates %>%
      mutate(B1_dep = W1_arr)%>%
      relocate(B1_dep, .before = W1_arr)
    
  }
  
  
  
  phen <- data.frame(Ring = Ring, Sp = Sp, Colony = col.name,
                     ID = geo_trip, Migrant = Migrant, 
                     Year1 = start_year, Year2 = end_year,
                     phen_dates
  )
  
  phen_KDEdates <- extract_phen_KDEdates(results_list)
  
  
  phen.df <- phen %>%
    bind_cols(phen_KDEdates)%>%
    mutate(
      !!!setNames(
        lapply(seq_along(WAs$Breed3), function(i) {
          WAs$Ecoregion[match(paste0("W", i), WAs$Breed3)]
        }),
        paste0("WA", seq_along(WAs$Breed3))
      )
    )%>%
    mutate(
      across(
        ends_with("KDE"), 
        ~ coalesce(., get(sub("_KDE$", "", cur_column())))
      )
    )
  
  
  
  phen.list[i] <-list(phen.df)
  post_env<-ls()
  post_env <- post_env[!(post_env %in% c('pre_env'))]
  
  if(i == length(ids)){break}else{rm(list=setdiff(post_env, pre_env))}
}

#### 9. Saving results ####
phen.all <- plyr::ldply(phen.list, data.frame) %>%
  dplyr::select(Ring:Year2, 
         B1_dep, 
         W1_arr, W1_dep, W2_arr, W2_dep, W3_arr, W3_dep, W4_arr, W4_dep, W5_arr, W5_dep,
         B2_arr,
         B1_dep_KDE, 
         W1_arr_KDE, W1_dep_KDE, W2_arr_KDE, W2_dep_KDE, W3_arr_KDE, W3_dep_KDE, W4_arr_KDE, W4_dep_KDE, W5_arr_KDE, W5_dep_KDE,
         B2_arr_KDE,
         WA1, WA2, WA3, WA4, WA5)


write_csv(phen.all, file = 'Results/Phenology/Carryover_effects/CALspp_phenology.csv')

track.all <- plyr::ldply(track.list, data.frame)
write_csv(track.all, file = 'Results/Phenology/Carryover_effects/CALspp_phenology_tracks_corrected_v2.csv')
