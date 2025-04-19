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
               'tidyverse', 'lubridate', 'readxl', 'readr', # data carpenter
               'sp', 'sf', 'spData', 'trip', # spatial data
               'migrateR', 'adehabitatHR', 'track2KBA', # migrate behaviour and kernels
               'segclust2d', # spatial segmentation and clustering
               'units', 
               'ggplot2', 'ggpubr', 'ggspatial',# plots
               'patchwork', 'cowplot', # composing plots
               'caTools'# runmean and runmin calcs
)

source(paste0(wd,'/RCode/Migratory patterns/06_FUNCTIONS_automated_phenology_v9.5b.R'))
source(paste0(wd,'/RCode/Migratory patterns/06B_FUNCTIONS_Manual_interactive_reclassification_v2.R'))
#### 1. Load files  and information ####
GLS.files <- read_csv('Raw Data/GLS.files3.csv', col_names = TRUE)
colonies <- read_excel('Raw Data/Colony_coordinates.xlsx')
equinox<-read_excel('Raw Data/Equinox_by_year.xlsx')
colony_days <- read_csv('Raw Data/colony_days_75_25_v2.csv')

phen.all <- read_csv('Results/Phenology/Carryover_effects/CALspp_phenology_results_v9.5.csv')
phen.all <- phen.all %>%
  distinct(ID, .keep_all = T)
seg_all <- read_csv('Results/Phenology/Carryover_effects/CALspp_segmentation_results1.csv')
seg_all2 <- read_csv('Results/Phenology/Carryover_effects/CALspp_segmentation_results2.csv')

seg_all <- seg_all[,-1]
seg_all2 <- seg_all2[,-1]

tracks_phen <- read_csv('Results/Phenology/Carryover_effects/CALspp_phenology_tracks_v9.5.csv')
# tracks_phen <- read_csv('Results/Phenology/Carryover_effects/CALspp_phenology_tracks_corrected.csv')
# colnames(tracks_phen)[1] <- 'ID'
# Which tracks I should check?
# 
# count(seg_all2$Confidence)
# seg_all %>%
#   dplyr::count(Metrics, sort = TRUE)
# 
# metrics_split <- strsplit(seg_all$Metrics, "-")
# 
# # Unlist the split data into a vector of individual metrics
# metrics_vector <- unlist(metrics_split)
# 
# # Convert to a data frame to count the occurrences of each metric
# metric_count <- data.frame(Metric = metrics_vector) %>%
#   dplyr::count(Metric, sort = TRUE)
# 
# # View the results
# print(metric_count)
# 
# # If a track has less than 0.33 segments it is mandatory to check
# # If it has between 0.33 and 0.5 is highly recommended to check
# table(seg_all2$ID, seg_all2$Confidence)
# 

# tracks_0 <- seg_all2 %>%
#   group_by(ID) %>%                         
#   filter(any(Confidence == 0)) %>%   
#   dplyr::summarise(Conf_mean = mean(Confidence, na.rm = T))%>%
#   left_join(.,auxi_mig, by = join_by(ID == geo_trip))
# tracks_33 <- seg_all2 %>%
#   group_by(ID) %>%                         
#   filter(any(Confidence <= 0.33)) %>%   
#   dplyr::summarise(count = n_distinct(ID))
# tracks_50 <- seg_all2 %>%
#   group_by(ID) %>%                         
#   filter(any(Confidence <= 0.50)) %>%   
#   dplyr::summarise(count = n_distinct(ID))

# Also you can check the mean confidence value of each track (< 0.80 is highly recommended to check)

#Backup
tracks_phen0 <- tracks_phen
phen.all0 <- phen.all

track.list <- list()
segs.list <- list()


# ids_check <- unique(c('19087002001', '2165003001', 'BR667001001', 'N805001001', 'N808001001',
#          'BD799003001', 'BD132003001', 'BL215002001', 'BM291002001', 'BM300001001',
#          '2208001001', 'S942002001', 'BH061001001', 'BD768001001', 'BD770001001', 
#          'BD772001001', 'BD776003001', '18B430002001', '19425002001', '2204003001',
#          'BM290001001', 'BX602002001', 'J550001001','Q665001001', 'T057001001', 
#          'BD756001001', 'BR932003001', 'N619003001', 'S928002001', 'BD892002001',
#          'BX571001001', 'BX064001001', 'BD780002001', 'BL218002001', 'BX643002001',
#          'P903001001', 'Q671001001', 'V396018001001', '1652004001', '1650001001',
#          '1716003001', 'BD800001001', 'BD808001001', 'N812001001', '1692002001',
#          'J554001001', 'BQ120003001', '33001001', 'BH061001001', 'BX329002001',
#          'BX333001001', 'BQ114003001', 'BR707002001', 'BR931001001', 'V396041001001'
#          ))

# ids_check <- unique(c('33001001','S928002001'))

ids <-unique(tracks_phen$ID) 
ids_length <- length(ids)
pre_env<-ls()
for (i in 1:ids_length) {
  # i = 600
  # i = which(ids == 'BD891002001')
  print(i)
  # x11(width = 10, height = 8)
  seg_df <- seg_all2 %>%
    filter(ID == ids[i])
  
  track <- tracks_phen %>%
    filter(ID == ids[i])
  # track <-track[-c(375,376),]
  if(any(is.na(track$Breed2))){
    track <- track %>%
      tidyr::fill(Breed, .direction = "down")%>%
      tidyr::fill(Breed2, .direction = "down")
  }
  utrack <- supervisedClass(seg_df = seg_df, 
                            track_df = track,
                            byseg = T, x11 = F)
  # Update Breed
  track2 <- track
  track2$Breed2 <- utrack$Breed2
  
  # any(is.na(utrack$Breed2))
  # 
  seg_df2 <- utrack %>%
    dplyr::select(Breed2, State:Check)%>%
    distinct()%>%
    drop_na()
  
  seg_df3 <- seg_df %>%
    add_column(Breed2 = seg_df2$Breed2)
  
  utrack2 <- supervisedClass(seg_df = seg_df3, 
                             track_df = track2,
                             byseg = F, x11 = F)
  # dev.off()
  # Update Breed
  track3<- track2
  track3$Breed2 <- utrack2$Breed2
  
  if(any(is.na(track3$Breed2))){
    track3 <- track3 %>%
      tidyr::fill(Breed, .direction = "down")%>%
      tidyr::fill(Breed2, .direction = "down")
  }

  # Store data reclassified
  track.list[[i]] <- list(track3)
  # names( track.list[[i]]) <- unique(track$ID)
  segs.list[[i]] <- list(seg_df3)
  
  cat('Track reclassified \n')
  post_env<-ls()
  post_env <- post_env[!(post_env %in% c('pre_env'))]
  
  if(i == ids_length){break}else{rm(list=setdiff(post_env, pre_env))}
  save.image("D:/Dropbox/Diego Vicente/TESIS/01_GLS_analysis/RCode/Migratory patterns/06B_Manual_interactive_reclassification_v2.RData")
}

track_all2 <- plyr::ldply(track.list, data.frame)
write_csv(track_all2, file = 'Results/Phenology/Carryover_effects/CALspp_phenology_tracks_corrected.csv')

seg_all3 <- plyr::ldply(segs.list, data.frame)
write.csv(seg_all2, file = 'Results/Phenology/Carryover_effects/CALspp_segmentation_results_corrected.csv')


