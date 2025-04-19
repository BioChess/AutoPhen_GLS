#### SET UP GEOLOCATOR ANALYSES ####

#### Updated: 14-05-2024

#### Created by Diego Vicente-Sastre
#### Departament Biologia Evolutiva, Ecologia i Ciéncies Ambientals, 
#### University of Barcelona (BEECA, UB)
#### Institut de Recerca de la Biodiversitat (IRBio)

################################################################################

clean_positions <- function (data = df_trip, dist_threshold = 40){
  if (!'dist' %in% names(data) || !is.numeric(data$dist)) {
    stop('Distance variable cannot be found or is not a numeric variable\n')
  }
  
  while(any(data$dist > dist_threshold, na.rm = T)) {
    # print = 1
    # Remove rows that exceed threshold
    data <- data[data$dist <= dist_threshold, ]
    
    # Calculate distances between consecutive points
    
    x1 <- data[-1, ]
    x2 <- data[-nrow(data), ]
    data$dist <- c(sqrt((x1$x - x2$x)^2 + (x1$y - x2$y)^2), NA)
    data$R2n <- (data$x - data$x[1])^2 + (data$y - data$y[1])^2
    data$dt <- c(unclass(x1$date) - unclass(x2$date), NA)
    data$dx <- c(x1$x - x2$x, NA)
    data$dy <- c(x1$y - x2$y, NA)
    abs.angle <- ifelse(data$dist < 1e-07, NA, atan2(data$dy,data$dx))
    
  }
  data <- drop_na(data)
  return(data)
}


################################################################################

remove_duplicates <- function(data, equinox_days = 20, 
                              Sep_eq = Sep_eq, Mar_eq = Mar_eq, dup_threshold = 4,
                              cdays = NULL){
  
  if(!'y' %in% names(data)){
    stop('Columns format must be ltraj names \n')
  }
  n <- nrow(data)
  
  data <- data %>%
    mutate(Equinox = ifelse(date >= (as.Date(Sep_eq) - days(equinox_days)) & 
                              date <= (as.Date(Sep_eq) + days(equinox_days)), 'Sep',
                            ifelse(date >= (as.Date(Mar_eq) - days(equinox_days)) & 
                                     date <= (as.Date(Mar_eq) + days(equinox_days)), 'Mar', 'No'))
    )
  if (!is.null(cdays)){
    data <- data %>%
      filter(!(as.Date(date) %in% cdays$date))
    cat('Days at the colony has been removed \n')
    
  }
  # Remove equinox positions when there is no latitudinal movement
  if(length(unique(data$y[data$Equinox == 'Sep'])) == 1){
    data <- data %>%
      filter(Equinox != 'Sep' | date == date[.$Equinox == 'Sep'][1] |
               date == date[.$Equinox == 'Sep'][length(data$date[data$Equinox == 'Sep'])] )
    
    cat('September equinox positions without latitudinal movement has been removed \n')
    
  }
  
  if(length(unique(data$y[data$Equinox == 'Mar'])) == 1){
    data <- data %>%
      filter(Equinox != 'Mar' | date == date[.$Equinox == 'Mar'][1] |
               date == date[.$Equinox == 'Mar'][length(data$date[data$Equinox == 'Mar'])] )
    cat('March equinox positions without latitudinal movement has been removed \n')
  }
  
  if(sum(duplicated(data$y)) > dup_threshold){
    data <- data %>%
      distinct(y, .keep_all = T)
    cat('Other duplicated latitudes has been removed \n')
  }
  
  if(any(duplicated(tail(data$y, 5)))){
    data <- data[1:(nrow(data) - 5),]
  } else if (any(duplicated(head(data$y, 5)))){
    data <- data[5:nrow(data),]
  }
  
  if(nrow(data) == n){
    cat('No duplicates found \n')
  } else {
    # Correct R2n values, according to the real distance to the colony
    data <- data %>%
      mutate(R2n = (x - col.lon)^2 + (y - col.lat)^2)
  }
  
  return(data)
}



################################################################################

trajmetrics <- function(trj){
  
  if (!requireNamespace("trajr", quietly = TRUE)) {
    stop("The 'trajr' package is required but is not installed. 
         Please install it before proceeding.")
  }
  
  if (!inherits(trj, "Trajectory")) {
    stop("Input 'trj' must be of class 'Trajectory' from the 'trajr' package.")
  }
  
  out <- tryCatch({
    # Rediscretize the trajectory based on the mean step length
    trj_r <- trajr::TrajRediscretize(trj, R = mean(trajr::TrajStepLengths(trj)))
    
    # Turning angles (mean)
    ang_vector <- Mod(trajr::TrajMeanVectorOfTurningAngles(trj_r))  # mean vector angle
    ang_mean <- Arg(trajr::TrajMeanVectorOfTurningAngles(trj_r))    # mean angle
    
    # Expected square displacement of a correlated random walk
    sq_dp <- trajr::TrajExpectedSquareDisplacement(trj)
    Emax <- trajr::TrajEmax(trj)
    
    # Straightness and sinuosity
    stght <- trajr::TrajStraightness(trj_r)
    snsty <- trajr::TrajSinuosity2(trj_r)
    
    # Directional change (cannot be calculated in trj_r)
    dc <- mean(TrajDirectionalChange(trj))
    sddc <- sd(TrajDirectionalChange(trj))
    
    # Create a data frame for output
    data.frame(
      sinuosity = snsty,
      straightness = stght,
      ang_vector = ang_vector,
      ang_mean = ang_mean,
      Emax = Emax,
      sq_disp = sq_dp,
      dc_mean = dc,
      dc_sd = sddc
    )
  }, error = function(e) {
    # Return an informative error message if something goes wrong
    cat("An error occurred while processing the trajectory data: ", e$message, "\n")
    return(NULL)
  })
  return(out)
}
################################################################################

################################################################################
# Determining optimun R2n limit

R2n_limit <- function(x = track, plot = F){
  if (!requireNamespace("segclust2d", quietly = TRUE)) {
    stop("The 'segclust2d' package is required but is not installed. 
         Please install it before proceeding.")
  }
  
  if (!is.data.frame(x)) {
    stop("'x' must be a data.frame.")
  }
  
  if (!is.logical(plot) || length(plot) != 1) {
    stop("'plot' must be a logical value (TRUE or FALSE).")
  }
  
  out <- tryCatch({
    suppressWarnings({
      x_segclust <- segclust(x,
                             Kmax = 3,
                             lmin = 5,
                             ncluster = 3,
                             seg.var = c('R2n'),
                             subsample = F)
      
      x_seg<-segment(x_segclust)
    
      mean(c(x$R2n[x_seg$begin[2]], x$R2n[x_seg$begin[3]]), na.rm = T)
    })
  }, error = function(e) {
    message("An error occurred while calculating changepoints: ", e$message)
    return(NA)
  })
  
  if (plot) {
    # plot(cpt_result)
    plot(x_segclust)
  }
  
  return(out)
}
################################################################################
combm <- function(v = add.metrics, minv, maxv) {
  
  combinations <- list()
  for (i in minv:maxv) {
    combs <- combn(v, i, simplify = FALSE)
    combinations <- append(combinations, combs)
  }
  
  return(combinations)
}

################################################################################

################################################################################

run_clust <- function(var_comb, df_seg, 
                      method = "complete",
                      num_clusters = 2, R2n = R2n) {
  
  if (!is.data.frame(df_seg)) stop("The 'df_seg' input must be a data frame.")
  if (!is.numeric(num_clusters) || num_clusters <= 0) stop("'num_clusters' must be a positive integer.")
  
  formula <- as.formula(paste("~", paste(var_comb, collapse = " + ")))

  recipe1 <- recipe(formula, data = df_seg) %>%
    step_normalize(all_predictors())
  
  recipes <- list(
    recipe_1 = recipe1,
    recipe_2 = recipe(~ R2n_mean, data = df_seg) %>%
      step_normalize(all_predictors())
    # recipe_2 = recipe(~ R2n_mean + sd.dx + sd.dy, data = segs2a) %>%
    #   step_normalize(all_predictors())
  )
 
  prepped_recipes <- lapply(recipes, prep, training = df_seg)

  baked_data <- lapply(prepped_recipes, bake, new_data = df_seg)

  dist_matrices <- lapply(baked_data, dist)

  hclust_res <- lapply(dist_matrices, function(dist_matrix) {
    hclust_fit <- hclust(dist_matrix, method = method)
    cutree(hclust_fit, k = num_clusters)
  })
  
  df_clust <- df_seg %>%
    mutate(Clust_SM = as.character(hclust_res$recipe_1),
           Clust_BnB = as.character(hclust_res$recipe_2))

  df_clust <- df_clust %>%
    mutate(
      S_M = ifelse(Clust_SM %in% Clust_SM[which(ang_vector == min(.$ang_vector, na.rm = T))], 'Stationary', 'Mig'),
      B_nB = ifelse(Clust_BnB %in% Clust_BnB[which(R2n_mean == min(.$R2n_mean, na.rm = T))], 'B', 'Non-B'),
      Breed = case_when(
        S_M == 'Stationary' & B_nB == 'B' ~ 'B',
        S_M == 'Stationary' & B_nB == 'Non-B' & (end - begin) <= 5 & straightness > 0.5 ~ 'mig',
        S_M == 'Stationary' & B_nB == 'Non-B' ~ 'W',
        S_M == 'Mig' & B_nB == 'Non-B' ~ 'mig',
        S_M == 'Mig' & B_nB == 'B' ~ case_when(
          (end - begin) <= 5 ~ 'B',
          row_number() == 1 ~ 'B',
          row_number() == nrow(df_clust) & R2n_mean <= R2n ~ 'B',
          row_number() != nrow(df_clust) & R2n_mean <= R2n & ang_vector < 0.75 ~ 'B',
          TRUE ~ 'mig'
        ),
        
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(
      B_W = ifelse(Breed %in% c('B', 'W'),
                   (straightness < 0.5) + (ang_vector < 0.5) + (sinuosity > 1),
                   0),
      mig = ifelse(Breed == 'mig',
                   (straightness > 0.5) + (ang_vector > 0.5) + (sinuosity < 1),
                   0)
    ) %>%
    mutate(
      Confidence = case_when(
        row_number() == 1 ~ 3,
        row_number() == nrow(.) ~ 3,
        # Zero confidence for impossible classifications
        Breed == 'B' & R2n_mean > R2n ~ 0,
        Breed == 'W' & R2n_mean < R2n ~ 0,
        Breed == 'mig' & straightness < 0.2 & ang_vector < 0.9 ~ 0,
        B_W == 3 | mig == 3 ~ 3,
        B_W == 2 | mig == 2 ~ 2,
        B_W == 1 | mig == 1 ~ 1,
        TRUE ~ 0
      )
    ) %>%
    dplyr::select(-B_W, -mig)
  
  # Calculate the sum of confidence for this recipe
  mean_confidence <- mean(df_clust$Confidence, na.rm = TRUE)/3
  
  out <-list(
    Confidence = mean_confidence,
    df_clust = df_clust,
    metrics = var_comb
  )
  return(out)
}

################################################################################

segment_class <- function (mode_segclust =  mode_segclust,
                           option.by = 'seg',
                           add.metrics = NULL,
                           col.lon = NULL, col.lat = NULL){
  
  required_packages <- c("trajr", "tidyverse", "segclust2d", "ClusterR", "tidymodels", "tidyclust", "furrr")
  invisible(lapply(required_packages, function(p) {
    if (!require(p, character.only = TRUE)) {
      stop(paste("Package", p, "is required but not loaded. Please install it and try again."))
    }
  }))
  
  # Track positions classified by clusters
  track <- segclust2d::augment(mode_segclust)
 
  track_trj <- TrajFromCoords(track[,c('x','y','date','R2n', 'dist', 'state')],
                              spatialUnits = 'km', timeUnits = 's')
  segs <- segment(mode_segclust)
  segs <- segs %>%
    dplyr::select(state, begin:state_ordered)
  
  trj.list <- list()
  if(option.by == 'seg'){
    for (i in 1:nrow(segs)){
      # i= 1
      # print(i)
      begin <- segs[i,]$begin
      end <- segs[i,]$end
      trj_stt <- track_trj[begin:end,]
      stt <- unique(trj_stt$state)
      
      metrics <- trajmetrics(trj_stt)
      metrics$R2n_mean <- mean(trj_stt$R2n)
      metrics$R2n_sd <- sd(trj_stt$R2n)
      metrics$state <- unique(trj_stt$state)
      metrics$begin <- begin
      metrics$end <- end
      
      metrics[,c(1:10)] <- round(metrics[,c(1:10)], digits = 3)
      
      trj.list[i] <- list(metrics)
    }
    trj_df <- plyr::ldply(trj.list, data.frame)
    segs2 <- segs %>%
      left_join(., trj_df, by = c('begin', 'end', 'state'))
    
  } else if(option.by == 'cluster'){
    for (i in unique(segs$state)){
      # i= 1
      trj_stt <- track_trj  %>%
        filter(state == i)
      
      metrics <- Trajmetrics(trj_stt)
      metrics$R2n_mean <- mean(trj_stt$R2n)
      metrics$R2n_sd <- sd(trj_stt$R2n)
      metrics$state <- unique(trj_stt$state)
      metrics[,c(1:10)] <- round(metrics[,c(1:10)], digits = 3)
      trj.list[i] <- list(metrics)
    }
    trj_df <- plyr::ldply(trj.list, data.frame)
    segs2 <- segs %>%
      left_join(., trj_df, by = 'state')
    
  }
  
  R2n <- R2n_limit(track, plot = F)
  if(max(R2n) <= median(segs2$R2n_mean)){
    R2n <- max(R2n)
  } else {
    R2n <- min(R2n)
  }

  if(!is.null(add.metrics)){
    comb_metrics <- combm(add.metrics, minv = 1, maxv =3)
    
    # Add the base variables to each combination
    comb_metrics2 <- lapply(comb_metrics, function(vars) {
      expanded_vars <- unlist(lapply(vars, function(v) {
        if (v %in% names(pair.m)) pair.m[[v]] else v
      }))
      unique(c(expanded_vars, 
               setdiff(c("straightness", "sinuosity", "ang_vector"), expanded_vars)))
    })
  } else {
    comb_metrics2 <- list()
    comb_metrics2[[1]] <- list("straightness", "sinuosity", "ang_vector")
  }

  
  
  if (!is.null(segs2)) {
    plan(multisession)  # or plan(multicore) if supported
    clustering_results <- future_map(comb_metrics2, function(var_combination) {
      run_clust(var_combination, df_seg = segs2, method = "complete", num_clusters = 2, R2n = R2n)
    },
    .options = furrr_options(seed = TRUE)
    )
    
    # Find the best result by the sum of confidence
    best <- tibble(
      Confidence = map_dbl(clustering_results, "Confidence"),
      Metrics = map_int(clustering_results, ~ length(.x$metrics))
    ) %>%
      mutate(pos = row_number()) %>%
      arrange(desc(Confidence), Metrics) %>%
      slice(1) %>%
      pull(pos)

    # Retrieve the best clustering result
    segs3 <- clustering_results[[best]]
    return(segs3)
  } else {
    stop("Error: clustering results were not computed correctly.")
  }
}
################################################################################

equinox_interference <- function(df_seg = segments1, track = df_trip2,
                                 Sep_eq, Mar_eq){
  # require(dplyr)
  eq_dates <- c(seq(Sep_eq - days(20), Sep_eq + days(20), by = 'day'),
                seq(Mar_eq - days(20), Mar_eq + days(20), by = 'day')
  )
  
  df_seg2 <- df_seg %>%
    dplyr::mutate(f_day = as.Date(track$date[begin]),
                  l_day = as.Date(track$date[end])) %>%
    rowwise() %>%
    dplyr::mutate(Eq = any(eq_dates >= f_day & eq_dates <= l_day)) %>%
    dplyr::mutate(Confidence = ifelse(Eq == T & Breed == 'mig' &
                                       Confidence < 3, 
                                      # 2*sd.x < sd.y,
                                      Confidence/2, Confidence))%>%
    ungroup()
  
  
  return(df_seg2)
  
}

################################################################################


# stght_th <- 0.5
# angv_th <- 0.5
# snsty_th <- 1


segment_correction <- function(df_seg = segments2){
  
  median_sq_disp <- median(df_seg$sq_disp, na.rm = T)
  mean_R2n_mean <- mean(df_seg$R2n_mean, na.rm = T)

  if('Eq' %in% colnames(df_seg)){
    # Equinox interference
    df_seg <- df_seg %>%
      dplyr::mutate(Breed = case_when(Breed == 'mig' & Eq == T & 
                                     sq_disp < median_sq_disp & 2*sd.x < sd.y ~ 
                                       ifelse(R2n_mean < mean_R2n_mean ,'B', 'W'), 
                                     Breed == 'B' & Eq == T & (sd.x > 1 | sd.y > 1) ~
                                       ifelse(sd.x > 4*sd.y, 'mig', Breed),
                                     T ~ Breed))
  }
  
  # Breed_0 <- df_seg$Breed
  
  # Wintering area classified as breeding area
  segs4 <- df_seg %>%
    dplyr::mutate(Breed = ifelse(Breed == 'B' & R2n_mean > mean_R2n_mean &
                                   straightness < 0.2, 'W', Breed))
  
  # Migrating area classified as breeding area
  segs5 <- segs4 %>%
    mutate(
      Breed = ifelse(Breed == 'B' & Confidence < 1.5 & 
                       (lead(Breed) == 'mig' | lag(Breed) == 'mig') &
                       straightness > 0.2, 
                     case_when(
                       sq_disp > median_sq_disp ~ 'mig',
                       is.nan(sq_disp) ~ Breed,
                       TRUE ~ Breed
                     ), 
                     Breed)
    )
  # Breeding area classified as wintering area
  segs6 <- segs5 %>%
    dplyr::mutate(Breed = ifelse(Breed == 'W' & R2n_mean < mean_R2n_mean , 'B', Breed))
  
  # Migrating area classified as wintering area
  segs7 <- segs6 %>%
    dplyr::mutate(Breed = ifelse(Breed == 'W' & sq_disp > median_sq_disp & 
                                   sinuosity < 1 & straightness > 0.2, 'mig', Breed))
  
  # Impsible situations
  out <- segs7 %>%
    mutate(
      Breed = case_when(
        Breed != 'B' & lead(Breed) == 'B' & lag(Breed) == 'B' ~ 'B',
        Breed == 'B' & lead(Breed) == 'W' & lag(Breed) == 'W' ~ 'W',
        TRUE ~ Breed
      )
    ) %>%
    mutate(
      Breed = ifelse(row_number() != 1 & row_number() != nrow(.) &
                             Breed == 'B' & lead(Breed) == 'mig' &
                             lag(Breed) == 'mig','W', Breed)) %>%
    mutate(
      Breed_0 = df_seg$Breed
    ) %>%
    mutate(
      Confidence = case_when(
        Breed_0 != Breed ~ 3,
      #   Breed_0 == 'B' & Breed2 == 'W' ~ 1,  
      #   Breed_0 != 'B' & Breed2 == 'B' ~ 1,                
        TRUE ~ Confidence                         
      )
    )
  
  return(out)
}




################################################################################
my_theme <- theme_bw() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(size = 20),
        legend.title = element_text(size = 18,face = "bold"),
        legend.text = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) 

plot_segmentation <- function(segclust = NULL, map = NULL, 
                              segclass = NULL, col.lon = NULL, col.lat = NULL,
                              breed = T, data = df_trip4,
                              world = world){
  
 
  
  if(!is.null(segclust)){
    segclust <- segclust +
      labs(title = 'Segmentation/Clustering')+
      xlab('Indice') + ylab('Value')+
      guides(color = FALSE, size = FALSE)+ #Remove legend
      guides(fill=guide_legend(title='State'))+ # New legend
      my_theme
  } else {
    segclust <- NULL
  }
  
  if(!is.null(map)){
    map <- map +
      labs(title = 'Segmentation/Clustering map')+
      xlab('Longitude') + ylab('Latitude')+
      guides(color = FALSE, size = FALSE)+ #Remove legend
      my_theme
    
  } else {
    map <- NULL
  }
  
  if(!is.null(segclass)){
    # Plot to understand the segmentation index
    segclass_plot <- ggplot(segclass, aes(x = mu.x, y = mu.y)) +
      geom_point(x = col.lon, y = col.lat, 
                 shape = 24, fill = 'orchid3', size = 5)+
      geom_point() +
      geom_errorbar(aes(ymin = (mu.y-sd.y), ymax = (mu.y+sd.y)),
                    width = 0.2, color = "blue") +
      geom_errorbarh(aes(xmin = (mu.x-sd.x), xmax = (mu.x+sd.x)),
                     height = 0.2, color = "red")+
      labs(title = 'Trip periods classification')+
      xlab('Longitude') + ylab('Latitude')+
      my_theme
  } else {
    segclass_plot <- NULL
  }
  
  if(!is.null(breed)){
    require(rnaturalearth)
    world <- ne_countries(scale = "large", returnclass = "sf")
    data <- data %>%
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
    breaks <- unique(data$Color)
    labels <- unique(data$lab)
    
    # Breeding vs non-breeding plot
    breed_plot <- ggplot()+
      # World
      {if(!is.null(world))
        geom_sf(data = world, fill=grey(0.9), color=grey(0.6), lwd = 0.2)
      } +
      
      coord_sf(xlim = c(min(data$x), max(data$x)), 
               ylim = c(min(data$y), max(data$y)), expand = TRUE) +
      # Trip
      geom_path(data= data, aes(x, y))+
      {if(length(unique(data$Shape)) > 1)
        geom_point(data= data, aes(x, y, fill = Color,
                                   shape = Shape), size = 3)
      } +
      {if(length(unique(data$Shape)) == 1)
        geom_point(data= data, aes(x, y, fill = Color),
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
    
  } else {
    breed_plot <- NULL
  }

  plot_list <- list(segclust = segclust_plot, map = map_seg, segclass = segclass_plot, breed = breed_plot)
  return(plot_list)
}

