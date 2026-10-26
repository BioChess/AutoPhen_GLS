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
    # Coordinates are projected and centred on the colony (0, 0), so the
    # net squared colony distance (NSCD) is x^2 + y^2 (km^2)
    data$R2n <- data$x^2 + data$y^2
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
                              edge_positions = 5, cdays = NULL){
  # Repeated values are checked on latitude ('lat', in degrees), because GLS
  # modelling can produce identical latitudes. 'x' and 'y' are projected
  # coordinates (km) centred on the colony.
  
  if(!all(c('y', 'lat') %in% names(data))){
    stop('Data must have ltraj columns (x, y) and latitude (lat) \n')
  }
  n <- nrow(data)
  
  # Equinox periods (only for GLS data; equinox_days = NULL for other data)
  if (is.null(equinox_days)) {
    data$Equinox <- 'No'
  } else {
    data <- data %>%
      mutate(Equinox = ifelse(date >= (as.Date(Sep_eq) - days(equinox_days)) & 
                                date <= (as.Date(Sep_eq) + days(equinox_days)), 'Sep',
                              ifelse(date >= (as.Date(Mar_eq) - days(equinox_days)) & 
                                       date <= (as.Date(Mar_eq) + days(equinox_days)), 'Mar', 'No'))
      )
  }
  if (!is.null(cdays)){
    data <- data %>%
      filter(!(as.Date(date) %in% cdays$date))
    cat('Days at the colony has been removed \n')
    
  }
  # Remove equinox positions when there is no latitudinal movement
  if(length(unique(data$lat[data$Equinox == 'Sep'])) == 1){
    data <- data %>%
      filter(Equinox != 'Sep' | date == date[.$Equinox == 'Sep'][1] |
               date == date[.$Equinox == 'Sep'][length(data$date[data$Equinox == 'Sep'])] )
    
    cat('September equinox positions without latitudinal movement has been removed \n')
    
  }
  
  if(length(unique(data$lat[data$Equinox == 'Mar'])) == 1){
    data <- data %>%
      filter(Equinox != 'Mar' | date == date[.$Equinox == 'Mar'][1] |
               date == date[.$Equinox == 'Mar'][length(data$date[data$Equinox == 'Mar'])] )
    cat('March equinox positions without latitudinal movement has been removed \n')
  }
  
  # Positions fixed at the colony (first and last positions of the track,
  # at x = 0, y = 0) are repeated on purpose: they are not checked
  at_colony <- abs(data$x) < 1e-6 & abs(data$y) < 1e-6
  lat_check <- ifelse(at_colony, NA, data$lat)
  dup_lat   <- duplicated(lat_check) & !is.na(lat_check)

  if(sum(dup_lat) > dup_threshold){
    data <- data[!dup_lat, ]
    at_colony <- at_colony[!dup_lat]
    cat('Other duplicated latitudes has been removed \n')
  }
  
  lat_free <- data$lat[!at_colony]
  if(any(duplicated(tail(lat_free, edge_positions)))){
    drop <- tail(which(!at_colony), edge_positions)
    data <- data[-drop, ]
  } else if (any(duplicated(head(lat_free, edge_positions)))){
    drop <- head(which(!at_colony), edge_positions - 1)
    data <- data[-drop, ]
  }
  
  if(nrow(data) == n){
    cat('No duplicates found \n')
  }
  
  # NSCD: squared distance to the colony (km^2); the colony is at (0, 0)
  data <- data %>%
    mutate(R2n = x^2 + y^2)
  
  return(data)
}



################################################################################

################################################################################
# agreement_score(): agreement between the class of each segment and its
# movement metrics (0-1). Breeding and wintering segments are expected to be
# tortuous and slow (straightness < 0.5, ang_vector < 0.5, sinuosity > 1,
# mean step length not much larger than the median of the track) and
# migration segments the opposite. The score is the proportion of these four
# criteria that agree with the class. It is 0 for impossible classes:
# breeding far from the colony (NSCD above the limit), wintering near the
# colony, or migration with very low straightness and ang_vector.
# The score is a heuristic: it describes how typical the movement of a
# segment is for its class, not the probability that the class is correct.

agreement_score <- function(df_seg, R2n, adv) {
  speed_ref <- median(df_seg$mu.dist, na.rm = TRUE)
  stationary <- (df_seg$straightness < adv$score_straightness) +
    (df_seg$ang_vector < adv$score_ang_vector) +
    (df_seg$sinuosity > adv$score_sinuosity) +
    (df_seg$mu.dist <= adv$score_speed * speed_ref)
  moving <- (df_seg$straightness > adv$score_straightness) +
    (df_seg$ang_vector > adv$score_ang_vector) +
    (df_seg$sinuosity < adv$score_sinuosity) +
    (df_seg$mu.dist > adv$score_speed * speed_ref)
  n_agree <- ifelse(df_seg$Breed %in% c('B', 'W'), stationary, moving)

  impossible <- with(df_seg,
    (Breed == 'B' & R2n_mean > R2n) |
    (Breed == 'W' & R2n_mean < R2n) |
    (Breed == 'mig' & straightness < adv$impossible_straightness &
       ang_vector < adv$impossible_ang_vector))

  score <- n_agree / 4
  score[dplyr::coalesce(impossible, FALSE)] <- 0
  round(score, 3)
}

################################################################################
# fast_w_to_mig(): wintering segments at the start or end of the wintering
# stage that are much faster than the main (longest) wintering segment are
# reclassified as migration. Only consecutive segments from the edge are
# changed, so the core of the wintering stage is never changed.

fast_w_to_mig <- function(df_seg, adv) {
  b <- df_seg$Breed
  w <- which(b == 'W')
  if (length(w) < 2) return(b)
  core <- w[which.max((df_seg$end - df_seg$begin)[w])]
  fast <- dplyr::coalesce(df_seg$mu.dist > adv$fast_w_ratio * df_seg$mu.dist[core] &
                            df_seg$mu.dist > adv$fast_w_min, FALSE)
  lo <- core; while (lo > 1 && b[lo - 1] == 'W') lo <- lo - 1
  hi <- core; while (hi < length(b) && b[hi + 1] == 'W') hi <- hi + 1
  i <- lo; while (i < core && fast[i]) { b[i] <- 'mig'; i <- i + 1 }
  i <- hi; while (i > core && fast[i]) { b[i] <- 'mig'; i <- i - 1 }
  b
}

################################################################################
# flag_segments(): segments to review ('Check') and the reason(s):
#   low agreement   agreement lower than adv$check_score
#   corrected       class changed by the correction rules
#   equinox         change of stage within an equinox window (GLS)
#   direct B-W      change between breeding and wintering without migration
#   far W           wintering segment more than adv$w_far_km away from the
#                   main (longest) wintering segment (stopover or loop?)

flag_segments <- function(df_seg, adv) {
  n <- nrow(df_seg)
  change <- (df_seg$Breed != dplyr::lag(df_seg$Breed)) |
    (df_seg$Breed != dplyr::lead(df_seg$Breed))
  change[is.na(change)] <- FALSE
  direct <- (df_seg$Breed == 'W' & dplyr::lag(df_seg$Breed) == 'B') |
    (df_seg$Breed == 'B' & dplyr::lag(df_seg$Breed) == 'W')
  direct <- dplyr::coalesce(direct, FALSE)

  far_W <- rep(FALSE, n)
  w <- which(df_seg$Breed == 'W')
  if (length(w) > 1) {
    core <- w[which.max((df_seg$end - df_seg$begin)[w])]
    d <- sqrt((df_seg$mu.x - df_seg$mu.x[core])^2 +
                (df_seg$mu.y - df_seg$mu.y[core])^2)
    far_W <- df_seg$Breed == 'W' & d > adv$w_far_km
  }

  reasons <- cbind(
    'low agreement' = dplyr::coalesce(df_seg$Agreement < adv$check_score, FALSE),
    'corrected'     = df_seg$Breed != df_seg$Breed_0,
    'equinox'       = df_seg$Eq & change,
    'direct B-W'    = direct,
    'far W'         = far_W)
  df_seg$Check_reason <- apply(reasons, 1, function(r) paste(colnames(reasons)[r], collapse = '; '))
  df_seg$Check <- ifelse(df_seg$Check_reason == '', 'Ok', 'Check')
  df_seg
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
    # Rediscretize the trajectory based on the mean step length. In very
    # short or compact segments (e.g. positions around the colony) this can
    # fail; then a shorter step length is used (half, then a quarter).
    R0 <- mean(trajr::TrajStepLengths(trj))
    trj_r <- NULL
    for (f in c(1, 0.5, 0.25)) {
      trj_r <- tryCatch(trajr::TrajRediscretize(trj, R = R0 * f),
                        error = function(e) NULL)
      if (!is.null(trj_r) && nrow(trj_r) > 2) break
    }
    if (is.null(trj_r)) stop('the segment could not be rediscretised')
    
    # Turning angles (mean)
    ang_vector <- Mod(trajr::TrajMeanVectorOfTurningAngles(trj_r))  # mean vector angle
    ang_mean <- Arg(trajr::TrajMeanVectorOfTurningAngles(trj_r))    # mean angle
    
    # Expected square displacement of a correlated random walk
    sq_dp <- trajr::TrajExpectedSquareDisplacement(trj)
    Emax <- trajr::TrajEmax(trj)
    
    # Straightness and sinuosity
    stght <- trajr::TrajStraightness(trj_r)
    # Sinuosity (Benhamou 2004) has units of 1/sqrt(length), so its value
    # depends on the units and on the step length. We multiply it by the
    # square root of the mean step length of the rediscretised trajectory to
    # obtain a dimensionless index that depends only on the turning angles
    # (and on the variation of step lengths). With constant step length:
    #   sinuosity = 2 / sqrt((1 + c) / (1 - c)), c = mean cosine of turning angles
    # sinuosity = 1 when c = 0.6; higher values = more tortuous paths.
    snsty <- trajr::TrajSinuosity2(trj_r) *
      sqrt(mean(trajr::TrajStepLengths(trj_r)))
    
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
    # If the metrics cannot be calculated, they are NA (the classification
    # rules treat NA conditions as FALSE)
    cat("Trajectory metrics not calculated for one segment: ", e$message, "\n")
    return(data.frame(sinuosity = NA_real_, straightness = NA_real_,
                      ang_vector = NA_real_, ang_mean = NA_real_,
                      Emax = NA_real_, sq_disp = NA_real_,
                      dc_mean = NA_real_, dc_sd = NA_real_))
  })
  return(out)
}
################################################################################

################################################################################
# Determining optimun R2n limit

R2n_limit <- function(x = track, plot = F, Kmax = 3, lmin = 5, ncluster = 3){
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
  
  # Consecutive (nearly) identical NSCD values (e.g. positions fixed at the
  # same distance from the colony) stop segclust2d. Only the first position
  # of each run is kept (change in distance to the colony < 0.5 km). This does
  # not change the distance values used for the threshold.
  d_col <- sqrt(x$R2n)
  x <- x[c(TRUE, abs(diff(d_col)) >= 0.5), ]
  
  out <- tryCatch({
    suppressWarnings({
      x_segclust <- segclust(x,
                             Kmax = Kmax,
                             lmin = lmin,
                             ncluster = ncluster,
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
                      num_clusters = 2, R2n = R2n,
                      adv, lmin = 5) {
  
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
        S_M == 'Stationary' & B_nB == 'Non-B' & (end - begin) <= lmin &
          straightness > adv$short_straightness ~ 'mig',
        S_M == 'Stationary' & B_nB == 'Non-B' ~ 'W',
        S_M == 'Mig' & B_nB == 'Non-B' ~ 'mig',
        S_M == 'Mig' & B_nB == 'B' ~ case_when(
          (end - begin) <= lmin ~ 'B',
          row_number() == 1 ~ 'B',
          row_number() == nrow(df_clust) & R2n_mean <= R2n ~ 'B',
          row_number() != nrow(df_clust) & R2n_mean <= R2n &
            ang_vector < adv$breeding_ang_vector ~ 'B',
          TRUE ~ 'mig'
        ),
        
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(R2n_limit = R2n)
  df_clust$Agreement <- agreement_score(df_clust, R2n = R2n, adv = adv)

  # Mean agreement of this combination of metrics
  mean_agreement <- mean(df_clust$Agreement, na.rm = TRUE)
  
  out <-list(
    Agreement = mean_agreement,
    df_clust = df_clust,
    metrics = var_comb
  )
  return(out)
}

################################################################################

segment_class <- function (mode_segclust =  mode_segclust,
                           option.by = 'seg',
                           add.metrics = NULL,
                           col.lon = NULL, col.lat = NULL,
                           adv, lmin = 5){
  
  required_packages <- c("trajr", "tidyverse", "segclust2d", "recipes", "furrr")
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
    
  } else {
    stop("option.by must be 'seg'.")
  }
  
  R2n <- R2n_limit(track, plot = F, Kmax = adv$nscd_Kmax, lmin = lmin,
                   ncluster = adv$nscd_ncluster)
  if(max(R2n) <= median(segs2$R2n_mean)){
    R2n <- max(R2n)
  } else {
    R2n <- min(R2n)
  }

  if(!is.null(add.metrics)){
    comb_metrics <- combm(add.metrics, minv = 1, maxv = adv$max_extra_metrics)
    
    # Add the base variables to each combination
    comb_metrics2 <- lapply(comb_metrics, function(vars) {
      expanded_vars <- unlist(lapply(vars, function(v) {
        if (v %in% names(adv$paired_metrics)) adv$paired_metrics[[v]] else v
      }))
      unique(c(expanded_vars, 
               setdiff(adv$core_metrics, expanded_vars)))
    })
  } else {
    comb_metrics2 <- list()
    comb_metrics2[[1]] <- as.list(adv$core_metrics)
  }

  
  
  if (!is.null(segs2)) {
    plan(multisession)  # or plan(multicore) if supported
    clustering_results <- future_map(comb_metrics2, function(var_combination) {
      run_clust(var_combination, df_seg = segs2, method = adv$clust_method,
                num_clusters = adv$clust_k, R2n = R2n, adv = adv, lmin = lmin)
    },
    .options = furrr_options(seed = TRUE)
    )
    
    # Best combination of metrics: highest mean agreement (and fewest
    # metrics if tied)
    best <- tibble(
      Agreement = map_dbl(clustering_results, "Agreement"),
      Metrics = map_int(clustering_results, ~ length(.x$metrics))
    ) %>%
      mutate(pos = row_number()) %>%
      arrange(desc(Agreement), Metrics) %>%
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
                                 Sep_eq, Mar_eq, equinox_days = 20){
  # Flags the segments that overlap an equinox window (Eq = TRUE). The
  # agreement score is not changed: segments at a change of stage within an
  # equinox window are flagged for review (see flag_segments()).
  if (is.null(equinox_days)) {
    # No equinox effect (data other than GLS)
    return(df_seg %>%
             dplyr::mutate(f_day = as.Date(track$date[begin]),
                           l_day = as.Date(track$date[end]),
                           Eq = FALSE))
  }
  eq_dates <- c(seq(Sep_eq - days(equinox_days), Sep_eq + days(equinox_days), by = 'day'),
                seq(Mar_eq - days(equinox_days), Mar_eq + days(equinox_days), by = 'day')
  )
  
  df_seg2 <- df_seg %>%
    dplyr::mutate(f_day = as.Date(track$date[begin]),
                  l_day = as.Date(track$date[end])) %>%
    rowwise() %>%
    dplyr::mutate(Eq = any(eq_dates >= f_day & eq_dates <= l_day)) %>%
    ungroup()
  
  
  return(df_seg2)
  
}

################################################################################

segment_correction <- function(df_seg = segments2, adv){
  # Missing metrics (NA) must not change the class: conditions with NA are
  # treated as FALSE (no correction)
  ok <- function(cond) dplyr::coalesce(cond, FALSE)
  
  median_sq_disp <- median(df_seg$sq_disp, na.rm = T)
  mean_R2n_mean <- mean(df_seg$R2n_mean, na.rm = T)

  if('Eq' %in% colnames(df_seg)){
    # Equinox interference
    df_seg <- df_seg %>%
      dplyr::mutate(Breed = case_when(Breed == 'mig' & Eq == T & 
                                     sq_disp < median_sq_disp &
                                       adv$eq_sd_ratio_mig * sd.x < sd.y ~ 
                                       ifelse(R2n_mean < mean_R2n_mean ,'B', 'W'), 
                                     Breed == 'B' & Eq == T &
                                       (sd.x > adv$eq_sd_limit | sd.y > adv$eq_sd_limit) ~
                                       ifelse(sd.x > adv$eq_sd_ratio_breed * sd.y, 'mig', Breed),
                                     T ~ Breed))
  }
  
  # Breed_0 <- df_seg$Breed
  
  # Wintering area classified as breeding area
  segs4 <- df_seg %>%
    dplyr::mutate(Breed = ifelse(ok(Breed == 'B' & R2n_mean > mean_R2n_mean &
                                   straightness < adv$corr_straightness), 'W', Breed))
  
  # Migrating area classified as breeding area
  segs5 <- segs4 %>%
    mutate(
      Breed = ifelse(ok(Breed == 'B' & Agreement < adv$corr_score & 
                       (lead(Breed) == 'mig' | lag(Breed) == 'mig') &
                       straightness > adv$corr_straightness), 
                     case_when(
                       ok(sq_disp > median_sq_disp) ~ 'mig',
                       is.nan(sq_disp) ~ Breed,
                       TRUE ~ Breed
                     ), 
                     Breed)
    )
  # Breeding area classified as wintering area
  segs6 <- segs5 %>%
    dplyr::mutate(Breed = ifelse(ok(Breed == 'W' & R2n_mean < mean_R2n_mean), 'B', Breed))
  
  # Migrating area classified as wintering area
  segs7 <- segs6 %>%
    dplyr::mutate(Breed = ifelse(ok(Breed == 'W' & sq_disp > median_sq_disp & 
                                   sinuosity < adv$corr_sinuosity &
                                   straightness > adv$corr_straightness), 'mig', Breed))
  
  # Fast wintering segments at the start or end of the wintering stage
  # (e.g. a loop before the return migration) → migration
  segs7$Breed <- fast_w_to_mig(segs7, adv)

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
      Breed = ifelse(ok(row_number() != 1 & row_number() != nrow(.) &
                             Breed == 'B' & lead(Breed) == 'mig' &
                             lag(Breed) == 'mig'),'W', Breed)) %>%
    mutate(
      Breed_0 = df_seg$Breed,
      Agreement_0 = df_seg$Agreement
    )
  # Agreement of the final classes (also for the corrected segments)
  out$Agreement <- agreement_score(out, R2n = out$R2n_limit[1], adv = adv)
  
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
    world <- ne_countries(scale = "medium", returnclass = "sf")
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

  plot_list <- list(segclust = segclust, map = map, segclass = segclass_plot, breed = breed_plot)
  return(plot_list)
}



################################################################################
# classify_track(): full step 1 for ONE track
################################################################################
# Runs, for one track, the same steps shown in 01_Segmentation_classification.qmd:
# projection, cleaning, segmentation/clustering, classification of segments,
# equinox flags, corrections and classification of each position.
#
# Arguments:
#   track  data frame with columns Date_Time (POSIXct, UTC), lon, lat
#   info   one row of the metadata (geo_trip, Ring, Sp, Colony, Colony.lon,
#          Colony.lat)
#   Kmax, lmin, ncluster, eq_window_days, adv   parameters (parameters.R)
#   make_plots  TRUE to return the maps of the segmentation
#
# Returns a list with: segments (one row per segment), summary (one row per
# track), track (one row per position) and plot (or NULL).

classify_track <- function(track, info, Kmax, lmin, ncluster, eq_window_days,
                           adv, use_equinox = TRUE, make_plots = FALSE) {

  geo_trip <- info$geo_trip
  col.lon  <- info$Colony.lon
  col.lat  <- info$Colony.lat
  track    <- track %>% arrange(Date_Time)

  start_year <- year(track$Date_Time[1])
  end_year   <- year(track$Date_Time[nrow(track)])

  # --- Projection centred on the colony (km) and trajectory -----------------
  proj_colony <- paste0('+proj=aeqd +lat_0=', col.lat, ' +lon_0=', col.lon,
                        ' +units=km +datum=WGS84')
  xy_km <- sf::sf_project(from = sf::st_crs(4326), to = proj_colony,
                          pts = as.matrix(track[, c('lon', 'lat')]))
  ltraj <- adehabitatLT::as.ltraj(xy = xy_km, date = track$Date_Time,
                                  id = geo_trip)
  track_df <- ltraj[[1]]
  pos <- match(track_df$date, track$Date_Time)
  track_df <- track_df %>%
    mutate(lon = track$lon[pos], lat = track$lat[pos], R2n = x^2 + y^2)

  # --- Equinoxes (only GLS data) ---------------------------------------------
  # First September and March equinoxes after the start of the track (this
  # also works for tracks that do not start in summer, e.g. winter breeders)
  next_date <- function(start, month_day) {
    d <- as.Date(paste0(year(start), month_day))
    if (d < as.Date(start)) d <- as.Date(paste0(year(start) + 1, month_day))
    d
  }
  Sep_eq <- next_date(track$Date_Time[1], '-09-21')
  Mar_eq <- next_date(track$Date_Time[1], '-03-21')
  eq_days <- if (use_equinox) eq_window_days else NULL

  # --- Cleaning --------------------------------------------------------------
  track_df2 <- clean_positions(data = track_df, dist_threshold = adv$max_step)
  track_df3 <- remove_duplicates(data = track_df2, equinox_days = eq_days,
                                 Sep_eq = Sep_eq, Mar_eq = Mar_eq,
                                 dup_threshold = adv$max_repeated,
                                 edge_positions = adv$edge_positions)

  # --- Segmentation and clustering -------------------------------------------
  mode_segclust <- segclust(track_df3,
                            Kmax = Kmax, lmin = lmin, ncluster = ncluster,
                            seg.var = c('x', 'y'), order.var = 'x',
                            diag.var = c('dx', 'dy', 'R2n', 'x', 'y', 'dist'),
                            subsample = FALSE)
  bic_max <- BIC(mode_segclust) %>% filter(BIC == max(BIC))

  # --- Classification of segments --------------------------------------------
  add.metrics <- c(names(adv$paired_metrics), adv$extra_metrics)
  segments1 <- segment_class(mode_segclust = mode_segclust, option.by = 'seg',
                             add.metrics = add.metrics,
                             col.lon = col.lon, col.lat = col.lat,
                             adv = adv, lmin = lmin)
  segments2 <- equinox_interference(df_seg = segments1$df_clust,
                                    track = track_df3, Sep_eq, Mar_eq,
                                    equinox_days = eq_days)
  segments3 <- segment_correction(df_seg = segments2, adv = adv) %>%
    flag_segments(adv = adv)
  segments3 <- cbind(ID = geo_trip, Sp = info$Sp, segments3)

  seg_df <- data.frame(Ring = info$Ring, Sp = info$Sp, Colony = info$Colony,
                       ID = geo_trip,
                       n_seg = bic_max$nseg, n_cluster = bic_max$ncluster,
                       Breed_class = paste(segments3$Breed, collapse = '-'),
                       n_seg_B = sum(segments3$Breed == 'B', na.rm = TRUE),
                       n_seg_mig = sum(segments3$Breed == 'mig', na.rm = TRUE),
                       n_seg_W = sum(segments3$Breed == 'W', na.rm = TRUE),
                       Agreement_0 = round(mean(segments3$Agreement_0, na.rm = TRUE), 3),
                       Agreement = round(mean(segments3$Agreement, na.rm = TRUE), 3),
                       Metrics = paste(segments1$metrics, collapse = '-'),
                       n_check = sum(segments3$Check == 'Check'),
                       Direct_BW = any((segments3$Breed == 'W' & lag(segments3$Breed) == 'B') |
                                         (segments3$Breed == 'B' & lag(segments3$Breed) == 'W'),
                                       na.rm = TRUE))

  # --- Class of each position ------------------------------------------------
  track_df5 <- segclust2d::augment(mode_segclust) %>%
    mutate(Breed = segments3$Breed[findInterval(row_number(), segments3$begin)]) %>%
    distinct() %>%
    add_column(ID = geo_trip, .before = 1) %>%
    mutate(Breed2 = case_when(
      Breed == 'mig' & cumsum(replace_na(Breed == 'W', FALSE)) > 0 &
        rev(cumsum(rev(replace_na(Breed == "W", FALSE)))) == 0 ~ 'mig2',
      Breed == "mig" & rev(cumsum(rev(replace_na(Breed == "W", FALSE)))) > 0 ~ "mig1",
      TRUE ~ Breed)) %>%
    mutate(x_km = x, y_km = y, x = lon, y = lat) %>%
    dplyr::select(-lon, -lat) %>%
    mutate(Equinox = case_when(
      !use_equinox ~ 'No',
      date >= (Sep_eq - days(eq_window_days)) &
        date <= (Sep_eq + days(eq_window_days)) ~ 'Sep',
      date >= (Mar_eq - days(eq_window_days)) &
        date <= (Mar_eq + days(eq_window_days)) ~ 'Mar',
      TRUE ~ 'No'))

  # --- Optional plots ----------------------------------------------------------
  plot_out <- NULL
  if (make_plots) {
    map_seg <- segmap(mode_segclust, pointsize = 3, coord.names = c("x", "y"),
                      interactive = FALSE)
    plot_list <- plot_segmentation(segclust = NULL, map = map_seg,
                                   segclass = NULL, col.lon = 0, col.lat = 0,
                                   breed = TRUE, data = track_df5, world = NULL)
    plot_out <- plot_list$map | plot_list$breed
  }

  list(segments = segments3, summary = seg_df, track = track_df5,
       plot = plot_out)
}
