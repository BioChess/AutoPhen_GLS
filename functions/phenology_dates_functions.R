################################################################################

check_phen <- function(data, column, min_length = 5) {
  # Short sequences (fewer than 'min_length' positions) that are surrounded by
  # the same class are given that class. For example, B-B-mig-B-B becomes
  # B-B-B-B-B. The first and last sequences are never changed.
  x <- data[[column]]
  rle_ph <- rle(x)
  n_seq  <- length(rle_ph$lengths)
  ends   <- cumsum(rle_ph$lengths)
  starts <- ends - rle_ph$lengths + 1

  if (n_seq >= 3) {
    for (l in 2:(n_seq - 1)) {
      if (rle_ph$lengths[l] >= min_length) next

      prev_val <- rle_ph$values[l - 1]
      next_val <- rle_ph$values[l + 1]

      # Missing classes (NA) are not used to fill other sequences
      if (!is.na(prev_val) && !is.na(next_val) && prev_val == next_val) {
        x[starts[l]:ends[l]] <- prev_val
      }
    }
  }

  data[[column]] <- x
  return(data)
}


################################################################################


phen_periods <- function(breed_vector) {
  # Compute run-length encoding (RLE) of the input vector
  rle_breed <- rle(breed_vector)
  
  # Initialize counters for "B" and "W"
  b_count <- 0
  w_count <- 0
  
  # Modify values to create "B1", "B2", etc.
  values_new <- sapply(rle_breed$values, function(x) {
    if (x == "B") {
      b_count <<- b_count + 1
      return(paste0("B", b_count))
    } else if (x == "W") {
      w_count <<- w_count + 1
      return(paste0("W", w_count))
    } else {
      return(x)  # Keep "mig1" and "mig2" unchanged
    }
  })
  
  # Reconstruct the modified vector using inverse RLE
  breed3_vector <- inverse.rle(list(lengths = rle_breed$lengths, values = values_new))
  
  return(breed3_vector)  # Return the modified vector
}

################################################################################

subset_to_kde <- function(data, column) {
  # Apply run-length encoding (RLE) to the specified column
  rle_ph <- rle(data[[column]])
  
  # Create a data frame mapping indices to stages
  stages_df <- data.frame(
    start = cumsum(c(1, head(rle_ph$lengths, -1))),  # Start indices
    end = cumsum(rle_ph$lengths),  # End indices
    stage = rle_ph$values
  )
  
  # Identify unique B/W stages
  bw_stages <- stages_df$stage[grep("^B|^W", stages_df$stage)]
  
  # Function to subset each period with full previous and next stage
  get_subset <- function(label) {
    idx <- which(stages_df$stage == label)  # Find the index of the stage
    if (length(idx) == 0) return(NULL)  # Return NULL if not found
    
    # Identify previous and next stages
    prev_idx <- ifelse(idx > 1, idx - 1, NA)
    next_idx <- ifelse(idx < nrow(stages_df), idx + 1, NA)
    
    # Extract full ranges for previous, current, and next stages
    included_indices <- c(
      if (!is.na(prev_idx)) seq(stages_df$start[prev_idx], stages_df$end[prev_idx]),
      seq(stages_df$start[idx], stages_df$end[idx]),
      if (!is.na(next_idx)) seq(stages_df$start[next_idx], stages_df$end[next_idx])
    )
    
    # Remove NA values and return subset
    included_indices <- na.omit(included_indices)
    return(data[included_indices, ])
  }
  
  # Apply the function to all B/W stages
  subsets <- lapply(bw_stages, get_subset)
  
  # Name each subset with the corresponding stage
  names(subsets) <- bw_stages
  
  return(subsets)
}


################################################################################

calculate_kde <- function(ph_seg, stage, proj = 4326, level = 90,
                          min_positions = 5, bandwidth = 'href',
                          fallback_grid = 1000, fallback_extent = 4,
                          exclude_equinox = TRUE) {
  # Only the positions of this stage (ph_seg also contains the previous and
  # next stages, which are needed to find the dates but not for the KDE)
  stage_points <- ph_seg[ph_seg$Breed3 == stage, ]

  # GLS: positions around the equinoxes have large latitude errors. They are
  # not used for the KDE if enough positions remain.
  if (exclude_equinox && 'Equinox' %in% names(stage_points)) {
    no_eq <- stage_points[stage_points$Equinox %in% 'No', ]
    if (nrow(no_eq) > min_positions) stage_points <- no_eq
  }

  if (nrow(stage_points) <= min_positions) {
    # Not enough points for KDE
    return(NULL)
  }

  # Equal-area projection centred on the positions of the stage
  p_tracks <- suppressMessages(
    projectTracks(dataGroup = stage_points, projType = 'azim', custom = TRUE))
  hVals <- findScale(tracks = p_tracks, scaleARS = FALSE, sumTrips = NULL)

  # Calculate KDE
  KDE <- estSpaceUse(tracks = p_tracks, scale = hVals[[bandwidth]], levelUD = level,
                     polyOut = TRUE)
  KDE90_sf <- KDE$UDPolygons

  # Alternative KDE calculation if estSpaceUse() fails: the grid extent is
  # increased until the KDE fits in the grid
  if (is.null(KDE90_sf)) {
    for (ext in fallback_extent * c(1, 2, 4, 8)) {
      KDE90 <- tryCatch({
        kuds <- kernelUD(p_tracks, h = (hVals[[bandwidth]] * 1000),
                         grid = fallback_grid, extent = ext, same4all = FALSE)
        getverticeshr(kuds, percent = level)
      }, error = function(e) NULL)
      if (!is.null(KDE90)) break
    }
    if (is.null(KDE90)) return(NULL)
    KDE90_sf <- sf::st_transform(sf::st_as_sf(KDE90), crs = proj)
  }

  return(KDE90_sf)
}

# Function to find first and last points inside KDE
find_dates <- function(ph_seg, KDE90_sf, proj = proj, stage) {
  if (is.null(KDE90_sf)) {
    # If KDE is NULL, return NA for points and dates
    return(list(arr_point = data.frame(Longitude = NA, Latitude = NA), arr_date = NA,
                dep_point = data.frame(Longitude = NA, Latitude = NA), dep_date = NA
    ))
  }
  
  ph_seg_sf <- ph_seg%>%
    st_as_sf(coords = c("Longitude", "Latitude"),
             crs = proj, agr = "constant", remove = FALSE)
  
  # Filter ph_seg using KDE polygons. Only the positions of this stage and of
  # the migrations next to it are used: positions of a previous or next
  # breeding/wintering stage can also fall inside the KDE (e.g. when the
  # wintering area is close to the colony), but they are not part of this
  # stage.
  KDE90_points <- st_filter(ph_seg_sf, KDE90_sf)
  KDE90_points <- st_drop_geometry(KDE90_points)
  KDE90_points <- KDE90_points[KDE90_points$Breed3 %in% c(stage, 'mig1', 'mig2'), ]
  
  # Find departure and arrival points
  if (nrow(KDE90_points) == 0) {
    return(list(arr_point = data.frame(Longitude = NA, Latitude = NA), arr_date = NA,
                dep_point = data.frame(Longitude = NA, Latitude = NA), dep_date = NA
                ))
  }
  
  # Initialize results
  dep_point <- data.frame(Longitude = NA, Latitude = NA)
  dep_date <- NA
  arr_point <- data.frame(Longitude = NA, Latitude = NA)
  arr_date <- NA
  
  # Calculate departure or arrival based on period type
  if (nrow(KDE90_points) > 0) {
    if (stage == "B1") {
      # Calculate departure only
      dep_index <- min(which(ph_seg$DateTime == max(KDE90_points$DateTime)))
      dep_point <- ph_seg[dep_index + 1, c("Longitude", "Latitude")]
      dep_date <- as.Date(ph_seg$DateTime[dep_index + 1])
    } else if (stage == "B2") {
      # Calculate arrival only
      arr_index <- min(which(ph_seg$DateTime == min(KDE90_points$DateTime)))
      arr_point <- ph_seg[arr_index, c("Longitude", "Latitude")]
      arr_date <- as.Date(ph_seg$DateTime[arr_index])
    } else if (grepl("^W", stage)) {
      # Calculate both arrival and departure for W periods
      dep_index <- min(which(ph_seg$DateTime == max(KDE90_points$DateTime)))
      dep_point <- ph_seg[dep_index + 1, c("Longitude", "Latitude")]
      dep_date <- as.Date(ph_seg$DateTime[dep_index + 1])
      
      arr_index <- min(which(ph_seg$DateTime == min(KDE90_points$DateTime)))
      arr_point <- ph_seg[arr_index, c("Longitude", "Latitude")]
      arr_date <- as.Date(ph_seg$DateTime[arr_index])
    }
  }
  
  return(list(arr_point = arr_point, arr_date = arr_date, dep_point = dep_point, dep_date = dep_date))
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


phen_map <- function(track, col.lat = NULL, col.lon = NULL, results_list, world) {
  # Ensure Longitude and Latitude are numeric, removing NAs
  track <- track %>%
    mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
    filter(!is.na(Longitude) & !is.na(Latitude))
  
  # Debugging: Check if Longitude and Latitude are numeric
  if (!is.numeric(track$Longitude) || !is.numeric(track$Latitude)) {
    stop("Error: Longitude and Latitude must be numeric.")
  }
  
  # Check if df_format is empty after removing NA values
  if (nrow(track) == 0) {
    stop("Error: df_format contains no valid coordinates after removing NAs.")
  }
  
  df_sf <- track%>%
    st_as_sf(coords = c("Longitude", "Latitude"),
             crs = 4326, agr = "constant", remove = FALSE)
  # Define coordinate sets for mapping
  coordsets <- sf::st_bbox(df_sf)
  coordsets[c("xmin", "ymin")] <- coordsets[c("xmin", "ymin")] - 10
  coordsets[c("xmax", "ymax")] <- coordsets[c("xmax", "ymax")] + 10
  
  # Initialize ggplot
  map_trip <- ggplot() +
    # Add the world landmask
    geom_sf(data = world, fill = grey(0.75), color = grey(0.75), lwd = 0.2) +
    
    # Add the track
    geom_path(data = df_sf, aes(x = Longitude, y = Latitude), color = 'black') +
    
    # Add colony position
    geom_point( 
               aes(x = col.lon, y = col.lat), 
               fill = 'orchid3', color = "black", 
               shape = 23, size = 5) +
    
    # Longitude and Latitude axis settings
    scale_x_continuous(name = "Longitude", breaks = seq(-180, 180, by = 15)) +
    scale_y_continuous(name = "Latitude", breaks = seq(-90, 90, by = 15)) +
    
    # Axis labels and theme
    xlab('Longitude') + ylab('Latitude') +
    my_theme +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_text(size = 8))
  
  # Add KDE and points dynamically
  for (k in names(results_list)) {
    result <- results_list[[k]]
    
    # Skip if there is no KDE data
    if (is.null(result$KDE90_sf)) next
    
    # Define aesthetics based on period type
    kde_color <- ifelse(k == "B1", "orange",
                        ifelse(k == "B2", "orange4", "steelblue"))
    dep_color <- ifelse(k == "B1", "steelblue4",
                        ifelse(k == "B2", "steelblue2", "orange2"))
    arr_color <- ifelse(grepl("^W", k), "orange4", dep_color)
    
    # Add KDE polygon
    map_trip <- map_trip +
      geom_sf(data = result$KDE90_sf, col = kde_color, fill = kde_color, alpha = 0.5)
    
    # Add arrival and departure points safely
    if (!is.null(result$arr_point) && nrow(result$arr_point) > 0) {
      result$arr_point <- result$arr_point %>%
        mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
        filter(!is.na(Longitude) & !is.na(Latitude))
      
      map_trip <- map_trip +
        geom_point(data = result$arr_point, 
                   aes(x = Longitude, y = Latitude), 
                   fill = arr_color, color = arr_color, 
                   shape = 21, size = 4)
    }
    
    if (!is.null(result$dep_point) && nrow(result$dep_point) > 0) {
      result$dep_point <- result$dep_point %>%
        mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
        filter(!is.na(Longitude) & !is.na(Latitude))
      
      map_trip <- map_trip +
        geom_point(data = result$dep_point, 
                   aes(x = Longitude, y = Latitude), 
                   fill = dep_color, color = dep_color, 
                   shape = 21, size = 4)
    }
  }
  
  # Set coordinate limits
  map_trip <- map_trip +
    coord_sf(
      xlim = c(coordsets["xmin"], coordsets["xmax"]), 
      ylim = c(coordsets["ymin"], coordsets["ymax"]), 
      expand = FALSE
    )
  
  return(map_trip)
}

################################################################################

extract_phen_dates <- function(results_list) {
  # Initialize an empty list to store the extracted dates
  kde_dates <- list()
  
  # Extract all unique keys (e.g., "B1", "W1", "W2", etc.)
  phases <- names(results_list)
  
  # Iterate over each phase and extract relevant dates
  for (phase in phases) {
    
    # Handle BA phases
    if (phase == 'B1') {
      kde_dates[[paste0(phase, "_dep_KDE")]] <- results_list[[phase]][["dep_date"]]
    }
    if (phase == 'B2') {
      kde_dates[[paste0(phase, "_arr_KDE")]] <- results_list[[phase]][["arr_date"]]
    }
    if (startsWith(phase, "W")) {
      # Handle WA phases
      if ("arr_date" %in% names(results_list[[phase]])) {
        kde_dates[[paste0(phase, "_arr_KDE")]] <- results_list[[phase]][["arr_date"]]
      }
      if ("dep_date" %in% names(results_list[[phase]])) {
        kde_dates[[paste0(phase, "_dep_KDE")]] <- results_list[[phase]][["dep_date"]]
      }
    }
  }
  
  kde_dates_df <- as.data.frame(kde_dates)
  
  if(!('B2_arr_KDE' %in% colnames(kde_dates_df))){
    kde_dates_df$B2_arr_KDE <- as.Date(NA)
  }
  
  
  kde_dates_df <- kde_dates_df %>%
    mutate(across(everything(), as.Date))
  
  return(kde_dates_df)  # Return the list of extracted dates
}


################################################################################
# phenology_track(): full step 3 for ONE track
################################################################################
# Runs, for one track, the same steps shown in 03_Phenology_dates.qmd: checks
# the sequence of stages, delimits breeding and wintering areas (KDE) and
# extracts arrival and departure dates.
#
# Arguments:
#   track.ph  positions of one track classified by script 01 (or 02), with
#             columns ID, date, x (longitude), y (latitude), Breed2
#   info      one row of the metadata (geo_trip, Ring, Sp, Colony, Colony.lon,
#             Colony.lat)
#   min_stage_length, kde_level, adv   parameters (parameters.R)
#   world     world map (sf) for the plot; make_plot TRUE to return the map
#
# Returns a list with: dates (one row with the phenological dates) and map.

phenology_track <- function(track.ph, info, min_stage_length, kde_level, adv,
                            use_equinox = TRUE, world = NULL,
                            make_plot = FALSE, areas = NULL,
                            area_col = NULL) {

  track.ph <- track.ph %>% arrange(date)
  if (!'Equinox' %in% names(track.ph)) track.ph$Equinox <- 'No'

  # --- Sequence of stages (B1, mig1, W1, ..., mig2, B2) ----------------------
  track.ph2 <- check_phen(track.ph, "Breed2", min_length = min_stage_length)

  # Movements between two wintering stages (without breeding between them)
  # are part of the same wintering stage: W-mig-W becomes one W
  r <- rle(track.ph2$Breed2)
  is_W <- r$values == 'W'
  for (k in seq_along(r$values)) {
    if (r$values[k] %in% c('mig1', 'mig2', 'mig')) {
      prev_W <- which(is_W[seq_len(k - 1)])
      next_W <- which(is_W[-seq_len(k)]) + k
      if (length(prev_W) > 0 && length(next_W) > 0) {
        between <- r$values[max(prev_W):min(next_W)]
        if (!any(between == 'B')) r$values[k] <- 'W'
      }
    }
  }
  track.ph2$Breed2 <- inverse.rle(r)
  track.ph2$Breed3 <- phen_periods(track.ph2$Breed2)

  BW_periods <- unique(track.ph2$Breed3[grep("^B|^W", track.ph2$Breed3)])
  B_periods  <- grep("^B", BW_periods, value = TRUE)
  # If the track does not start with breeding, the first breeding stage is B2
  # (all breeding stages are renumbered: B1 -> B2, B2 -> B3, ...)
  if (length(B_periods) > 0 && B_periods[1] == 'B1' &&
      track.ph2$Breed2[1] != 'B') {
    is_B <- grepl('^B[0-9]+$', track.ph2$Breed3)
    track.ph2$Breed3[is_B] <- paste0('B', as.integer(sub('B', '', track.ph2$Breed3[is_B])) + 1)
  }

  track_f <- formatFields(dataGroup = track.ph2, fieldID = "ID",
                          fieldDateTime = "date", fieldLon = "x",
                          fieldLat = "y")

  Migrant <- ifelse(all(track_f$Breed3 == 'B1'), 'Resident', 'Migrant')

  # --- KDE and dates for each breeding/wintering stage ------------------------
  ph_seg.list  <- subset_to_kde(track_f, "Breed3")
  results_list <- list()
  for (k in names(ph_seg.list)) {
    ph_seg <- ph_seg.list[[k]]
    stg <- ph_seg[!(ph_seg$Breed3 %in% c("mig1", "mig2")), ]
    if (is.null(ph_seg) || nrow(stg) == 0) {
      results_list[[k]] <- list(KDE90_sf = NULL,
                                arr_point = data.frame(Longitude = NA, Latitude = NA),
                                arr_date = NA,
                                dep_point = data.frame(Longitude = NA, Latitude = NA),
                                dep_date = NA)
      next
    }
    KDE90_sf <- calculate_kde(ph_seg, stage = k, proj = 4326, level = kde_level,
                              min_positions = adv$kde_min_positions,
                              bandwidth = adv$kde_bandwidth,
                              fallback_grid = adv$kde_fallback_grid,
                              fallback_extent = adv$kde_fallback_extent,
                              exclude_equinox = use_equinox)
    phen_dates <- find_dates(ph_seg, KDE90_sf, proj = 4326, stage = k)
    results_list[[k]] <- c(list(KDE90_sf = KDE90_sf), phen_dates)
  }
  if (Migrant == 'Resident' && 'B1' %in% names(results_list)) {
    results_list[["B1"]][["dep_date"]] <- as.Date(NA)
  }

  # --- Output ----------------------------------------------------------------
  dates <- data.frame(Ring = info$Ring, Sp = info$Sp, Colony = info$Colony,
                      ID = info$geo_trip, Migrant = Migrant,
                      Year1 = year(track.ph$date[1]),
                      Year2 = year(track.ph$date[nrow(track.ph)]),
                      extract_phen_dates(results_list))

  # --- Dates to review --------------------------------------------------------
  # (1) Dates within an equinox window (GLS): latitude is uncertain, so the
  #     date may be wrong (e.g. arrival at the colony during the March
  #     equinox).
  date_cols <- grep('_KDE$', names(dates), value = TRUE)
  eq_days <- unique(as.Date(track.ph$date[track.ph$Equinox != 'No']))
  in_eq <- date_cols[sapply(date_cols, function(d) {
    x <- as.Date(dates[[d]]); !is.na(x) && x %in% eq_days })]
  dates$Dates_in_equinox <- paste(sub('_KDE$', '', in_eq), collapse = '; ')

  # (2) Overlap between the KDE of consecutive breeding and wintering stages:
  #     arrival and departure are poorly defined when the areas overlap.
  stg <- names(results_list)[grepl('^B|^W', names(results_list))]
  stg <- stg[order(sapply(stg, function(k) {
    i <- which(track.ph2$Breed3 == k); if (length(i)) min(i) else Inf }))]
  overlaps <- character(0)
  for (k in seq_len(max(0, length(stg) - 1))) {
    a <- results_list[[stg[k]]]$KDE90_sf; b <- results_list[[stg[k + 1]]]$KDE90_sf
    if (!is.null(a) && !is.null(b)) {
      ov <- tryCatch(any(lengths(sf::st_intersects(
                       sf::st_union(sf::st_transform(a, 4326)),
                       sf::st_union(sf::st_transform(b, 4326)))) > 0),
                     error = function(e) FALSE)
      if (ov) overlaps <- c(overlaps, paste0(stg[k], '-', stg[k + 1]))
    }
  }
  dates$KDE_overlap <- paste(overlaps, collapse = '; ')
  dates$Check <- ifelse(dates$Dates_in_equinox != '' | dates$KDE_overlap != '',
                        'Check', 'Ok')

  map <- NULL
  if (make_plot) {
    map <- phen_map(track = track_f, col.lat = info$Colony.lat,
                    col.lon = info$Colony.lon, results_list, world)
  }
  # --- Optional: name of each wintering area ----------------------------------
  if (!is.null(areas)) {
    wa <- name_wintering_areas(track.ph2, areas, area_col)
    for (k in seq_len(nrow(wa))) dates[[paste0(wa$Breed3[k], '_area')]] <- wa$Area[k]
  }

  list(dates = dates, map = map)
}


################################################################################
# name_wintering_areas(): OPTIONAL. Gives a name to each wintering stage from a
# map of areas (e.g. marine ecoregions): the area that contains most positions
# of the stage. If no position falls inside any area, the name is NA (no
# default area is assigned).
#
# @track    : positions with columns x (longitude), y (latitude) and Breed3
# @areas    : sf polygons with the areas (any CRS)
# @area_col : column of 'areas' with the name of each area

name_wintering_areas <- function(track, areas, area_col) {
  w <- track[grepl('^W[0-9]+$', track$Breed3), ]
  stages <- unique(w$Breed3)
  if (length(stages) == 0) return(data.frame(Breed3 = character(), Area = character()))
  pts <- sf::st_as_sf(w, coords = c('x', 'y'), crs = 4326, remove = FALSE)
  areas <- sf::st_transform(areas[, area_col], 4326)
  joined <- sf::st_drop_geometry(sf::st_join(pts, areas, left = FALSE))
  out <- data.frame(Breed3 = stages, Area = NA_character_)
  for (k in seq_along(stages)) {
    a <- joined[[area_col]][joined$Breed3 == stages[k]]
    if (length(a) > 0) out$Area[k] <- names(which.max(table(a)))
  }
  out
}
