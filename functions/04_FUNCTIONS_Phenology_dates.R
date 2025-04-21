################################################################################

check_phen <- function(data, column, min_length = 5) {
  # Apply run-length encoding (RLE) to the specified column
  rle_ph <- rle(data[[column]])  
  
  # Initialize position in the original vector
  pos <- 1  
  
  # Iterate through each sequence in RLE
  for (l in seq_along(rle_ph$lengths)) {
    # Skip the first and last sequences
    if (l == 1 || l == length(rle_ph$lengths)) {
      pos <- pos + rle_ph$lengths[l]
      next
    }
    
    # Skip sequences longer than or equal to the minimum length
    if (rle_ph$lengths[l] >= min_length) {
      pos <- pos + rle_ph$lengths[l]
      next
    }
    
    # Check if the previous and next sequences have the same value
    if (rle_ph$values[l - 1] == rle_ph$values[l + 1]) {
      # Replace the short sequence with the value of the surrounding sequences
      data[[column]][pos:(pos + rle_ph$lengths[l] - 1)] <- rle_ph$values[l - 1]
    }
    
    # Update position in the vector
    pos_actual <- pos_actual + rle_ph$lengths[l]
  }
  
  return(data)  # Return modified data
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

calculate_kde <- function(ph_seg, proj= 4326) {
  # Filter non-migratory points
  non_migratory <- ph_seg[!(ph_seg$Breed3 %in% c("mig1", "mig2")), ]
  
  if (nrow(non_migratory) <= 5) {
    # Not enough points for KDE
    return(NULL)
  }
  
  # Project tracks for KDE calculation
  p_tracks <- projectTracks(dataGroup = non_migratory, projType = 'azim', custom = FALSE)
  hVals <- findScale(tracks = p_tracks, scaleARS = FALSE, sumTrips = NULL)
  
  # Calculate KDE
  KDE <- estSpaceUse(tracks = p_tracks, scale = hVals$href, levelUD = 90, polyOut = TRUE)
  KDE90_sf <- KDE$UDPolygons
   # mapKDE(KDE$UDPolygons)
  # Handle alternative KDE calculation if KDE90_sf is NULL
  if (is.null(KDE90_sf)) {
    kuds <- kernelUD(p_tracks, h = (hVals$href * 1000), grid = 1000, extent = 4, same4all = FALSE)
    KDE90 <- getverticeshr(kuds, percent = 90)
    KDE90_sf <- st_as_sf(spTransform(KDE90, CRS(proj)))
  }
  
  return(KDE90_sf)
}

# Function to find first and last points inside KDE
find_dates <- function(ph_seg, KDE90_sf, proj = proj) {
  if (is.null(KDE90_sf)) {
    # If KDE is NULL, return NA for points and dates
    return(list(arr_point = data.frame(Longitude = NA, Latitude = NA), arr_date = NA,
                dep_point = data.frame(Longitude = NA, Latitude = NA), dep_date = NA
    ))
  }
  
  ph_seg_sf <- ph_seg%>%
    st_as_sf(coords = c("Longitude", "Latitude"),
             crs = proj, agr = "constant", remove = FALSE)
  
  # Filter ph_seg using KDE polygons
  KDE90_points <- st_filter(ph_seg_sf, KDE90_sf)
  KDE90_points <- st_drop_geometry(KDE90_points)
  
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
    if (k == "B1") {
      # Calculate departure only
      dep_index <- min(which(ph_seg$DateTime == max(KDE90_points$DateTime)))
      dep_point <- ph_seg[dep_index + 1, c("Longitude", "Latitude")]
      dep_date <- as.Date(ph_seg$DateTime[dep_index + 1])
    } else if (k == "B2") {
      # Calculate arrival only
      arr_index <- min(which(ph_seg$DateTime == min(KDE90_points$DateTime)))
      arr_point <- ph_seg[arr_index, c("Longitude", "Latitude")]
      arr_date <- as.Date(ph_seg$DateTime[arr_index])
    } else if (grepl("^W", k)) {
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
