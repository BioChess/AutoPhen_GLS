# =============================================================================
# GLS SUPERVISED CALIBRATION ANALYSES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Calibration of GLS light loggers using rooftop calibration
#              method. Extracts zenith angles and error distribution parameters
#              (alpha shape and scale) for each calibration period.
# =============================================================================


# =============================================================================
# 1. INSTALL PACKAGES (run only the first time) ----
# =============================================================================


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(devtools, remotes, install = TRUE)

# Packages from GitHub (run only once)
# pacman::p_load_gh(
#   "SWotherspoon/SGAT",
#   "SWotherspoon/BAStag",
#   "bmcclintock/momentuHMM",
#   "SLisovski/TwGeos",
#   "SLisovski/GeoLight"
# )


# =============================================================================
# 2. LOAD PACKAGES AND FUNCTIONS ----
# =============================================================================

pacman::p_load(
  'devtools', 'remotes',       # installing external packages
  'readxl',                    # opening Excels and csv
  'tidyverse',                 # managing data
  'lubridate',                 # dealing with dates
  'evaluate',                  # check errors
  'maps', 'rworldmap',         # creating maps
  'GeoLight', 'BAStag', 'TwGeos', # geolocators
  'MASS'                        # fitdistr, needed by thresholdCalibration
)

# Detach purrr to avoid conflicts with base map()
if ("purrr" %in% .packages()) {
  detach("package:purrr", unload = TRUE)
}

# Load custom functions
source("functions/01_FUNCTIONS_GLS_analysis.R")


# =============================================================================
# 3. READ METADATA AND DEFINE THRESHOLDS ----
# =============================================================================

mdata     <- read_csv('input/GLS.files.csv',       col_names = TRUE)
info_file <- read_csv('input/GLS.info.csv',        col_names = TRUE)
calib.info <- read_csv('input/Info.calibrations.csv')

# Light threshold: value above which it is considered daytime (typically 1.5-2)
threshold <- 2

# Map extent (used as visual aid during calibration)
xlim <- c(-20, 0)
ylim <- c(15, 40)
map  <- rworldmap::getMap(resolution = "coarse")

# Method for thresholdCalibration ("gamma" or "log-normal")
method <- "gamma"


# =============================================================================
# 4. MAIN LOOP: CALIBRATION PARAMETERS PER GLS ----
# =============================================================================
# For each GLS, up to 3 calibration periods are processed:
#   - Step 1: pre-deployment calibration 1
#   - Step 2: pre-deployment calibration 2 (if available)
#   - Step 3: post-recovery calibration
# =============================================================================

gls.dir        <- paste0(getwd(), '/input/GLS/')
calib.list_SGAT <- list()  # stores one row per GLS

for (i in seq_len(nrow(mdata))) {
# i = 1
  # --- 4.1. Read GLS information -------------------------------------------

  info      <- mdata[i, ]
  lightfile <- if_else(!is.na(info$File.adjlux), info$File.adjlux, info$File.lux)
  lightfile_c <- paste0(gls.dir, lightfile)

  if (!file.exists(lightfile_c) || is.na(lightfile_c)) {
    message('Skipping geo (', i, '): light file not found')
    next
  }

  geoID     <- info$geo_alfanum
  Producer  <- info$Producer
  col.name  <- info$Colony
  sp        <- info$Sp
  year_rec  <- info$Year_rec
  start_date <- as.Date(info$start_date, tz = 'GMT')
  end_date   <- as.Date(info$end_date,   tz = 'GMT')

  cat('Processing geo (', i, ')', geoID, 'for year', year_rec, '\n')

  # Calibration info for this GLS
  calib <- calib.info %>% filter(geo_alfanum == geoID)

  if (nrow(calib) == 0) {
    message("Skipping geo (", i, ") ", geoID, ": no calibration info found")
    next
  }
  
  # Update Calibrations column: "Yes" if at least one calibration period has coordinates
  has_calib_period <- any(!is.na(calib[, c("Lon1", "Lon2", "Lon3")]))
  mdata[i, "Calibrations"] <- ifelse(has_calib_period, "Yes", "No")
  
  # --- 4.2. Read light data -------------------------------------------------

  lu <- read_lightGLS(
    producer   = Producer,
    ID         = geoID,
    file       = lightfile_c
    # start_date = start_date  # uncomment to subset to deployment period
    # end_date   = end_date
  )

  if (is.null(lu) || nrow(lu) == 0) {
    message('Skipping geo (', i, '): no light data')
    next
  }

  # --- 4.3. Loop through calibration periods --------------------------------

  twl_error.list <- list()

  for (step in 1:3) {

    calib_start <- as.POSIXct(calib[[paste0("calib_start", step)]], tz = "GMT")
    calib_end   <- as.POSIXct(calib[[paste0("calib_end",   step)]] + 1, tz = "GMT")
    n_days <- round(as.numeric(difftime(calib_end, calib_start, units = 'days')))
    loc    <- calib[[paste0("loc_calib", step)]]
    lon    <- calib[[paste0("Lon", step)]]
    lat    <- calib[[paste0("Lat", step)]]

    # Column suffix: PD = pre-deployment, PR = post-recovery
    step_label <- if (step %in% c(1, 2)) paste0("PD", step) else "PR3"
    # File name label
    step_name  <- if (step %in% c(1, 2)) "predep" else "postrec"

    calib_SGAT <- c(NA_real_, NA_real_, NA_real_, NA_real_)

    # length() == 0 when the column exists but has no value for this GLS
    # is.na()  == TRUE when the value exists but is NA
    if (length(loc) > 0 && length(calib_end) > 0 &&
        !is.na(loc)  && !is.na(calib_end)) {

      # Open an external graphics window: the RStudio plot pane does not support
      # mouse/keyboard event handling required by preprocessLight() inside calib_period()
      if (.Platform$OS.type == "windows") {
        windows()
      } else if (Sys.info()[["sysname"]] == "Darwin") {
        quartz()
      } else {
        x11()
      }

      # Interactive calibration: review light image, then close the window to continue
      twl <- calib_period(
        data      = lu,
        start     = calib_start,
        end       = calib_end,
        lon       = lon,
        lat       = lat,
        threshold = threshold
      )

      # Close external device and return focus to RStudio
      if (dev.cur() > 1) dev.off()

      if (!is.data.frame(twl) || nrow(twl) == 0) {
        message('No valid twilights for step ', step, ' of geo ', geoID)

      } else {

        # Save .trn file for this calibration period
        trn <- twl %>%
          mutate(
            date_time  = Twilight,
            type       = ifelse(Rise == TRUE, 'Sunrise', 'Sunset'),
            confidence = Marker
          ) %>%
          dplyr::select(date_time, type, confidence)

        if (!dir.exists('output/calibrations'))
          dir.create('output/calibrations', recursive = TRUE, showWarnings = FALSE)

        file_trn <- paste(geoID, sp, year_rec, col.name, step_name, '.trn', sep = '_')
        write.table(trn,
                    file      = paste0('output/calibrations/', file_trn),
                    sep       = ',',
                    col.names = FALSE,
                    row.names = FALSE,
                    quote     = FALSE)
        cat('  Trn saved for calibration period:', step_name, '\n')

        # Remove deleted twilights and ensure alternating Rise/FALSE pattern
        # before SGAT calibration to avoid NA in twilight deviations
        twl_clean <- twl %>%
          filter(Deleted == FALSE) %>%
          arrange(Twilight)

        # Drop isolated twilights that break the sunset/sunrise alternation
        valid <- rle(twl_clean$Rise)$lengths
        if (any(valid > 1)) {
          keep <- !duplicated(twl_clean$Rise) | c(diff(twl_clean$Rise) != 0, TRUE)
          twl_clean <- twl_clean[keep, ]
        }

        # Ensure even number of rows (paired events)
        if (nrow(twl_clean) %% 2 != 0) twl_clean <- twl_clean[-nrow(twl_clean), ]

        # SGAT threshold calibration
        calib_SGAT <- try(
          thresholdCalibration(twl_clean$Twilight, twl_clean$Rise, lon, lat, method = method),
          silent = TRUE
        )
        if (inherits(calib_SGAT, 'try-error')) {
          message('  Error in SGAT calibration at step ', step, ' of geo ', geoID)
          calib_SGAT <- c(NA_real_, NA_real_, NA_real_, NA_real_)
        }
      }

    } else {
      message('  No calibration information for step ', step, ' of geo ', geoID)
    }

    # Build one row for this calibration step
    twl_error <- bind_cols(
      geoID      = geoID,
      Producer   = Producer,
      Year       = year_rec,
      Colony     = col.name,
      n_days     = n_days,
      zenith     = calib_SGAT[1],
      zenith0    = calib_SGAT[2],
      alpha.mean = calib_SGAT[3],
      alpha.sd   = calib_SGAT[4],
      Loc        = loc
    )
    # Rename metric columns with period-specific suffix
    colnames(twl_error)[5:10] <- paste0(colnames(twl_error)[5:10], '_', step_label)

    twl_error.list[[step]] <- twl_error
  }

  # Collapse the 3 steps into one row per GLS
  twl_error_GLS <- bind_rows(twl_error.list) %>%
    group_by(geoID, Producer, Year, Colony) %>%
    summarise(across(everything(), ~ na.omit(.x)[1]), .groups = 'drop')

  calib.list_SGAT[[i]] <- twl_error_GLS
  
  # Update Calib.param column: "Yes" if thresholdCalibration succeeded for at least one step
  has_calib_param <- any(sapply(twl_error.list, function(x) {
    !all(is.na(x[, grep("zenith_", colnames(x), value = TRUE)]))
  }))
  mdata[i, "Calib.param"] <- ifelse(has_calib_param, "Yes", "No")
  cat('Geo (', i, ')', geoID, 'done\n\n')
}

write_csv(mdata, file = 'input/GLS.files.csv')
cat('Metadata updated with Calibrations and Calib.param columns\n')

# =============================================================================
# 5. JOIN ALL GLS CALIBRATION RESULTS ----
# =============================================================================

calib.all_SGAT <- bind_rows(calib.list_SGAT)

str(calib.all_SGAT)

# Convert parameter columns to numeric (already numeric from loop, just enforce type)
calib.all_SGAT <- calib.all_SGAT %>%
  mutate(across(c(3, 5:8, 10:13, 15:18), as.numeric))

# =============================================================================
# 6. CHECK AND REMOVE OUTLIERS ----
# =============================================================================

vars <- c(
  "zenith_PD1",     "zenith_PD2",     "zenith_PR3",
  "zenith0_PD1",    "zenith0_PD2",    "zenith0_PR3",
  "alpha.mean_PD1", "alpha.mean_PD2", "alpha.mean_PR3",
  "alpha.sd_PD1",   "alpha.sd_PD2",   "alpha.sd_PR3"
)

prd <- c("Migrate Technology", "Biotrack")

calib_list <- list()

for (p in prd) {
  cal.df <- calib.all_SGAT %>% filter(Producer == p)
  cal.df <- calib_out(cal.df, vars)
  calib_list[[p]] <- cal.df
}

calib.all <- bind_rows(calib_list) %>% arrange(geoID)


# =============================================================================
# 7. CALCULATE MEAN CALIBRATION PARAMETERS PER GLS ----
# =============================================================================

calib.all_SGAT2 <- calib.all %>%
  group_by(Year, Colony, Producer) %>%
  mutate(
    mean_zenith = ifelse(
      is.na(zenith_PD1) | is.na(zenith_PD2) | is.na(zenith_PR3),
      mean(c(zenith_PD1, zenith_PD2, zenith_PR3), na.rm = TRUE),
      zenith_PD1
    ),
    mean_zenith0 = ifelse(
      is.na(zenith0_PD1) | is.na(zenith0_PD2) | is.na(zenith0_PR3),
      mean(c(zenith0_PD1, zenith0_PD2, zenith0_PR3), na.rm = TRUE),
      zenith0_PD1
    ),
    mean_alpha.mean = ifelse(
      is.na(alpha.mean_PD1) | is.na(alpha.mean_PD2) | is.na(alpha.mean_PR3),
      mean(c(alpha.mean_PD1, alpha.mean_PD2, alpha.mean_PR3), na.rm = TRUE),
      alpha.mean_PD1
    ),
    mean_alpha.sd = ifelse(
      is.na(alpha.sd_PD1) | is.na(alpha.sd_PD2) | is.na(alpha.sd_PR3),
      mean(c(alpha.sd_PD1, alpha.sd_PD2, alpha.sd_PR3), na.rm = TRUE),
      alpha.sd_PD1
    )
  ) %>%
  filter(!duplicated(geoID))

# Diagnostic plots
par(mfrow = c(3, 2))
boxplot(calib.all_SGAT2[, 23:24], main = "Zenith angles")
boxplot(calib.all_SGAT2[, 25:26], main = "Alpha parameters")
hist(calib.all_SGAT2$mean_zenith,      main = "Mean zenith",      xlab = "")
hist(calib.all_SGAT2$mean_zenith0,     main = "Mean zenith0",     xlab = "")
hist(calib.all_SGAT2$mean_alpha.mean,  main = "Mean alpha.mean",  xlab = "")
hist(calib.all_SGAT2$mean_alpha.sd,    main = "Mean alpha.sd",    xlab = "")

write_csv(calib.all_SGAT2, file = 'output/Supervised_calibrations_SGAT.csv')
cat('Calibration parameters saved to output/Supervised_calibrations_SGAT.csv\n')


# =============================================================================
# 8. SUMMARY PER COLONY, YEAR AND PRODUCER ----
# =============================================================================
# Useful for GLS that were not calibrated or have invalid calibration periods

calib_SGAT_summary <- calib.all_SGAT2 %>%
  group_by(Colony, Year, Producer) %>%
  dplyr::summarise(
    n_geos          = n_distinct(geoID),
    mean_zenith     = mean(mean_zenith,     na.rm = TRUE),
    mean_zenith0    = mean(mean_zenith0,    na.rm = TRUE),
    mean_alpha.mean = mean(mean_alpha.mean, na.rm = TRUE),
    mean_alpha.sd   = mean(mean_alpha.sd,   na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Year) %>%
  distinct()

write_csv(calib_SGAT_summary, file = 'output/Summary_SGAT_calibrations.csv')
cat('Summary saved to output/Summary_SGAT_calibrations.csv\n')
