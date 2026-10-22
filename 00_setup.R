# =============================================================================
# 00. SETUP: INSTALL ALL REQUIRED PACKAGES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Installs every package needed to run the AutoPhen_GLS workflow.
#              Run this script ONCE, before scripts 01-06.
#
#              Two options:
#                A) (Recommended) Restore the exact package versions recorded
#                   in renv.lock. This is the most reliable way to reproduce
#                   the results of the paper.
#                B) Install the current versions from CRAN and the pinned
#                   versions of the GitHub-only packages. Use this option if
#                   option A fails on your system.
#
# Requirements: R >= 4.6 (tested with R 4.6.1).
#               On Windows, Rtools is needed to install packages from GitHub
#               (https://cran.r-project.org/bin/windows/Rtools/).
# =============================================================================


# =============================================================================
# 1. CHECK R VERSION ----
# =============================================================================

if (getRversion() < "4.6.0") {
  warning("This workflow was tested with R 4.6.1. You are using R ",
          getRversion(), ". Some packages may not install or work as expected.")
}


# =============================================================================
# 2. OPTION A (RECOMMENDED): RESTORE EXACT VERSIONS WITH renv ----
# =============================================================================
# renv reads renv.lock and installs the same package versions we used, inside
# a project-specific library. It does not change the packages installed in
# your main R library.

use_renv <- TRUE   # set to FALSE to use option B

if (use_renv) {
  if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

  if (file.exists("renv.lock")) {
    renv::restore(prompt = FALSE)
    message("Packages restored from renv.lock. You can now run scripts 01-06.")
  } else {
    message("renv.lock not found. Using option B instead.")
    use_renv <- FALSE
  }
}


# =============================================================================
# 3. OPTION B: INSTALL PACKAGES MANUALLY ----
# =============================================================================

if (!use_renv) {

  # --- 3.1. Packages from CRAN ------------------------------------------------
  cran_pkgs <- c(
    # Installing packages
    "remotes", "pacman",
    # Data handling
    "tidyverse", "plyr", "lubridate", "readxl", "evaluate", "units",
    # GLS processing (steps 1-3)
    "MASS", "suncalc", "roll", "mice", "DescTools", "CircStats",
    "baytrends", "chron", "geosphere",
    # Spatial data and maps
    "sf", "sp", "spData", "maps", "rworldmap",
    "rnaturalearth", "rnaturalearthdata",
    # Segmentation, classification and trajectory metrics (step 4)
    "adehabitatLT", "adehabitatHR", "segclust2d", "trajr", "furrr", "caTools",
    # Kernel density estimation (step 6)
    "track2KBA",
    # Interactive reclassification (step 5)
    "shiny", "leaflet", "DT",
    # Plots
    "ggplot2", "ggpubr", "ggspatial", "cowplot", "patchwork",
    # Simulations
    "marcher"
  )

  missing_pkgs <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

  # --- 3.2. Packages only available on GitHub ---------------------------------
  # These packages are not on CRAN. We pin each one to a specific commit so
  # that everybody installs exactly the same code, even if the packages
  # change (or stop being maintained) in the future.
  github_pkgs <- c(
    "SLisovski/GeoLight@17afeabdd047",   # GeoLight 2.1.0  (2024-05-29)
    "SLisovski/TwGeos@028f4992c306",     # TwGeos 0.1.2    (2024-01-16)
    "SWotherspoon/BAStag@f49584a64f98",  # BAStag 0.1.3    (2021-03-29)
    "SWotherspoon/SGAT@d0894847fb9b"     # SGAT 0.1.3      (2022-05-25)
  )

  for (pkg in github_pkgs) {
    remotes::install_github(pkg, upgrade = "never", dependencies = TRUE)
  }

  message("Packages installed. You can now run scripts 01-06.")
}


# =============================================================================
# 4. FOR MAINTAINERS: CREATE OR UPDATE renv.lock ----
# =============================================================================
# Run these lines (only the maintainer, not users) after installing and
# testing all packages with option B. They create renv.lock, which records the
# exact version (and GitHub commit) of every package used by the workflow.
#
# renv::init(bare = TRUE)          # creates the renv/ folder, no packages yet
# renv::install(c(cran_pkgs, github_pkgs))
# renv::snapshot(type = "all")     # writes renv.lock
#
# After this, commit renv.lock, .Rprofile and renv/activate.R to GitHub.
