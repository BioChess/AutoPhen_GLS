# =============================================================================
# 00. SETUP: INSTALL ALL REQUIRED PACKAGES
# =============================================================================
# Author: Diego Vicente-Sastre
# Description: Installs every package needed to run the AutoPhen_GLS workflow.
#              Run this script ONCE, before the other scripts.
#
#              Two options:
#                A) (Recommended) Restore the exact package versions recorded
#                   in renv.lock. This is the most reliable way to reproduce
#                   the results of the paper.
#                B) Install the current versions from CRAN (and, optionally,
#                   the pinned versions of the GitHub-only packages needed for
#                   GLS pre-processing). Use this option if option A fails.
#
# Requirements: R >= 4.6 (tested with R 4.6.1).
#               On Windows, Rtools is needed to install packages from GitHub
#               (https://cran.r-project.org/bin/windows/Rtools/).
# =============================================================================


# =============================================================================
# 1. SETTINGS ----
# =============================================================================

# Use renv.lock to install the exact package versions (option A)?
use_renv <- TRUE

# Install also the packages for the optional GLS pre-processing scripts
# (folder gls_preprocessing/)? Not needed for the main workflow.
install_gls_preprocessing <- FALSE

if (getRversion() < "4.6.0") {
  warning("This workflow was tested with R 4.6.1. You are using R ",
          getRversion(), ". Some packages may not install or work as expected.")
}


# =============================================================================
# 2. PACKAGE LISTS ----
# =============================================================================

# --- 2.1. Main workflow (steps 1-3) -------------------------------------------
core_pkgs <- c(
  # Installing packages
  "remotes", "pacman",
  # Data handling
  "tidyverse", "plyr", "lubridate", "readxl", "units",
  # Spatial data and maps
  "sf", "spData", "rnaturalearth", "rnaturalearthdata",
  # Segmentation, classification and trajectory metrics (step 1)
  "adehabitatLT", "adehabitatHR", "segclust2d", "trajr", "recipes", "furrr",
  "caTools",
  # Interactive reclassification (step 2)
  "shiny", "leaflet", "DT",
  # Kernel density estimation (step 3)
  "track2KBA",
  # Plots
  "ggplot2", "ggpubr", "ggspatial", "cowplot", "patchwork",
  # Simulations
  "marcher"
)

# --- 2.2. Optional GLS pre-processing (folder gls_preprocessing/) -------------
gls_cran_pkgs <- c(
  "evaluate", "MASS", "suncalc", "roll", "mice", "DescTools", "CircStats",
  "baytrends", "chron", "geosphere", "sp", "maps", "rworldmap"
)

# Packages only available on GitHub. Each one is pinned to a specific commit
# so that everybody installs exactly the same code, even if the packages
# change (or stop being maintained) in the future.
gls_github_pkgs <- c(
  "SLisovski/GeoLight@17afeabdd047",   # GeoLight 2.1.0  (2024-05-29)
  "SLisovski/TwGeos@028f4992c306",     # TwGeos 0.1.2    (2024-01-16)
  "SWotherspoon/BAStag@f49584a64f98",  # BAStag 0.1.3    (2021-03-29)
  "SWotherspoon/SGAT@d0894847fb9b"     # SGAT 0.1.3      (2022-05-25)
)


# =============================================================================
# 3. OPTION A (RECOMMENDED): RESTORE EXACT VERSIONS WITH renv ----
# =============================================================================
# renv reads renv.lock and installs the same package versions we used, inside
# a project-specific library. It does not change the packages installed in
# your main R library.

if (use_renv) {
  if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

  if (file.exists("renv.lock")) {
    renv::restore(prompt = FALSE)
    message("Packages restored from renv.lock.")
  } else {
    message("renv.lock not found. Using option B instead.")
    use_renv <- FALSE
  }
}


# =============================================================================
# 4. OPTION B: INSTALL PACKAGES MANUALLY ----
# =============================================================================

if (!use_renv) {

  cran_pkgs <- core_pkgs
  if (install_gls_preprocessing) cran_pkgs <- c(cran_pkgs, gls_cran_pkgs)

  missing_pkgs <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

  if (install_gls_preprocessing) {
    for (pkg in gls_github_pkgs) {
      remotes::install_github(pkg, upgrade = "never", dependencies = TRUE)
    }
  }

  message("Packages installed.")
}


# =============================================================================
# 5. FOR MAINTAINERS: CREATE OR UPDATE renv.lock ----
# =============================================================================
# Run these lines (only the maintainer, not users) after installing and
# testing all packages with option B (with install_gls_preprocessing <- TRUE).
# They create renv.lock, which records the exact version (and GitHub commit)
# of every package used in the repository.
#
# renv::init(bare = TRUE)          # creates the renv/ folder, no packages yet
# renv::install(c(core_pkgs, gls_cran_pkgs, gls_github_pkgs))
# renv::snapshot(type = "all")     # writes renv.lock
#
# After this, commit renv.lock, .Rprofile and renv/activate.R to GitHub.

