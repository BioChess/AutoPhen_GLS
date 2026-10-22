# AutoPhen_GLS

**A semi-automated R workflow to extract migration phenology dates from animal tracking data, developed for light-level geolocator (GLS) data.**

AutoPhen_GLS takes the positions of a migratory track (one annual cycle) and returns four key phenological dates:

| Date | Code | Definition |
|---|---|---|
| Departure from the breeding area | `B1_dep` | Exit from the 90% KDE of the pre-migration breeding positions |
| Arrival at the wintering area | `W1_arr` | First position inside the 90% KDE of the wintering positions |
| Departure from the wintering area | `W1_dep` | Exit from the 90% KDE of the wintering positions |
| Arrival at the breeding area | `B2_arr` | First position inside the 90% KDE of the post-migration breeding positions |

If an individual uses more than one wintering area, dates are given for each one (`W2_arr`, `W2_dep`, …).

The workflow combines path segmentation and clustering (`segclust2d`), trajectory metrics (`trajr`, `adehabitatLT`) and kernel density estimation (`track2KBA`), together with our own classification rules and an interactive tool for supervised correction (`shiny` + `leaflet`).

This code accompanies the manuscript:

> *Automated detection of migration phenology: a segmentation and clustering method for tracking data.* Movement Ecology (under revision).

---

## Status

This repository is currently being updated. The version used for the submitted manuscript is available at <https://doi.org/10.5281/zenodo.18088716>.

---

## Scope and intended use

- The workflow starts from **positions** (longitude, latitude, date-time). It was developed and tested with positions estimated from **GLS data of seabirds (Procellariiformes)**.
- It assumes that each track **starts and ends at a known breeding colony** (central-place breeders) and covers one annual cycle (breeding → migration → non-breeding → migration → breeding).
- It is a **framework that users are expected to adapt** to their species. Default parameter values were chosen for our study species; please read the comments next to each parameter before applying it to other systems (e.g. short-distance migrants, species with stopovers, Southern Hemisphere breeders, or GPS data).
- The reclassification of segments (step 2) is **optional and supervised**.

---

## Repository structure

```
AutoPhen_GLS/
├── 00_setup.R                          Install packages
├── 01_Segmentation_classification.qmd  Step 1: segmentation, clustering and classification
├── 02_Supervised_reclassification.R    Step 2 (optional): interactive reclassification
├── 03_Phenology_dates.qmd              Step 3: breeding/wintering areas and phenological dates
├── functions/                          Functions used in steps 1-3
├── input/                              Example positions and metadata
├── output/                             Results
├── simulations/                        Simulated tracks with known phenology (validation)
├── analyses/                           Code to reproduce the analyses, tables and figures of the paper
└── gls_preprocessing/                  Optional: example of how GLS positions were obtained from raw light data
```

---

## Workflow

| Step | Script | What it does | Main input | Main output |
|---|---|---|---|---|
| 0 | `00_setup.R` | Installs all required packages | — | R library ready |
| 1 | `01_Segmentation_classification` | Segmentation and clustering of the track (`segclust2d`), classification of segments into Breeding / Migration / Wintering, flags for segments that need review | `input/GLS_positions_example.csv`, `input/GLS_metadata_example.csv` | `output/Phenology_tracks.csv`, `output/Segmentation_results.csv` |
| 2 | `02_Supervised_reclassification` | *Optional.* Interactive map to review and correct flagged segments or positions | Outputs of step 1 | `output/*_corrected.csv` |
| 3 | `03_Phenology_dates` | Breeding and wintering areas as 90% KDE and extraction of phenological dates | Outputs of step 1 or 2 | `output/Phenology.csv` |

Functions:

| File | Used in |
|---|---|
| `functions/segmentation_classification_functions.R` | Step 1 |
| `functions/supervised_reclassification_functions.R` | Step 2 |
| `functions/phenology_dates_functions.R` | Step 3 |

---

## Installation

### Requirements

- **R ≥ 4.6** (tested with R 4.6.1).
- **RStudio** is recommended (open `AutoPhen_GLS.Rproj` so that all paths are relative to the project folder).
- On Windows, **Rtools** (matching your R version) is needed to install packages from GitHub.

### Steps

1. Download or clone the repository:

   ```bash
   git clone https://github.com/BioChess/AutoPhen_GLS.git
   ```

2. Open `AutoPhen_GLS.Rproj` in RStudio.

3. Run `00_setup.R`. This restores the exact package versions recorded in `renv.lock`.

If you find a problem with a new R version, please open an [issue](https://github.com/BioChess/AutoPhen_GLS/issues).

---

## Quick start with the example data

The folder `input/` contains positions of seven Cory's shearwater (*Calonectris borealis*) tracks (Veneguera, Gran Canaria, 2002–2003), estimated from GLS data.

1. Run `00_setup.R` (only once).
2. Run `01_Segmentation_classification`, then (optionally) `02_Supervised_reclassification`, then `03_Phenology_dates`.
3. Results are written to `output/`.

---

## Input data format

| File | Content |
|---|---|
| `input/GLS_positions_example.csv` | One row per position: track ID (`geo_trip`), date-time (`Date_Time`, UTC), longitude (`lon`) and latitude (`lat`) in decimal degrees (WGS84) |
| `input/GLS_metadata_example.csv` | One row per track: track ID, ring, species, colony name and colony coordinates (`Colony.lon`, `Colony.lat`) |

---

## Parameters

All settings are in `parameters.R`, which is loaded by scripts 01 and 03. It has two parts:

- **User parameters**: a few settings that may need to be adapted to your data or species (type of data, positions per day, equinox window, maximum number of segments, range of clusters, minimum duration of segments and stages, KDE level). Each one is explained in the file.
- **Advanced parameters**: thresholds of the classification rules. They are part of the method and are listed to make it transparent, but most users should not change them.

Phenological dates are defined as: **arrival** = first position inside the KDE of the breeding or wintering area; **departure** = first position outside the KDE after the last position inside it.

---

## GLS pre-processing (optional)

The folder `gls_preprocessing/` shows how the example positions were obtained from raw GLS light data (calibration, twilight detection and SGAT modelling). It is **not required** to run the workflow: positions can be obtained with any GLS software (e.g. `SGAT`, `probGLS`, `FLightR`) or come from other tracking devices. See `gls_preprocessing/README.md`.

---

## Citation

If you use this workflow, please cite the paper (reference will be added on publication) and the Zenodo archive of the version you used.

---

## License

Code is released under the GNU General Public License v3.0 (see `LICENSE`). Example data in `input/`, `output/`, `simulations/` and `gls_preprocessing/` are released under CC BY 4.0.

## Contact

Questions, bugs and suggestions: please open an [issue](https://github.com/BioChess/AutoPhen_GLS/issues).
