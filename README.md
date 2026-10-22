# AutoPhen_GLS

**A semi-automated R workflow to extract migration phenology dates from light-level geolocator (GLS) tracking data.**

AutoPhen_GLS takes raw light data recorded by GLS loggers and returns, for each annual track, four key phenological dates:

| Date | Code | Definition |
|---|---|---|
| Departure from the breeding area | `B1_dep` | Exit from the 90% KDE of the pre-migration breeding positions |
| Arrival at the wintering area | `W1_arr` | First position inside the 90% KDE of the wintering positions |
| Departure from the wintering area | `W1_dep` | Exit from the 90% KDE of the wintering positions |
| Arrival at the breeding area | `B2_arr` | First position inside the 90% KDE of the post-migration breeding positions |

If an individual uses more than one wintering area, dates are given for each one (`W2_arr`, `W2_dep`, …).

> Note: the exact definition of departure dates is being revised (see *Status* below).

The workflow combines existing R packages for GLS processing (`GeoLight`, `TwGeos`, `SGAT`), path segmentation and clustering (`segclust2d`), trajectory metrics (`trajr`, `adehabitatLT`) and kernel density estimation (`track2KBA`), together with our own classification rules and an interactive tool for supervised correction (`shiny` + `leaflet`).

This code accompanies the manuscript:

> *Automated detection of migration phenology: a segmentation and clustering method for tracking data.* Movement Ecology (under revision).

---

## Status

This repository is currently being updated. The version used for the submitted manuscript is available at <https://doi.org/10.5281/zenodo.18088716>.

---

## Scope and intended use

- The workflow is designed for **GLS data** and has been tested on **seabirds (Procellariiformes)** breeding in the North Atlantic and Mediterranean.
- It assumes that each track **starts and ends at a known breeding colony** (central-place breeders) and covers one annual cycle (breeding → migration → non-breeding → migration → breeding).
- It is a **framework that users are expected to adapt** to their species. Default parameter values were chosen for our study species; please read the comments next to each parameter before applying it to other systems (e.g. short-distance migrants, species with stopovers, Southern Hemisphere breeders, or GPS data).
- Two steps are **supervised**: the calibration of the loggers (script 01) and the optional reclassification of segments (script 05).

---

## Workflow

Scripts are numbered and must be run in order. Each `.R` script has a `.qmd` (Quarto) version with the same code plus explanations; you can use either.

| Step | Script | What it does | Main input | Main output |
|---|---|---|---|---|
| 0 | `00_setup.R` | Installs all required packages (exact versions via `renv`) | `renv.lock` | R library ready |
| 1 | `01_Supervised_GLS_calibrations` | Supervised calibration of each logger (zenith angle and twilight error parameters) | `input/GLS/` raw light files, `input/GLS.files.csv`, `input/GLS.info.csv`, `input/Info.calibrations.csv` | `output/Supervised_calibrations_SGAT.csv`, `output/Summary_SGAT_calibrations.csv` |
| 2 | `02_Trn_automated_calculation` | Automated twilight detection and cleaning (`twilightCalc`, `twilight_cleanup`, `move_twilights`) | Raw light files + calibrations | `output/trn/*.trn`, `input/GLS.files2.csv` |
| 3 | `03_SGAT_locations_from_trn` | Positions from twilights with the SGAT threshold model; optimisation of `tol` around the equinoxes; optional anchoring of colony days | `.trn` files, `input/colony_dates.csv` | `output/Raw_SGAT_tracks_tol*.csv`, `input/GLS.files3.csv` |
| 4 | `04_Automated_migratory_classification` | Segmentation and clustering (`segclust2d`), segment classification into Breeding / Migration / Wintering, equinox flags | SGAT tracks | `output/Phenology_tracks.csv`, `output/Segmentation_results.csv` |
| 5 | `05_Supervised_reclassification` | *Optional.* Interactive map to review and correct flagged segments or positions | Outputs of step 4 | `output/*_corrected.csv` |
| 6 | `06_Phenology_dates` | Breeding and wintering areas as 90% KDE and extraction of phenological dates | Outputs of step 4 or 5 | `output/Phenology.csv` |

Custom functions are in `functions/`:

| File | Used in |
|---|---|
| `01_FUNCTIONS_GLS_analysis.R` | Steps 1–3 |
| `FUNCTION_twilight_cleanup.R`, `FUNCTION_move_twilights.R` | Step 2 (adapted from SEATRACK, see *Credits*) |
| `02_FUNCTIONS_automated_phenology.R` | Step 4 |
| `03_FUNCTIONS_Supervised_reclassification.R` | Step 5 |
| `04_FUNCTIONS_Phenology_dates.R` | Step 6 |

The simulation study (validation with tracks of known phenology) is in `Simulations/`.

---

## Installation

### Requirements

- **R ≥ 4.6** (tested with R 4.6.1). Earlier versions (≥ 4.2.3) may work but are not tested.
- **RStudio** is recommended (open `AutoPhen_GLS.Rproj` so that all paths are relative to the project folder).
- On Windows, **Rtools** (matching your R version) is needed to install packages from GitHub.

### Steps

1. Download or clone the repository:

   ```bash
   git clone https://github.com/BioChess/AutoPhen_GLS.git
   ```

2. Open `AutoPhen_GLS.Rproj` in RStudio.

3. Run `00_setup.R`. This restores the exact package versions recorded in `renv.lock`, including the GitHub-only packages `SGAT`, `BAStag`, `TwGeos` and `GeoLight`.

Some packages used here are no longer under active development (e.g. `SGAT`, `BAStag`). To keep the workflow running in the future, their versions are pinned to specific GitHub commits in `renv.lock`. If you find a problem with a new R version, please open an [issue](https://github.com/BioChess/AutoPhen_GLS/issues).

---

## Quick start with the example data

The repository includes raw data from seven Cory's shearwater (*Calonectris borealis*) GLS tracks (Veneguera, Gran Canaria, 2002–2003) in `input/GLS/`, so the whole workflow can be run end to end:

1. Run `00_setup.R` (only once).
2. Run scripts `01` to `06` in order. Step 1 is interactive (calibration periods are selected by the user); the files produced for the example data are already in `output/`, so you can start from step 2 if you prefer.
3. Results are written to `output/`.

---

## Input data format

| File | Content |
|---|---|
| `input/GLS.files.csv` | One row per logger and trip: ring, colony, colony coordinates, species, logger ID, producer, model, deployment and recovery dates, raw file names |
| `input/GLS.info.csv` | Logger information |
| `input/Info.calibrations.csv` | Calibration periods (dates and coordinates) for each logger |
| `input/colony_dates.csv` | Days when the bird was at the colony (from light and wet-dry data). Optional |
| `input/GLS/` | Raw light (`.lig`, `.lux`) and wet-dry (`.act`, `.deg`) files |

Metadata files are updated along the workflow (`GLS.files.csv` → `GLS.files2.csv` → `GLS.files3.csv`).

---

## Parameters

Default values and the reasoning behind them will be listed here (and in a single parameter file) in the revised version.

---

## Credits

- `FUNCTION_twilight_cleanup.R` and `FUNCTION_move_twilights.R` are adapted from the SEATRACK automated GLS procedure: Bråthen VS, Moe B, Amélineau F, et al. (2021) *An automated procedure (v2.0) to obtain positions from light-level geolocators in large-scale tracking of seabirds.* NINA Report 1893. <https://brage.nina.no/nina-xmlui/handle/11250/2735757>
- GLS processing relies on `GeoLight` (Lisovski & Hahn 2012), `TwGeos` (Lisovski et al.), `BAStag` and `SGAT` (Wotherspoon et al.).

---

## Citation

If you use this workflow, please cite the paper (reference will be added on publication) and the Zenodo archive of the version you used.

---

## License

Code is released under the GNU General Public License v3.0 (see `LICENSE`). Example data in `input/` and `output/` are released under CC BY 4.0.

## Contact

Questions, bugs and suggestions: please open an [issue](https://github.com/BioChess/AutoPhen_GLS/issues).
