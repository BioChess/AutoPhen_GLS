# GLS pre-processing (optional)

These scripts show how the example positions in `input/GLS_positions_example.csv` were obtained from raw light data recorded by light-level geolocators (GLS). They are provided as an example only: **the main workflow starts from positions**, which can be obtained with any GLS software or tracking device.

Run the scripts in order, from the project root (open `AutoPhen_GLS.Rproj`):

| Step | Script | What it does | Main output |
|---|---|---|---|
| 1 | `01_Supervised_GLS_calibrations` | Supervised calibration of each logger (zenith angle and twilight error parameters) | `output/Supervised_calibrations_SGAT.csv`, `output/Summary_SGAT_calibrations.csv` |
| 2 | `02_Trn_automated_calculation` | Automated twilight detection (`twilightCalc`, `GeoLight`) and cleaning | `output/trn/*.trn`, `input/GLS.files2.csv` |
| 3 | `03_SGAT_locations_from_trn` | Positions from twilights with the SGAT threshold model; optimisation of `tol` around the equinoxes; optional anchoring of positions on colony days | `output/Raw_SGAT_tracks_tol*.csv`, `input/GLS.files3.csv` |

All paths are relative to the project root (e.g. `gls_preprocessing/input/`).

## Extra packages

These scripts need four packages that are only available on GitHub (`GeoLight`, `TwGeos`, `BAStag`, `SGAT`). They are installed with `00_setup.R` when `install_gls_preprocessing <- TRUE`.

## Twilight cleaning functions

Step 2 uses `twilight_cleanup()` and `move_twilights()` from the SEATRACK automated GLS procedure:

> Bråthen VS, Moe B, Amélineau F, et al. (2021) *An automated procedure (v2.0) to obtain positions from light-level geolocators in large-scale tracking of seabirds - a method description for the SEATRACK project.* NINA Report 1893. <https://brage.nina.no/nina-xmlui/handle/11250/2735757>

These functions are not distributed in this repository. To run step 2, obtain them from the SEATRACK project and save them in `gls_preprocessing/functions/`. The output of step 2 for the example data is already included (`output/trn/`), so steps 1 and 3 can be run without them.

## Input files

| File | Content |
|---|---|
| `input/GLS.files.csv` | One row per logger and trip: ring, colony, colony coordinates, species, logger ID, producer, model, deployment and recovery dates, raw file names |
| `input/GLS.info.csv` | Logger information |
| `input/Info.calibrations.csv` | Calibration periods (dates and coordinates) for each logger |
| `input/colony_dates.csv` | Days when the bird was at the colony (from light and wet-dry data). Optional |
| `input/GLS/` | Raw light (`.lig`) and wet-dry (`.act`) files |
