# Sensitivity tests

Each variant changes **one** parameter of the workflow (or one pair of related
thresholds). All other parameters keep the values of `parameters.R`. The
results are compared with the reference: the classification and dates of the
Calonectris tracks after the supervised review (and, optionally, the dates of a
second analyst).

The data and the results are not in the repository (`data/` and
`output/sensitivity/` are in `.gitignore`).

## Files

| File | What it does |
|---|---|
| `sensitivity_variants.R` | List of variants: parameter, value, value in the main analysis, and step that changes |
| `sensitivity_functions.R` | Functions to run a variant and to compare it with the reference |
| `01_prepare_cdays_positions.R` | Positions with other treatments of the colony days (only for the `cdays_` variants) |
| `02_run_variants.R` | Runs the variants chosen in `to_run` (one after the other) |
| `03_evaluate_variants.R` | Tables: results per track, per variant and per date |
| `04_sensitivity_figures.R` | Figures: one panel per parameter, and a summary of all tests |

## How to use

All scripts are run from the project root (open `AutoPhen_GLS.Rproj`).

1. Put the reference files in `data/reference/`:
   - `CAL_tracks_supervised.csv`: class of each position after the review (`ID`, `date`, `Breed2`)
   - `CAL_phen_supervised.csv`: dates after the review (`ID`, `B1_dep`, `W1_arr`, `W1_dep`, ..., `B2_arr`)
   - `Phen_other_analysts.csv` (optional): dates estimated visually by other analysts, before and without this workflow (`geo_trip`, `Analyst`, `B1_dep`, `W1_arr`, ...)
2. Only for the colony-day variants: run `01_prepare_cdays_positions.R`.
3. In `02_run_variants.R`, run `to_run <- 'default'` first. Then run the other
   variants, one or several at a time (e.g. `to_run <- c('kmax_25', 'kmax_35')`).
   Variants that are already done are skipped.
4. Run `03_evaluate_variants.R` and `04_sensitivity_figures.R`. They use all
   the variants that are done, so they can be run at any moment.

## Notes

- Variants with step `'01'` run scripts 01 and 03 again (about 40 s per track
  without parallel sessions). Variants with step `'03'` only change the dates
  and use the classification of `default` (a few minutes).
- The result of `segclust()` of each track is saved in
  `output/sensitivity/_segclust_cache/`. Variants that only change the
  classification thresholds use it and do not repeat the segmentation. For
  these variants, the run time is not reported.
- `segclust2d` does not accept segments shorter than 5 positions (2.5 days with
  2 positions per day), so only longer minimum durations are tested.
- With `n_workers > 1` several tracks are processed at the same time. The
  results are the same as with one session.
