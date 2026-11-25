# Comparison with baseline methods

This folder compares the workflow with three simple alternatives, on the
simulated tracks (known truth) and on the *Calonectris* GLS tracks with a
supervised reference (reviewer comment R2-M7).

| Method | What it does |
| --- | --- |
| `segclust_simple_<v>` | The same segclust2d segmentation as the workflow, with a simple rule per segment: breeding if its mean NSCD is below the NSCD limit of the workflow (calculated for each track); otherwise migration if its net speed is higher than v km/day, and wintering if not. No double clustering, no agreement score, no correction rules. The speed threshold is the same for all tracks and species, so three values are tested (50, 100 and 200 km/day). |
| `migrateR` | Net squared displacement models (Bunnefeld et al. 2011; Spitz et al. 2017) on the distance to the colony. If the best model (AIC) is `migrant` or `mixmig`, the dates are the start and end of the outward and return movements (5% and 95% of the change). |
| `hmm` | Hidden Markov model with three states (momentuHMM) on step length and distance to the colony. States are named by their means (breeding = closest to the colony; migration = longest steps). |

All methods use the same positions as the workflow (projection centred on the
colony and the same cleaning) and the same final rule (short stages between
two stages of the same class are joined), so only the classification changes.
As in the workflow, each stationary non-breeding stage is a wintering area
and the movement between two wintering areas is migration.
The segclust2d results are read from the cache of the workflow runs.

## Steps

1. Run the workflow on the simulations (`simulations/02_run_workflow.R`) and
   the default variant of the sensitivity tests (`analyses/sensitivity/02_run_variants.R`).
   Their results and segclust2d caches are used here.
2. `01_run_baselines.R`: runs the three methods (settings at the top).
3. `02_evaluate_baselines.R`: compares all methods with the truth / reference
   (tables and figure in `output/baselines/_results/`).

Packages: `momentuHMM` (CRAN) and `migrateR`
(`remotes::install_github('dbspitz/migrateR/migrateR')`). The dates of
migrateR are calculated from the model coefficients, as in
`migrateR::mvmt2dt()`, because that function fails in recent R versions.
