# mc-dormancy-strategies

Simulation code for a metacommunity model of dormancy strategy, dispersal, and
stability, and the figures of the paper "Dormancy strategy, dispersal, and
environmental predictability regulate metacommunity diversity and stability
across scales". The model asks how a bank of dormant propagules, and the way its
emergence is controlled, changes diversity, stability, and synchrony across
scales along a gradient of dispersal rate, in environments that differ in
temporal predictability and spatial synchrony.

Nathan Wisnoski, Mississippi State University.

## The model

The model is a discrete-time metacommunity of organisms with a seasonal active
stage and a bank of dormant propagules, with 10 species in 100 patches placed on
a 100 x 100 torus. In each time step, for each species in each patch:

1. **Emergence.** A fraction of the propagule pool becomes active; the rest stay
   dormant. Emergence is a fixed fraction (bet-hedging dormancy), or follows a
   cue: the predictable part of the coming conditions (anticipatory) or the
   conditions of the previous step (responsive). All three share one mean
   emergence fraction; a treatment with no dormancy is the control.
2. **Reproduction.** Active individuals produce propagules at a rate set by the
   match between the patch conditions and the species' niche optimum (a
   Gaussian response), reduced by Beverton-Holt competition among the active
   individuals of every species. Propagule numbers are Poisson draws.
3. **Survival in the bank.** Dormant propagules survive with a fixed
   probability.
4. **Dispersal.** A fraction of the new propagules emigrates. Emigrants leave
   the home patch and land in other patches under an exponential kernel.
   Dormant propagules do not disperse.

**Environment.** Each patch has its own long-term mean, and conditions fluctuate
around it through a periodic layer and an autocorrelated noise layer at constant
total temporal variance. Two properties are varied independently: the fraction of
temporal variance that is periodic (predictable or unpredictable), and whether
patches fluctuate together or independently (synchronous or asynchronous).

**Matched comparison.** A cued seed bank realizes a different emergence fraction
than its nominal mean, so bet-hedging is also run at a ladder of emergence
fractions, and each cued run is compared with the bet-hedger interpolated to the
fraction it realized.

The model, its parameters, and the analyses are described in the Methods and
Supporting Information of the paper. The table below gives each parameter of the
main sweep with its name in the code.

| Symbol | Meaning | Value | Name in the code |
| --- | --- | --- | --- |
| M | patches | 100 | `patches` |
| – | landscape | 100 x 100 grid, torus | `x_dim`, `y_dim` |
| S | species | 10 | `species` |
| z_i | niche optima | evenly spaced on [0.25, 0.75] | `niche_optima` |
| σ | niche breadth | 0.1 | `niche_breadth` |
| f_max | propagules per active individual at the niche optimum | 5 | `f_max` |
| α_ii, α_ij | intra- and interspecific competition | 0.05, 0.05 x Uniform(0, 1) | `intra`, `max_inter`, `scale` |
| θ̄_x | patch means | evenly spaced on [0.25, 0.75], correlation length 10 | `patch_mean_span`, `patch_mean_range` |
| σ_θ | temporal standard deviation of conditions | 0.13 | set by `total_var` and `scale_half` |
| p | predictable fraction of temporal variance | 0.9 or 0.1 | `pred_frac` |
| ω | angular frequency of the cycle | 2π/10 | `w` |
| a | noise color | -0.6, 0, 0.6 | `noise_colors` |
| l | noise correlation length | 10 (asynchronous) or 5000 (synchronous) | `spat_range` |
| ḡ | mean emergence fraction | 0.1; 1 without dormancy; bet-hedging also 0.05–0.25 | `g_mean` |
| β | responsiveness to the cue | 1.5 (cued), 0 (bet-hedging) | `beta` |
| s | survival in the seed bank | 0.95 | `survival` |
| m | dispersal rate | 20 values, 10^-5 to 1 | `dispersal_rates` |
| k | kernel exponent | 0.1 | `kernel_decay` |
| – | assembly, burn-in, and recorded steps | 100, 300, 300 | `initialization`, `burn_in`, `timesteps` |
| – | replicates | 10 | `nreps` |

A replicate is one draw of everything random about the setting: patch
positions, competition coefficients, patch means, and the realized environment
of each regime. Every treatment is run on every replicate, so contrasts are
paired within replicate, regime, and dispersal rate. Each run's random seed is
set by its replicate, regime, and dispersal rate and is stored in the output;
the sensitivity sweep uses the same seeds, so its baseline repeats the main
sweep run for run.

## Repository layout

```text
simulations/     model functions and simulation drivers
analysis/        figures and checks
docs/            tables written by the checks, and the software versions used
figures/         figure output; figures/main/ and ms_supplement/ are the paper's
sim_output/      simulation output (not tracked)
```

## Requirements

R 4.4 or later, with the packages `here`, `dplyr`, `tidyr`, `readr`, `tibble`,
`ggplot2`, `scales`, `patchwork`, `vegan`, `foreach`, and `doParallel`; and a
LaTeX installation with `latexmk` for Fig. 1.
`analysis/session_info.R` records the versions used, in `docs/session_info.txt`.

## Reproducing the paper

Run all scripts from anywhere inside the repository; paths are resolved from the
repository root. Every setting is defined at the top of the script that uses
it, and each script runs its full analysis as written. Landscapes,
environments, and individual simulations are seeded, so each script reproduces
its output exactly.

### 1. Run the simulations

| Sweep | Command | Output |
| --- | --- | --- |
| Main sweep | `Rscript simulations/metacom_regimes.R` | `sim_output/regimes_<timestamp>/` |
| Sensitivity sweep | `Rscript simulations/metacom_sensitivity.R` | `sim_output/sensitivity_<timestamp>/` |

The main sweep takes about half an hour on 19 cores, and the sensitivity sweep
about 7 minutes. Set `test_run <- TRUE` at the top of either script for a small
version that finishes in seconds. Each sweep writes its results as CSV files
together with `run_settings.csv`; a simulation that failed keeps its row in the
results, flagged in the `ran` column.

### 2. Make the figures

Run these in order once both sweeps have finished. The first three scripts read
the most recent matching directory in `sim_output/`; to analyze a specific run,
set `sweep_dir` or `sensitivity_dir` at the top of the script.

| Step | Command | Output |
| --- | --- | --- |
| Diversity, stability, and synchrony | `Rscript analysis/seedbank_stability.R` | `figures/main/`, `figures/supplement/` |
| Stability split by route | `Rscript analysis/stability_pathways.R` | `figures/supplement/pathways_*` |
| Sensitivity to parameters held fixed | `Rscript analysis/sensitivity.R` | `figures/supplement/sensitivity.pdf` |
| Environment panel of Fig. 1 | `Rscript analysis/env_time_series.R` | `figures/main/panels/` |
| Fig. 1 | `latexmk -pdf -jobname=model -outdir=figures/main -auxdir=analysis/build analysis/model_figure.tex` | `figures/main/model.pdf` |
| Supplementary figures, numbered | `Rscript analysis/ms_supplement.R` | `figures/ms_supplement/` |

The paper's figures are those in `figures/main/` and `figures/ms_supplement/`;
`docs/ms_supplement_figures.csv` lists what each supplementary figure shows.

## Additional checks

These scripts support statements in the Methods but draw no figure in the
paper. `seedbank_checks.R` reads both sweeps and
`check_bad_years.R` the main sweep; the others need no simulation output.

| Check | Command | Output |
| --- | --- | --- |
| Lost runs, paired contrasts, and low survival in the bank | `Rscript analysis/seedbank_checks.R` | `docs/seedbank_*.csv`, `figures/supplement/sensitivity_low_survival.pdf` |
| The environmental regimes vary as intended | `Rscript analysis/test_env_scenarios.R` | `docs/env_regime_diagnostics.csv` |
| How often conditions fall below replacement | `Rscript analysis/check_bad_years.R` | `docs/env_bad_years.csv` |
| Species' niches against the conditions they meet | `Rscript analysis/niche_overlap.R` | `figures/supplement/niche_overlap.pdf` |
| R and package versions, and their citations | `Rscript analysis/session_info.R` | `docs/session_info.*`, `docs/package_citations.bib` |

## Citation

Citation metadata are in `CITATION.cff`; GitHub shows them under "Cite this
repository". Each release is archived on Zenodo with its own DOI.

## License

GPL-3.0. See `LICENSE`.
