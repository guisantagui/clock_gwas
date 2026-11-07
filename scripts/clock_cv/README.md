# Brain clock cross-fitting (out-of-fold age gaps for the GWAS)

Retrains the brain transcriptomic clock so that every genotyped individual,
and every case, gets an age gap from a model that never saw that individual.

Inputs (from `scripts/parse_counts/parse_counts.R` and the HPC notebook):
- `results/counts_parsed/counts_merged.csv`: genes × samples, columns = `make.names(specimenID)`
- `results/counts_parsed/sample_metadata_merged.csv`: aligned sample metadata
- `results/counts_parsed/individual_metadata_merged_rnaseq_filt.csv`: `has_WGS` / `has_snpArray`
- `../brain_clock/data/LINCS/lincs_genes.txt`: gene universe (`--gene_universe all` uses every gene)
- `region_map.csv` (here): tissue label → harmonised region (edit as needed)

## Design
| role | who | in training | gets OOF gap |
|---|---|---|---|
| `always_train` | controls without genomics (BrainSeq pII, NABEC/ROSMAP/TBI w/o genotypes) | every fold | no |
| `cv` | genotyped individuals (cases + controls) and cases without genomics | when not in the held-out fold (controls only) | yes |
| `external` | `brainSeq_pI` (as in the paper) | never | averaged over fold models |

- Folds are grouped by individual and stratified by substudy × diagnosis × age bin.
  Pass `--related` (two columns of individualIDs, e.g. KING duplicates/1st degree)
  to force related donors into the same fold.
- Training labels: controls with a numeric age ≥ 20, not top-coded, RIN ≥ 6 (NA allowed: NABEC has no RIN).
  Top-coded (`90+`, `89+`), missing-age, under-20 and low-RIN samples are still predicted.
- Cerebellum is excluded by default (as in the paper).
- Inside each outer fold, everything is fit on the training controls and applied to the held-out
  samples one at a time, without their ages:
  1. gene filter (≥ 80 % zeros; mean count < 1 in every 5-year age bin)
  2. log2(counts + 1), quantile normalisation to the training reference
  3. `SmartSVA::smartsva.cpp` (as in the paper; `--sva_method sva` for sva::sva) on training (`mod = ~ age + region`, `mod0 = ~ region`, `n.sv` from `num.sv` leek, `--n_sv rmt` for EstDimRMT + 1, or a fixed number),
     then a per-sample SV projection, the `fast` branch of `fsva` without its batch re-scaling
     (see `sva_fit()` in `cv_utils.R`). `fsva(method = "exact")` needs one SVD of the full training matrix
     per held-out sample, which is not feasible here.
  4. elastic net (`glmnet`), alpha and lambda tuned by inner CV grouped by individual; inner folds
     reuse the outer fold's preprocessing.
- `--fold 0` trains the released clock on all eligible controls; `clock.rds` holds everything
  needed to apply it (genes, QN target, SVA projection, glmnet fit).

## Run
```bash
bash scripts/clock_cv/launch_cv.sh            # K=5, N_REP=3, 3 folds in parallel
```
or step by step:
```bash
Rscript scripts/clock_cv/01_make_folds.R --K 5 --n_rep 3
Rscript scripts/clock_cv/02_run_fold.R --rep 1 --fold 1 --cores 4   # repeat for every rep/fold, plus --fold 0
bash scripts/life_tables/get_us_life_tables.sh          # US life tables (CDC/NCHS) -> data/life_tables/
Rscript scripts/clock_cv/03_collect_gaps.R --life_table data/life_tables/us_life_tables_e89_e90.csv
```
For a quick check, add `--debug_n_genes 600 --debug_n_train 300 --n_sv 3 --outdir results/clock_cv_debug` to 02.

## Outputs (`results/clock_cv/`)
- `folds.csv`, `folds_summary.txt`: fold assignment, counts by role/fold, 90+ among genotyped, genotyped N per region
- `rep<r>/fold<k>/`: `predictions.csv`, `coefficients.csv`, `tuning.csv`, `clock.rds`, `log.txt`
- `sample_gaps.csv`: OOF predicted age (mean over repeats), `gap_raw`, `gap_resid`
  (`predicted_age ~ ns(age, 4) + substudy_batch + region + outer fold`, fit on OOF controls per repeat, then averaged over repeats), `gap_resid_ageonly`,
  `gap_ind_mean`, `gap_regional_dev` (sample gap − individual mean; individuals with ≥ 2 regions)
- `individual_gaps.csv`: mean gap, and BLUP from `gap_resid ~ region + (1 | individual)` if `lme4` is installed
- `oof_performance_controls.csv`, `gap_by_substudy_diagnosis.csv`

Top-coded ages ("90+", "89+") get floor + e(floor) from the US life tables by sex, averaged over
2006–2019 (`--life_table_years`; 2020–2022 are lowered by COVID-19): 90+ = 94.81 (F) / 94.10 (M),
89+ = 94.18 / 93.41. Without `--life_table` (or `--topcode_age_female/male`) they get no gap. `age_topcoded` flags them
for the 90+ covariate and the sensitivity analysis.
