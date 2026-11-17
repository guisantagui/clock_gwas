# Project status: brain clock age-acceleration GWAS

_Last updated: 2026-11-17._ Short, living summary of where the project stands.
Dated write-ups with numbers and conclusions go in [`reports/`](reports/README.md).

## Goal
Genetic determinants of deviation from healthy brain ageing. Phenotype = out-of-sample
transcriptomic age gap (predicted − chronological, residualised on controls), globally and per
brain region; plus disease associations of (regional) acceleration.

## Where we are

| Workstream | State | Details |
|---|---|---|
| Metadata / ID harmonisation | Done for duplicates and NABEC names (Nov 12) | Mayo AMP-AD corrections: **not confirmed applied** |
| Clock cross-fitting (`scripts/clock_cv/`) | v2 run on HPC (Nov 12, all genes); results in `results/clock_cv/` | [2026-11-17 report](reports/2026-11-17_clock_cv_v2_all_genes.md) |
| Top-coded ages (90+) | Life-table imputation implemented (`dc74ca4`) | **Not applied in v2**: 03 was run without `--life_table` |
| Genotypes (`scripts/genotypes/`) | Inventory, AMP-AD liftover, GTEx/NABEC prep, sample map | `data/genotypes/sample_map.tsv` has AMP-AD, GTEx, LBP; **NABEC not in it yet** |
| Relatedness (KING) | Not started | Needed for `--related` in fold assignment |
| Joint imputation / merge | Not started | Plan: TOPMed panel, platform + ancestry PCs |
| GWAS (REGENIE) | Not started | `scripts/gwas/prepare_inptus.sh` is a stub |

## Blocking issues (fix before v3)
1. **NABEC has no out-of-fold gaps in v2.** After the NABEC ID fix (RNA-seq `UKY-1066` vs VCF `KEN-1066`),
   every NABEC individual has `has_WGS = FALSE`, so all 236 went to `always_train`. Fix: run
   `scripts/genotypes/make_sample_map.R --nabec_info ... --vcf_samples ...` to write
   `data/genotypes/genotyped_individuals.txt`, then pass `GENOTYPED=` to `launch_cv.sh`.
2. **Top-coded individuals have no gap in v2** (1920 samples). Rerun `03_collect_gaps.R` with
   `--life_table data/life_tables/us_life_tables_e89_e90.csv` (cheap, no refitting) and make sure the
   life-table file is on the HPC for the next launch.
3. **Mayo metadata corrections** from AMP-AD: confirm whether applied.

## Next steps
1. Fix 1–3 above, then launch v3 (`bash scripts/clock_cv/launch_cv.sh` with `GENOTYPED` and `LIFE_TABLE`).
2. Add NABEC to the sample map; run KING on the merged genotypes → `--related` for folds.
3. Choose spline df (3 vs 4) for the residualisation by AIC on controls.
4. Write REGENIE phenotype/covariate files (`gap_blup`, per-cohort INT, n_regions, 90+ indicator).
5. Sanity checks: old global-SVA pipeline vs per-fold SVA case–control gaps; LBP missing ages (274 samples).

## Decisions taken
- Individual-level cross-fitting: training pool = all controls; folds over genotyped individuals
  (cases and controls), grouped by individual, stratified by study × diagnosis × age bin.
- Preprocessing (gene filter, QN, SmartSVA with `mod = ~age + region`) fit on training controls per
  outer fold; held-out samples projected without their ages.
- Gene universe = all genes (not LINCS).
- Residualisation `predicted_age ~ ns(age) + substudy_batch + region + outer fold`, fit on OOF controls.
- GWAS phenotype: `gap_blup` (from `gap_resid ~ region + (1 | individual)`), inverse-normal within cohort,
  n_regions covariate; `gap_mean` as sensitivity.
- 90+: excluded from training labels; gap from sex-specific life-table expectation (90+ = 94.81 F /
  94.10 M); 90+ indicator in downstream models; sensitivity analysis excluding them.
- BrainSeq treated as no-genomics (SNP arrays not accessible). `brainSeq_pI` kept as external validation.

## Conventions
- **Git/GitHub is done manually by the user**: Claude does not commit, push, pull, merge or open PRs;
  it only edits files in the working tree and says what changed.
- `results/`, `data/`, `GWAS/` are git-ignored and contain controlled-access IDs: never commit them,
  and keep reports to aggregate numbers.
- Reports: `reports/YYYY-MM-DD_<topic>.md`, one per analysis/run, listed in `reports/README.md`.
