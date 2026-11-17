# Clock cross-fitting v2: all genes, NABEC IDs fixed

- **Run:** HPC, 2026-11-12, `scripts/clock_cv/launch_cv.sh` (K = 5, 3 repeats + final model).
- **Code:** commit `e24a75b` (`--gene_universe all`).
- **Outputs:** `results/clock_cv/` (v2); compared against `results/clock_cv_v1/` (v1: Nov 6–7, LINCS genes,
  old metadata).
- **Settings:** cerebellum excluded; brainSeq_pI external; train controls with age ≥ 20, RIN ≥ 6, not top-coded;
  SmartSVA (`~ age + region` vs `~ region`), n.sv by leek (3 SVs per fold); elastic net, alpha grid
  0.05–1, 10 inner folds, `lambda.min`.

## Design as run

| Role | Individuals | Samples |
|---|---|---|
| `always_train` (controls w/o genomics) | 507 (BrainSeq pII 258, **NABEC 236**, TBI 8, ROSMAP 5) | 925 |
| `cv` (get OOF gaps) | 2598 (GTEx 391, LBP 410, Mayo 276, MSBB 315, ROSMAP 1163, TBI 43) | 7696 |
| `external` (brainSeq_pI) | 22 | 22 |

Per fold: ~2700 training samples / 870 individuals; ~26,000 of 34,018 genes pass the expression filter;
alpha = 0.05 selected (≈ 1800 genes).

## Performance (out-of-fold, controls)

| Set | n samples | MAE (y) | r | bias |
|---|---|---|---|---|
| All | 2866 | 4.51 | 0.88 | −0.02 |
| GTEx | 2453 | 4.51 | 0.82 | −0.03 |
| LBP | 89 | 4.44 | 0.88 | 0.05 |
| MSBB | 74 | 4.73 | 0.69 | −0.16 |
| Mayo | 28 | 5.81 | 0.56 | 0.12 |
| ROSMAP | 200 | 4.22 | 0.48 | −0.16 |
| brainSeq_pI (external) | 22 | 4.98 | 0.96 | 3.22 |

- Per fold: MAE 4.05–5.31, r 0.84–0.90. Coefficient correlations between fold models are 0.79–0.92.
- Weakest regions: posterior cingulate (r 0.35, n 55), temporal cortex (0.59), inferior frontal gyrus (0.65).
- Low r in ROSMAP, Mayo and MSBB reflects their narrow control age range (mostly 70–90), not a large error:
  the MAEs are similar to GTEx.

## Case vs control gaps (`gap_resid`, controls centred at 0)

| Substudy | Cases n samples | Mean gap (y) | SD |
|---|---|---|---|
| LBP | 133 | +5.51 | 5.44 |
| Mayo | 200 | +2.68 | 4.35 |
| MSBB | 721 | +1.39 | 5.41 |
| ROSMAP | 1203 | +1.38 | 4.46 |
| GTEx | 238 | +0.95 | 3.91 |

Cases have older transcriptomic age in every cohort, consistent with the published clock.

## v2 vs v1
- Overall OOF MAE 4.51 vs 4.63, but not like-for-like: v1 included NABEC controls (MAE 6.33). Restricted to
  GTEx, MAE is essentially unchanged (4.51 vs 4.46). Using all genes instead of LINCS does **not** improve accuracy.
- Individual gaps agree moderately between versions (n = 2620 shared: r = 0.85 for `gap_mean` and `gap_blup`).
  This is the expected sensitivity to the gene set and folds, and is worth keeping in mind for GWAS power.
- Case gaps keep the same direction and similar size (LBP 4.8 → 5.5, Mayo 3.1 → 2.7, ROSMAP 0.8 → 1.4).

## Problems found
1. **NABEC got no OOF gaps.** After the ID fix, all NABEC individuals have `has_WGS = FALSE` in
   `individual_metadata_merged_rnaseq_filt.csv` (RNA-seq IDs `UKY-…` vs VCF IDs `KEN-…`), so all 236 were
   put in `always_train`. In v1, 181 NABEC individuals were `cv`. Their current gaps would be in-sample and
   cannot be used. Fix: generate `data/genotypes/genotyped_individuals.txt` with
   `make_sample_map.R --nabec_info`, and launch with `GENOTYPED=`.
2. **Top-coded ages have no gap.** All 1920 top-coded samples have `gap_resid = NA`: the collect step ran
   without `--life_table`. Top-coded genotyped individuals: ROSMAP 547, MSBB 111, Mayo 46, LBP 19 (723 in total).
3. **Missing ages:** 274 LBP samples and 2 Mayo samples have no age, so they have no gap. 161 TBI samples
   are under the training age range.

Genotyped individuals with a usable gap in v2: GTEx 391, LBP 222, Mayo 228, MSBB 185, ROSMAP 541
(**1567**), rising to about 2100 with 90+ imputation and about 2400 with NABEC.

## Conclusions
- The cross-fitting clock performs about as well as the published one out of fold (MAE ≈ 4.5 y), and
  reproduces the case > control gap. The pipeline is ready for production runs.
- v2 cannot be used for the GWAS as it stands. A v3 run needs NABEC's genotyped status from the VCFs and the
  life-table imputation. Problem 2 alone only needs `03_collect_gaps.R` rerun.
- The gene universe choice (LINCS vs all) barely changes accuracy. "All genes" is kept as decided.
