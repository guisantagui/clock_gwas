#!/bin/bash
################################################################################
# Brain clock cross-fitting: run the whole pipeline.                           #
#                                                                              #
# Local:  bash scripts/clock_cv/launch_cv.sh                                   #
# SLURM:  sbatch --array=0-15 scripts/clock_cv/launch_cv.sh                    #
#         (task id t -> rep = t / (K+1) + 1, fold = t % (K+1); fold 0 is the   #
#          final model; run 01 before and 03 after the array.)                 #
# Run from the repository root. Every input can be overridden with an         #
# environment variable, e.g.  COUNTS=/path/counts.csv bash launch_cv.sh       #
# Extra arguments are passed to 02_run_fold.R.                                 #
################################################################################
#SBATCH --job-name=clock_cv
#SBATCH --time=0-04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=clock_cv_%A_%a.out

set -euo pipefail

# Inputs
COUNTS=${COUNTS:-results/counts_parsed/counts_merged.csv}                    # genes x samples
SAMP_MDAT=${SAMP_MDAT:-results/counts_parsed/sample_metadata_merged.csv}     # aligned with counts
IND_MDAT=${IND_MDAT:-results/counts_parsed/individual_metadata_merged_rnaseq_filt.csv}
LINCS=${LINCS:-../brain_clock/data/LINCS/lincs_genes.txt}                     # gene universe
GENOTYPED=${GENOTYPED:-data/genotypes/genotyped_individuals.txt}             # make_sample_map.R (optional)
LIFE_TABLE=${LIFE_TABLE:-data/life_tables/us_life_tables_e89_e90.csv}        # 90+ ages (optional)

# Settings
K=${K:-5}
N_REP=${N_REP:-3}
PARALLEL=${PARALLEL:-3}   # local: folds run at the same time
CORES=${CORES:-4}         # cores per fold for cv.glmnet
OUTDIR=${OUTDIR:-results/clock_cv}
SCRIPTS=scripts/clock_cv

for f in "$COUNTS" "$SAMP_MDAT" "$IND_MDAT" "$LINCS"; do
        [[ -f $f ]] || { echo "Missing input: $f"; exit 1; }
done
mkdir -p "$OUTDIR"
MDAT_ARGS=(--samp_mdat "$SAMP_MDAT" --ind_mdat "$IND_MDAT")
OPT_ARGS=()
if [[ -f $GENOTYPED ]]; then OPT_ARGS+=(--genotyped "$GENOTYPED"); echo "Genotyped list: $GENOTYPED"
else echo "No $GENOTYPED: using has_WGS / has_snpArray from $IND_MDAT"; fi
LT_ARGS=()
if [[ -f $LIFE_TABLE ]]; then LT_ARGS+=(--life_table "$LIFE_TABLE"); fi
FOLD_ARGS=("${MDAT_ARGS[@]}" --folds "$OUTDIR/folds.csv" --counts "$COUNTS" --lincs_genes "$LINCS" --cores "$CORES" --outdir "$OUTDIR" "$@")

run_fold() {  # rep fold
        Rscript $SCRIPTS/02_run_fold.R --rep "$1" --fold "$2" "${FOLD_ARGS[@]}"
}

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        t=$SLURM_ARRAY_TASK_ID
        run_fold $(( t / (K + 1) + 1 )) $(( t % (K + 1) ))
        exit 0
fi

Rscript $SCRIPTS/01_make_folds.R "${MDAT_ARGS[@]}" --K $K --n_rep $N_REP --outdir "$OUTDIR" \
        ${OPT_ARGS[@]+"${OPT_ARGS[@]}"}

# Final clock first: it also builds the counts cache, so parallel folds don't all parse the CSV.
run_fold 1 0

# Outer folds, PARALLEL at a time.
running=0
for rep in $(seq 1 $N_REP); do
        for fold in $(seq 1 $K); do
                ( run_fold $rep $fold > "$OUTDIR/rep${rep}_fold${fold}.out" 2>&1 \
                        && echo "rep $rep fold $fold done" || echo "rep $rep fold $fold FAILED (see $OUTDIR/rep${rep}_fold${fold}.out)" ) &
                running=$((running + 1))
                if (( running >= PARALLEL )); then wait; running=0; fi
        done
done
wait

Rscript $SCRIPTS/03_collect_gaps.R "${MDAT_ARGS[@]}" --outdir "$OUTDIR" --folds "$OUTDIR/folds.csv" \
        ${OPT_ARGS[@]+"${OPT_ARGS[@]}"} ${LT_ARGS[@]+"${LT_ARGS[@]}"}
