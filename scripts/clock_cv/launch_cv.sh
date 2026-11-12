#!/bin/bash
################################################################################
# Brain clock cross-fitting: run the whole pipeline.                           #
#                                                                              #
# Local:  bash scripts/clock_cv/launch_cv.sh                                   #
# SLURM:  sbatch --array=0-15 scripts/clock_cv/launch_cv.sh                    #
#         (task id t -> rep = t / (K+1) + 1, fold = t % (K+1); fold 0 is the   #
#          final model; run 01 before and 03 after the array.)                 #
# Run from the repository root. Extra arguments are passed to 02_run_fold.R.   #
################################################################################
#SBATCH --job-name=clock_cv
#SBATCH --time=0-04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=clock_cv_%A_%a.out

set -euo pipefail

K=${K:-5}
N_REP=${N_REP:-3}
PARALLEL=${PARALLEL:-3}   # local: folds run at the same time
CORES=${CORES:-4}         # cores per fold for cv.glmnet
OUTDIR=${OUTDIR:-results/clock_cv}
GENOTYPED=${GENOTYPED:-data/genotypes/genotyped_individuals.txt}   # from make_sample_map.R
GENO_ARG=(); [[ -f $GENOTYPED ]] && GENO_ARG=(--genotyped "$GENOTYPED")
LIFE_TABLE=${LIFE_TABLE:-data/life_tables/us_life_tables_e89_e90.csv}
LT_ARG=(); [[ -f $LIFE_TABLE ]] && LT_ARG=(--life_table "$LIFE_TABLE")
SCRIPTS=scripts/clock_cv

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        t=$SLURM_ARRAY_TASK_ID
        rep=$(( t / (K + 1) + 1 ))
        fold=$(( t % (K + 1) ))
        Rscript $SCRIPTS/02_run_fold.R --fold $fold --rep $rep --cores $CORES --outdir $OUTDIR "$@"
        exit 0
fi

Rscript $SCRIPTS/01_make_folds.R --K $K --n_rep $N_REP --outdir $OUTDIR ${GENO_ARG[@]+"${GENO_ARG[@]}"}

# Build the counts cache once, so that parallel folds don't all parse the CSV.
Rscript $SCRIPTS/02_run_fold.R --fold 0 --rep 1 --cores $CORES --outdir $OUTDIR "$@"

jobs=()
for rep in $(seq 1 $N_REP); do
        for fold in $(seq 1 $K); do
                jobs+=("$rep $fold")
        done
done
printf '%s\n' "${jobs[@]}" | xargs -P $PARALLEL -L 1 bash -c \
        'Rscript '"$SCRIPTS"'/02_run_fold.R --rep $0 --fold $1 --cores '"$CORES"' --outdir '"$OUTDIR"' '"$*"' > '"$OUTDIR"'/rep$0_fold$1.out 2>&1 && echo "rep $0 fold $1 done"'

Rscript $SCRIPTS/03_collect_gaps.R --outdir $OUTDIR --folds $OUTDIR/folds.csv ${GENO_ARG[@]+"${GENO_ARG[@]}"} ${LT_ARG[@]+"${LT_ARG[@]}"}
