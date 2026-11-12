#!/bin/bash

wd=$(pwd)
script_dir=$HOME/projects/clock_gwas/scripts/genotypes
if [ "$wd" == "$script_dir" ]; then
    cd ../..
    echo "Moving to ${HOME}/projects/clock_gwas/"
fi

sbatch --array=1-21,23 scripts/genotypes/02_ampad_liftover.sh