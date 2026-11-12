#!/bin/bash

wd=$(pwd)
script_dir=$HOME/projects/clock_wwas_scripts/genotypes
if [ "$wd" == "$script_dir" ]; then
    cd ../..
fi

source scripts/genotypes/config.sh

tar -xzOf ${G}/NABEC/genomics/phg001711.v1.NABEC_WGS.sample-info.MULTI.tar.gz sample/sample-info.csv | head -5
cut -d- -f1 $OUT/nabec/samples.txt | sort | uniq -c