#!/bin/bash

wd=$(pwd)
script_dir=$HOME/projects/clock_wwas_scripts/genotypes
if [ "$wd" == "$script_dir" ]; then
    cd ../..
fi

source scripts/genotypes/config.sh

activate_env

tar -xzOf ${G}/NABEC/genomics/phg001711.v1.NABEC_WGS.sample-info.MULTI.tar.gz sample/sample-info.csv | head -5
tar -xzOf ${G}/NABEC/genomics/phg001711.v1.NABEC_WGS.sample-info.MULTI.tar.gz sample/sample-info.csv | grep -v "^#" | head -5
tar -xzOf ${G}/NABEC/genomics/phg001711.v1.NABEC_WGS.sample-info.MULTI.tar.gz sample/sample-info.csv | grep -v "^#" | grep -E "KEN-10(66|69)"
cut -d- -f1 $OUT/nabec/samples.txt | sort | uniq -c

mkdir -p $OUT/vcf_samples
bcftools query -l $OUT/ampad/ampad_chr22.bcf > $OUT/vcf_samples/ampad.txt
cp $OUT/nabec/samples.txt $OUT/vcf_samples/nabec.txt