 #!/bin/bash
 source scripts/genotypes/config.sh
 

 # Check transition/transversion ratio of each file (AMP_AD and NABEC). bcftools stats reads every variant
 # in the BCF file and prints a summary report: counts of SNPs and indels, allele frequencies, quality
 # distributions and so on. Each section of the report starts with a tag; the line starting with TSTV holds
 # the transition/transversion counts. awk '$1=="TSTV"{print "ts/tv:", $5}' keeps only that line and prints
 # its 5th column, which is the ratio.
 #     - Transition: base substitition within the same chemical type.
 #     - Transversion: base substitition across different chemical type.
 # This ratio should be close to 2-2.1. We see it a bit higher, because chromosome 22 is CG-rich and has many
 # CpG sites, where transitions are especially common--> higher ratio is expected.
 # Important result: both agree.
 bcftools stats $OUT/ampad/ampad_chr22.bcf | awk '$1=="TSTV"{print "ts/tv:", $5}'
 bcftools stats $OUT/nabec/nabec_chr22.bcf | awk '$1=="TSTV"{print "ts/tv:", $5}'
 bcftools index -n $OUT/nabec/nabec_chr22.bcf; wc -l < $OUT/nabec/samples.txt; head -3 $OUT/nabec/samples.txt
 #bcftools index -n $OUT/nabec/ampad_chr22.bcf; wc -l < $OUT/nabec/samples.txt; head -3 $OUT/nabec/samples.txt