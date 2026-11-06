#!/bin/bash

counts="/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/counts_parsed/counts_merged.csv"

Rscript normalize_dataset.R "$counts"
