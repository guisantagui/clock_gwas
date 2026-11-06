# Postmortem human brain expression data curation and parsing
This directory contains the code used to curate and parse the postmortem human brain expression datasets used in our study into unified metadata and counts matrices. This step was semi-manual and involved harmonizing clinical metadata, identifying control samples and unifying categories and units across studies.

## Integrated studies
We integrated the expressiond ata from the following studies:
- [**The RNAseq Harmonization Study (rnaSeqReprocessing)**](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn9702085)
- [**The Living Brain Project (LBP)**](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn26337520)
- [**The Aging, Dementia and Traumatic Brain Injury Study(TBI)**](https://aging.brain-map.org)
- [**The Genotype-Tissue Expression (GTEx) project (GTEx)**](https://www.gtexportal.org/home/downloads/adultgtex/bulk_tissue_expression)
- [**BrainSeq Phase 1**](http://eqtl.brainseq.org/phase1/)
- [**BrainSeq Phase 2**](https://eqtl.brainseq.org/phase2/)
- Pseudobulk data of brain samples from [**ageAnno**](https://relab.xidian.edu.cn/AgeAnno/#/)

## Running the pipeline
To execute the full pipeline, run:
```bash
sbatch dataset_parsing.sh
```

## Parsing workflow
This curation/parsing pipeline can be divided in the following steps:

### 1. Parsing of AMP-AD datasets
#### 1.1. Parsing of AMP-AD individual metadata
- Explore diagnosis variables to determine which individuals are healthy, which suffered from neurodegeneration and which needed to be removed due to having ambiguous or non-existent age information or diagnosis out of the scope of our study.
- Homogenized information such as ethnicity and sex accross the different metadata datasets.

**Script**: `amp_ad_metadata_parsing.R`.

**Required files**, whose access need to be requested in synapse.org and can be accessed from the following locations:
- From [The Living Brain Project (LBP) metadata directory](https://www.synapse.org/Synapse:syn27127033):
    - `LBP_individual_metadata.csv` (syn58589260).
    - Put in `./data/metadata/LBP`.
- From [The RNAseq Harmonization Study (RNAseq Harmonization) metadata directory](https://www.synapse.org/Synapse:syn27000096):
    - `RNAseq_Harmonization_Mayo_combined_metadata.csv` (syn27000373)
    - `RNAseq_Harmonization_MSBB_combined_metadata.csv` (syn27000243)
    - `RNAseq_Harmonization_ROSMAP_combined_metadata.csv` (syn27034471).
    - Put in `./data/metadata/RNAseq_Harmonization`.
- From [The MayoRNAseq metadata directory](https://www.synapse.org/Synapse:syn23634010):
    - `MayoRNAseq_individual_metadata.csv` (syn23277389)
    - Put in `./data/metadata/rnaSeqReprocessing`.
    - The reason for including this file is that there are several samples that in the RNAseq Harminozation version have the Thal and Braak scores missing but in this version are not, so we included this metadata file to complete them.

**Output**: unified CSV metadata files for each dataset, and a metadata file comining all the selected samples from the RNAseq Harmonization study (`RNAseq_Harmonization_ind_all.csv`), all written to `./results/metadata_parsed/`.

#### 1.2. Homogenizing RNAseq harmonization metadata and counts matrices
- Merge count matrices of the RNAseq Harmonization study.
- Add clinical information of ROSMAP study is added to the RNAseq Harmonization combined metadata file.

**Script**:  `rnaHarm_parseCounts.R`.

**Required files**: generated in the previous step, and other files obtained from synapse.org. The whole set of files needed by this script and their required location  is the following:
- The metadata files generated in the previous step (saved in `results/metadata_parsed/`).
- `ROSMAP_clinical.csv` (syn3191087), obtained from the [ROSMAP study metadata subdirectory](https://www.synapse.org/Synapse:syn3157322).
- Count matrices:
    - `Mayo_gene_all_counts_matrix_clean.txt` (syn21544635), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Mayo subdirectory](https://www.synapse.org/Synapse:syn20825471). This file needs to be placed in `./data/expression/RNAseqHarm/`.
    - `MSBB_gene_all_counts_matrix_clean.txt` (syn21544666), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) MSSM subdirectory](https://www.synapse.org/Synapse:syn20957610). This file needs to be placed in `./data/expression/RNAseqHarm/`.
    - `ROSMAP_batch1_gene_all_counts_matrix_clean.txt` (syn22283382), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 1 subdirectory](https://www.synapse.org/Synapse:syn22279877). This file needs to be placed in `./data/expression/RNAseqHarm/`.
    - `ROSMAP_batch2_gene_all_counts_matrix_clean.txt` (syn22301601), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 2 subdirectory](https://www.synapse.org/Synapse:syn22296752). This file needs to be placed in `./data/expression/RNAseqHarm/`.
    - `ROSMAP_batch3_gene_all_counts_matrix_clean.txt` (syn22314230), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 3 subdirectory](https://www.synapse.org/Synapse:syn22300974). This file needs to be placed in `./data/expression/RNAseqHarm/`.
    - `ROSMAP_batch4_gene_all_counts_matrix_clean.txt` (syn25817661), obtained from [The RNAseq_Harmonization Gene Expression (Raw Gene Counts) Rosmap batch 4 subdirectory](https://www.synapse.org/Synapse:syn25810549). This file needs to be placed in `./data/expression/RNAseqHarm/`.

**Output**:
- A counts matrix of all the RNAseq_Harmonization study datasets: `RNAseqHarm_allCounts.csv`. Saved in `results/parsed/`.
- A metadata file including clinical variables: `RNAseq_Harmonization_ind_all_ROSMAPBtch.csv`. Saved in `results/metadata_parsed`.

#### 1.3. Homogenizing LBP metadata.
- Combine individual, biospecimen and assay metadata of the LBP study into a single file, following the same structure as the metadata of the RNAseq harmonization study.

**Script**: `parse_LBP_metdat.R`.

**Required files**:
- `LBP_ind_all.csv`, which should have been generated by `amp_ad_metadata_parsing.R` in step 1.1, and placed into `./results/metadata_parsed/`.
- `LBP_individual_metadata.csv` (syn58589260), `LBP_assay_RNAseq_metadata.csv` (syn64555267) and `LBP_biospecimen_metadata.csv` (syn64504156), obtained from [LBP's Metadata subdirectory](https://www.synapse.org/Synapse:syn27127033). These files need to be placed into `./data/metadata/LBP/`.

**Output**: a metadata file of the LBP samples following the same structure of `RNAseq_Harmonization_ind_all.csv`: `LBP_metadata_unified.csv`. Saved in `results/parsed/`.

### 2. Parsing of GTEx datasets.
- Merge the GTEx brain tissue counts files into a single file.
- Explore metadata files to select the samples that fit the requirements of our study:
    - having an unambiguous age
    - not having suffered psychiatric or drug abuse disorders
- Divide the samples in healthy controls or affected by neurodegeneration
- Harmonize metadata the file structure to based on `RNAseq_Harmonization_ind_all.csv`.

**Script**: `GTEx_parsing.R`.
**Required files**:
- GTEx brain bulk counts files. Open access, obtained from [The GTEx portal bulk tissue expression section](https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression). These files need to be placed in `./data/gtex/counts/`. In particular, the files we used were:
    - `gene_reads_2017-06-05_v8_brain_amygdala.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_anterior_cingulate_cortex_ba24.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_caudate_basal_ganglia.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_cerebellar_hemisphere.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_cerebellum.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_cortex.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_frontal_cortex_ba9.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_hippocampus.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_hypothalamus.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_nucleus_accumbens_basal_ganglia.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_putamen_basal_ganglia.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_spinal_cord_cervical_c-1.gct.gz`
    - `gene_reads_2017-06-05_v8_brain_substantia_nigra.gct.gz`
- The GTEx metadata files:
    - `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`
    - `GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt`
    - available upon request at dbGaP (protected data access information can be found [here](https://www.gtexportal.org/home/protectedDataAccess)).
    - These files need to be placed in `data/gtex/`.

**Output**:
- A counts matrix integrating all the selected GTEx datasets: `GTEx_allCounts.csv`.
- The corresponding metadata file: `GTEx_metadata_unified.csv`.
- Both written at `results/parsed/`.

### 3. Merge GTEx, LBP and RNAseq Harmonization datasets
- Merge parsed counts and metadata matrices coming from GTEx, LBP and RNAseq Harmonization generated in the previous steps.

**Script**: `merge_GTEx_LBP_rnaSeqHarm.R`

**Required files**:
- `LBP_FlagshipPaper_featureCounts.csv` (syn64567831). Obtained from the [LBP Gene Expression (RNAseq) counts subdirecotry](https://www.synapse.org/Synapse:syn52132890). Must be placed in `data/expression/LBP/`.
- `LBP_metadata_unified.csv`. Should have been previously generated by `rnaHarm_parseCounts.R` and placed in `results/parsed/` in step 1.2.
- RNAseq Harmonization counts matrix file `RNAseqHarm_allCounts.csv`. Should have been previously generated by `rnaHarm_parseCounts.R`, and placed in `results/parsed/` in step 1.2.
- RNAseq Harmonization metadata file `RNAseq_Harmonization_ind_all_ROSMAPBtch.csv`, Should have been created by `rnaHarm_parseCounts.R` and placed in `results/metadata_parsed/` in step 1.2.
- GTEx counts matrix file `GTEx_allCounts.csv` and its corresponding metadata file `GTEx_metadata_unified.csv`, which should have been previously generated by `GTEx_parsing.R` and placed in `results/parsed/` in step 2.

**Output**:
- Combined counts matrix: `combined_counts.csv`
- Combined metadata matrix: `combined_metDat.csv`.
- Both saved in `results/parsed/`.

### 4. Obtention of TBI data and merge with the combined dataset
- Download data from [the Aging, Dementia and TBI study](https://aging.brain-map.org).
- Merge files into a single counts matrix (keeping only samples that didn't suffered a TBI, and have non-ambiguous ages).
- Merge TBI counts matrix with the `combined_counts.csv` matrix previously generated in step 3.
- Merge TBI metadata with `combined_metDat.csv` matrix, previously generated in step 3.

**Script**: `tbi_get.R`

**Required files**:
- `tbi_data_files.csv`. Freely available [here](https://aging.brain-map.org/data/tbi_data_files.csv). Must be placed in `data/dataset_parsing/`.
- `DonorInformation.csv`. Freely available [here](https://aging.brain-map.org/api/v2/data/query.csv?criteria=model::ApiTbiDonorDetail,rma::options[num_rows$eqall]). Must be placed in `data/dataset_parsing/`.
- `combined_counts.csv` and `combined_metDat.csv`. Should have been previously generated by `merge_GTEx_LBP_rnaSeqHarm.R` in `results/parsed/`.

**Output**:
- `TBI_metadata_unified.csv`
- `combined_metDat_wTBI.csv`
- `combined_counts_wTBI.csv`
- All saved into `results/parsed/`.

### 5. Obtention of BrainSeq Phase 2 data
- Parse and harmonize data from [BrainSeq Phase 2 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7000204/#S10) into counts and metadata matrices.
- Filter samples to keep only control individuals without ambiguous ages and of 18 years old or older.

**Script**: `brainSeq_pII_parsing.R`

**Required files**: `rse_gene_unfiltered.Rdata`. Freely available [here](https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_gene_unfiltered.Rdata). Must be placed in `data/brainSeq_pII/`.

**Output**: 
- `brainseq_pII_counts.csv`
- `brainseq_pII_metadata.csv`
- Both written into `results/parsed/`.

### 6. Obtention of BrainSeq Phase 1 data
- Parse and harmonize data from [BrainSeq Phase 1 study](https://doi.org/10.1038/s41593-018-0197-y) into counts and metadata matrices. This dataset has much less samples than the Phase 2. Because of this, we used entirely this dataset as independent validation dataset. In [BrainSeq Phase 1](http://eqtl.brainseq.org/phase1/) there are some samples that come from individuals also present in BrainSeq Phase 2. We removed those samples, and retained only the controls of 20 years old or older, finally keeping 22 samples.

**Script**: `brainSeq_pI_parsing.R`

**Required files**: `rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda`. Freely available [here](https://s3.us-east-2.amazonaws.com/jaffe-nat-neuro-2018/rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda). Must be placed in `./data/brainSeq_pI/`.

**Output**:
- `brainseq_pI_counts.csv`
- `brainseq_pI_metadata.csv`
- Both written into `results/parsed/`.

### 7. Pre-processing, integration and annotation of AgeAnno brain scRNA-seq data
- Pre-process brain scRNA-seq Seurat object from [AgeAnno](https://relab.xidian.edu.cn/AgeAnno/#/) was pre-processed to remove low quality samples and cells.
- Integrate samples.
- Annotate cell types.
- Generate a matrix of pseudo-bulk counts per cell type within each sample.

**Script**: `ageAnno_preproc.R`

**Required files**:
- `brain.rds`: Seurat Object of AgeAnno's brain samples. Available [here](https://drive.google.com/drive/folders/160i9KmFJ0tEYP2QBT5IjbJqhqvuu5rBW). Must be placed in `./data/ageAnno/`.
- `age_info.csv`: age of death of each one of the donors. Already included in `data/ageAnno`.
- `CellCycleGenes_Human.csv`: cell cycle genes. Obtained from [HBC Training](https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html). This file is already included in `./data/utility_files/`.
- `scRNAmarker.txt`: cell markers used to annotate each one of the cell types in the brain in [AgeAnno's original publication](https://doi.org/10.1093/nar/gkac847). Can be obtained from [AgeAnno's GitHub page](https://github.com/vikkihuangkexin/AgeAnno?tab=readme-ov-file). This file is already included in `./data/utility_files/`.

**Output** files will be written in `results/parsed/ageAnno/`, and are:
- Multiple plots of the pre-processing steps.
- Seurat objects for each sample and for the integrated objects
- Plots of the integration process
- Files used in downstream analysis:
    - `ageanno_brain_pb_counts.csv`: the pseudo-bulk counts of each cell type per sample.
    - `ageanno_brain_pb_metdat.csv`: the parsed metadata matching pseudo-bulk rows.

### 8. Integrate all clinical samples into a single file and pre-process
#### 8.1 Add brainSeq Phase II, AgeAnno and perturbation compilation to merged dataset
- Merge counts and metadata files generated in steps 4, 5, 6 and 7 with `combined_metDat_wTBI.csv` and `combined_counts_wTBI.csv`, which were generated in step 4.

**Script**: `merge_clinical.R`.

**Required files**: If previous steps ran correctly, the required files should be already placed in the corresponding directories.

**Output**:
- `merged_counts.csv`
- `merged_metdat.csv`
- Both saved at `results/parsed/merged/`.
#### 8.2 Filter clinical dataset
- Filter merged counts and metadata datasets are by:
    - Keeping only genes present in LINCS L1000 level 3
    - Removing genes with zero or near zero expression
    - Removing samples with RIN < 6

**Script**: `filt_dataset.R`.

**Output**: Unless specified, `merged_counts.csv` and `merged_metdat.csv` will be overwritten.