library(dplyr)
################################################################################
# Brain clock GWAS: merge all counts files, and create a global metadata       #
# file that corresponds to it row by row.                                      #
################################################################################

ind_mdat_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/individual_metadata_combined.csv"
samp_mdat_file <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/specimen_metadata_combined.csv"

gtex_counts_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/GTEx"
rnaseq_harm_counts_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/RNASeq_Harm"
lbp_counts_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/LBP"

metDat1_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/GTEx/phs000424.v11.pht002743.v10.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz"

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/counts_parsed/"

plotUtils::create_dir_if_not(outDir)

# Functions
################################################################################

createGTExMat <- function(fileDir){
        brainFiles <- list.files(fileDir, full.names = T)
        brainFiles <- brainFiles[!grepl("gtf|txt", brainFiles)]
        brainFiles <- brainFiles[grepl("gct", brainFiles)]
        
        tissues <- gsub("gene_tpm_2017-06-05_v8_",
                        "",
                        gsub(".gct.gz",
                             "",
                             basename(brainFiles)))
        
        
        brainDF <- data.frame()
        sampInfo <- data.frame(matrix(nrow = 0, ncol = 2 + ncol(metDat1) - 1,
                                      dimnames = list(NULL,
                                                      c("sample",
                                                        "tissue",
                                                        colnames(metDat1)[2:ncol(metDat1)]))))
        pb <- txtProgressBar(min = 0,
                             max = length(brainFiles),
                             style = 3)
        
        for(i in seq_along(brainFiles)){
                setTxtProgressBar(pb, i)
                f <- brainFiles[i]
                tiss <- tissues[i]
                tpms <- read.delim(file=f, skip=2)
                if(ncol(brainDF) == 0){
                        brainDF <- tpms
                }else{
                        tpms <- tpms[match(brainDF$Name, tpms$Name), ]
                        tpms <- tpms[, !colnames(tpms) %in% c("id",
                                                              "Name",
                                                              "Description")]
                        brainDF <- cbind.data.frame(brainDF, tpms)
                }
                samps <- colnames(tpms)[!colnames(tpms) %in% c("id",
                                                               "Name",
                                                               "Description")]
                subjs <- sapply(samps,
                                function(x) paste(strsplit(x,
                                                           split = ".",
                                                           fixed = T)[[1]][1],
                                                  strsplit(x,
                                                           split = ".",
                                                           fixed = T)[[1]][2],
                                                  sep = "."))
                toBind <- data.frame(sample = samps,
                                     subj = subjs,
                                     tissue = rep(tiss, length(samps)))
                toBindMetDt <- metDat1[match(gsub(".", "-", toBind$sample, fixed = T),
                                             metDat1$SAMPID), ]
                toBindMetDt <- toBindMetDt[, colnames(toBindMetDt) != "SAMPID"]
                toBind <- cbind.data.frame(toBind, toBindMetDt)
                
                sampInfo <- rbind.data.frame(sampInfo, toBind)
        }
        close(pb)
        out <- list(sample_info = sampInfo, combinedDF = brainDF)
        return(out)
}

getENSEMBL_ID_info <- function(ensembl_list,
                               mart,
                               removeNAs = T,
                               filt = "ensembl_gene_id_version"){
        #ensembl_list <- gsub("\\_.*", "", ensembl_list)
        print("Mapping the ensembl_gene_id_version...")
        if(filt == "ensembl_gene_id_version"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "ensembl_gene_id_version",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }else if(filt == "ensembl_gene_id"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }
        gene_coords <- getBM(attributes = atts,
                             filters = filt,
                             values = gsub("_PAR_Y", "", ensembl_list),
                             mart = mart)
        
        gene_coords$size <- gene_coords$end_position - gene_coords$start_position
        if(filt == "ensembl_gene_id_version"){
                mapped_matchVec <- match(gsub("_PAR_Y", "", ensembl_list),
                                         gene_coords$ensembl_gene_id_version)
                outDF <- data.frame(ensembl_gene_id_version_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    ensembl_gene_id_version_new = gene_coords$ensembl_gene_id_version[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                unMapped <- outDF$ensembl_gene_id_version_orig[is.na(outDF$ensembl_gene_id_version_new)]
                
                unMapped_toMap <- gsub("\\..*", "", unMapped)
                
                print("Mapping the ensembl_gene_id_version IDs that couldn't be mapped to any of those to the respective ensembl_gene_id...")
                gene_coords_unmapped <- getBM(attributes = atts,
                                              filters = "ensembl_gene_id",
                                              values = unMapped_toMap,
                                              mart = mart)
                
                gene_coords_unmapped$size <- gene_coords_unmapped$end_position - gene_coords_unmapped$start_position
                
                gene_coords_unmapped <- gene_coords_unmapped[match(unMapped_toMap, gene_coords_unmapped$ensembl_gene_id), ]
                
                
                isUnmap <- is.na(outDF$ensembl_gene_id_version_new)
                outDF$ensembl_gene_id[isUnmap] <- gene_coords_unmapped$ensembl_gene_id
                outDF$ensembl_gene_id_version_new[isUnmap] <- gene_coords_unmapped$ensembl_gene_id_version
                outDF$symbol[isUnmap] <- gene_coords_unmapped$hgnc_symbol
                outDF$chromosome[isUnmap] <- gene_coords_unmapped$chromosome_name
                outDF$start[isUnmap] <- gene_coords_unmapped$start_position
                outDF$end[isUnmap] <- gene_coords_unmapped$end_position
                outDF$size[isUnmap] <- gene_coords_unmapped$size
                naVec <- is.na(outDF$ensembl_gene_id_version_new)
        }else if(filt == "ensembl_gene_id"){
                mapped_matchVec <- match(gsub("_PAR_Y", "", ensembl_list),
                                         gene_coords$ensembl_gene_id)
                outDF <- data.frame(ensembl_gene_id_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                naVec <- is.na(outDF$ensembl_gene_id)
        }
        if(removeNAs){
                nNAs <- sum(naVec)
                propNAs <- round(nNAs/nrow(outDF) * 100, digits = 2)
                outDF <- outDF[!naVec, ]
                print(sprintf("%s ensembl_gene_ids were removed due to inability to map to ENSEMBL database.", as.character(nNAs)))
                print(sprintf("This represents the %s%% of the total number of ensembl_gene_ids submitted (n = %s)",
                              as.character(propNAs),
                              as.character(length(ensembl_list))))
        }
        return(outDF)
}

sum_dups <- function(m) {
        rowsum(m, group = rownames(m), reorder = FALSE)
}

# Load metadata files
################################################################################
ind_mdat <- read.csv(ind_mdat_file, row.names = 1)
samp_mdat <- read.csv(samp_mdat_file, row.names = 1)


# Parse GTEx counts
################################################################################

metDat1 <- read.delim(metDat1_f, sep = '\t', comment.char = "#")

countOut <- createGTExMat(gtex_counts_dir)

gtex_counts <- countOut$combinedDF

gtex_geneInfo <- gtex_counts[, c("Name", "Description")]

################################################################################
# Load RNAseqHarmonization counts                                              #
################################################################################

rnaseq_harm_count_files <- list.files(rnaseq_harm_counts_dir, full.names = TRUE)

rnaseq_harm_count_files <- rnaseq_harm_count_files[
        grep("matrix_clean.txt", rnaseq_harm_count_files, fixed = TRUE)
]


rnaseqHarm_counts <- sapply(
        rnaseq_harm_count_files,
        function(x) plotUtils::read_table_fast(x, sep = "\t", row.names = 1)
)

names(rnaseqHarm_counts) <- make.unique(
        gsub("\\_.*", "", basename(rnaseq_harm_count_files))
)

#rnaseqHarm_counts <- do.call(cbind.data.frame, rnaseqHarm_counts)

#colnames(rnaseqHarm_counts) <- sapply(
#        colnames(rnaseqHarm_counts),
#        function(x) strsplit(x, split = ".", fixed = TRUE)[[1]][4]
#)

# ROSMAP
################################################################################

rosmap_counts <- cbind.data.frame(
        rnaseqHarm_counts$ROSMAP,
        rnaseqHarm_counts$ROSMAP.1,
        rnaseqHarm_counts$ROSMAP.2,
        rnaseqHarm_counts$ROSMAP.3
)

rosmap_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/RNAseq_Harmonization/RNAseq_Harmonization_ROSMAP_combined_metadata.csv")
rosmap_mdat[!make.names(rosmap_mdat$specimenID) %in% colnames(rosmap_counts), ]

rosmap_mdat[rosmap_mdat$individualID == "R9938478", ]

colnames(rosmap_counts)[!colnames(rosmap_counts) %in% make.names(samp_mdat$specimenID)]
rosmap_specs <- make.names(samp_mdat$specimenID[samp_mdat$substudy == "ROSMAP"])
rosmap_not_in_counts <- rosmap_specs[!rosmap_specs %in% colnames(rosmap_counts)]
rosmap_mdat[rosmap_mdat$specimenID %in% rosmap_not_in_counts, ]

# There are 2 samples in the expression data not in the metadata. Remove them
sampls_not_in_mdat_rosmap <- colnames(rosmap_counts)[!colnames(rosmap_counts) %in% make.names(samp_mdat$specimenID)]
rosmap_counts <- rosmap_counts[, !colnames(rosmap_counts) %in% sampls_not_in_mdat_rosmap]
# There are 16 specimenIDs in the integrated metadata that are not in counts
# file, as they are from other assays than RNAseq. Remove them.
samp_mdat <- samp_mdat[!samp_mdat$specimenID %in% rosmap_not_in_counts, ]


# MSBB
################################################################################

length(colnames(rnaseqHarm_counts$MSBB)[!colnames(rnaseqHarm_counts$MSBB) %in% make.names(samp_mdat$specimenID)])

msbb_in_counts_not_in_mdat <- colnames(rnaseqHarm_counts$MSBB)[
        !colnames(rnaseqHarm_counts$MSBB) %in% make.names(samp_mdat$specimenID)
]

# There are 24 specimens in the RNAseq counts matrix from MSBB that are not
# in the sample metadata file, which was built with the
# RNAseq_Harmonization_MSBB_combined_metadata.csv file. Check
# MSBB_biospecimen_metadata.csv (syn21893059)
# 23 of these 24 specimens that are not in the metadata are samples with individualID
# == "". We removed them from metadata, so do the same with counts

MSBB_biosp_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/MSBB/MSBB_biospecimen_metadata.csv")

MSBB_comb_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/RNAseq_Harmonization/RNAseq_Harmonization_MSBB_combined_metadata.csv")

MSBB_comb_mdat[match(msbb_in_counts_not_in_mdat, MSBB_comb_mdat$specimenID), ]

colnames(rnaseqHarm_counts$MSBB)[!colnames(rnaseqHarm_counts$MSBB) %in% make.names(MSBB_comb_mdat$specimenID)]

# There's just one that's in the counts file, but not in the metadata file:
# hB_RNA_10892
colnames(rnaseqHarm_counts$MSBB)[!colnames(rnaseqHarm_counts$MSBB) %in% make.names(MSBB_biosp_mdat$specimenID)]

MSBB_comb_mdat$specimenID[MSBB_comb_mdat$individualID == ""] %in% msbb_in_counts_not_in_mdat

msbb_counts <- rnaseqHarm_counts$MSBB[, !colnames(rnaseqHarm_counts$MSBB) %in% msbb_in_counts_not_in_mdat]

# Mayo
################################################################################

length(colnames(rnaseqHarm_counts$Mayo)[!colnames(rnaseqHarm_counts$Mayo) %in% make.names(samp_mdat$specimenID)])
samp_mdat$specimenID[samp_mdat$substudy == "Mayo"][!make.names(samp_mdat$specimenID[samp_mdat$substudy == "Mayo"]) %in% colnames(rnaseqHarm_counts$Mayo)]

#View(samp_mdat[samp_mdat$substudy == "Mayo", ])

# All the samples in the Mayo counts file are in the metadata file.
# 46 samples in the metadata file are not in the counts file

mayo_not_in_counts <- samp_mdat$specimenID[samp_mdat$substudy == "Mayo"][!make.names(samp_mdat$specimenID[samp_mdat$substudy == "Mayo"]) %in% colnames(rnaseqHarm_counts$Mayo)]

mayo_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/RNAseq_Harmonization/RNAseq_Harmonization_Mayo_combined_metadata.csv")

mayo_mdat[match(mayo_not_in_counts, mayo_mdat$specimenID), ]

samp_mdat[match(mayo_not_in_counts, samp_mdat$specimenID), ]

samp_mdat[match(mayo_not_in_counts, samp_mdat$specimenID), ]$specimenID


# Load counts files from  syn51401682
mayo_ctrls <- read.table("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/Mayo/MayoRNAseq_RNAseq_PSPPACtrl_geneCounts.tsv.gz", sep = "\t", header = TRUE)

make.names(mayo_not_in_counts) %in% colnames(mayo_ctrls)


sum(colnames(mayo_ctrls) %in% make.names(samp_mdat$specimenID))

intersect(make.names(colnames(mayo_ctrls)), make.names(colnames(rnaseqHarm_counts$Mayo)))
#mayo_ctrls[, "X11333_TCX"]

grep("18202", colnames(rnaseqHarm_counts$Mayo))

# Check Mayo data from https://www.synapse.org/Synapse:syn20818650
################################################################################

# Check syn8466815 from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn9702085

mayo_tcx <- read.table("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/RNASeq_Harm/MAYO_CBE_TCX_Counts.tsv", sep = "\t", header = TRUE)

make.names(mayo_not_in_counts) %in% colnames(mayo_tcx)

# No overlap

mayo_sampswap <- read.table("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/RNASeq_Harm/Mayo_sampleSwapSeqAtMSBB_gene_all_counts_matrix_clean.txt", sep = "\t", header = TRUE)

make.names(mayo_not_in_counts) %in% mayo_sampswap


# Now check syn51401791
mayo_tcx_all <- read.table("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/Mayo/MayoRNAseq_RNAseq_TCX_geneCounts-278.tsv.gz", se= "\t", header = TRUE)

dim(mayo_tcx_all)

dim(mayo_tcx)


length(intersect(colnames(mayo_tcx_all), colnames(rnaseqHarm_counts$Mayo)))/ncol(rnaseqHarm_counts$Mayo)

# Compare harmonized metadata and original metadatas to see how many individuals
# the harmonized one is discarding.
mayo_indiv_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/RNAseq_Harmonization/MayoRNAseq_individual_metadata.csv")

unique_tcx_inds <- unique(gsub("\\_.*", "", colnames(mayo_tcx)[6:ncol(mayo_tcx)]))

mayo_indiv_mdat[match(unique_tcx_inds, make.names(mayo_indiv_mdat$individualID)), ]

mayo_comb_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/RNAseq_Harmonization/RNAseq_Harmonization_Mayo_combined_metadata.csv")

# 265 individuals in mayo with rnaseq data
mayo_biosp_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/RNAseq_Harmonization/MayoRNAseq_biospecimen_metadata.csv")

mayo_inds_in_biosp <- intersect(
        mayo_indiv_mdat$individualID, mayo_biosp_mdat$individualID
)

mayo_biosp_mdat_filt <- mayo_biosp_mdat[mayo_biosp_mdat$individualID %in% mayo_inds_in_biosp & mayo_biosp_mdat$assay == "rnaSeq", ]

sum(!mayo_indiv_mdat$individualID %in% mayo_comb_mdat$individualID)
sum(!unique(mayo_biosp_mdat_filt$individualID) %in% mayo_comb_mdat$individualID)


# Considerable overlap (50%)

# 62 individuals from MSBB individual metadata not in the harmonization study
msbb_indiv_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/MSBB/MSBB_individual_metadata.csv")
sum(!msbb_indiv_mdat$individualID %in% MSBB_comb_mdat$individualID)

rosmap_clin_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/ROSMAP/ROSMAP_clinical.csv")

rosmap_comb_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/RNAseq_Harmonization/RNAseq_Harmonization_ROSMAP_combined_metadata.csv")
sum(!rosmap_clin_mdat$individualID %in% rosmap_comb_mdat$individualID)
sum(!rosmap_comb_mdat$individualID %in% rosmap_clin_mdat$individualID)

rosmap_biosp <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/ROSMAP/ROSMAP_biospecimen_metadata.csv")


rosmap_inds_in_biosp <- intersect(rosmap_clin_mdat$individualID, rosmap_biosp$individualID)

rosmap_biosp_filt <- rosmap_biosp[rosmap_biosp$individualID %in% rosmap_inds_in_biosp, ]

# 257 individuals in ROSMAP with RNAseq data but not in RNASeq Harmonization
# study
sum(!unique(rosmap_biosp_filt$individualID[rosmap_biosp_filt$assay == "rnaSeq"]) %in% rosmap_comb_mdat$individualID)

################################################################################
# There are 257 + 62 + 265 individuals with RNAseq data that are not in the    #
# RNAseq harmonization data. We might either include them for the GWAS or      #  
# keep them out for validation. For now let's just use the ones we selected    #
# (the) harmonized ones, and keep the rest in the mind for using as validation.#
################################################################################

samp_mdat <- samp_mdat[!samp_mdat$specimenID %in% mayo_not_in_counts, ]

mayo_counts <- rnaseqHarm_counts$Mayo

rnaseqHarm_counts_matrix <- cbind.data.frame(
        rosmap_counts,
        msbb_counts,
        mayo_counts
) 

rnaseqHarm_counts_matrix <- rnaseqHarm_counts_matrix[
        5:nrow(rnaseqHarm_counts_matrix),
] %>% as.matrix()

################################################################################
# LBP                                                                          #
################################################################################

lbp_counts <- read.csv(file.path(lbp_counts_dir, "LBP_FlagshipPaper_featureCounts.csv"), row.names = 1)

all(colnames(lbp_counts) %in% samp_mdat$specimenID)

LBP_not_in_counts <- samp_mdat$specimenID[samp_mdat$substudy == "LBP"][!samp_mdat$specimenID[samp_mdat$substudy == "LBP"] %in% colnames(lbp_counts)]

table(samp_mdat[match(LBP_not_in_counts, samp_mdat$specimenID), ]$assay)

# The 542 specimens that are not in the counts file are either snrnaSeq or
# wholeGenomeSeq, so just remove them from metadata

samp_mdat <- samp_mdat[!samp_mdat$specimenID %in% LBP_not_in_counts, ]

sum(unique(samp_mdat$individualID) %in% ind_mdat$individualID[ind_mdat$has_RNAseq & ind_mdat$has_WGS])

table(samp_mdat$substudy)

# GTEX
################################################################################



sum(unique(samp_mdat$individualID) %in% ind_mdat$individualID[ind_mdat$has_RNAseq & ind_mdat$has_WGS & !ind_mdat$used_in_clock_training])


gtex_geneInfo

gtex_counts <- gtex_counts[, 3:ncol(gtex_counts)]
rownames(gtex_counts) <- gtex_geneInfo$Name


all(colnames(gtex_counts) %in% make.names(samp_mdat$specimenID))

all(make.names(samp_mdat$specimenID[samp_mdat$substudy == "GTEx"]) %in% colnames(gtex_counts))

# All GTEx corresponds well to metadata

# brainSeq
################################################################################

brainSeqDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/brainSeq"
bspI_counts_f <- file.path(brainSeqDir, "brainseq_pI_counts.csv")
bspII_counts_f <- file.path(brainSeqDir, "brainseq_pII_counts.csv")
bspI_mdat_f <- file.path(brainSeqDir, "brainseq_pI_metadata.csv")
bspII_mdat_f <- file.path(brainSeqDir, "brainseq_pII_metadata.csv")

bspI_counts <- read.csv(bspI_counts_f, row.names = 1)
bspII_counts <- read.csv(bspII_counts_f, row.names = 1)
bspI_mdat <- read.csv(bspI_mdat_f, row.names = 1)
bspII_mdat <- read.csv(bspII_mdat_f, row.names = 1)

bs_genes <- intersect(rownames(bspI_counts), rownames(bspII_counts))
bs_counts <- cbind.data.frame(bspI_counts[bs_genes, ], bspII_counts[bs_genes, ])
bs_mdat <- rbind.data.frame(bspI_mdat, bspII_mdat)

colnames(bs_counts) == bs_mdat$specimenID

# TBI
################################################################################

tbi_counts_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/counts_parsed/tbi_counts.csv"
tbi_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/tbi_mdat.csv"

tbi_counts <- read.csv(tbi_counts_f, row.names = 1)
tbi_mdat <- read.csv(tbi_mdat_f, row.names = 1)

dim(tbi_counts)
dim(tbi_mdat)

all(colnames(tbi_counts) == tbi_mdat$specimenID)

# NABEC
################################################################################

nabec_counts <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/counts_parsed/nabec_counts_parsed.csv", row.names = 1)
nabec_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/nabec_metdat_parsed.csv", row.names = 1)

dim(nabec_counts)
dim(nabec_mdat)

# Merge counts
################################################################################

rnaseqHarm_counts_matrix


dup_genes_gtex <- gsub(
        "\\..*",
        "",
        gtex_geneInfo$Name
)[which(duplicated(gsub("\\..*", "", gtex_geneInfo$Name)))]

gtex_geneInfo[gsub("\\..*", "", gtex_geneInfo$Name) %in% dup_genes_gtex, ]



geneInfo_GTEx_f <- file.path(outDir, "geneInfo_GTEx.tsv")
geneInfo_LBP_f <- file.path(outDir, "geneInfo_LBP.tsv")
geneInfo_rnaseqHarm_f <- file.path(outDir, "geneInfo_rnaseqHarm.tsv")

if (!file.exists(geneInfo_GTEx_f)){
        library(biomaRt)
        human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        geneInfo_GTEx <- getENSEMBL_ID_info(
                gtex_geneInfo$Name,
                mart = human,
                filt = "ensembl_gene_id_version"
        )
        write.table(
                geneInfo_GTEx,
                sep = '\t',
                col.names = TRUE,
                file = geneInfo_GTEx_f,
                row.names = FALSE
        )
}else{
        geneInfo_GTEx <- read.table(
                geneInfo_GTEx_f, sep = '\t', header = TRUE
        )
}

if (!file.exists(geneInfo_LBP_f)){
        library(biomaRt)
        human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        geneInfo_LBP <- getENSEMBL_ID_info(
                rownames(lbp_counts),
                mart = human,
                filt = "ensembl_gene_id_version"
        )
        write.table(
                geneInfo_LBP,
                sep = '\t',
                col.names = TRUE,
                file = geneInfo_LBP_f,
                row.names = FALSE
        )
}else{
        geneInfo_LBP <- read.table(
                geneInfo_LBP_f, sep = '\t', header = TRUE
        )
}

if (!file.exists(geneInfo_rnaseqHarm_f)){
        library(biomaRt)
        human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        geneInfo_rnaSeqHarm <- getENSEMBL_ID_info(
                rownames(rnaseqHarm_counts_matrix),
                mart = human,
                filt = "ensembl_gene_id_version"
        )
        write.table(
                geneInfo_rnaSeqHarm,
                sep = '\t',
                col.names = TRUE,
                file = geneInfo_rnaseqHarm_f,
                row.names = FALSE
        )
}else{
        geneInfo_rnaSeqHarm <- read.table(
                geneInfo_rnaseqHarm_f, sep = '\t', header = TRUE
        )
}

clock_coefs <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/models/secnd_round/mod_alpha1_coefs.csv")

all(clock_coefs$ensembl_gene_id %in% geneInfo_GTEx$ensembl_gene_id)

all(clock_coefs$ensembl_gene_id %in% geneInfo_LBP$ensembl_gene_id)

all(clock_coefs$ensembl_gene_id %in% geneInfo_rnaSeqHarm$ensembl_gene_id)

head(geneInfo_GTEx)

any(duplicated(geneInfo_GTEx$ensembl_gene_id))


rnaseqHarm_nuNames <- geneInfo_rnaSeqHarm$ensembl_gene_id[
        match(
                rownames(rnaseqHarm_counts_matrix),
                geneInfo_rnaSeqHarm$ensembl_gene_id_version_orig
        )
]

rnaseqHarm_counts_matrix <- rnaseqHarm_counts_matrix[!is.na(rnaseqHarm_nuNames), ]
rownames(rnaseqHarm_counts_matrix) <- rnaseqHarm_nuNames[!is.na(rnaseqHarm_nuNames)]

rnaseqHarm_counts_matrix <- sum_dups(rnaseqHarm_counts_matrix)


GTEx_nuNames <- geneInfo_GTEx$ensembl_gene_id[
        match(
                rownames(gtex_counts),
                geneInfo_GTEx$ensembl_gene_id_version_orig
        )
]

gtex_counts <- as.matrix(gtex_counts)

gtex_counts <- gtex_counts[!is.na(GTEx_nuNames), ]
rownames(gtex_counts) <- GTEx_nuNames[!is.na(GTEx_nuNames)]

gtex_counts <- sum_dups(gtex_counts)


LBP_nuNames <- geneInfo_LBP$ensembl_gene_id[
        match(
                rownames(lbp_counts),
                geneInfo_LBP$ensembl_gene_id_version_orig
        )
]

lbp_counts <- as.matrix(lbp_counts)

lbp_counts <- lbp_counts[!is.na(LBP_nuNames), ]
rownames(lbp_counts) <- LBP_nuNames[!is.na(LBP_nuNames)]

lbp_counts <- sum_dups(lbp_counts)

commGenes <- rownames(rnaseqHarm_counts_matrix) %>%
        intersect(rownames(gtex_counts)) %>%
        intersect(rownames(lbp_counts)) %>%
        intersect(rownames(bs_counts)) %>%
        intersect(rownames(tbi_counts)) %>%
        intersect(rownames(nabec_counts))

all(clock_coefs$ensembl_gene_id %in% commGenes)

counts_merged <- cbind(
        gtex_counts[commGenes, ],
        rnaseqHarm_counts_matrix[commGenes, ],
        lbp_counts[commGenes, ],
        bs_counts[commGenes, ],
        tbi_counts[commGenes, ],
        nabec_counts[commGenes, ]
)


dim(counts_merged)

dim(samp_mdat)

samp_mdat_rnaseq <- samp_mdat[samp_mdat$assay == "rnaSeq", ]

head(samp_mdat_rnaseq)

# Add the rest of the metadata columns
bs_mdat$cohort <- bs_mdat$substudy
bs_mdat$diagnosis <- "control"
bs_mdat <- cbind.data.frame(
        bs_mdat,
        data.frame(
                matrix(
                        nrow = nrow(bs_mdat),
                        ncol = length(colnames(samp_mdat_rnaseq)[!colnames(samp_mdat_rnaseq) %in% colnames(bs_mdat)]),
                        dimnames = list(rownames(bs_mdat), colnames(samp_mdat_rnaseq)[!colnames(samp_mdat_rnaseq) %in% colnames(bs_mdat)])
                )
        )
)

tbi_mdat$cohort <- tbi_mdat$TBI

tbi_mdat$diagnosis <- ifelse(tbi_mdat$diagn_4BrainClck == "Control", "control", "traumatic brain injury")


colnames(samp_mdat_rnaseq)[!colnames(samp_mdat_rnaseq) %in% colnames(tbi_mdat)]

tbi_mdat <- cbind.data.frame(
        tbi_mdat,
        data.frame(
                matrix(
                        nrow = nrow(tbi_mdat),
                        ncol = length(colnames(samp_mdat_rnaseq)[!colnames(samp_mdat_rnaseq) %in% colnames(tbi_mdat)]),
                        dimnames = list(rownames(tbi_mdat), colnames(samp_mdat_rnaseq)[!colnames(samp_mdat_rnaseq) %in% colnames(tbi_mdat)])
                )
        )
)

colnames(samp_mdat_rnaseq)[!colnames(samp_mdat_rnaseq) %in% colnames(tbi_mdat)]

colnames(samp_mdat_rnaseq)[!colnames(samp_mdat_rnaseq) %in% colnames(nabec_mdat)]

samp_mdat_rnaseq <- rbind.data.frame(
        samp_mdat_rnaseq,
        bs_mdat[, colnames(samp_mdat_rnaseq)],
        tbi_mdat[, colnames(samp_mdat_rnaseq)],
        nabec_mdat[colnames(samp_mdat_rnaseq)]
)

samp_mdat_rnaseq <- samp_mdat_rnaseq[
        match(colnames(counts_merged),
              make.names(samp_mdat_rnaseq$specimenID)), 
]

table(ind_mdat[!ind_mdat$individualID %in% samp_mdat_rnaseq$individualID, ]$substudy)

inds_not_in_counts <- ind_mdat[!ind_mdat$individualID %in% samp_mdat_rnaseq$individualID, ]

sum(unique(samp_mdat_rnaseq$individualID) %in% ind_mdat$individualID[ind_mdat$has_RNAseq & ind_mdat$has_WGS & !ind_mdat$used_in_clock_training])

# Check where the individuals not in the sample metadata, and therefore not
# in the expression data, are coming from:
mayo_comb_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/RNAseq_Harmonization/RNAseq_Harmonization_Mayo_combined_metadata.csv")

inds_not_in_counts_mayo <- inds_not_in_counts$individualID[inds_not_in_counts$substudy == "Mayo"]

mayo_comb_mdat[mayo_comb_mdat$individualID %in% inds_not_in_counts_mayo, ]
# For mayo there are the TCX that don't appear in the harmonized counts.


inds_not_in_counts_lbp <- inds_not_in_counts$individualID[inds_not_in_counts$substudy == "LBP"]

lbp_bsp_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/AMP_AD/LBP/LBP_biospecimen_metadata.csv")

table(lbp_bsp_mdat[lbp_bsp_mdat$individualID %in% inds_not_in_counts_lbp, ]$assay)
# For LBP they are individuals that only have snRNAseq and/or WGS

inds_not_in_counts_gtex <- inds_not_in_counts$individualID[inds_not_in_counts$substudy == "GTEx"]

gtex_samp_mdat <- read.delim(
        "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/GTEx/phs000424.v11.pht002743.v10.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz",
        sep = '\t',
        comment.char = "#"
)

gtex_samp_mdat$individualID <- sapply(
        gtex_samp_mdat$SAMPID,
        function(x) paste(
                strsplit(x, split = "-")[[1]][1],
                strsplit(x, split = "-")[[1]][2],
                sep = "-"
        )
)


head(gtex_samp_mdat[gtex_samp_mdat$individualID %in% inds_not_in_counts_gtex, ])

gtex_samp_mdat_brain <- gtex_samp_mdat[gtex_samp_mdat$SMTS == "Brain", ]

table(
        gtex_samp_mdat[gtex_samp_mdat$individualID %in% inds_not_in_counts_gtex, ]$SMTS,
        gtex_samp_mdat[gtex_samp_mdat$individualID %in% inds_not_in_counts_gtex, ]$SMMTRLTP
)

table(
        gtex_samp_mdat_brain[gtex_samp_mdat_brain$individualID %in% inds_not_in_counts_gtex, ]$SMTSD,
        gtex_samp_mdat_brain[gtex_samp_mdat_brain$individualID %in% inds_not_in_counts_gtex, ]$SMGEBTCHT
)

table(
        gtex_samp_mdat_brain[gtex_samp_mdat_brain$individualID %in% inds_not_in_counts_gtex, ]$SMTSD,
        gtex_samp_mdat_brain[gtex_samp_mdat_brain$individualID %in% inds_not_in_counts_gtex, ]$SMAFRZE
)


table(ind_mdat[!ind_mdat$individualID %in% samp_mdat$individualID, ]$substudy)

sum(make.names(gtex_samp_mdat_brain$SAMPID[gtex_samp_mdat_brain$SMMTRLTP == "RNA:Total RNA"]) %in% colnames(countOut$combinedDF)[3:ncol(countOut$combinedDF)])
sum(!make.names(gtex_samp_mdat_brain$SAMPID[gtex_samp_mdat_brain$SMMTRLTP == "RNA:Total RNA"]) %in% colnames(countOut$combinedDF)[3:ncol(countOut$combinedDF)])

gtex_samps_not_in_counts <- gtex_samp_mdat_brain$SAMPID[gtex_samp_mdat_brain$SMMTRLTP == "RNA:Total RNA"][!make.names(gtex_samp_mdat_brain$SAMPID[gtex_samp_mdat_brain$SMMTRLTP == "RNA:Total RNA"]) %in% colnames(countOut$combinedDF)[3:ncol(countOut$combinedDF)]]

gtex_samps_not_in_counts

unique(
        sapply(
                gtex_samps_not_in_counts,
                function(x) paste(
                        strsplit(x, split = "-")[[1]][1],
                        strsplit(x, split = "-")[[1]][2],
                        sep = "-"
                )
        )
)

# For GTEx these are samples that were excluded from the final freeze

ind_mdat_filt <- ind_mdat[ind_mdat$individualID %in% samp_mdat_rnaseq$individualID, ]


(ind_mdat_filt$individualID[ind_mdat_filt$substudy == "GTEx"] %in% unique(gtex_samp_mdat$individualID[gtex_samp_mdat$SMAFRZE %in% c("WGS", "DEEPWGS")])) == ind_mdat_filt$has_WGS[ind_mdat_filt$substudy == "GTEx"]

inds_not_in_counts_gtex %in% sapply(colnames(countOut$combinedDF)[3:ncol(countOut$combinedDF)], function(x) paste(strsplit(x, split = ".", fixed = TRUE)[[1]][1:2], collapse = "-"))

outDir

samp_mdat_rnaseq$substudy_batch <- samp_mdat_rnaseq$substudy

samp_mdat_rnaseq$substudy_batch[samp_mdat_rnaseq$substudy == "ROSMAP" & make.names(samp_mdat_rnaseq$specimenID) %in% colnames(rnaseqHarm_counts$ROSMAP)] <- "ROSMAP_1"
samp_mdat_rnaseq$substudy_batch[samp_mdat_rnaseq$substudy == "ROSMAP" & make.names(samp_mdat_rnaseq$specimenID) %in% colnames(rnaseqHarm_counts$ROSMAP.1)] <- "ROSMAP_2"
samp_mdat_rnaseq$substudy_batch[samp_mdat_rnaseq$substudy == "ROSMAP" & make.names(samp_mdat_rnaseq$specimenID) %in% colnames(rnaseqHarm_counts$ROSMAP.2)] <- "ROSMAP_3"
samp_mdat_rnaseq$substudy_batch[samp_mdat_rnaseq$substudy == "ROSMAP" & make.names(samp_mdat_rnaseq$specimenID) %in% colnames(rnaseqHarm_counts$ROSMAP.3)] <- "ROSMAP_4"


write.csv(counts_merged, file = file.path(outDir, "counts_merged.csv"))

write.csv(samp_mdat_rnaseq, file = file.path(outDir, "sample_metadata_merged.csv"))


# Add BS, NABEC and TBI individual metadata
nabec_ind_mdat <- read.csv(
        "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/nabec_individual_metdat_parsed.csv",
        row.names = 1
)

colnames(ind_mdat_filt)[!colnames(ind_mdat_filt) %in% colnames(nabec_ind_mdat)]


bs_ind_mdat <- read.csv(
        "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/brainSeq_individual_metdat_parsed.csv",
        row.names = 1
)

colnames(ind_mdat_filt)[!colnames(ind_mdat_filt) %in% colnames(bs_ind_mdat)]

colnames(bs_ind_mdat)[!colnames(bs_ind_mdat) %in% colnames(tbi_mdat)]

tbi_ind_mdat <- tbi_mdat

colnames(bs_ind_mdat)[!colnames(bs_ind_mdat) %in% colnames(tbi_ind_mdat)]
tbi_ind_mdat$has_RNAseq <- TRUE
tbi_ind_mdat$has_WGS <- FALSE
tbi_ind_mdat$has_snpArray <- FALSE
tbi_ind_mdat$used_in_clock_training <- FALSE
tbi_ind_mdat <- unique(tbi_ind_mdat[, colnames(bs_ind_mdat)])

length(unique(nabec_ind_mdat$individualID)) == nrow(nabec_ind_mdat)
length(unique(tbi_ind_mdat$individualID)) == nrow(tbi_ind_mdat)
length(unique(bs_ind_mdat$individualID)) == nrow(bs_ind_mdat)

bs_dup_inds <- bs_ind_mdat$individualID[duplicated(bs_ind_mdat$individualID)]

bs_ind_mdat$ageDeath <- round(bs_ind_mdat$ageDeath)
bs_ind_mdat <- unique(bs_ind_mdat)

length(unique(bs_ind_mdat$individualID)) == nrow(bs_ind_mdat)

ind_mdat_filt <- rbind.data.frame(
        ind_mdat_filt,
        bs_ind_mdat[bs_ind_mdat$individualID %in% samp_mdat_rnaseq$individualID, colnames(ind_mdat_filt)],
        nabec_ind_mdat[nabec_ind_mdat$individualID %in% samp_mdat_rnaseq$individualID, colnames(ind_mdat_filt)],
        tbi_ind_mdat[tbi_ind_mdat$individualID %in% samp_mdat_rnaseq$individualID, colnames(ind_mdat_filt)]
)

length(unique(ind_mdat_filt$individualID))

dim(ind_mdat_filt)

ind_mdat_filt[
        ind_mdat_filt$has_RNAseq & (ind_mdat_filt$has_WGS | ind_mdat_filt$has_snpArray),
]

write.csv(ind_mdat_filt, file = file.path(outDir, "individual_metadata_merged_rnaseq_filt.csv"))

# Generate a version with only samples not used during training
counts_merged_not_train <- counts_merged[, 
                                         samp_mdat_rnaseq$individualID[match(colnames(counts_merged), make.names(samp_mdat_rnaseq$specimenID))] %in% ind_mdat_filt$individualID[!ind_mdat_filt$used_in_clock_training]
]

dim(counts_merged_not_train)


samp_mdat_rnaseq_not_train <- samp_mdat_rnaseq[
        match(
                colnames(counts_merged_not_train),
                make.names(samp_mdat_rnaseq$specimenID)
        ), 
]

write.csv(counts_merged_not_train, file = file.path(outDir, "counts_merged_not_train.csv"))

write.csv(samp_mdat_rnaseq_not_train, file = file.path(outDir, "sample_metadata_merged_not_train.csv"))
