################################################################################
# Brain clock: Merge GTEx, LBP and RNAseq Harmonization datasets               #
################################################################################

if (!require("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
}
if (!require("biomaRt")){
        BiocManager::install("biomaRt", update = F)
}
if (!require(plotUtils)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}

library(biomaRt)
library(plotUtils)

# Directory stuff
################################################################################
setwd("../..")

outDir <- "./results/parsed/"

create_dir_if_not(outDir)

# Functions
################################################################################

# This function accepts a vector of ensembl IDs and returns info about 
# the updated version, symbols, start pos, end pos and size
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

sum_dups <- function(m){
        dup_genes <- rownames(m)[duplicated(rownames(m))]
        for(i in seq_along(dup_genes)){
                dg <- dup_genes[i]
                nu_vec <- colSums(m[rownames(m) == dg, ])
                matches <- which(rownames(m) == dg)
                m[matches[1], ] <- nu_vec
                m <- m[-matches[-1], ]
        }
        return(m)
}

# Load the data
################################################################################
LBP_counts_file <- "./data/expression/LBP/LBP_FlagshipPaper_featureCounts.csv"
GTEx_counts_file <- "./results/parsed/GTEx_allCounts.csv"
RNAseqHarm_counts_file <- "./results/parsed/RNAseqHarm_allCounts.csv"

LBP_metDat_file <- "./results/parsed/LBP_metadata_unified.csv"
GTEx_metDat_file <- "./results/parsed/GTEx_metadata_unified.csv"
RNAseqHarm_metDat_file <- "./results/metadata_parsed/RNAseq_Harmonization_ind_all_ROSMAPBtch.csv"

LBP_geneInfo_file <- "./results/parsed/LBP_geneInfo.csv"
GTEx_geneInfo_file <- "./results/parsed/GTEx_geneInfo.csv"
RNAseqHarm_geneInfo_file <- "./results/parsed/RNAseqHarm_geneInfo.csv"

LBP_counts <- read_table_fast(LBP_counts_file, row.names = 1)
GTEx_counts <- read_table_fast(GTEx_counts_file, row.names = 1)
RNAseqHarm_counts <- read_table_fast(RNAseqHarm_counts_file, row.names = 1)

LBP_metDat <- read.csv(LBP_metDat_file, row.names = 1)
GTEx_metDat <- read.csv(GTEx_metDat_file, row.names = 1)
RNAseqHarm_metDat <- read.csv(RNAseqHarm_metDat_file, row.names = 1)

# Create gene info dataframes (for homogenizing gene IDs).
if (!(file.exists(LBP_geneInfo_file) & file.exists(GTEx_geneInfo_file) & file.exists(RNAseqHarm_geneInfo_file))){
        human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
        LBP_geneInfo <- getENSEMBL_ID_info(rownames(LBP_counts), mart = human)
        GTEx_geneInfo <- getENSEMBL_ID_info(GTEx_counts$Name,
                                            mart = human,
                                            filt = "ensembl_gene_id_version")
        RNAseqHarm_geneInfo <- getENSEMBL_ID_info(rownames(RNAseqHarm_counts),
                                                  mart = human)

        # Add "PAR_Y" to the gene_id of the genes that map to chromosome Y
        GTEx_geneInfo$ensembl_gene_id[grepl("PAR_Y", GTEx_geneInfo$ensembl_gene_id_version_orig)] <- paste0(GTEx_geneInfo$ensembl_gene_id[grepl("PAR_Y", GTEx_geneInfo$ensembl_gene_id_version_orig)], "_PAR_Y")
        # Add "PAR_Y" to the gene_id of the genes that map to chromosome Y
        LBP_geneInfo$ensembl_gene_id[grepl("PAR_Y", LBP_geneInfo$ensembl_gene_id_version_orig)] <- paste0(LBP_geneInfo$ensembl_gene_id[grepl("PAR_Y", LBP_geneInfo$ensembl_gene_id_version_orig)], "_PAR_Y")
        RNAseqHarm_geneInfo$ensembl_gene_id[grepl("PAR_Y", RNAseqHarm_geneInfo$ensembl_gene_id_version_orig)] <- paste0(RNAseqHarm_geneInfo$ensembl_gene_id[grepl("PAR_Y", RNAseqHarm_geneInfo$ensembl_gene_id_version_orig)], "_PAR_Y")

        write.csv(LBP_geneInfo, paste0(outDir, "LBP_geneInfo.csv"))
        write.csv(GTEx_geneInfo, paste0(outDir, "GTEx_geneInfo.csv"))
        write.csv(RNAseqHarm_geneInfo, paste0(outDir, "RNAseqHarm_geneInfo.csv"))
}else{
        LBP_geneInfo <- read.csv(LBP_geneInfo_file, row.names = 1)
        GTEx_geneInfo <- read.csv(GTEx_geneInfo_file, row.names = 1)
        RNAseqHarm_geneInfo <- read.csv(RNAseqHarm_geneInfo_file, row.names = 1)
}

# Merge counts DFs
################################################################################

LBP_counts <- LBP_counts[rownames(LBP_counts) %in% LBP_geneInfo$ensembl_gene_id_version_orig, ]

any(duplicated(LBP_geneInfo$ensembl_gene_id[match(rownames(LBP_counts), LBP_geneInfo$ensembl_gene_id_version_orig)]))

rownames(LBP_counts) <- LBP_geneInfo$ensembl_gene_id[match(rownames(LBP_counts), LBP_geneInfo$ensembl_gene_id_version_orig)]

# Remove the samples that are not in the cleaned-up metadata
LBP_counts <- LBP_counts[, make.names(colnames(LBP_counts)) %in% LBP_metDat$specimenID]

rownames(GTEx_counts) <- GTEx_counts$Name

GTEx_counts <- GTEx_counts[, !colnames(GTEx_counts) %in% c("id", "Name", "Description")]

GTEx_counts <- GTEx_counts[rownames(GTEx_counts) %in% GTEx_geneInfo$ensembl_gene_id_version_orig, ]

any(duplicated(GTEx_geneInfo$ensembl_gene_id[match(rownames(GTEx_counts), GTEx_geneInfo$ensembl_gene_id_version_orig)]))

rownames(GTEx_counts) <- GTEx_geneInfo$ensembl_gene_id[match(rownames(GTEx_counts), GTEx_geneInfo$ensembl_gene_id_version_orig)]


RNAseqHarm_counts[1:10, 1:10]

RNAseqHarm_counts <- RNAseqHarm_counts[rownames(RNAseqHarm_counts) %in% RNAseqHarm_geneInfo$ensembl_gene_id_version_orig, ]

any(duplicated(RNAseqHarm_geneInfo$ensembl_gene_id[match(rownames(RNAseqHarm_counts), RNAseqHarm_geneInfo$ensembl_gene_id_version_orig)]))

rownames(RNAseqHarm_counts) <- RNAseqHarm_geneInfo$ensembl_gene_id[match(rownames(RNAseqHarm_counts),
                                                                         RNAseqHarm_geneInfo$ensembl_gene_id_version_orig)]

commonGenes <- intersect(rownames(RNAseqHarm_counts),
                         intersect(rownames(GTEx_counts),
                                   rownames(LBP_counts)))

combinedCounts <- rbind.data.frame(t(LBP_counts[commonGenes, ]), 
                                   t(GTEx_counts[commonGenes, ]),
                                   t(RNAseqHarm_counts[commonGenes, ]))

# Remove the genes that have zero counts across all the samples
combinedCounts <- combinedCounts[, colSums(combinedCounts) != 0]
# Save merged counts file
write_table_fast(combinedCounts, f = paste0(outDir, "combined_counts.csv"))

# Combine metadata dataframes into a single file
################################################################################
colnames(RNAseqHarm_metDat) <- gsub("ROSMAP_batch",
                                    "batch_seq",
                                    colnames(RNAseqHarm_metDat))

RNAseqHarm_metDat$batch_rna <- rep(NA, nrow(RNAseqHarm_metDat))
RNAseqHarm_metDat$batch_lib <- rep(NA, nrow(RNAseqHarm_metDat))

RNAseqHarm_metDat <- RNAseqHarm_metDat[, c(colnames(GTEx_metDat),
                                           "mmse30_first_ad_dx",
                                           "mmse30_lv")]

GTEx_metDat$mmse30_first_ad_dx <- rep(NA, nrow(GTEx_metDat))
GTEx_metDat$mmse30_lv <- rep(NA, nrow(GTEx_metDat))

LBP_metDat$mmse30_first_ad_dx <- rep(NA, nrow(LBP_metDat))
LBP_metDat$mmse30_lv <- rep(NA, nrow(LBP_metDat))

combined_metDat <- rbind.data.frame(GTEx_metDat,
                                    LBP_metDat,
                                    RNAseqHarm_metDat)

# Convert PMI units of metadata to the same (hours)
pmi_hours <- combined_metDat$pmi[grepl("hour",
                                       combined_metDat$pmi)]

pmi_hours_conv <- c()
for(i in seq_along(pmi_hours)){
        pmi <- pmi_hours[i]
        pmi_split <- strsplit(pmi, split = ",")[[1]]
        hours <- gsub("[^0-9.-]", "", pmi_split[1])
        mins <- gsub("[^0-9.-]", "", pmi_split[2])
        converted <- as.numeric(hours) + as.numeric(mins)/60
        pmi_hours_conv <- c(pmi_hours_conv, converted)
}

pmi_colon <- combined_metDat$pmi[grepl(":", combined_metDat$pmi)]

pmi_colon_conv <- c()
for(i in seq_along(pmi_colon)){
        pmi <- pmi_colon[i]
        pmi_split <- strsplit(pmi, split = ":", fixed = T)[[1]]
        hours <- gsub("[^0-9.-]", "", pmi_split[1])
        mins <- gsub("[^0-9.-]", "", pmi_split[2])
        converted <- as.numeric(hours) + as.numeric(mins)/60
        pmi_colon_conv <- c(pmi_colon_conv, converted)
}

combined_metDat$pmi[grepl(":",
                          combined_metDat$pmi)] <- as.character(pmi_colon_conv)
combined_metDat$pmi[grepl("hour",
                          combined_metDat$pmi)] <- as.character(pmi_hours_conv)
combined_metDat$pmi <- as.numeric(combined_metDat$pmi)

# MSBB pmi is in minutes, so convert to hours
combined_metDat$pmi[combined_metDat$substudy == "MSBB"] <- combined_metDat$pmi[combined_metDat$substudy == "MSBB"]/60

combined_metDat$tissue <- tolower(combined_metDat$tissue)

# Keep only RNAseq data
combined_metDat <- combined_metDat[combined_metDat$assay == "rnaSeq", ]

write.csv(combined_metDat, file = paste0(outDir, "combined_metDat.csv"))
#write.csv(geneInfo_combined, file = paste0(outDir, "combined_geneInfo.csv"))