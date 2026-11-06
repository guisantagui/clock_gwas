################################################################################
# Create merged metadata file for WGS data.
################################################################################

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/"

wgs_gtex_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/gtex_biosp_wgs_mdat_parsed.csv"
wgs_ampad_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/wgs_amp_ad_mdat_comb.csv"
indiv_mdat_combined_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed/individual_metadata_combined.csv"

wgs_gtex_mdat <- read.csv(wgs_gtex_mdat_f, row.names = 1)
wgs_ampad_mdat <- read.csv(wgs_ampad_mdat_f, row.names = 1)
indiv_mdat_combined <- read.csv(indiv_mdat_combined_f, row.names = 1)

head(wgs_gtex_mdat)
head(wgs_ampad_mdat)

wgs_ampad_mdat$assay <- gsub("wholeGenomeSeq", "WGS", wgs_ampad_mdat$assay)

common_colnames <- intersect(colnames(wgs_gtex_mdat), colnames(wgs_ampad_mdat))

wgs_mdat_comb <- rbind.data.frame(
        wgs_gtex_mdat[, common_colnames],
        wgs_ampad_mdat[, common_colnames]
)

match(wgs_mdat_comb$individualID, indiv_mdat_combined$individualID)

table(wgs_ampad_mdat$assay)
head(indiv_mdat_combined)


wgs_gtex_mdat[!make.names(wgs_gtex_mdat$individualID) %in% make.names(indiv_mdat_combined$individualID), ]
wgs_ampad_mdat[!make.names(wgs_ampad_mdat$individualID) %in% make.names(indiv_mdat_combined$individualID), ]

wgs_mdat_comb[!make.names(wgs_mdat_comb$individualID) %in% make.names(indiv_mdat_combined$individualID), ]


indiv_mdat_combined$individualID[indiv_mdat_combined$substudy != "GTEx" & indiv_mdat_combined$has_WGS & !indiv_mdat_combined$individualID %in% wgs_ampad_mdat$individualID]

# All the individual IDs in the combined individual metadata are in the wgs
# data. There are some samples from AMP-AD that have WGS data but not RNAseq, 
# but we can remove them.

wgs_mdat_comb <- wgs_mdat_comb[make.names(wgs_mdat_comb$individualID) %in% make.names(indiv_mdat_combined$individualID), ]

wgs_mdat_comb <- cbind.data.frame(
        wgs_mdat_comb,
        indiv_mdat_combined[match(wgs_mdat_comb$individualID, indiv_mdat_combined$individualID), colnames(indiv_mdat_combined) != "individualID"]
)


write.csv(wgs_mdat_comb, file = file.path(outDir, "wgs_metadata_mapping.csv"))
