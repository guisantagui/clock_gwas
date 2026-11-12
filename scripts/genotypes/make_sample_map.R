################################################################################
# Map the sample IDs of each WGS callset to the individualIDs used by the     #
# clock (individual_metadata_merged_rnaseq_filt.csv), and flag the            #
# individuals that have an age gap.                                            #
#                                                                              #
# AMP-AD joint VCF: WGS specimenID -> individualID through the WGS            #
#                   biospecimen metadata (ROSMAP, MSBB, Mayo).                 #
# LBP:              CRAM/VCF sample name = specimenID (5721-EC-*) ->           #
#                   individualID through LBP_biospecimen_metadata.csv.        #
# GTEx:             VCF sample = donor ID (GTEX-xxxx) = individualID.         #
# NABEC:            dbGaP sample-info.csv from phg001711 (optional; give its  #
#                   path with --nabec_info once copied from Iris).            #
#                                                                              #
# Output (contains controlled-access IDs: keep it in data/, never commit):    #
#   data/genotypes/sample_map.tsv  cohort, vcf_id, individualID, has_gap, sex  #
#   data/genotypes/keep_<cohort>.txt  VCF IDs of individuals with a gap        #
################################################################################
suppressPackageStartupMessages(library(argparser))
p <- arg_parser("Build the WGS sample map.")
p <- add_argument(p, "--meta", default = "data/metadata/AMP_AD", help = "AMP-AD metadata dir.")
p <- add_argument(p, "--ind", default = "results/counts_parsed/individual_metadata_merged_rnaseq_filt.csv",
                  help = "Individual metadata.")
p <- add_argument(p, "--gaps", default = "results/clock_cv/individual_gaps.csv", help = "Individual gaps.")
p <- add_argument(p, "--gtex_samples", default = NA, help = "bcftools query -l of the GTEx VCF (optional).")
p <- add_argument(p, "--nabec_info", default = NA, help = "NABEC sample-info.csv from phg001711 (optional).")
p <- add_argument(p, "--outdir", default = "data/genotypes", help = "Output directory.")
a <- parse_args(p)
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

ind <- read.csv(a$ind, row.names = 1, stringsAsFactors = FALSE)
gaps <- read.csv(a$gaps, stringsAsFactors = FALSE)
gap_ids <- gaps$individualID[!is.na(gaps$gap_mean)]

rd <- function(f) read.csv(file.path(a$meta, f), stringsAsFactors = FALSE, colClasses = "character")

# AMP-AD joint VCF (ROSMAP + MSBB + Mayo)
bio <- do.call(rbind, lapply(c("ROSMAP", "MSBB", "MayoRNAseq"), function(co) {
        b <- rd(sprintf("WGS_Harmonization/%s_biospecimen_metadata.csv", co))
        w <- rd(sprintf("WGS_Harmonization/%s_assay_wholeGenomeSeq_metadata.csv", co))
        data.frame(cohort = "ampad", vcf_id = w$specimenID,
                   individualID = b$individualID[match(w$specimenID, b$specimenID)])
}))

# LBP
lb <- rd("LBP/LBP_biospecimen_metadata.csv")
lb <- lb[lb$assay == "wholeGenomeSeq", ]
lbp <- data.frame(cohort = "lbp", vcf_id = lb$specimenID, individualID = lb$individualID)

# GTEx: donor IDs are used directly
gtex_ids <- if (!is.na(a$gtex_samples)) readLines(a$gtex_samples) else
        ind$individualID[ind$substudy == "GTEx"]
gtex <- data.frame(cohort = "gtex", vcf_id = gtex_ids, individualID = gtex_ids)

maps <- list(bio, lbp, gtex)

# NABEC
if (!is.na(a$nabec_info)) {
        # dbGaP sample-info.csv: '#' lines describe the columns, no header row.
        # Column 1 = SubjectID (as in the RNA-seq metadata, e.g. UKY-1066),
        # column 2 = IID, the sample ID used in the VCF (e.g. KEN-1066).
        ni <- read.csv(a$nabec_info, header = FALSE, comment.char = "#",
                       stringsAsFactors = FALSE, colClasses = "character")
        ni <- data.frame(individualID = trimws(ni[[1]]), vcf_id = trimws(ni[[2]]))
        nab_ids <- ind$individualID[ind$substudy == "NABEC"]
        message(sprintf("NABEC: %d WGS samples, %d match RNA-seq individualIDs",
                        nrow(ni), sum(ni$individualID %in% nab_ids)))
        flag_wrong <- ni$individualID %in% nab_ids &
                !(ind$has_WGS | ind$has_snpArray)[match(ni$individualID, ind$individualID)]
        if (any(flag_wrong)) {
                message(sum(flag_wrong), " NABEC individuals have WGS but has_WGS is FALSE in ", a$ind,
                        " (e.g. ", paste(head(ni$individualID[flag_wrong], 3), collapse = ", "), ")")
        }
        maps[[4]] <- data.frame(cohort = "nabec", vcf_id = ni$vcf_id, individualID = ni$individualID)
}

m <- do.call(rbind, maps)
m$in_clock_metadata <- m$individualID %in% ind$individualID
m$has_gap <- m$individualID %in% gap_ids
m$sex <- ind$sex[match(m$individualID, ind$individualID)]
write.table(m, file.path(a$outdir, "sample_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (co in unique(m$cohort)) {
        writeLines(unique(m$vcf_id[m$cohort == co & m$has_gap]),
                   file.path(a$outdir, sprintf("keep_%s.txt", co)))
}
print(with(m, table(cohort, has_gap)))
if (any(duplicated(m$individualID[m$has_gap]))) {
        message("Individuals with more than one WGS sample (keep one after KING):")
        print(m[m$individualID %in% m$individualID[m$has_gap][duplicated(m$individualID[m$has_gap])], ])
}
