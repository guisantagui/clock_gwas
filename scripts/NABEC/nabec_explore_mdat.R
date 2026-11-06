subject_pheno <- read.table("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/NABEC/phs001300.v6.pht006722.v2.p1.c1.NABEC_Subject_Phenotypes.GRU.txt")
dim(subject_pheno)

mRNA_sampAtts <- read.table(
        "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/NABEC/phs001353.v6.pht006729.v1.p1.c1.NABEC_Frontal_Cortex_mRNA_Seq_Sample_Attributes.GRU.txt",
        sep = "\t", header = TRUE
)

mRNA_sampAtts$subjectID <- sapply(
        mRNA_sampAtts$SAMPLE_ID,
        function(x) paste(strsplit(x, split = "-")[[1]][1:2], collapse = "-")
)

totalRNA_sampAtts <- read.table(
        "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/NABEC/phs001354.v6.pht006731.v1.p1.c1.NABEC_Frontal_Cortex_Total_RNA_Seq_Sample_Attributes.GRU.txt",
        sep = "\t", header = TRUE
)

totalRNA_sampAtts$subjectID <- sapply(
        totalRNA_sampAtts$SAMPLE_ID,
        function(x) paste(strsplit(x, split = "-")[[1]][1:2], collapse = "-")
)

wgs_sampAtts <- read.table(
        "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/NABEC/phs002636.v5.pht011803.v1.p1.c1.NABEC_WGS_Sample_Attributes.GRU.txt",
        sep = "\t", header = TRUE
)

sum(wgs_sampAtts$SAMPLE_ID %in% mRNA_sampAtts$subjectID)
sum(wgs_sampAtts$SAMPLE_ID %in% totalRNA_sampAtts$subjectID | wgs_sampAtts$SAMPLE_ID %in% mRNA_sampAtts$subjectID)

length(intersect(wgs_sampAtts$SAMPLE_ID, totalRNA_sampAtts$subjectID))

totalRNA_sampAtts

counts_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/NABEC/counts.csv"
counts <- as.data.frame(data.table::fread(counts_f))


sra_runTable_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/NABEC/NABEC_SraRunTable.csv"

totalRNA_sampAtts$subjectID %in% subject_pheno$dbGaP_Subject_ID
counts[1:10, 1:10]

sra_runTable <- read.csv(sra_runTable_f)

keep_RNAstudies <- c(
        "NABEC: mRNA Sequencing of Human Cerebral Frontal Cortex",
        "NABEC: Total RNA Sequencing of Human Cerebral Frontal Cortex"
)
sra_runTable <- sra_runTable[sra_runTable$study_name %in% keep_RNAstudies & sra_runTable$Assay.Type == "RNA-Seq", ]


colnames(counts) <- sra_runTable$biospecimen_repository_sample_id[match(colnames(counts), sra_runTable$Run)]
rownames(counts) <- counts[, 1]
counts <- counts[, 2:ncol(counts)]
counts[1:10, 1:10]

totalRNA_sampAtts$analyte_type_fine <- "totalRNA"
mRNA_sampAtts$analyte_type_fine <- "mRNA"

rna_sampAtts <- rbind.data.frame(
        totalRNA_sampAtts,
        mRNA_sampAtts
)

rna_sampAtts[!rna_sampAtts$subjectID %in% subject_pheno$dbGaP_Subject_ID, ]

rna_sampAtts$age <- subject_pheno$AGE[
        match(
                make.names(rna_sampAtts$subjectID),
                make.names(subject_pheno$dbGaP_Subject_ID)
        )
]        

unique(rna_sampAtts[!rna_sampAtts$subjectID %in% subject_pheno$dbGaP_Subject_ID, ]$subjectID)
   
colnames(counts) %in% totalRNA_sampAtts$SAMPLE_ID
colnames(counts) %in% mRNA_sampAtts$SAMPLE_ID



mdat_4bas <- rna_sampAtts[match(colnames(counts), rna_sampAtts$SAMPLE_ID), ]
counts <- counts[, !is.na(mdat_4bas$SAMPLE_ID) & !is.na(mdat_4bas$age)]
mdat_4bas <- mdat_4bas[!is.na(mdat_4bas$SAMPLE_ID) & !is.na(mdat_4bas$age), ]
rownames(mdat_4bas) <- mdat_4bas$SAMPLE_ID

write.csv(
        mdat_4bas,
        file = "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/NABEC/rna_mdat_merged.csv"
)

write.csv(
        counts,
        file = "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/NABEC/counts_filt.csv"
)

rownames(counts) <- gsub("\\..*", "", rownames(counts))

library(brainAgeShiftR)

bas <- create_brainAgeShiftObj(counts = counts, metadata = mdat_4bas)

bas <- normalizeCounts(bas, useTrainMeans = TRUE)

bas <- predictAge(bas)
bas$metadata

plot(bas$metadata$age, bas$metadata$predicted_age)

library(ggplot2)

ggplot(
        data = bas$metadata,
        mapping = aes(x = age, y = predicted_age, colour = analyte_type_fine)
) + geom_point()

cor(bas$metadata$age, bas$metadata$predicted_age, method = "spearman")

?create_brainAgeShiftObj

sra_runTable$biospecimen_repository_sample_id %in% totalRNA_sampAtts$SAMPLE_ID
