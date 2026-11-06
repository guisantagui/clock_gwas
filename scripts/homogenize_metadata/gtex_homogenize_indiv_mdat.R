################################################################################
# GWAS of Brain clock: go through the individual metadata files of each one of #
# the studies to determine which individuals are controls and which not, and   #
# generate a unified individual metadata file that we can use to assess what   # 
# individuals have both RNAseq and genomics data                               #
################################################################################

wd <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/"
setwd(wd)

sample_mdat_f <- "./data/metadata/GTEx/phs000424.v11.pht002743.v10.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz"
indivd_mdat_f <- "./data/metadata/GTEx/phs000424.v11.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz"
outDir <- "./results/metadata_parsed/"


sample_mdat <- read.delim(sample_mdat_f, sep = '\t', comment.char = "#")
indivd_mdat <- read.delim(indivd_mdat_f, sep = '\t', comment.char = "#")

dim(sample_mdat)
head(sample_mdat)


table(sample_mdat$SMMTRLTP)
sample_mdat[sample_mdat$SMMTRLTP == "DNA:DNA Genomic", ]

head(indivd_mdat)
#sample_mdat <- sample_mdat[sample_mdat$SMTS == "Brain", ]

sample_mdat$individualID <- sapply(
        sample_mdat$SAMPID,
        function(x) paste(
                strsplit(x, split = "-")[[1]][1],
                strsplit(x, split = "-")[[1]][2],
                sep = "-"
        )
)

inds_w_brainExp <- unique(sample_mdat$individualID[sample_mdat$SMTS == "Brain"])

#length(unique(sample_mdat$individualID))
indivd_mdat <- indivd_mdat[indivd_mdat$SUBJID %in% inds_w_brainExp, ]
indivd_mdat$has_RNAseq <- indivd_mdat$SUBJID %in% sample_mdat$individualID[sample_mdat$SMMTRLTP ==  "RNA:Total RNA" & sample_mdat$SMTS == "Brain"]
indivd_mdat$has_WGS <- indivd_mdat$SUBJID %in% sample_mdat$individualID[sample_mdat$SMMTRLTP %in%  c("DNA:DNA Genomic", "DNA:DNA Somatic")]
indivd_mdat$has_snpArray <- FALSE
sum(indivd_mdat$has_RNAseq & (indivd_mdat$has_WGS | indivd_mdat$has_snpArray))
sum(!(indivd_mdat$has_RNAseq & (indivd_mdat$has_WGS | indivd_mdat$has_snpArray)))

indivs_without_genomics <- indivd_mdat$SUBJID[!(indivd_mdat$has_RNAseq & (indivd_mdat$has_WGS | indivd_mdat$has_snpArray))]
dim(indivd_mdat)
#View(indivd_mdat)

# Some individuals appear to not have genomics data. Let's see why
sample_mdat[grepl(indivs_without_genomics[2], sample_mdat$SAMPID), c("SAMPID", "SMTS", "SMMTRLTP")]


clock_train_test_f <- "./data/metadata/brain_clock/train_test.csv"

clock_train_test <- read.csv(clock_train_test_f, row.names = 1)

clock_training_inds <- unique(
        clock_train_test$individualID[clock_train_test$train_test_split == "train"]
)

indivd_mdat$used_in_clock_training <- make.names(
        indivd_mdat$SUBJID
) %in% make.names(clock_training_inds)
sum(indivd_mdat$used_in_clock_training)

indivd_mdat$SEX <- ifelse(indivd_mdat$SEX == 1, "male", "female")

race_encoding <- data.frame(
        value = 1:4,
        race = c("Asian", "Black", "White", "American Indian or Alaska Native")
)

indivd_mdat$RACE <- race_encoding$race[
        match(indivd_mdat$RACE, race_encoding$value)
]

table(indivd_mdat$ETHNCTY)

indivd_mdat$pmi <- indivd_mdat$TRISCHD / 60

#head(indivd_mdat)
#head(amp_ad_ind_mdat_merged)

indivd_mdat$substudy <- "GTEx"
indivd_mdat$apoeGenotype <- NA
indivd_mdat[, c("SUBJID", "COHORT", "substudy", "SEX", "RACE", "AGE")]

#amp_ad_ind_mdat_merged$used_in_clock_training

indivd_mdat$diagnosis <- NA

indivsDrugSchiz <- indivd_mdat$SUBJID[indivd_mdat$MHSCHZ != 0 | indivd_mdat$MHCOCAINE5 != 0 | indivd_mdat$MHSDRGABS]
ctrls <- indivd_mdat$SUBJID[
        indivd_mdat$MHALZDMT == 0 &
        indivd_mdat$MHALZHMR == 0 &
        indivd_mdat$MHALS == 0 &
        indivd_mdat$MHPRKNSN == 0 &
        indivd_mdat$MHDMNTIA == 0 &
        indivd_mdat$MHCVD == 0
]
ctrls <- ctrls[!is.na(ctrls)]
indivd_mdat$diagnosis[indivd_mdat$SUBJID %in% ctrls] <- "control"
table(indivd_mdat$MHDMNTIA[indivd_mdat$MHALZHMR == 1])
table(indivd_mdat$MHALZDMT[indivd_mdat$MHALZHMR == 1])
indivd_mdat$diagnosis[indivd_mdat$MHALZDMT == 1] <- "Alzheimer Disease"
indivd_mdat$diagnosis[indivd_mdat$MHDMNTIA == 1] <- "other dementia"
indivd_mdat$diagnosis[indivd_mdat$MHALS == 1] <- "Amyotropic Lateral Sclerosis"
indivd_mdat$diagnosis[indivd_mdat$MHPRKNSN == 1] <- "Parkinson's disease"
indivd_mdat$diagnosis[indivd_mdat$MHSCHZ == 1] <- "Schizophrenia"
indivd_mdat$diagnosis[indivd_mdat$MHCVD == 1] <- "Cerebrovascular disease"
indivd_mdat[is.na(indivd_mdat$diagnosis), c("MHALZDMT", "MHALZHMR", "MHALS", "MHPRKNSN", "MHSCHZ", "MHCVD")]
indivd_mdat$diagnosis[is.na(indivd_mdat$diagnosis)] <- "Unknown"


#indivd_mdat$
#head(indivd_mdat)

indivd_mdat[is.na(indivd_mdat$diagnosis), 71:94]

which(colnames(indivd_mdat) == "MHALS")
colnames(indivd_mdat)[71:94]


indivd_mdat_parsed <- indivd_mdat[, c("SUBJID", "COHORT", "SEX", "RACE", "AGE", "pmi", "diagnosis", "has_RNAseq", "has_WGS", "has_snpArray", "used_in_clock_training")]

indivd_mdat_parsed$substudy <- "GTEx"

indivd_mdat_parsed$apoeGenotype <- NA
indivd_mdat_parsed$Braak <- NA
indivd_mdat_parsed$CERAD <- NA
indivd_mdat_parsed$CDR <- NA
indivd_mdat_parsed$plaqueMean <- NA
indivd_mdat_parsed$thal <- NA


ctrls_brn_clck <- indivd_mdat$SUBJID[
        indivd_mdat$MHALZDMT == 0 &
        indivd_mdat$MHALZHMR == 0 &
        indivd_mdat$MHALS == 0 &
        indivd_mdat$MHPRKNSN == 0 &
        indivd_mdat$MHDMNTIA == 0
]

ctrls_brn_clck <- ctrls_brn_clck[!is.na(ctrls_brn_clck)]

indivd_mdat_parsed$diagn_4BrainClck <- "Control"

indivd_mdat_parsed$diagn_4BrainClck[!indivd_mdat_parsed$SUBJID %in% ctrls_brn_clck] <- "Rest"

colnames(indivd_mdat_parsed) <- c(
        "individualID",
        "cohort",
        "sex",
        "race",
        "ageDeath",
        "pmi",
        "diagnosis",
        "has_RNAseq",
        "has_WGS",
        "has_snpArray",
        "used_in_clock_training",
        "substudy",
        "apoeGenotype",
        "Braak",
        "CERAD",
        "CDR",
        "plaqueMean",
        "thal",
        "diagn_4BrainClck"
)

indivd_mdat_parsed <- indivd_mdat_parsed[, c(
        "individualID",
        "cohort",
        "substudy",
        "sex",
        "race",
        "ageDeath",
        "apoeGenotype",
        "pmi",
        "diagnosis",
        "Braak",
        "CERAD",
        "CDR",
        "plaqueMean",
        "thal",
        "diagn_4BrainClck",
        "has_RNAseq",
        "has_WGS",
        "has_snpArray",
        "used_in_clock_training"
)]

write.csv(
        indivd_mdat_parsed,
        file = file.path(outDir, "gtex_indiv_mdat_parsed.csv")
)

# Parse sample-level metadata
spec_mdat_parsed <- data.frame(
        specimenID = sample_mdat$SAMPID,
        platform = sample_mdat$SMGEBTCH,
        RIN = sample_mdat$SMRIN,
        libraryPrep = rep(NA, nrow(sample_mdat)),
        libraryPreparationMethod = rep(NA, nrow(sample_mdat)),
        runType = rep(NA, nrow(sample_mdat)),
        readLength = sample_mdat$SMRDLGTH,
        individualID = sample_mdat$individualID,
        organ = tolower(sample_mdat$SMTS),
        tissue = sapply(sample_mdat$SMTSD,
                        function(x) ifelse(
                                length(strsplit(x, split = " - ")[[1]]) == 2,
                                strsplit(x, split = " - ")[[1]][2],
                                strsplit(x, split = " - ")[[1]][1])),
        assay = gsub("RNASEQ", "rnaSeq", sample_mdat$SMAFRZE)
)

spec_mdat_parsed$assay[spec_mdat_parsed$assay == ""] <- NA

# keep only brain samples
spec_mdat_parsed_brain_rnaseq <- spec_mdat_parsed[
        spec_mdat_parsed$organ == "brain" &
        spec_mdat_parsed$assay == "rnaSeq" &
        !is.na(spec_mdat_parsed$assay),
]

spec_mdat_parsed_wgs <- spec_mdat_parsed[
        (spec_mdat_parsed$assay == "WGS" | spec_mdat_parsed$assay == "DEEPWGS") & spec_mdat_parsed$individualID %in% spec_mdat_parsed_brain_rnaseq$individualID,
]

spec_mdat_parsed_wgs <- spec_mdat_parsed_wgs[!is.na(spec_mdat_parsed_wgs$specimenID), ]

spec_mdat_parsed_brain_rnaseq$exclude <- NA
spec_mdat_parsed_brain_rnaseq$excludeReason <- ""
spec_mdat_parsed_brain_rnaseq$sex <- indivd_mdat_parsed$sex[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$race <- indivd_mdat_parsed$race[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$apoeGenotype <- indivd_mdat_parsed$apoeGenotype[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$ageDeath <- indivd_mdat_parsed$ageDeath[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$pmi <- indivd_mdat_parsed$pmi[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$Braak <- indivd_mdat_parsed$Braak[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$CERAD <- indivd_mdat_parsed$CERAD[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$diagn_4BrainClck <- indivd_mdat_parsed$diagn_4BrainClck[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$diagnosis <- indivd_mdat_parsed$diagnosis[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$CDR <- indivd_mdat_parsed$CDR[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$plaqueMean <- indivd_mdat_parsed$plaqueMean[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$thal <- indivd_mdat_parsed$thal[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$substudy <- indivd_mdat_parsed$substudy[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]

spec_mdat_parsed_brain_rnaseq$cohort <- indivd_mdat_parsed$cohort[
        match(
                spec_mdat_parsed_brain_rnaseq$individualID,
                indivd_mdat_parsed$individualID
        )
]


write.csv(
        spec_mdat_parsed_brain_rnaseq,
        file = file.path(outDir, "gtex_biosp_brainRnaseq_mdat_parsed.csv")
)

table(spec_mdat_parsed$assay)
spec_mdat_parsed_wgs <- spec_mdat_parsed[
        spec_mdat_parsed$assay %in% c("DEEPWGS", "WGS") & spec_mdat_parsed$individualID %in% spec_mdat_parsed_brain_rnaseq$individualID,
]

write.csv(
        spec_mdat_parsed_wgs,
        file = file.path(outDir, "gtex_biosp_wgs_mdat_parsed.csv")
)
