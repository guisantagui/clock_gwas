################################################################################
# Brain clock: Parse GTEx counts and metadata files into a single counts file  #
# and a curated metadata file                                                  #
################################################################################

setwd("../..")

brainCountsDir <- "./data/gtex/counts/"
outDir <- "./results/parsed/"

metDatFile1 <- "./data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

metDat1 <- read.delim(metDatFile1)

metDatFile2 <- "./data/gtex/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"

metDat2 <- read.delim(metDatFile2)

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

countOut <- createGTExMat(brainCountsDir)

brainCounts <- countOut$combinedDF


gtex_metDat_unified <- data.frame(specimenID = metDat1$SAMPID[make.names(metDat1$SAMPID) %in% colnames(brainCounts)],
                                  platform = metDat1$SMGEBTCHT[make.names(metDat1$SAMPID) %in% colnames(brainCounts)],
                                  RIN = metDat1$SMRIN[make.names(metDat1$SAMPID) %in% colnames(brainCounts)],
                                  libraryPrep = rep(NA, sum(make.names(metDat1$SAMPID) %in% colnames(brainCounts))),
                                  libraryPreparationMethod =rep(NA, sum(make.names(metDat1$SAMPID) %in% colnames(brainCounts))),
                                  runType =rep(NA, sum(make.names(metDat1$SAMPID) %in% colnames(brainCounts))),
                                  readLength = metDat1$SMRDLGTH[make.names(metDat1$SAMPID) %in% colnames(brainCounts)])

gtex_metDat_unified$individualID <- sapply(gtex_metDat_unified$specimenID,
                                           function(x) paste(strsplit(x, split = "-")[[1]][1],
                                                             strsplit(x, split = "-")[[1]][2],
                                                             sep = "-"))

gtex_metDat_unified$organ <- tolower(metDat1$SMTS[make.names(metDat1$SAMPID) %in% colnames(brainCounts)])
gtex_metDat_unified$tissue <- metDat1$SMTSD[make.names(metDat1$SAMPID) %in% colnames(brainCounts)]

gtex_metDat_unified$tissue <- gsub("Brain - ", "", gtex_metDat_unified$tissue)
gtex_metDat_unified$assay <- metDat1$SMAFRZE[make.names(metDat1$SAMPID) %in% colnames(brainCounts)]
gtex_metDat_unified$assay <- gsub("RNASEQ", "rnaSeq", gtex_metDat_unified$assay)
gtex_metDat_unified$exclude <- rep(NA, nrow(gtex_metDat_unified))
gtex_metDat_unified$excludeReason <- rep("", nrow(gtex_metDat_unified))
gtex_metDat_unified$sex <- c("male",
                             "female")[metDat2$SEX[match(gtex_metDat_unified$individualID,
                                                         metDat2$SUBJID)]]

gtex_metDat_unified$race <- metDat2$RACE[match(gtex_metDat_unified$individualID,
                                               metDat2$SUBJID)]

# Not clear what the race numbers mean. The field is not included in metadata
# explanatory file

gtex_metDat_unified$ageDeath <- metDat2$AGE[match(gtex_metDat_unified$individualID,
                                                  metDat2$SUBJID)]

gtex_metDat_unified$apoeGenotype <- rep(NA, nrow(gtex_metDat_unified))
gtex_metDat_unified$pmi <- metDat2$TRISCH[match(gtex_metDat_unified$individualID,
                                                metDat2$SUBJID)]
gtex_metDat_unified$Braak <- rep(NA, nrow(gtex_metDat_unified))
gtex_metDat_unified$diagn_4BrainClck <- rep("Control", nrow(gtex_metDat_unified))
gtex_metDat_unified$substudy <- rep("GTEx", nrow(gtex_metDat_unified))
gtex_metDat_unified$batch_rna <- metDat1$SMNABTCH[make.names(metDat1$SAMPID) %in% colnames(brainCounts)]
gtex_metDat_unified$batch_lib <- rep(NA, nrow(gtex_metDat_unified))
gtex_metDat_unified$batch_seq <- metDat1$SMGEBTCH[make.names(metDat1$SAMPID) %in% colnames(brainCounts)]

# These categories code neurodegeneration, schizophrenia and drug abuse.
# We are going to remove the samples that have drug abuse or schizophrenia from
# the dataset and encode the ones that have 0 in all the variables as controls
sum(metDat2$MHALZDMT != 0)
sum(metDat2$MHALZHMR != 0)
metDat2$MHALS
metDat2$MHCOCAINE5
sum(metDat2$MHPRKNSN != 0)
metDat2$MHSCHZ
metDat2$MHSDRGABS
metDat2$MHDMNTIA

indivsDrugSchiz <- metDat2$SUBJID[metDat2$MHSCHZ != 0 | metDat2$MHCOCAINE5 != 0 | metDat2$MHSDRGABS]
ctrls <- metDat2$SUBJID[metDat2$MHALZDMT == 0 &
                                metDat2$MHALZHMR == 0 &
                                metDat2$MHALS == 0 &
                                metDat2$MHPRKNSN == 0 &
                                metDat2$MHDMNTIA == 0]

gtex_metDat_unified <- gtex_metDat_unified[!gtex_metDat_unified$individualID %in% indivsDrugSchiz, ]
gtex_metDat_unified$diagn_4BrainClck[!gtex_metDat_unified$individualID %in% ctrls] <- "Rest"


write.csv(gtex_metDat_unified,
          file = paste0(outDir, "GTEx_metadata_unified.csv"))

# Remove from counts dataframe the samples that are not in metadata and
# save as csv.

removed <- colnames(brainCounts)[!colnames(brainCounts) %in% c(colnames(brainCounts)[1:3], make.names(gtex_metDat_unified$specimenID))]


# Check if we're removing something we shouldnt
all(sapply(removed, function(x) paste(strsplit(x, split = ".", fixed = T)[[1]][1],
                                      strsplit(x, split = ".", fixed = T)[[1]][2],
                                      sep = ".")) %in% make.names(indivsDrugSchiz))
# All fine

brainCounts <- brainCounts[, colnames(brainCounts) %in% c(colnames(brainCounts)[1:3],
                                                          make.names(gtex_metDat_unified$specimenID))]

write.csv(brainCounts, file = paste0(outDir, "GTEx_allCounts.csv"))