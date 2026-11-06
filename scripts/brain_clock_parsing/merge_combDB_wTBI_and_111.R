################################################################################
# Brain clock:  Merge counts data from databases with the perturbation files   #
# Javier and Sascha generated.                                                 #
################################################################################

countsFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_counts_wTBI.csv"
counts_111_sas <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/111_sascha_complete_counts_df.csv"
counts_111_jav <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/javier_111_counts.csv"
metDatSasFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/Sascha_metadata.csv"
metDatJavFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/javier_metadata.csv"
metDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat_wTBI.csv"
outCountsFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_counts_wTBI_wPert111.csv"
outMetDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat_wTBI_wPert111.csv"

# Functions
################################################################################

# read.csv but faster
readCsvFst <- function(path){
        df <- data.frame(data.table::fread(path))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

# write.csv, but faster
writeCsvFst <- function(df, file, rowNames = T, colNames = T){
        if(rowNames){
                rn <- rownames(df)
                df <- data.table::data.table(df)
                df[, V1 := rn]
                data.table::setcolorder(df, c("V1", setdiff(names(df), "V1")))
        }
        data.table::fwrite(df, file, col.names = colNames)
}

# Load data
################################################################################
counts <- readCsvFst(countsFile)
counts_111_sas <- readCsvFst(counts_111_sas)
counts_111_jav <- readCsvFst(counts_111_jav)
metDatSas <- readCsvFst(metDatSasFile)
metDatJav <- readCsvFst(metDatJavFile)
metDat <- readCsvFst(metDatFile)

head(colnames(counts_111_jav))

# Merge 111 files
counts_111_jav <- counts_111_jav[match(rownames(counts_111_sas),
                                       rownames(counts_111_jav)), ]

counts_111 <- cbind.data.frame(counts_111_sas, counts_111_jav)
rownames(counts_111) <- gsub("\\..*", "", rownames(counts_111))

# Merge 111 metadata files
metDatSas$Intervention <- paste(metDatSas$Gene,
                                metDatSas$Method,
                                sep = "_")

colnames(metDatSas) <- gsub("Biosample.Name", "Celltype", colnames(metDatSas))

metDat_111 <- rbind.data.frame(metDatSas[, c("GSE.Series", "Celltype",
                                             "Intervention", "Experimental_GSM",
                                             "Control_GSM")],
                               metDatJav[, c("GSE.Series", "Celltype",
                                             "Intervention", "Experimental_GSM",
                                             "Control_GSM")])

# Filter cell types to keep only the ones we're interested in

# Cell types to be kept
cellTypes2KeepJav <- c("astrocyte",
                       "Cultured cortical neurons,",
                       "GABAergic neurons",
                       "hCMEC/D3",
                       "human brain",
                       "human neurons",
                       "monocyte-derived dendritic cells",
                       "monocytes-derived dendritic cells ",
                       "NSC ")

cellTypes2KeepSas <- c("human motor neurons",
                       "Motor Neuron",
                       "Neural precursor cells",
                       "neuroepithelial stem cells",
                       "Neuron",
                       "NHA",
                       "NPC",
                       "NPCs",
                       "NSC",
                       "NSCs",
                       "U5-NPC (G1 Phase)")

cellTypesKeep <- c(cellTypes2KeepJav,
                   cellTypes2KeepSas)
metDat_111 <- metDat_111[metDat_111$Celltype %in% cellTypesKeep, ]

# Iterate over GSE.Series to select each GSM from the table, in case it exists,
# and create new metadata file that has info for each gsm about gseSeries, cell
# type, perturbation and experimental group (if it's control or experimental).
# When a unique GSM maps to more than one column in counts_111 compute average,
# as it is due to different runs of same sample.
counts_111_filt <- data.frame(matrix(nrow = nrow(counts_111),
                                     ncol = 0,
                                     dimnames = list(rownames(counts_111), NULL)))

metDat_111_parsed <- data.frame(matrix(nrow = 0, ncol = 5,
                                       dimnames = list(NULL,
                                                       c("specimenID",
                                                         "GSE_series",
                                                         "cell_type",
                                                         "perturbation",
                                                         "exper_group"))))


#dupGSESeries <- metDat_111$GSE.Series[duplicated(metDat_111$GSE.Series)]

#metDat_111[metDat_111$GSE.Series == dupGSESeries[1], ]

for(i in seq_along(metDat_111$GSE.Series)){
        gse <- metDat_111$GSE.Series[i]
        xprs <- strsplit(metDat_111$Experimental_GSM[i], split = ";")[[1]]
        ctls <- strsplit(metDat_111$Control_GSM[i], split = ";")[[1]]
        cellType <- metDat_111$Celltype[i]
        pert <- metDat_111$Intervention[i]
        gsms_kept <- c()
        gsms <- c(xprs, ctls)
        expGroup <- c()
        for(j in seq_along(gsms)){
                #j <- 1
                gsm <- gsms[j]
                idx <- grep(gsm, colnames(counts_111))
                if(length(idx) > 0){
                        counts_111_gsm <- counts_111[, idx]
                        if(!is.null(dim(counts_111_gsm))){
                                counts_111_gsm <- apply(counts_111_gsm, 1, mean)
                        }
                        counts_111_gsm <- matrix(counts_111_gsm,
                                                 nrow = length(counts_111_gsm),
                                                 ncol = 1,
                                                 dimnames = list(names(counts_111_gsm),
                                                                 gsm))
                        counts_111_filt <- cbind.data.frame(counts_111_filt, counts_111_gsm)
                        gsms_kept <- c(gsms_kept, gsm)
                        if(gsm %in% xprs){
                                group <- "expr"
                        }else if(gsm %in% ctls){
                                group <- "ctrl"
                        }
                        expGroup <- c(expGroup, group)
                }
        }
        toBindMetDF <- data.frame(specimenID = gsms_kept,
                                  GSE_series = rep(gse, length(gsms_kept)),
                                  cell_type = rep(cellType, length(gsms_kept)),
                                  perturbation = rep(pert, length(gsms_kept)),
                                  exper_group = expGroup)
        metDat_111_parsed <- rbind.data.frame(metDat_111_parsed,
                                              toBindMetDF)
}

#dupGSMs <- unique(colnames(counts_111_filt)[duplicated(colnames(counts_111_filt))])
#dupGSM <- dupGSMs[1]
#isCtrlVec <- c()
#for(dupGSM in dupGSMs){
#        isCtrl <- any(grepl(dupGSM, metDat_111$Control_GSM))
#        print(sprintf("Is control: %s", isCtrl))
#        isCtrlVec <- c(isCtrlVec, isCtrl)
#        metDat_111[grepl(dupGSM,
#                         metDat_111$Experimental_GSM) | grepl(dupGSM, metDat_111$Control_GSM), ]
#}
#notCtrlGSMs <- dupGSMs[!isCtrlVec]



#dupGSM <- notCtrlGSMs[5]
#dupGSM
#metDat_111[grepl(dupGSM,
#                 metDat_111$Experimental_GSM) | grepl(dupGSM, metDat_111$Control_GSM), ]


#counts_111_filt[1:10, grep(dupGSM, colnames(counts_111_filt))]

# We are dealing with some duplicates due to same GSMs appearing in 
# several rows of the metadata. In Counts are identical columns, so let's
# just get rid of them

metDat_111_parsed$perturbation[metDat_111_parsed$exper_group == "ctrl"] <- "none"
metDat_111_parsed <- metDat_111_parsed[!duplicated(metDat_111_parsed$specimenID), ]
counts_111_filt <- counts_111_filt[, !duplicated(colnames(counts_111_filt))]


dim(metDat_111_parsed)
dim(counts_111_filt)
dim(counts_111)

# Merge with global counts file (ROSMAP etc)
counts <- t(counts)
commonGenes <- intersect(rownames(counts), rownames(counts_111_filt))

# Merge the datasets
counts_merged <- cbind.data.frame(counts[commonGenes, ],
                                  counts_111_filt[commonGenes, ])

# Create a common metadata file
head(metDat)
head(metDat_111_parsed)

cols2Add_metDat_111 <- c("platform", "RIN", "libraryPrep",
                         "libraryPreparationMethod",
                         "runType",
                         "readLength",
                         "individualID",
                         "exclude",
                         "excludeReason",
                         "sex",
                         "race",
                         "ageDeath",
                         "apoeGenotype",
                         "pmi",
                         "Braak",
                         "batch_lib",
                         "batch_seq",
                         "mmse30_first_ad_dx",
                         "organ",
                         "mmse30_lv")

metDat_111_parsed$substudy <- rep("perts_111",
                                  nrow(metDat_111_parsed))
metDat_111_parsed$assay <- rep("rnaSeq",
                               nrow(metDat_111_parsed))
metDat_111_parsed$diagn_4BrainClck <- rep("perturbations",
                                          nrow(metDat_111_parsed))

# Change GSE_series to batch_rna, as GSE series correspond to experiments
colnames(metDat_111_parsed) <- gsub("GSE_series", "batch_rna", colnames(metDat_111_parsed))

# Chance cel_type for tissue, for later unifying the category as tissue_cellType
# (once merged)
colnames(metDat_111_parsed) <- gsub("cell_type", "tissue", colnames(metDat_111_parsed))
metDat_111_parsed <- cbind.data.frame(metDat_111_parsed,
                                      data.frame(matrix(nrow = nrow(metDat_111_parsed),
                                                        ncol = length(cols2Add_metDat_111),
                                                        dimnames = list(rownames(metDat_111_parsed),
                                                                        cols2Add_metDat_111))))

metDat_111_parsed$RIN <- max(metDat$RIN[!is.na(metDat$RIN)])

cols2Add_metDat_orig <- colnames(metDat_111_parsed)[!colnames(metDat_111_parsed) %in% colnames(metDat)]
colnames(metDat)[!colnames(metDat) %in% colnames(metDat_111_parsed)]

metDat$perturbation <- rep("none", nrow(metDat))
metDat$exper_group <- rep(NA, nrow(metDat))

metDat_merged <- rbind.data.frame(metDat,
                                  metDat_111_parsed[, colnames(metDat)])

writeCsvFst(counts_merged, file = outCountsFile)
writeCsvFst(metDat_merged, file = outMetDatFile)

