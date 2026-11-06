if(!require("biomaRt", quietly = T)) BiocManager::install("biomaRt", update = F)
library(biomaRt)
################################################################################
# Brain clock:  Merge counts data from databases with the perturbation files   #
# Javier and Sascha generated, and with sc data from ageAnno.                  #
################################################################################

countsFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_counts_wTBI.csv"
counts_111_sas <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/111_sascha_complete_counts_df.csv"
counts_111_jav <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/javier_111_counts.csv"
metDatSasFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/Sascha_metadata.csv"
metDatJavFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/javier_metadata.csv"
metDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat_wTBI.csv"
metDatSCFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/aging/data/gkac847_supplemental_files/Supplementary Tables.xlsx"
scDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/cell_type_matrices"
scFiles <- list.files(scDir, full.names = T)
outCountsFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_counts_wTBI_wPert111_wSC.csv"
outMetDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat_wTBI_wPert111_wSC.csv"


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

# Get gene info given one of several identifiers
getENSEMBL_ID_info <- function(ensembl_list,
                               mart,
                               removeNAs = T,
                               filt = "ensembl_gene_id_version",
                               remapUnmap = T){
        #ensembl_list <- gsub("X", "", colnames(countMat))
        #mart <- human
        #filt <- "entrezgene_id"
        
        
        print(sprintf("Mapping the %s...", filt))
        if(filt == "ensembl_gene_id_version"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "ensembl_gene_id_version",
                          "entrezgene_id",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }else if(filt == "ensembl_gene_id" | filt == "entrezgene_id" | filt == "hgnc_symbol"){
                atts <- c("hgnc_symbol",
                          "ensembl_gene_id",
                          "entrezgene_id",
                          "chromosome_name",
                          "start_position",
                          "end_position")
        }
        
        gene_coords <- getBM(attributes = atts,
                             filters = filt,
                             values = ensembl_list,
                             mart = mart)
        
        gene_coords$size <- gene_coords$end_position - gene_coords$start_position
        if(filt == "ensembl_gene_id_version"){
                mapped_matchVec <- match(ensembl_list,
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
                mapped_matchVec <- match(ensembl_list,
                                         gene_coords$ensembl_gene_id)
                outDF <- data.frame(ensembl_gene_id_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                naVec <- is.na(outDF$ensembl_gene_id)
        }else if(filt == "entrezgene_id"){
                mapped_matchVec <- match(ensembl_list,
                                         gene_coords$entrezgene_id)
                outDF <- data.frame(entrezgene_id_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    entrezgene_id = gene_coords$entrezgene_id[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                naVec <- is.na(outDF$ensembl_gene_id)
                unMapped <- outDF$entrezgene_id_orig[is.na(outDF$entrezgene_id)]
                # Try out with org.Hs.eg.db for the unmapped ones
                unMapped_ensembl <- mapIds(org.Hs.eg.db,
                                           keys = unMapped,
                                           keytype = "ENTREZID",
                                           column = c("ENSEMBL"),
                                           fuzzy = T)
                unMapped_ensembl_noNAs <- unMapped_ensembl[!is.na(unMapped_ensembl)]
                gene_coords_unMap <- getBM(attributes = atts,
                                           filters = "ensembl_gene_id",
                                           values = unMapped_ensembl_noNAs,
                                           mart = mart)
                gene_coords_unMap$size <- gene_coords_unMap$end_position - gene_coords_unMap$start_position
                gene_coords_unMap$entrezgene_id_orig <- names(unMapped_ensembl_noNAs)[match(gene_coords_unMap$ensembl_gene_id, unMapped_ensembl_noNAs)]
                outDF[match(gene_coords_unMap$entrezgene_id_orig,
                            outDF$entrezgene_id_orig),
                      c("ensembl_gene_id",
                        "symbol",
                        "entrezgene_id",
                        "chromosome",
                        "start", "end", "size")] <- gene_coords_unMap[, c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "chromosome_name",
                                                                          "start_position", "end_position", "size")]
                # Repeat mapping of unMapped but this time with symbols
                unMapped <- outDF$entrezgene_id_orig[is.na(outDF$entrezgene_id)]
                unMapped_symbols <- mapIds(org.Hs.eg.db,
                                           keys = unMapped,
                                           keytype = "ENTREZID",
                                           column = c("SYMBOL"),
                                           fuzzy = T)
                unMapped_symbols_noNAs <- unMapped_symbols[!is.na(unMapped_symbols)]
                gene_coords_unMap <- getBM(attributes = atts,
                                           filters = "hgnc_symbol",
                                           values = unMapped_symbols_noNAs,
                                           mart = mart)
                gene_coords_unMap$size <- gene_coords_unMap$end_position - gene_coords_unMap$start_position
                gene_coords_unMap$entrezgene_id_orig <- names(unMapped_symbols_noNAs)[match(gene_coords_unMap$hgnc_symbol, unMapped_symbols_noNAs)]
                outDF[match(gene_coords_unMap$entrezgene_id_orig,
                            outDF$entrezgene_id_orig),
                      c("ensembl_gene_id",
                        "symbol",
                        "entrezgene_id",
                        "chromosome",
                        "start", "end", "size")] <- gene_coords_unMap[, c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "chromosome_name",
                                                                          "start_position", "end_position", "size")]
                
                
                # pass the original entrez ids to the entrezgene_id column if they don't have id there
                boolCompEntrezIDs <- !is.na(outDF$entrezgene_id_orig) & is.na(outDF$entrezgene_id)
                outDF$entrezgene_id[boolCompEntrezIDs] <- outDF$entrezgene_id_orig[boolCompEntrezIDs]
                naVec <- is.na(outDF$ensembl_gene_id)
                
        }else if(filt == "hgnc_symbol"){
                mapped_matchVec <- match(ensembl_list,
                                         gene_coords$hgnc_symbol)
                outDF <- data.frame(symbol_orig = ensembl_list,
                                    ensembl_gene_id = gene_coords$ensembl_gene_id[mapped_matchVec],
                                    symbol = gene_coords$hgnc_symbol[mapped_matchVec],
                                    entrezgene_id = gene_coords$entrezgene_id[mapped_matchVec],
                                    chromosome = gene_coords$chromosome_name[mapped_matchVec],
                                    start = gene_coords$start_position[mapped_matchVec],
                                    end = gene_coords$end_position[mapped_matchVec],
                                    size = gene_coords$size[mapped_matchVec])
                # If selected, another round of mapping is done on the genes 
                # that are unmapped using org.Hs.eg.db, to maximize number of
                # maps.
                if(remapUnmap){
                        #naVec <- is.na(outDF$ensembl_gene_id)
                        unMapped <- outDF$symbol_orig[is.na(outDF$entrezgene_id)]
                        unMapped_ensembl <- mapIds(org.Hs.eg.db,
                                                   keys = unMapped,
                                                   keytype = "SYMBOL",
                                                   column = c("ENSEMBL"),
                                                   fuzzy = T)
                        unMapped_ensembl_noNAs <- unMapped_ensembl[!is.na(unMapped_ensembl)]
                        gene_coords_unMap <- getBM(attributes = atts,
                                                   filters = "ensembl_gene_id",
                                                   values = unMapped_ensembl_noNAs,
                                                   mart = mart)
                        gene_coords_unMap$size <- gene_coords_unMap$end_position - gene_coords_unMap$start_position
                        gene_coords_unMap$symbol_orig <- names(unMapped_ensembl_noNAs)[match(gene_coords_unMap$ensembl_gene_id, unMapped_ensembl_noNAs)]
                        
                        outDF[match(gene_coords_unMap$symbol_orig,
                                    outDF$symbol_orig),
                              c("ensembl_gene_id",
                                "symbol",
                                "entrezgene_id",
                                "chromosome",
                                "start", "end", "size")] <- gene_coords_unMap[, c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "chromosome_name",
                                                                                  "start_position", "end_position", "size")]
                }
                naVec <- is.na(outDF$size)
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

# This function adds up columns or rows that ended up with same name 
# after transformation to EMSEMBL.
sumDupsMat <- function(m, dups_in){
        if(dups_in == "cols"){
                m <- t(m)
        }
        featsOrd <- unique(rownames(m))
        dups <- unique(rownames(m)[duplicated(rownames(m))])
        
        for(i in seq_along(dups)){
                dup <- dups[i]
                dupMat <- m[rownames(m) == dup, ]
                sumsMat <- matrix(colSums(dupMat),
                                  nrow = 1,
                                  ncol = ncol(dupMat),
                                  dimnames = list(dup,
                                                  colnames(dupMat)))
                m <- m[rownames(m) != dup, ]
                m <- rbind(m, sumsMat)
        }
        m <- m[match(featsOrd, rownames(m)), ]
        if(dups_in == "cols"){
                m <- t(m)
        }
        return(m)
}

# Load data
################################################################################
counts <- readCsvFst(countsFile)
counts_111_sas <- readCsvFst(counts_111_sas)
counts_111_jav <- readCsvFst(counts_111_jav)
metDatSas <- readCsvFst(metDatSasFile)
metDatJav <- readCsvFst(metDatJavFile)
metDat <- readCsvFst(metDatFile)
# Parse SC metadata
metDatSC <- as.data.frame(readxl::read_xlsx(metDatSCFile, sheet = 1))
metDatSC <- metDatSC[(which(metDatSC$`Table S1. Age information of the sample in AgeAnno` == "scRNA") + 1):(which(metDatSC$`Table S1. Age information of the sample in AgeAnno` == "scATAC") - 1), ]
colnames(metDatSC) <- unlist(as.vector(metDatSC[1, ]))
metDatSC <- metDatSC[2:(nrow(metDatSC) - 1), ]
metDatSC <- metDatSC[grepl("brain",
                           metDatSC$sample), colnames(metDatSC) != "organs"]
rownames(metDatSC) <- 1:nrow(metDatSC)

# Load SC data and merge them in a single matrix
counts_sc <- data.frame()
for(f in scFiles){
        df <- readCsvFst(f)
        if(nrow(counts_sc) == 0){
                counts_sc <- df
        }else{
                counts_sc <- cbind.data.frame(counts_sc,
                                              df[match(rownames(counts_sc),
                                                       rownames(df)), ])
        }
}

# The numbered mid1, mid2, mid3 etc of the counts_sc DF correspond to the
# metdat samples in the same order, so let's match them
metDatSC$samp <- tolower(metDatSC$`age group`)
for(grp in unique(metDatSC$samp)){
        metDatSC$samp[metDatSC$samp == grp] <- paste0(grp,
                                                      1:sum(metDatSC$samp == grp))
}
metDatSC$samp == unique(metDatSC$samp)[1]

metDatSC_pars <- data.frame(specimenID = colnames(counts_sc),
                            individualID = gsub(".*_", "", colnames(counts_sc)),
                            organ = rep("brain", ncol(counts_sc)),
                            tissue = gsub("\\_.*", "", colnames(counts_sc)),
                            assay = rep("scRNAseq", ncol(counts_sc)),
                            substudy = rep("ageAnno", ncol(counts_sc)),
                            ageDeath = metDatSC$`age (year)`[match(gsub(".*_",
                                                                        "",
                                                                        colnames(counts_sc)),
                                                                   metDatSC$samp)])
# Remove the + from the age column and convert to numeric
metDatSC_pars$ageDeath <- as.numeric(gsub("+", "", metDatSC_pars$ageDeath, fixed = T))


dim(counts_sc)
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


# Convert gene symbols to ENSEMBL ids in the SC df
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

symbs <- rownames(counts_sc)
geneInfoSC <- getENSEMBL_ID_info(symbs, human, removeNAs = T,
                                 filt = "hgnc_symbol",
                                 remapUnmap = F)

newNames <- geneInfoSC$ensembl_gene_id[match(rownames(counts_sc),
                                             geneInfoSC$symbol_orig)]

counts_sc <- counts_sc[!is.na(newNames), ]
counts_sc <- as.matrix(counts_sc)
rownames(counts_sc) <- newNames[!is.na(newNames)]

# Sum duplicated rows after conversion
counts_sc <- sumDupsMat(counts_sc, "rows")

# Merge with global counts file (ROSMAP etc)
counts <- t(counts)
commonGenes <- intersect(intersect(rownames(counts),
                                   rownames(counts_111_filt)),
                         rownames(counts_sc))

# Merge the datasets
counts_merged <- cbind.data.frame(counts[commonGenes, ],
                                  counts_111_filt[commonGenes, ],
                                  counts_sc[commonGenes, ])

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

# Add the sc metadata to the merged metadata
metDatSC_pars$diagn_4BrainClck <- "scRNAseq"
metDatSC_pars$RIN <- 10

cols2Add <- colnames(metDat_merged)[!colnames(metDat_merged) %in% colnames(metDatSC_pars)]

metDatSC_pars <- cbind.data.frame(metDatSC_pars,
                                  data.frame(matrix(nrow = nrow(metDatSC_pars),
                                                    ncol = length(cols2Add),
                                                    dimnames = list(rownames(metDatSC_pars),
                                                                    cols2Add))))

metDat_merged <- rbind.data.frame(metDat_merged,
                                  metDatSC_pars[, colnames(metDat_merged)])

writeCsvFst(counts_merged, file = outCountsFile)
writeCsvFst(metDat_merged, file = outMetDatFile)

