################################################################################
# Brain clock: Preprocess counts of integrated dataset but LINCS, doing        #
# log2(counts + 1) and quantile normalization for matching format of LINCS.    #
# Remove lowly expressed genes both in general and by age ranks.               #
################################################################################

countsFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_counts_wTBI_wPert111.csv"
lincsGeneInfo <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/LINCS/parsed/lincs_geneInfo.csv"
metDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat_wTBI_wPert111.csv"
propZeroRem <- .8
rinFilt <- 6
outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/"

outName <- sprintf("%s%s_log2_quantNorm_preproc.csv", outDir,
                   gsub(".csv", "", basename(countsFile)))

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

# Function for performing quantile normalization
quantNorm <- function(m, axis = 2){
        m_rank <- apply(m, axis, function(x) floor(rank(x)))
        m_sort <- apply(m, axis, sort)
        means <- apply(m_sort,
                       1,
                       mean)
        names(means) <- rank(means)
        m_norm <- apply(m_rank,
                        axis,
                        function(x) means[match(as.character(x),
                                                names(means))])
        dimnames(m_norm) <- dimnames(m)
        return(m_norm)
}

quantNorm <- function(m, axis = 2){
        m_rank <- apply(m, axis, rank, ties.method = "first")
        m_sort <- apply(m, axis, sort)
        means <- apply(m_sort,
                       1,
                       mean)
        names(means) <- rank(means)
        m_norm <- apply(m_rank,
                        axis,
                        function(x) means[match(as.character(x),
                                                names(means))])
        dimnames(m_norm) <- dimnames(m)
        return(m_norm)
}

quantNorm <- function(m, axis = 2) {
        if (axis == 1) {
                m <- t(m)
        }
        
        # Get ranks of the matrix
        m_rank <- apply(m, 2, rank, ties.method = "average")
        
        # Sort the matrix
        m_sort <- apply(m, 2, sort)
        
        # Calculate row means of the sorted matrix
        means <- rowMeans(m_sort)
        
        # Create a normalized matrix with the same dimensions
        m_norm <- matrix(0, nrow = nrow(m), ncol = ncol(m))
        
        for (i in 1:ncol(m)) {
                m_norm[, i] <- means[rank(m[, i], ties.method = "average")]
        }
        
        if (axis == 1) {
                m_norm <- t(m_norm)
        }
        
        dimnames(m_norm) <- dimnames(m)
        return(m_norm)
}

# Check proportion of zeros
printZeroProp <- function(df){
        zeroProp <- sum(df == 0)/(nrow(df)*ncol(df))
        
        print(sprintf("Proportion of zeros: %s",
                      as.character(round(zeroProp,
                                         digits = 4))))
}

# Finds genes that are expressed at less than an average of 1 per sample 
# in each one of the ages
findLowExpGenes_perAge <- function(df, metDat){
        lowExpGenes_perAge <- list()
        ageVec <- unique(metDat$ageDeath)
        ageVec <- ageVec[!is.na(ageVec)]
        for(age in ageVec){
                sampsAge <- metDat$specimenID[metDat$ageDeath == age]
                sampsAge <- sampsAge[!is.na(sampsAge)]
                sampsAge <- make.names(sampsAge)
                sampsAge <- sampsAge[sampsAge %in% rownames(df)]
                df_age <- df[make.names(sampsAge), ]
                #make.names(sampsAge) %in% row.names(comb_TPMs)
                #any(is.na(tpms_age))
                #tpms_age[1:100, 1:10]
                #if(length(sampsAge) == 1){
                #        df_age <- matrix(df_age, ncol = 1, nrow = length(df_age),
                #                         dimnames = list(names(df_age),
                #                                         sampsAge))
                #        df_age <- data.frame(df_age)
                #}
                if(length(sampsAge) > 1){
                        lowExpGenes <- colnames(df_age)[colSums(df_age) < nrow(df_age)]
                        lowExpGenes_perAge[[paste("age",
                                                  as.character(age),
                                                  sep = "_")]] <- lowExpGenes
                }
        }
        lowExpGenes_perAge <- lowExpGenes_perAge[lapply(lowExpGenes_perAge,
                                                        length) != 0]
        return(lowExpGenes_perAge)
}

# Look for the intersection across all elements in a list.
recurs_intersect <- function(lst){
        if (length(lst) == 1) {
                return(lst[[1]])
        }
        
        # Recursive case: find the intersection of the first element with the
        # intersection of the rest of the list
        return(intersect(lst[[1]], recurs_intersect(lst[-1])))
}

# Identifies rows or columns with a high proportion of zeros
idHighPropZero <- function(df, axis = 2, propZeroThrshld){
        highPropZero <- dimnames(df)[[axis]][apply(df,
                                                   axis,
                                                   function(x) sum(x == 0)/length(x)) > propZeroThrshld]
        return(highPropZero)
}

# Removes rows and columns with a high proportion of zeros.
remHighPropZeros <- function(df, propZeroThrshld){
        cols2Rem <- idHighPropZero(df,
                                   axis = 2,
                                   propZeroThrshld = propZeroThrshld)
        rows2Rem <- idHighPropZero(df,
                                   axis = 1,
                                   propZeroThrshld = propZeroThrshld)
        df_filt <- df[!rownames(df) %in% rows2Rem,
                      !colnames(df) %in% cols2Rem]
        return(df_filt)
}

# Load data
################################################################################
counts <- readCsvFst(countsFile)
lincsGeneInfo <- readCsvFst(lincsGeneInfo)
metDat <- readCsvFst(metDatFile)

# Preprocess data
################################################################################

# Keep only genes in LINCS, to make quantil normalization as close as possible
commGenes <- intersect(lincsGeneInfo$ensembl_gene_id, rownames(counts))

# Log transform and quantile normalize
counts <- counts[commGenes, ]
counts_log2 <- log2(counts + 1)
counts_log2_qNorm <- quantNorm(counts_log2)
plot(density(as.matrix(counts_log2_qNorm)))

counts_log2_qNorm <- t(counts_log2_qNorm)

# Remove lowly expressed genes
printZeroProp(counts_log2_qNorm)

# Identify genes that are lowly expressed in each age. (less than average)
# of one count per sample
lowExpGenes_age <- findLowExpGenes_perAge(counts_log2_qNorm, metDat = metDat)
lowExpGenes_common <- recurs_intersect(lowExpGenes_age)

# Remove the genes in the intersection of all ages from the dataset
counts_log2_qNorm <- counts_log2_qNorm[, !colnames(counts_log2_qNorm) %in% lowExpGenes_common]
plot(density(as.matrix(counts_log2_qNorm)))

# Check variables and samples that have high proportion of zeros
counts_log2_qNorm <- remHighPropZeros(counts_log2_qNorm,
                                      propZeroThrshld = propZeroRem)

# Remove samples with a RIN below of the user-specified threshold
# Filter out samples with RIN < 6
sampsGoodRIN <- make.names(metDat$specimenID[metDat$RIN >= rinFilt])
counts_log2_qNorm <- counts_log2_qNorm[rownames(counts_log2_qNorm) %in% sampsGoodRIN, ]
plot(density(as.matrix(counts_log2_qNorm)))
mean(counts_log2_qNorm)
sd(counts_log2_qNorm)
writeCsvFst(counts_log2_qNorm, file = outName)
print(sprintf("%s saved at %s", basename(outName), dirname(outName)))