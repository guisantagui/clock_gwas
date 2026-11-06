################################################################################
# Brain clock: intersects merged dataset with LINCS genes, remove lowly        #
# expressed genes and samples with RIN lower than 6.                           #
################################################################################
if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}
if(!require("plotUtils", quietly = T)){
        devtools::install_github('guisantagui/plotUtils', upgrade = "never")
}
library(plotUtils)

# Functions
################################################################################

findLowExpGenes_perAgeBin <- function(df, metDat, bin_width = 5){
    lowExpGenes_perBin <- list()
    
    # Bin ages into 5-year intervals (e.g., 60–64, 65–69, etc.)
    metDat$ageBin <- cut(metDat$ageDeath,
                         breaks = seq(floor(min(metDat$ageDeath, na.rm = TRUE)),
                                      ceiling(max(metDat$ageDeath, na.rm = TRUE)),
                                      by = bin_width),
                         include.lowest = TRUE, right = FALSE)

    ageBins <- unique(metDat$ageBin)
    ageBins <- ageBins[!is.na(ageBins)]

    for (bin in ageBins) {
        sampsBin <- metDat$specimenID[metDat$ageBin == bin]
        sampsBin <- sampsBin[!is.na(sampsBin)]
        sampsBin <- make.names(sampsBin)
        sampsBin <- sampsBin[sampsBin %in% rownames(df)]
        df_bin <- df[make.names(sampsBin), ]
        
        if (length(sampsBin) > 1) {
            lowExpGenes <- colnames(df_bin)[colSums(df_bin) < nrow(df_bin)]
            lowExpGenes_perBin[[as.character(bin)]] <- lowExpGenes
        }
    }

    # Remove empty entries
    lowExpGenes_perBin <- lowExpGenes_perBin[lapply(lowExpGenes_perBin, length) != 0]
    return(lowExpGenes_perBin)
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

# Directory stuff
################################################################################
counts_f <- "../../results/parsed/merged/merged_counts.csv"
lincsGeneInfo <- "../../data/perturbation/lincs/lincs_genes.csv"
overwrite_merged <- T
metdat_f <- "../../results/parsed/merged/merged_metdat.csv"

propZeroRem <- .8
rinFilt <- 6
outDir <- "../../results/parsed/merged/"

create_dir_if_not(outDir)

if (overwrite_merged){
        outNameCounts <- counts_f
        outNameMetdat <- metdat_f
}else{
        outNameCounts <- sprintf("%s%s_filt.csv", outDir,
                                 gsub(".csv", "", basename(counts_f)))
        outNameMetdat <- sprintf("%s%s_filt.csv", outDir,
                                 gsub(".csv", "", basename(metdat_f)))
}


# Load data
################################################################################

counts <- read_table_fast(counts_f, row.names = 1)
lincs_gi <- read.csv(lincsGeneInfo, row.names = 1)
metdat <- read_table_fast(metdat_f, row.names = 1)

# Filtering
################################################################################

# Keep only genes in LINCS
commGenes <- intersect(lincs_gi$gene, colnames(counts))
counts <- counts[, commGenes]
# Ensure that only what is available in the metadata is in counts
counts <- counts[make.names(rownames(counts)) %in% make.names(metdat$specimenID), ]
# Remove genes with a high proportion of zeros:
counts <- counts[, apply(counts,
                         2, function(x) sum(x == 0)/length(x)) < propZeroRem]

# Find genes that at each 5 year interval age, express less than an average of one count per
# gene
lowExp_genes_perAgeBin <- findLowExpGenes_perAgeBin(counts, metdat)
to_rem <- recurs_intersect(lowExp_genes_perAgeBin)
counts <- counts[, !colnames(counts) %in% to_rem]

# Filter out samples with low RIN
counts <- counts[make.names(rownames(counts)) %in% make.names(metdat$specimenID[!is.na(metdat$RIN) & metdat$RIN >= rinFilt]), ]

# keep in the metadata only the samples that were retained in the counts
metdat <- metdat[match(make.names(rownames(counts)),
                       make.names(metdat$specimenID)), ]

# Save
################################################################################

# Genes in the final dataset
all_genes <- colnames(counts)
all_genes <- data.frame(ensembl_gene_id = all_genes)
write_table_fast(all_genes, f = sprintf("%sall_genes.csv", outDir))

# Processed counts file
write_table_fast(counts, outNameCounts)

# Filtered metadata file
write_table_fast(metdat,
                 f = outNameMetdat)