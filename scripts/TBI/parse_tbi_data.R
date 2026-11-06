library(data.table)

hgnc_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/metabognn_dev/data/utilities/hgnc_complete_set.txt"
hgnc <- read.csv(hgnc_f, sep = '\t')
tbi_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/tbi_raw"
tbi_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/TBI_metadata_unified.csv"
tbi_data_files <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/tbi_data_files.csv")
tbi_mdat_unparsed <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/DonorInformation.csv")


tbi_mdat <- read.csv(tbi_mdat_f, row.names = 1)

tbi_files <- list.files(tbi_dir, full.names = TRUE)

tbi_files <- tbi_files[grep(".txt", tbi_files, fixed = TRUE)]

countMat <- list()
for(i in seq_along(tbi_files)){
        print(sprintf("%i/%s", i, length(tbi_files)))
        tbi_f <- tbi_files[i]
        df <- data.frame(fread(tbi_f))
        df <- df[, c("gene_id", "expected_count")]
        samp_name <- gsub(".txt", "", basename(tbi_f))
        countMat[[samp_name]] <- df
}
sum(countMat$H15.09.110.TCx.01$gene_id %in% hgnc$entrez_id)

ids_per_file <- lapply(countMat, `[[`, "gene_id")   # run before the Reduce
table(lengths(ids_per_file))                        # should be one value, ~50K
ref <- ids_per_file[[1]]
odd <- names(ids_per_file)[!sapply(ids_per_file, function(x) setequal(x, ref))]
odd                                                 # files with a different gene set
sum(duplicated(ref))                                # duplicated IDs inflate/garble merge

# "H14.09.081.FWM.01" is shorter. Examine the file
data.frame(fread(tbi_files[grep("H14.09.081.FWM.01", tbi_files)]))
# Its truncated. Download again
tbi_data_files[tbi_data_files$donor_name == "H14.09.081", ]
# Its the 2nd one, as theres a typo in the specimen names: its labeled as PCx
# (parietal neocortex) in specimen name, but both in structure acronym and
# structure name says FWM.
baseURL <- "https://aging.brain-map.org"
truncated_url <- tbi_data_files[tbi_data_files$donor_name == "H14.09.081", ]$gene_level_fpkm_file_link[2]
trunc_filename <- paste(tbi_dir, "H14.09.081.FWM.01_2.txt", sep = "/")

comm <- sprintf("curl -o %s %s%s", trunc_filename, baseURL, truncated_url)
#system(comm)


countMat <- Reduce(
        function(x, y) merge(x, y, by = "gene_id", all = FALSE),
        lapply(names(countMat), function(n) {
                df <- countMat[[n]]
                colnames(df)[2] <- n
                df
        })
)

rownames(countMat) <- countMat$gene_id
countMat$gene_id <- NULL

countMat <- as.matrix(countMat)
countMat[1:10, 1:10]



nu_rownames <- hgnc$ensembl_gene_id[match(rownames(countMat), hgnc$entrez_id)]
countMat <- countMat[!is.na(nu_rownames), ]
rownames(countMat) <- nu_rownames[!is.na(nu_rownames)]

countMat <- rowsum(countMat, group = rownames(countMat))

# Keep only the samples in tbi_mdat, which was already filtered in the brain
# clock to remove ambiguities in age

countMat <- countMat[, match(tbi_mdat$specimenID, colnames(countMat))]
dim(countMat)

count_outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/counts_parsed"
mdat_outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/metadata_parsed"

write.csv(countMat, file.path(count_outDir, "tbi_counts.csv"))
write.csv(tbi_mdat, file.path(mdat_outDir, "tbi_mdat.csv"))
