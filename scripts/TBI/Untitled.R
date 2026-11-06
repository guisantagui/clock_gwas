library(data.table)

hgnc_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/metabognn_dev/data/utilities/hgnc_complete_set.txt"
hgnc <- read.csv(hgnc_f, sep = '\t')
tbi_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/tbi_raw"

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

countMat <- Reduce(
        function(x, y) merge(x, y, by = "gene_id", all = FALSE),
        countMat
)
