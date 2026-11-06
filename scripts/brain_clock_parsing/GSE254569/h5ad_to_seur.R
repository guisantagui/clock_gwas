if (!require("anndata", quietly = T)){
    install.packages("anndata",
                     repos='http://cran.us.r-project.org')
}
library(anndata)
library(Seurat)
h5ad <- "/work/projects/age_sorter/data/counts/GSE254569/GSE254569/GSE254569_adata_RNA.h5ad"
seur <- "/work/projects/age_sorter/data/counts/GSE254569/GSE254569/GSE254569_adata_RNA.rds"

sceasy::convertFormat(h5ad, from="anndata", to="seurat", outFile= seur)

data <- read_h5ad(h5ad)
counts_matrix <- reticulate::py_to_r(data$X)
obs_df <- reticulate::py_to_r(data$obs)
data <- CreateSeuratObject(counts = t(counts_matrix), meta.data = obs_df,min.features = 700, min.cells = 3)
#data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 700, min.cells = 3)
saveRDS(data, sprintf("%s/GSE254569_seur.rds"))

write_table_fast <- function (df, f, row.names = T, col.names = T, sep = ",", sep2 = c("", 
    "|", ""), quote = "auto") 
{
    if (row.names) {
        rn <- rownames(df)
        df <- data.table(df)
        df[, `:=`(V1, rn)]
        setcolorder(df, c("V1", setdiff(names(df), "V1")))
    }
    else {
        df <- data.table(df)
    }
    fwrite(df, f, col.names = col.names, sep = sep, sep2 = sep2, 
        quote = quote)
}