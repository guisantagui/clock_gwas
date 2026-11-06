import anndata as ad
import pandas as pd
import scipy.io
import os
import numpy as np

input_h5ad = "/work/projects/age_sorter/data/counts/GSE254569/GSE254569/GSE254569_adata_RNA.h5ad"

# Load the h5ad file
adata = ad.read_h5ad(input_h5ad)

# Define output directory (inside the input file's directory)
output_dir = os.path.join(os.path.dirname(input_h5ad), "extracted_data")
os.makedirs(output_dir, exist_ok=True)  # Create directory if it doesn’t exist

# Extract raw counts (keep sparse)
if "counts" in adata.layers:
    count_matrix = adata.layers["counts"]
else:
    count_matrix = adata.X

# Ensure matrix is sparse
if not scipy.sparse.issparse(count_matrix):
    count_matrix = scipy.sparse.csr_matrix(count_matrix)

# Save sparse count matrix in MTX format (for Seurat)
scipy.io.mmwrite(os.path.join(output_dir, "matrix.mtx"), count_matrix)

# Save cell metadata (obs)
adata.obs.to_csv(os.path.join(output_dir, "metadata.csv"))

# Save gene metadata (var)
adata.var.to_csv(os.path.join(output_dir, "features.csv"))

# Save cell barcodes separately for Seurat
adata.obs.index.to_series().to_csv(os.path.join(output_dir, "barcodes.tsv"), sep="\t", index=False, header=False)

# Save gene names separately for Seurat
adata.var.index.to_series().to_csv(os.path.join(output_dir, "features.tsv"), sep="\t", index=False, header=False)

print(f"✅ Files saved in: {output_dir}")