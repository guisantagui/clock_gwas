################################################################################
# Brain Clock: process, integrate, annotate and obtain the pseudo-bulk counts  # 
# of AgeAnno brain scRNAseq data                                               #
################################################################################

if(!require("remotes", quietly = T)) install.packages("remotes")
if(!require("BiocManager", quietly = T)) install.packages("BiocManager")
if(!require("Seurat", quietly = T)) remotes::install_github("satijalab/seurat",
                                                            "seurat5",
                                                            quiet = TRUE,
                                                            upgrade = "never")
library(Seurat)
if(!require("plyr", quietly = T)){
        install.packages("plyr",
                         repos='http://cran.us.r-project.org')
}
library(plyr)
library(ggplot2)
if(!require("ggExtra", quietly = T)){
        install.packages("ggExtra",
                         repos='http://cran.us.r-project.org')
}
library(ggExtra)
if(!require("cowplot", quietly = T)){
        install.packages("cowplot",
                         repos='http://cran.us.r-project.org')
}
library(cowplot)
if(!require("biomaRt", quietly = T)){
        install.packages("biomaRt",
                         repos='http://cran.us.r-project.org')
}
library(biomaRt)
if(!require(dplyr)) install.packages("dplyr",
                                     repos = 'http://cran.us.r-project.org')
library(dplyr)
if(!require("DoubletFinder", quietly = T)){
        remotes::install_github("chris-mcginnis-ucsf/DoubletFinder",
                                quiet = T,
                                upgrade = "never")
}
library(DoubletFinder)
if(!require("monocle3", quietly = T)){
        remotes::install_github("cole-trapnell-lab/monocle3",
                                quiet = T,
                                upgrade = "never")
}
library(monocle3)
if(!require("SeuratWrappers", quietly = T)){
        devtools::install_github('satijalab/seurat-wrappers', upgrade = "never")
}
library(SeuratWrappers)
if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}
if(!require("plotUtils", quietly = T)){
        devtools::install_github('guisantagui/plotUtils', upgrade = "never")
}
library(plotUtils)
library(tidyr)

# Functions
################################################################################

filterCells <- function(seur,mad.coeff = 3,pass = -1,org = "HUMAN",plot.path = "./"){
        seur@meta.data$Barcodes <- rownames(seur@meta.data)
        #Calculate percent.mito and percent.ribo metadata columns if they are not there
        if(!any(colnames(seur@meta.data) == "percent.mito")){
                if(org == "HUMAN"){
                        seur[["percent.mito"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "^MT-")
                } else if(org == "MOUSE"){
                        seur[["percent.mito"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "^mt-")
                } else {
                        stop("The specified organism is not supported")
                }
        }
        if(!any(colnames(seur@meta.data) == "percent.ribo")){
                if(org == "HUMAN"){
                        seur[["percent.ribo"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "(^RPL|^RPS|^MRP)")
                } else if(org == "MOUSE"){
                        seur[["percent.ribo"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "(^Rpl|^Rps|^Mrp)")
                } else {
                        stop("The specified organism is not supported")
                }
        }
        # Filtering cells based on percentage of mitochondrial transcripts
        cell.QC.stat <- seur@meta.data
        
        max.mito.thr <- median(cell.QC.stat$percent.mito) + mad.coeff*mad(cell.QC.stat$percent.mito)
        min.mito.thr <- median(cell.QC.stat$percent.mito) - mad.coeff*mad(cell.QC.stat$percent.mito)
        
        p1 <- ggplot(cell.QC.stat, aes(x=nFeature_RNA, y=percent.mito)) +
                geom_point() +
                geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
                geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
                annotate(geom = "text",
                         label = paste0(as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[2]),
                                        " cells removed\n",
                                        as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[1]),
                                        " cells remain"),
                         x = 6000,
                         y = -10)
        
        p <- ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100) 
        ggsave(paste0(plot.path,"Mitofilter_Marginal_Pass",pass,".png"),plot = p)
        
        cell.QC.stat <- cell.QC.stat %>%
                dplyr::filter(percent.mito <= max.mito.thr) %>% 
                dplyr::filter(percent.mito >= min.mito.thr)
        
        # Filtering cells based on number of genes and transcripts detected
        # Set low and hight thresholds on the number of detected genes
        min.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) - mad.coeff*mad(log10(cell.QC.stat$nFeature_RNA))
        max.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) + mad.coeff*mad(log10(cell.QC.stat$nFeature_RNA))
        
        # Set hight threshold on the number of transcripts
        max.count.thr <- median(log10(cell.QC.stat$nCount_RNA)) + mad.coeff*mad(log10(cell.QC.stat$nCount_RNA))
        
        p1 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point() +
                geom_smooth(method="lm") +
                geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
                geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
                geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2)
        
        p <- ggMarginal(p1, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountFilter_1_Pass",pass,".png"),plot = p)
        
        # Filter cells base on both metrics
        cell.QC.stat <- cell.QC.stat %>% 
                dplyr::filter(log10(nFeature_RNA) > min.features.thr) %>%
                dplyr::filter(log10(nFeature_RNA) < max.features.thr) %>%
                dplyr::filter(log10(nCount_RNA) < max.count.thr)
        
        lm.model <- lm(data = cell.QC.stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
        
        p2 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point() +
                geom_smooth(method="lm") +
                geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
                geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
                geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2) +
                geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
                annotate(geom = "text", label = paste0(dim(cell.QC.stat)[1], " QC passed cells"), x = 4, y = 3.8)
        
        p <- ggMarginal(p2, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountOutlier_2_Pass",pass,".png"),plot = p)
        
        # Cells to exclude lie below an intercept offset of -0.09
        cell.QC.stat$validCells <- log10(cell.QC.stat$nFeature_RNA) > (log10(cell.QC.stat$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09))
        
        p3 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point(aes(colour = validCells), fill = "black",pch=21) +
                scale_color_manual(breaks = c("TRUE", "FALSE"), 
                                   values = c("black","firebrick1")) +
                geom_smooth(method="lm") +
                geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") + 
                theme(legend.position="none") +
                annotate(geom = "text", label = paste0(as.numeric(table(cell.QC.stat$validCells)[2]), " QC passed cells\n",
                                                       as.numeric(table(cell.QC.stat$validCells)[1]), " QC filtered"), x = 4, y = 3.8)
        
        p <- ggMarginal(p3, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountOutlier_3_Pass",pass,".png"),plot = p)
        
        # Remove invalid cells
        cell.QC.stat <- cell.QC.stat %>% dplyr::filter(validCells)
        
        seur <- subset(seur, subset = Barcodes %in% cell.QC.stat$Barcodes)
        return(seur)
}

findNumPCs <- function(seur){
        pct <- seur[["pca"]]@stdev / sum(seur[["pca"]]@stdev) * 100
        cumu <- cumsum(pct)
        
        co1 <- which(cumu > 90 & pct < 5)[1]
        
        co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05),
                    decreasing = T)[1] + 1
        
        return(min(co1,co2))
}

sum_dups <- function(m){
        dup_genes <- rownames(m)[duplicated(rownames(m))]
        for(i in seq_along(dup_genes)){
                dg <- dup_genes[i]
                nu_vec <- colSums(m[rownames(m) == dg, ])
                matches <- which(rownames(m) == dg)
                m[matches[1], ] <- nu_vec
                m <- m[-matches[-1], ]
        }
        return(m)
}

# User parameters
################################################################################
setwd("../..")

brain_f <- "./data/ageAnno/brain.rds"
age_info_f <- "./data/ageAnno/age_info.csv"
cellCycGenesFile <- "./data/utility_files/CellCycleGenes_Human.csv"
aa_marker_f <- "./data/utility_files/scRNAmarker.txt"
out_dir <- "./results/parsed/ageAnno/"
create_dir_if_not(out_dir)

plot_dir <- sprintf("%splots/", out_dir)
seur_dir <- sprintf("%sseur/", out_dir)
create_dir_if_not(plot_dir)
create_dir_if_not(seur_dir)

# Load data
################################################################################
brain <- readRDS(brain_f)

aa_marker <- read_table_fast(aa_marker_f)
head(aa_marker)
colnames(aa_marker) <- make.names(unlist(as.vector(aa_marker[2, ])))
aa_marker <- aa_marker[3:nrow(aa_marker), ]
aa_marker <- aa_marker[aa_marker$Tissue == "Brain", ]

mdat <- read.csv(mdat_f, row.names = 1)
age_info <- mdat[, c("orig.ident", "age_death")] %>% distinct(orig.ident, age_death)
rownames(age_info) <- 1:nrow(age_info)
brain@meta.data <- cbind.data.frame(brain@meta.data,
                                    mdat[match(rownames(brain@meta.data),
                                               rownames(mdat)),
                                         c("cell_type",
                                           "cluster_id")])
brain@meta.data$age_death <- age_info$age_death[match(brain@meta.data$orig.ident,
                                                      age_info$orig.ident)]

cellcyclegenes <- read.csv(cellCycGenesFile)

mart <- useDataset("hsapiens_gene_ensembl",
                   useMart("ensembl"))
ensToSymb <- getBM(filters= "ensembl_gene_id",
                   attributes= c("ensembl_gene_id","hgnc_symbol"),
                   values=cellcyclegenes$geneID,mart= mart)
cellcyclegenes <- merge(cellcyclegenes,ensToSymb,by.x = 2, by.y = 1, all.x = T)

# Process data
################################################################################

# Remove samples with NAs in the ages, as these correspond to old samples with
# 90+ age (ambiguous age)
brain <- brain[, rownames(brain@meta.data)[!is.na(brain@meta.data$age_death)]]
# Remove samples less than 18 yo
brain <- brain[, rownames(brain@meta.data)[brain@meta.data$age_death >= 18]]

# Get ENSEMBL IDs from the symbols
ensemblIDs <- getBM(filters = "hgnc_symbol",
                    attributes = c("ensembl_gene_id","hgnc_symbol",
                                   "start_position",
                                   "end_position",
                                   "chromosome_name",
                                   "strand"),
                    values = rownames(brain),
                    mart = mart)

ensemblIDs <- ensemblIDs[nchar(ensemblIDs$chromosome_name) <= 3, ]

mod_coefs$ensembl_gene_id %in% ensemblIDs$ensembl_gene_id

viol_pre_nfeats <- VlnPlot(brain, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.0000)  + 
        geom_hline(yintercept = 300, linetype = "dashed", color = "red")
viol_pre_ncount <- VlnPlot(brain, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0.0000)
viol_pre_mitprc <- VlnPlot(brain, features = "percent.mt", group.by = "orig.ident", pt.size = 0.0000)

ggsave(filename = sprintf("%sviol_pre_nfeats.pdf",
                          out_dir),
       viol_pre_nfeats)
ggsave(filename = sprintf("%sviol_pre_ncount.pdf",
                          out_dir),
       viol_pre_ncount)
ggsave(filename = sprintf("%sviol_pre_mitprc.pdf",
                          out_dir),
       viol_pre_mitprc)

# Remove mid3, mid4, mid5 and old18, as they look off in the distributions.
samps_to_rem <- c("mid3", "mid4", "mid5", "old18")

keep_cells <- rownames(brain@meta.data)[!brain@meta.data$orig.ident %in% samps_to_rem]
brain <- brain[, keep_cells]

# Split object per original identity (sample), so we can preprocess each sample
# separately
brain_split <- Seurat::SplitObject(brain, split.by = "orig.ident")
brain_proc <- list()

for (i in seq_along(brain_split)){
        samp <- names(brain_split)[i]
        out_seur_file <- sprintf("%s%s.rds", seur_dir, samp)
        if (!file.exists(out_seur_file)){
                print(sprintf("Processing %s...", samp))
                seur_samp <- brain_split[[samp]]
                minCells <- 100
                minFeats <- 300
                nfeatures <- Matrix::colSums(x = seur_samp@assays$RNA$counts > 0)
                nCellsPassed <- length(which(x = nfeatures >= minFeats))
                cellsPassed <- names(nfeatures)[nfeatures >= minFeats]
                if (length(cellsPassed) > minCells){
                        seur_samp <- seur_samp[, cellsPassed]
                        plot_samp_dir <- sprintf("%s%s/", plot_dir, samp)
                        create_dir_if_not(plot_samp_dir)
                        seur_samp <- filterCells(seur_samp,
                                                 mad.coeff = 3,
                                                 pass = -1,
                                                 org = "HUMAN",
                                                 plot.path = plot_samp_dir)
                        # Log-normalize data
                        seur_samp <- NormalizeData(seur_samp)
                        # Find variable features --> Identify features that are outliers on a mean
                        # variability plot
                        seur_samp <- FindVariableFeatures(seur_samp,
                                                          selection.method = "vst",
                                                          nfeatures = 3000)
                        # Scale data --> Scales and centers features in the dataset
                        seur_samp <- ScaleData(seur_samp,
                                               features = VariableFeatures(object = seur_samp))
                        # Run PCA
                        if(ncol(seur_samp) < 100){
                                nPCs <- ncol(seur_samp) - 1
                        }else{
                                nPCs <- 100
                        }
                        seur_samp <- RunPCA(seur_samp,
                                            features = VariableFeatures(object = seur_samp),
                                            npcs = nPCs)
                        g2m_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "G2/M"]
                        g2m_genes <- intersect(g2m_genes, rownames(seur_samp))
                        s_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "S"]
                        s_genes <- intersect(s_genes, rownames(seur_samp))
                        seur_samp <- CellCycleScoring(seur_samp,
                                                      s.features = s_genes,
                                                      g2m.features = g2m_genes,
                                                      set.ident = F)
                        
                        ggsave(paste0(plot_samp_dir,
                                      "PCA_CellCycleGenes.png"),
                               DimPlot(seur_samp,
                                       group.by = "Phase"))
                        seur_samp$CC.Difference <- seur_samp$S.Score - seur_samp$G2M.Score
                        
                        seur_samp <- SCTransform(seur_samp,
                                                 vars.to.regress = "CC.Difference",
                                                 vst.flavor = "v2",
                                                 verbose = T)
                        
                        seur_samp <- RunPCA(seur_samp,
                                            assay = "SCT",
                                            npcs = nPCs,
                                            verbose = T)
                        
                        numPCs <- findNumPCs(seur_samp)
                        elbPlot <- ElbowPlot(seur_samp,
                                             ndims = nPCs) +
                                geom_vline(aes(xintercept = numPCs),
                                           colour = "red",
                                           linetype = "dashed")
                        
                        ggsave(paste0(plot_samp_dir,
                                      "ElbowPlot.png"),
                               elbPlot)
                        # This code works with DoubletFinder v2.0.6
                        sweep.list <- paramSweep(seur_samp,
                                                 PCs = 1:numPCs,
                                                 sct = T)
                        
                        sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
                        bcmvn <- find.pK(sweep.stats)
                        ## Estimate expected percentage of doublets from 10X Genomics 
                        # estimates from 3' Gene Expression v3.1 assay##
                        estDoublets <- c(0.4,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8)
                        numCellsRec <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
                        scatter.smooth(numCellsRec, estDoublets) #Looks linear
                        lm_doublets <- lm(estDoublets ~ numCellsRec)
                        summary(lm_doublets) #Perfect linear relationship r2 = 1
                        
                        nExp <- round(ncol(seur_samp) * (unname(predict(lm_doublets,
                                                                        data.frame(numCellsRec = ncol(seur_samp))))/100))
                        
                        pK <- as.numeric(levels(bcmvn$pK)[bcmvn$BCmetric == max(bcmvn$BCmetric)])
                        
                        seur_samp <- doubletFinder(seu = seur_samp,
                                                   PCs = 1:numPCs,
                                                   pN = 0.25, #default
                                                   pK = pK,
                                                   nExp = nExp,
                                                   reuse.pANN = NULL,
                                                   sct = TRUE,
                                                   annotations = NULL)
                        seur_samp$Doublets <- "Singlet"
                        
                        seur_samp$Doublets[seur_samp[[paste0("pANN_0.25_",
                                                             pK,
                                                             "_",
                                                             nExp)]] >= 0.5] <- "Doublet"
                        
                        seur_samp$Doublets <- factor(seur_samp$Doublets,
                                                     levels = c("Doublet",
                                                                "Singlet"))
                        
                        p <- ggplot(seur_samp@meta.data,
                                    aes(x=log10(nCount_RNA),
                                        y=log10(nFeature_RNA))) +
                                geom_point(aes(colour = Doublets), fill = "black",pch=21) +
                                scale_color_manual(breaks = c("Singlet", "Doublet"), 
                                                   values = c("black","firebrick1")) +
                                geom_smooth(method="lm") +
                                theme(legend.position="none") +
                                annotate(geom = "text",
                                         label = paste0(as.numeric(table(seur_samp@meta.data$Doublets)[2]),
                                                        " Singlets\n",
                                                        as.numeric(table(seur_samp@meta.data$Doublets)[1]),
                                                        " Doublets"),
                                         x = 4,
                                         y = 3.8)
                        
                        ggsave(paste0(plot_samp_dir,"Doublets.png"),plot = p)
                        seur_samp <- subset(seur_samp,
                                            subset = Doublets == "Singlet")
                        seur_samp <- RunUMAP(seur_samp,
                                             dims = 1:numPCs,
                                             n.neighbors = 20)
                        umap_seur_clusts <- DimPlot(seur_samp,
                                                    reduction = "umap",
                                                    group.by = "cell_type")
                        ggsave(filename = sprintf("%sumap_ref_cell_types.pdf",
                                                  plot_samp_dir),
                               plot = umap_seur_clusts)
                        saveRDS(seur_samp, file = out_seur_file)
                }else{
                        print(sprintf("%s doesn't have enough good quality cells. Skipping.",
                                      samp))
                }
        }else{
                print(sprintf("%s has already been processed. Loading.",
                              samp))
                seur_samp <- readRDS(out_seur_file)
        }
        brain_proc[[samp]] <- seur_samp
}
brain_proc

# Merge and inspect again
################################################################################

# Merge processed objects
brain_merged <- merge(x = brain_proc$youth4,
                      y = brain_proc[2:length(brain_proc_vec)])
DefaultAssay(brain_merged) <- "RNA"

viol_pos_nfeats <- VlnPlot(brain_merged, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.0000)  + 
        geom_hline(yintercept = 300, linetype = "dashed", color = "red")
viol_pos_ncount <- VlnPlot(brain_merged, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0.0000)
viol_pos_mitprc <- VlnPlot(brain_merged, features = "percent.mt", group.by = "orig.ident", pt.size = 0.0000)

ggsave(filename = sprintf("%sviol_pos_nfeats.pdf",
                          out_dir),
       viol_pre_nfeats)
ggsave(filename = sprintf("%sviol_pos_ncount.pdf",
                          out_dir),
       viol_pre_ncount)
ggsave(filename = sprintf("%sviol_pos_mitprc.pdf",
                          out_dir),
       viol_pre_mitprc)

# Let's remove samples old13, old14, old15 and old21, as they have really high 
# mitochondrial percentage.
high_mito_samps <- c("old13", "old14", "old15", "old21")
brain_merged <- brain_merged[, rownames(brain_merged@meta.data)[!brain_merged@meta.data$orig.ident %in% high_mito_samps]]

# Integrate
################################################################################
brain_merged <- NormalizeData(brain_merged)

g2m_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "G2/M"]
g2m_genes <- intersect(g2m_genes, rownames(brain_merged))
s_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "S"]
s_genes <- intersect(s_genes, rownames(brain_merged))

brain_merged <- CellCycleScoring(brain_merged,
                                 s.features = s_genes,
                                 g2m.features = g2m_genes)
# Remove SCT slot
brain_merged[["SCT"]] <- NULL
brain_merged[["RNA"]] <- split(brain_merged[["RNA"]], f = brain_merged$orig.ident)
brain_merged <- FindVariableFeatures(brain_merged, verbose = FALSE)

###Subset representative for each dataset
brain_merged <- SketchData(object = brain_merged,
                           ncells = 5000,
                           method = "LeverageScore",
                           sketched.assay = "sketch")

### Integration on the sketched cells accross samples
DefaultAssay(brain_merged) <- "sketch"
brain_merged <- FindVariableFeatures(brain_merged, verbose = F)

brain_merged <- ScaleData(brain_merged,
                          vars.to.regress=c("S.Score","G2M.Score","percent.mito"))
brain_merged <- RunPCA(brain_merged, verbose = F)

brain_merged <- IntegrateLayers(brain_merged,
                                method = RPCAIntegration,
                                orig = "pca",
                                new.reduction = "integrated.rpca")

saveRDS(brain_merged, file=sprintf("%sbrain_sketch_integ.rds", out_dir))

# Annotate cell types
################################################################################

# Cluster integrated data
brain_merged <- FindNeighbors(brain_merged,
                              reduction = "integrated.rpca", dims = 1:50)
brain_merged <- FindClusters(brain_merged, resolution = seq(0,1, by=0.1))

brain_merged <- RunUMAP(brain_merged,
                        reduction = "integrated.rpca",
                        dims = 1:50,
                        return.model = T,
                        verbose = F)


umap_orig.ident <- DimPlot(brain_merged,
                           reduction = "umap",
                           group.by = "orig.ident")
ggsave(filename = sprintf("%sumap_orig.ident.pdf",
                          out_dir),
       plot = umap_orig.ident,
       height = 6.8,
       width = 8)
umap_cell_type <- DimPlot(brain_merged,
                          reduction = "umap",
                          group.by = "cell_type",
                          label = T)
umap_seur_clusts_0.1 <- DimPlot(brain_merged,
                                reduction ="umap",
                                group.by = "sketch_snn_res.0.1",
                                label = T)

ggsave(filename = sprintf("%sumap_seur_clusts_0.1.pdf",
                          out_dir),
       plot = umap_seur_clusts_0.1,
       height = 6.8,
       width = 8)

clusts <- unique(brain_merged$sketch_snn_res.0.2)
unsupMarkers <- list()
for(i in seq_along(clusts)){
        clust <- clusts[i]
        rest <- clusts[clusts != clust]
        clustMarks <- FindMarkers(brain_merged,
                                  assay = "RNA",
                                  slot = "data",
                                  ident.1 = clust,
                                  ident.2 = rest)
        unsupMarkers[[clust]] <- clustMarks
}
saveRDS(unsupMarkers, file = sprintf("%smarkers_res0.2.rds",
                                     out_dir))

# Add EPAS1 to the endothelial markers, to double check,
# as VTN is not expressed in cluster 14.
to_bind <- data.frame(Tissue = "Brain",
                      Cell.type = "Endothelial",
                      Marker.gene = "EPAS1")

aa_marker <- rbind.data.frame(to_bind,
                              aa_marker)

marker_dotplot <- DotPlot(object = brain_merged,
                          features = aa_marker$Marker.gene,
                          group.by = 'sketch_snn_res.0.1') + RotatedAxis() +
        geom_vline(xintercept = c(2.5,
                                  11.5,
                                  16.5,
                                  21.5,
                                  25.5),
                   linetype = "dashed",
                   color = "gray") + 
        annotate("text", x = 1.25, y = 7.5, label = "Endothelial",
                 angle = 90, vjust = 1.2, size = 6, alpha = 0.5) +
        annotate("text", x = 7, y = 7.5, label = "Astrocytes",
                 angle = 90, vjust = 1.2, size = 6, alpha = 0.5) +
        annotate("text", x = 14, y = 7.5, label = "OPCs",
                 angle = 90, vjust = 1.2, size = 6, alpha = 0.5) +
        annotate("text", x = 19, y = 7.5, label = "Oligodendrocytes",
                 angle = 90, vjust = 1.2, size = 6, alpha = 0.5) +
        annotate("text", x = 23.5, y = 7.5, label = "Excitatory neurons",
                 angle = 90, vjust = 1.2, size = 6, alpha = 0.5) +
        annotate("text", x = 26.5, y = 7.5, label = "Inhibitory neurons",
                 angle = 90, vjust = 1.2, size = 6, alpha = 0.5)

ggsave(filename = sprintf("%smarker_dotplot.pdf", out_dir),
       plot = marker_dotplot,
       width = 10, height = 6)

# Cluster 14: endothelial?
# Cluster 6: Astrocytes
# Cluster 7: OPCs
# Cluster 1: Oligodendrocytes
# Clusters 0, 3, 5, 8, 9, 10, 12 and 13: Excitatory neurons
# Clusters 2, 4, 11 and 15: Inhibitory neurons

manual_annot <- data.frame(cluster = factor(0:15),
                           cell_type = c("Excitatory neurons",
                                         "Oligodendrocytes",
                                         "Inhibitory neurons",
                                         "Excitatory neurons",
                                         "Inhibitory neurons",
                                         "Excitatory neurons",
                                         "Astrocytes",
                                         "OPCs",
                                         "Excitatory neurons",
                                         "Excitatory neurons",
                                         "Excitatory neurons",
                                         "Inhibitory neurons",
                                         "Excitatory neurons",
                                         "Excitatory neurons",
                                         "Endothelial",
                                         "Inhibitory neurons"))

unique(brain_merged@meta.data$sketch_snn_res.0.1)

brain_merged@meta.data$cell_type_manual <- manual_annot$cell_type[match(brain_merged@meta.data$sketch_snn_res.0.1,
                                                                        manual_annot$cluster)]

umap_man_ann <- DimPlot(brain_merged,
                        reduction = "umap",
                        group.by = "cell_type_manual",
                        label = T)
ggsave(filename = sprintf("%sumap_man_ann.pdf",
                          out_dir),
       plot = umap_man_ann,
       height = 6.8,
       width = 8)

umap_man_ann_splitSamp <- DimPlot(brain_merged,
                                  reduction="umap",
                                  group.by = "cell_type_manual",
                                  label = T,
                                  split.by = "orig.ident",
                                  ncol = 4)
ggsave(filename = sprintf("%sumap_man_ann_splitSamp.pdf",
                          out_dir),
       plot = umap_man_ann_splitSamp)

# Obtain pseudo-bulk counts
################################################################################

counts <- LayerData(brain_merged, assay = "RNA", layer = "counts")
meta <- brain_merged@meta.data
groups <- paste(make.names(meta$cell_type_manual),
                meta$orig.ident,
                sep = "_")
pseudobulk <- t(rowsum(t(counts), group = groups))

nu_names <- ensemblIDs$ensembl_gene_id[match(make.names(rownames(pseudobulk)),
                                             make.names(ensemblIDs$hgnc_symbol))]

pseudobulk <- pseudobulk[!is.na(nu_names), ]
nu_names <- nu_names[!is.na(nu_names)]
rownames(pseudobulk) <- nu_names

# Add up duplicates that might happen when changing IDs
pseudobulk <- sum_dups(pseudobulk)

write_table_fast(pseudobulk, f = sprintf("%sageanno_brain_pb_counts.csv", out_dir))

# Parse metadata
################################################################################

meta_parsed <- data.frame(specimenID = colnames(pseudobulk),
                          platform = "",
                          RIN = 10,
                          individualID = gsub(".*_", "", colnames(pseudobulk)),
                          organ = "brain",
                          tissue = gsub("\\_.*", "", colnames(pseudobulk)),
                          assay = "scRNAseq",
                          substudy = "ageAnno",
                          diagn_4BrainClck = "scRNAseq")
meta_parsed$ageDeath <- age_info$age_death[match(meta_parsed$individualID,
                                                 age_info$orig.ident)]
write_table_fast(meta_parsed,
                 f = sprintf("%sageanno_brain_pb_metdat.csv", out_dir))