if (!require("preprocessCore", quietly = T)){
        BiocManager::install("preprocessCore", configure.args=c(preprocessCore = "--disable-threading"),
                             update = F)
}
library(preprocessCore)
if(!require(bigstatsr, quietly = T)){
        install.packages("bigstatsr",
                         repos='http://cran.us.r-project.org')
}
library(bigstatsr)
if(!require("limma", quietly = T)){
        BiocManager::install("limma", update = F)
}
library(limma)

if(!require("GEOquery", quietly = T)){
        BiocManager::install("GEOquery", update = F)
}

if(!require("SmartSVA", quietly = T)){
        install.packages("SmartSVA", repos='http://cran.us.r-project.org')
}

if(!require("isva", quietly = T)){
        install.packages("isva", repos='http://cran.us.r-project.org')
}

if(!require("qvalue", quietly = T)){
        BiocManager::install("qvalue")
}

library(SmartSVA)
if(!require(sva, quietly = T)) BiocManager::install("sva", update = F)

library(plotUtils)
library(ggplot2)

library(dplyr)
# Functions
################################################################################

big_PCA <- function(df, nComps = 20, stand = F){
        if(stand){
                print("Standarizing the input dataset...")
                # Remove genes with SD == 0
                df <- df[, apply(df, 2, sd) != 0]
                scaled_DF <- apply(df, 2, function(x) (x - mean(x))/sd(x))
                scal = T
                centr = T
        }else{
                scaled_DF <- df
                scal = F
                centr = F
        }
        
        scaled_DF_fbm <- FBM(nrow(scaled_DF),
                             ncol(scaled_DF),
                             init = scaled_DF)
        
        print("Computing SVD...")
        scaled_DF_svd <- big_SVD(scaled_DF_fbm,
                                 fun.scaling = big_scale(center = F,
                                                         scale = F),
                                 k = nComps)
        
        # Multiply U and D to get the scores
        print("Generating PCA output...")
        print("Computing scores matrix...")
        scores <- scaled_DF_svd$u %*% diag(scaled_DF_svd$d)
        rownames(scores) <- rownames(scaled_DF)
        colnames(scores) <- paste0("PC", as.character(1:ncol(scores)))
        rotation <- scaled_DF_svd$v
        rownames(rotation) <- colnames(scaled_DF)
        colnames(rotation) <- paste0("PC", as.character(1:ncol(rotation)))
        print("Computing eigenvalues...")
        eigenVals <- scaled_DF_svd$d^2
        
        print("Computing total variance of the initial dataset...")
        totVar <- sum(diag(t(scaled_DF) %*% scaled_DF))
        varExp <- eigenVals/totVar
        
        cumulative_sum <- function(v, index = 1, current_sum = 0) {
                if (index > length(v)){
                        return(c())
                }
                
                current_sum <- current_sum + v[index]
                
                return(c(current_sum, cumulative_sum(v, index + 1, current_sum)))
        }
        
        print("Computing PC cumulative variance...")
        cumVar <- cumulative_sum(varExp)
        pcSD <- apply(scores, 2, sd)
        
        summMat <- rbind(pcSD, varExp, cumVar)
        rownames(summMat) <- c("Standard deviation",
                               "Proportion of Variance",
                               "Cumulative Proportion")
        outList <- list(sdev = pcSD,
                        rotation = rotation,
                        scale = scal,
                        center = centr,
                        x = scores,
                        summary = summMat)
        print("Done")
        return(outList)
}

pcBiplot <- function(PC, x="PC1", y="PC2", varPlotFilt = NULL, biPlot = F,
                     colBy = "tissue", colROSMAPBatch = F, colLINCSCell = F){
        #metDat <- metDat
        if(colROSMAPBatch){
                metDat$substudy[metDat$substudy == "ROSMAP"] <- paste("ROSMAP",
                                                                      metDat$batch_seq[metDat$substudy == "ROSMAP"],
                                                                      sep = "_")
        }
        if(colLINCSCell){
                metDat$substudy[metDat$substudy == "LINCS"] <- paste("LINCS",
                                                                     metDat$tissue[metDat$substudy == "LINCS"],
                                                                     sep = "_")
        }
        data <- data.frame(obsnames=row.names(PC$x), PC$x)
        data <- data[, c("obsnames", x, y)]
        data <- cbind.data.frame(data, metDat[match(rownames(data),
                                                    make.names(metDat$specimenID)), colnames(metDat) != "specimenID"])
        #data$tissue <- metDat$tissue[match(rownames(data),
        #                                   make.names(metDat$specimenID))]
        
        #data$substudy <- metDat$substudy[match(rownames(data),
        #                                       make.names(metDat$specimenID))]
        #data$pmi <- metDat$pmi[match(rownames(data),
        #                             make.names(metDat$specimenID))]
        #data$individualID <- metDat$individualID[match(rownames(data),
        #                                               make.names(metDat$specimenID))]
        #data$age_death <- metDat$ageDeath[match(rownames(data),
        #                                        make.names(metDat$specimenID))]
        #data$diagnosis <- metDat$diagn_4BrainClck[match(rownames(data),
        #                                          make.names(metDat$specimenID))]
        #data$exper_group <- metDat$exper_group[match(rownames(data),
        #                                       make.names(metDat$specimenID))]
        propVar <- PC$summary[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        plot <- ggplot(data, aes_string(x = x, 
                                        y = y, 
                                        color = colBy)) + 
                #scale_discrete_manual("Time",
                #                      aesthetics = "colour",
                #                      values = c("cyan",
                #                                 "dodgerblue",
                #                                 "dodgerblue4")) +
                geom_hline(yintercept = 0, alpha = 0.6) +
                geom_vline(xintercept = 0, alpha = 0.6) +
                geom_point() + 
                xlab(sprintf("%s (%s%%)", x, propX)) +
                ylab(sprintf("%s (%s%%)", y, propY)) +
                #geom_text_repel() +
                #theme_minimal()
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if(biPlot){
                datapc <- data.frame(varnames=rownames(PC$rotation), 
                                     PC$rotation)
                mult <- min(
                        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
                        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
                )
                datapc <- transform(datapc,
                                    v1 = .7 * mult * (get(x)),
                                    v2 = .7 * mult * (get(y))
                )
                datapc$x0 <- rep(0, nrow(datapc))
                datapc$y0 <- rep(0, nrow(datapc))
                if(!is.null(varPlotFilt)){
                        datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
                }
                plot <- plot +
                        geom_text_repel(data=datapc, 
                                        aes(x=v1, y=v2, label=varnames), 
                                        color = "black", 
                                        size = 3) + 
                        geom_segment(data = datapc, aes(x=x0, 
                                                        y=y0, 
                                                        xend=v1, 
                                                        yend=v2, 
                                                        label = varnames),
                                     arrow = arrow(length=unit(0.2,"cm"),
                                                   type = "closed",
                                                   angle = 20), 
                                     alpha=0.75, 
                                     color="black", 
                                     size = 0.5)
        }
        return(plot)
}

# Plots Multi PC score/biplot
doPCAMultiPlot <- function(PC,
                           nComps,
                           varPlotFilt = NULL,
                           biPlot = F,
                           colBy = "tissue",
                           colROSMAPBatch = F){
        plotList <- list()
        for(j in 2:(nComps)){
                for(i in 1:(nComps - 1)){
                        if(j > i){
                                scPlot <- pcBiplot(PC,
                                                   biPlot = F,
                                                   x = sprintf("PC%s", i),
                                                   y = sprintf("PC%s", j),
                                                   varPlotFilt = varPlotFilt,
                                                   colBy = colBy,
                                                   colROSMAPBatch = colROSMAPBatch)
                                if(j < nComps){
                                        scPlot <- scPlot +
                                                theme(axis.title.x = element_blank(),
                                                      axis.text.x = element_blank())
                                }
                                if(i > 1){
                                        scPlot <- scPlot +
                                                theme(axis.title.y = element_blank(),
                                                      axis.text.y = element_blank())
                                }
                        }else{
                                scPlot <- NA
                        }
                        plotList[[sprintf("PC%s_PC%s", i, j)]] <- scPlot
                }
        }
        multPlot <- ggarrange(plotlist = plotList,
                              common.legend = T,
                              ncol = nComps - 1,
                              nrow = nComps - 1,
                              widths = c(1, rep(.8, nComps-2)),
                              heights = c(rep(.8, nComps-2), 1))
        return(multPlot)
}

# Directory stuff
################################################################################

counts_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/counts_parsed/counts_merged.csv"
mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/counts_parsed/sample_metadata_merged.csv"
out_dir <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/preprocessing"
coefs_f <- "/Users/guillem.santamaria/Documents/postdoc/manuscripts/brain_clock/submission_advanced_science/round_2/supp_tabs_material/coefs_wInter_wBrain.csv"

brainseq_counts_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/brainseq_pII_counts.csv"
brainseq_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/brainseq_pII_metadata.csv"

GSE216281_counts_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/GSE216281/counts.csv"
GSE216281_mdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/counts/GSE216281/401_2023_2597_MOESM2_ESM.xlsx"


sysdat_f <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/scripts/brainAgeShiftR/inst/extdata/sysdata.rda"

create_dir_if_not(out_dir)


# Load data
################################################################################
library(GEOquery)

counts <- read_table_fast(counts_f, row.names = 1)
mdat <- read_table_fast(mdat_f, row.names = 1)
coefs <- read.csv(coefs_f, row.names = 1)
ind_mdat <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/results/counts_parsed/individual_metadata_merged_rnaseq_filt.csv")

load(sysdat_f)

counts <- counts[bg_genes, !mdat$tissue %in% c("Cerebellar Hemisphere", "cerebellum", "Cerebellum")]
mdat <- mdat[!mdat$tissue %in% c("Cerebellar Hemisphere", "cerebellum", "Cerebellum"), ]

brainseq_counts <- read.csv(brainseq_counts_f, row.names = 1)
brainseq_mdat <- read.csv(brainseq_mdat_f, row.names = 1)

GSE216281_counts <- read.csv(GSE216281_counts_f, row.names = 1)
GSE216281_mdat <- as.data.frame(readxl::read_xlsx(GSE216281_mdat_f, sheet = 2))

colnames(GSE216281_mdat) <- unlist(as.vector(GSE216281_mdat[1, ]))

GSE216281_mdat <- GSE216281_mdat[2:nrow(GSE216281_mdat), ]


GSE216281_id_mapping_f <- file.path(out_dir, "GSE216281_id_mapping.csv")
if(!file.exists(GSE216281_id_mapping_f)){
        GSE216281_geo <- GEOquery::getGEO("GSE216281", GSEMatrix = FALSE)
        
        GSE216281_sampMdat <- GEOquery::GSMList(GSE216281_geo)
        
        library(rentrez)
        GSE216281_id_mapping <- list()
        for(i in seq_along(GSE216281_sampMdat)){
                print(i/length(GSE216281_sampMdat))
                samp_mdat <- Meta(GSE216281_sampMdat[[i]])
                series_id <- samp_mdat$series_id
                geo_accession <- samp_mdat$geo_accession
                srx <- gsub(".*=", "", samp_mdat$relation[2])
                sra_search <- entrez_search(
                        db = "sra",
                        term = srx
                )
                xml <- entrez_fetch(
                        db = "sra",
                        id = sra_search$ids,
                        rettype = "runinfo",
                        parsed = FALSE
                )
                xml <- strsplit(xml, split = ",")[[1]]
                srr <- gsub("\\..*", "", basename(xml[grep("SRR", xml)][2]))
                tissue <- samp_mdat$source_name_ch1
                descript <- samp_mdat$description
                paper_id <- gsub("Paper ID: ", "", descript[1])
                xlsx_id <- descript[2]
                tobind <- data.frame(accession = series_id,
                                     geo_accession = geo_accession,
                                     srx = srx,
                                     srr = srr,
                                     paper_id = paper_id,
                                     xslx_id = xlsx_id,
                                     tissue = tissue)
                GSE216281_id_mapping[[i]] <- tobind
        }
        GSE216281_id_mapping <- do.call(rbind.data.frame, GSE216281_id_mapping)
        
        
        
        run_mat <- matrix(GSE216281_id_mapping$srr)
        write.table(
                run_mat, file = file.path(out_dir, "accessions_sra.txt"),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t"
        )
        write.csv(GSE216281_id_mapping, file = GSE216281_id_mapping_f)
}else{
        GSE216281_id_mapping <- read.csv(GSE216281_id_mapping_f, row.names = 1)
}


GSE216281_mdat <- cbind.data.frame(
        GSE216281_mdat,
        GSE216281_id_mapping[match(GSE216281_mdat$`Sample ID`, GSE216281_id_mapping$paper_id), ]
)

# THere are a couple of IDs that don't match. Just remove them
GSE216281_mdat$`Sample ID`[!GSE216281_mdat$`Sample ID` %in% GSE216281_id_mapping$paper_id]
GSE216281_id_mapping[!GSE216281_id_mapping$paper_id %in% GSE216281_mdat$`Sample ID`, ]

GSE216281_mdat <- GSE216281_mdat[!is.na(GSE216281_mdat$tissue), ]
GSE216281_counts <- GSE216281_counts[, match(GSE216281_mdat$srr, colnames(GSE216281_counts))]
colnames(GSE216281_counts) <- GSE216281_mdat$paper_id

rownames(GSE216281_counts) <- gsub("\\..*", "", rownames(GSE216281_counts))

plot(density(as.matrix(GSE216281_counts)))


head(GSE216281_mdat)

GSE216281_counts[bg_genes, ]
missing_genes <- bg_genes[!bg_genes %in% rownames(GSE216281_counts)]

#library(biomaRt)

#mart <- useEnsembl(
#        biomart = "genes",
#        dataset = "hsapiens_gene_ensembl",
#        mirror = "useast"
#)


#res <- getBM(
#        attributes = c(
#                "ensembl_gene_id",
#                "external_gene_name",
#                "gene_biotype",
#                "ensembl_gene_id_version"
#        ),
#        filters = "ensembl_gene_id",
#        values = missing_genes,
#        mart = mart
#)

#res



comm_genes <- intersect(
        intersect(rownames(brainseq_counts), rownames(counts)),
        rownames(GSE216281_counts)
)

brainseq_mdat$ageDeath <- as.character(round(brainseq_mdat$ageDeath))
brainseq_mdat$CERAD <- NA
brainseq_mdat$cohort <- "brainSeq_pII"
brainseq_mdat$cohort <- "brainSeq_pII"
brainseq_mdat$diagnosis <- "control"
brainseq_mdat$substudy_batch <- brainseq_mdat$substudy

colnames(GSE216281_mdat) <- gsub(
        "Sample ID", "specimenID", colnames(GSE216281_mdat)
)

GSE216281_mdat$individualID <- GSE216281_mdat$specimenID

colnames(GSE216281_mdat) <- gsub(
        "Braak stage", "Braak", colnames(GSE216281_mdat)
)

GSE216281_mdat$Braak <- as.numeric(GSE216281_mdat$Braak)

colnames(GSE216281_mdat) <- gsub(
        "Sex", "sex", colnames(GSE216281_mdat)
)

GSE216281_mdat$sex <- ifelse(GSE216281_mdat$sex == "M", "male", "female")

GSE216281_mdat$tissue <- "Frontal cortex"

colnames(GSE216281_mdat) <- gsub(
        "Age at death", "ageDeath", colnames(GSE216281_mdat)
)

GSE216281_mdat$diagn_4BrainClck <- "Rest"

colnames(GSE216281_mdat) <- gsub(
        "Diagnosis", "diagnosis", colnames(GSE216281_mdat)
)

diagn_tab <- data.frame(
        diagn_key = c(
                "Ctrl", "iLBD", "PD", "PDD"
        ),
        diagnosis = c(
                "control", "iLBD", "Parkinson's disease", "Parkinson's disease"
        )
)

GSE216281_mdat$diagnosis <- diagn_tab$diagnosis[
        match(GSE216281_mdat$diagnosis, diagn_tab$diagn_key)
]

table(GSE216281_mdat$diagnosis)

GSE216281_mdat$substudy <- "GSE216281"
GSE216281_mdat$substudy_batch <- GSE216281_mdat$substudy
colnames(GSE216281_mdat) <- gsub(
        "PMD (min)", "pmi", colnames(GSE216281_mdat), fixed = TRUE
)
GSE216281_mdat$pmi <- as.numeric(GSE216281_mdat$pmi) / 60

GSE216281_mdat$RIN <- round(as.numeric(GSE216281_mdat$RIN), digits = 1)

comm_cols_mdat <- intersect(
        intersect(colnames(mdat), colnames(brainseq_mdat)),
        colnames(GSE216281_mdat)
)

mdat <- rbind.data.frame(mdat[, comm_cols_mdat], brainseq_mdat[, comm_cols_mdat], GSE216281_mdat[, comm_cols_mdat])
counts <- cbind.data.frame(counts[comm_genes, ], brainseq_counts[comm_genes, ], GSE216281_counts[comm_genes, ])


# Normalize
################################################################################
counts_log2 <- log2(counts + 1)

#counts_log2_qnorm <- normalize.quantiles(as.matrix(counts_log2),
#                                         keep.names = TRUE)

#write_table_fast(counts_log2_qnorm, f = file.path(out_dir, "counts_log2_qnorm.csv"))

# Perform PCA
################################################################################

counts_log2_pca <- big_PCA(t(counts_log2), nComps = 20, stand = TRUE)

#counts_qnorm_pca <- big_PCA(t(counts_log2_qnorm), nComps = 20, stand = TRUE)

rownames(mdat) <- make.names(mdat$specimenID)

mdat$age_numeric <- mdat$ageDeath
mdat$age_numeric[mdat$age_numeric == ""] <- NA
mdat$age_numeric <- gsub("+", "", mdat$age_numeric, fixed = TRUE)
mdat$age_numeric <- as.numeric(mdat$age_numeric)

metDat <- mdat

pca_substudy <- pcBiplot(counts_log2_pca, colBy = "substudy_batch", biPlot = F,
                         colROSMAPBatch = F)


ggsave(
        filename = file.path(out_dir, "counts_log2.pdf"),
        plot = pca_substudy,
        height = 10, width = 10
)

pca_age <- pcBiplot(counts_log2_pca, colBy = "age_numeric", biPlot = F,
                         colROSMAPBatch = F)


# Run SVA
################################################################################

# Keep only samples with age informnation, and re-normalize again

counts <- counts[, !is.na(mdat$age_numeric)]
#mdat <- mdat[!is.na(mdat$age_numeric), ]

#counts_log2 <- log2(counts + 1)

counts_log2 <- counts_log2[, !is.na(mdat$age_numeric)]
mdat <- mdat[!is.na(mdat$age_numeric), ]

counts_log2_qnorm <- normalize.quantiles(as.matrix(counts_log2),
                                         keep.names = TRUE)

all(make.names(mdat$specimenID) == colnames(counts_log2))

svobj_f <- file.path(out_dir, "svobj.rds")
if (!file.exists(svobj_f)){
        mod_onlyAge = model.matrix(~age_numeric, data=mdat)
        mod0_onlyAge = model.matrix(~1, data=mdat)
        
        n.sv = num.sv(counts_log2_qnorm, mod_onlyAge, method="leek")
        
        svobj = smartsva.cpp(
                as.matrix(counts_log2_qnorm),
                mod_onlyAge,
                mod0_onlyAge,
                n.sv=n.sv)
        
        saveRDS(svobj, file = file.path(out_dir, "svobj.rds"))
}else{
        svobj <- readRDS(svobj_f)
}


#counts_log2_qnorm_adj <- removeBatchEffect(counts_log2_qnorm, covariates = svobj$sv)

#counts_qnorm_adj_pca <- big_PCA(t(counts_log2_qnorm_adj), nComps = 20, stand = TRUE)

#metDat <- mdat
#pca_substudy_adj <- pcBiplot(
#        counts_qnorm_adj_pca, colBy = "substudy", biPlot = F,
#        colROSMAPBatch = F)

# Predict age with brainAgeShiftR
################################################################################
if (!require("brainAgeShiftR", quietly = TRUE)){
        devtools::install_github("guisantagui/brainAgeShiftR")
}
library(brainAgeShiftR)

# Shift adjusted data to positive space
#offset <- min(counts_log2_qnorm_adj)
#offset <- 20
#expression_4clock <- counts_log2_qnorm_adj# - offset

baso <- create_brainAgeShiftObj(
        counts,
        metadata = mdat
)

baso <- normalizeCounts(baso, useTrainMeans = TRUE, skip_log2 = FALSE)

baso <- predictAge(baso)

plot(
        baso$metadata$age_numeric[baso$metadata$substudy == "GSE216281"],
        baso$metadata$predicted_age[baso$metadata$substudy == "GSE216281"]
)

cor(
        baso$metadata$age_numeric[baso$metadata$substudy == "GSE216281"],
        baso$metadata$predicted_age[baso$metadata$substudy == "GSE216281"]
)


plot(
        baso$metadata$age_numeric[baso$metadata$substudy == "ROSMAP"],
        baso$metadata$predicted_age[baso$metadata$substudy == "ROSMAP"]
)

cor(
        baso$metadata$age_numeric[baso$metadata$substudy == "ROSMAP"],
        baso$metadata$predicted_age[baso$metadata$substudy == "ROSMAP"]
)


head(baso$metadata)

plot(baso$metadata$age_numeric,
     baso$metadata$predicted_age)


res_pred <- baso$metadata

res_pred_ctrls <- res_pred[res_pred$diagn_4BrainClck == "Control", ]

res_pred_filt <- res_pred[res_pred$RIN >= 6, ]

res_pred_filt <- res_pred_filt[(res_pred_filt$diagn_4BrainClck == "Control" | (res_pred_filt$diagn_4BrainClck == "Rest" & res_pred_filt$Braak >= 4 & !is.na(res_pred_filt$Braak))), ]

plot(res_pred_ctrls$age_numeric,
     res_pred_ctrls$predicted_age)

r2 <- 1 - sum((res_pred_ctrls$age_numeric[res_pred_ctrls$RIN >= 6] - res_pred_ctrls$predicted_age[res_pred_ctrls$RIN >= 6])^2)/sum((res_pred_ctrls$age_numeric[res_pred_ctrls$RIN >= 6] - mean(res_pred_ctrls$age_numeric[res_pred_ctrls$RIN >= 6]))^2)

ggplot(res_pred_filt, mapping = aes(x = age_numeric, y = predicted_age, col = diagn_4BrainClck)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")


ggplot(res_pred[res_pred$RIN >= 7, ], mapping = aes(x = age_numeric, y = predicted_age, col = diagn_4BrainClck)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")

cor(res_pred_ctrls$age_numeric[res_pred_ctrls$RIN >= 6],
    res_pred_ctrls$predicted_age[res_pred_ctrls$RIN >= 6])

# Compute acceleration
################################################################################

res_pred$accel_raw <- res_pred$predicted_age - res_pred$age_numeric

accel_model <- lm(accel_raw ~ age_numeric + substudy_batch, data = res_pred)
res_pred$acceleration <- residuals(accel_model)

# Try different df in the splines and select the best one with AIC
accel_mod <- lm(
        accel_raw ~ splines::ns(age_numeric, df=3) + substudy_batch,
        data = res_pred[res_pred$diagnosis == "control", ]
)

accel_mod_1 <- lm(
        accel_raw ~ splines::ns(age_numeric, df=4) + substudy_batch,
        data = res_pred[res_pred$diagnosis == "control", ]
)

accel_mod_2 <- lm(
        accel_raw ~ splines::ns(age_numeric, df=5) + substudy_batch,
        data = res_pred[res_pred$diagnosis == "control", ]
)

AIC(accel_model, accel_mod, accel_mod_1, accel_mod_2)
# df = 4 drops AIC by 9.5 points. Use it


#res_pred$acceleration_nonlin <- residuals(lm(accel_raw ~ splines::ns(age_numeric, df=3) + substudy_batch, data = res_pred))
res_pred$acceleration_nonlin <- res_pred$accel_raw - predict(
        accel_mod_1, 
        newdata = res_pred
)

accel_by_cohort_bxplot_lin <- ggplot(
        res_pred,
        mapping = aes(
                x = substudy_batch,
                y = acceleration
        )
) +
        geom_boxplot() + 
        theme_bw() +
        theme(
                axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 1,
                        size = 12
                )
        ) +
        xlab("Substudy") +
        ylab("Acceleration (linear)")

accel_by_cohort_bxplot_nonlin <- ggplot(
        res_pred, mapping = aes(x = substudy_batch, y = acceleration_nonlin)
) +
        geom_boxplot() + 
        theme_bw() +
        theme(
                axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 1,
                        size = 12
                )
        ) +
        xlab("Substudy") +
        ylab("Acceleration (nonlinear)")

ggsave(
        filename = file.path(out_dir, "accel_by_cohort_bxplot_lin.pdf"),
        plot = accel_by_cohort_bxplot_lin, height = 5, width = 6
)

ggsave(
        filename = file.path(out_dir, "accel_by_cohort_bxplot_nonlin.pdf"),
        plot = accel_by_cohort_bxplot_nonlin, height = 5, width = 6
)

ggplot(res_pred, mapping = aes(x = age_numeric, y = acceleration_nonlin, col = diagnosis)) +
        geom_point()

ggplot(res_pred, mapping = aes(x = age_numeric, y = acceleration, col = diagnosis)) +
        geom_point()

mean(res_pred$acceleration_nonlin[res_pred$diagnosis == "control"])
mean(res_pred$acceleration_nonlin)

dim(res_pred)

res_pred_braak <- res_pred[!is.na(res_pred$Braak), ]

summary(lm(acceleration_nonlin ~ Braak, data = res_pred_braak))


res_pred_braak <- res_pred_braak %>%
        group_by(individualID, Braak, diagnosis) %>%
        summarise(
                mean_acceleration = mean(acceleration, na.rm = TRUE),
                mean_acceleration_nonlin = mean(acceleration_nonlin, na.rm = TRUE),
                .groups = "drop"
        )

braak_acceleration <- ggplot(res_pred_braak, mapping = aes(x = Braak, y = mean_acceleration_nonlin)) +
        geom_point() +
        geom_smooth(method = "lm")

summary(lm(mean_acceleration_nonlin ~ Braak, data = res_pred_braak))

ggsave(
        filename = file.path(out_dir, "braak_acceleration.pdf"),
        plot = braak_acceleration, height = 3, width = 4
)

res_pred$diagnosis <- factor(
        res_pred$diagnosis,
        levels = sort(unique(res_pred$diagnosis),
                      decreasing = TRUE)
)

diseases <- table(res_pred$diagnosis)

diseases <- names(diseases[diseases >= 3])

diseases <- diseases[diseases != "control"]

# simplify tissues

gtex_mdat <- read.delim(
        "/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas//data/metadata/GTEx/phs000424.v11.pht002743.v10.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz",
        sep = '\t', comment.char = "#")

table(gtex_mdat[gtex_mdat$SMTSD == "Brain - Cortex", ]$SMUBRTRM)

table(gtex_mdat[gtex_mdat$SMTSD == "Brain - Caudate (basal ganglia)", ]$SMUBRTRM)
table(gtex_mdat[gtex_mdat$SMTSD == "Brain - Putamen (basal ganglia)", ]$SMUBRTRM)
table(gtex_mdat[gtex_mdat$SMTSD == "Brain - Nucleus accumbens (basal ganglia)", ]$SMUBRTRM)
table(gtex_mdat[gtex_mdat$SMTSD == "Brain - Amygdala", ]$SMUBRTRM)
table(gtex_mdat[gtex_mdat$SMTSD == "Brain - Frontal Cortex (BA9)", ]$SMUBRTRM)
table(gtex_mdat[gtex_mdat$SMTSD == "Brain - Hypothalamus", ]$SMUBRTRM)
table(res_pred$tissue[res_pred$substudy == "GTEx"])

table(res_pred$substudy[res_pred$tissue == "frontal cortex"])

table(res_pred$tissue[res_pred$substudy == "ROSMAP"])

tissue_dict <- data.frame(tissue = unique(res_pred$tissue))

tissue_dict$tissue_simplified <- c(
        "Limbic",
        "Cortex_ACC_PCC",
        "Basal_ganglia",
        "Frontal_cortex",
        "Frontal_cortex",
        "Limbic",
        "Diencephalon",
        "Basal_ganglia",
        "Basal_ganglia",
        "Spinal_cord",
        "Nigral",
        "Frontal_cortex",
        "Cortex_ACC_PCC",
        "Basal_ganglia",
        "Frontal_cortex", 
        "Temporal_cortex",
        "Temporal_cortex",
        "Limbic",
        "Frontal_cortex",
        "Frontal_cortex",
        "Frontal_cortex",
        "Cortex_mixed",
        "Cortex_mixed",
        "Frontal_cortex",
        "Limbic",
        "Frontal_cortex"
)

tissue_dict$tissue_corrected <- c(
        "Amygdala",
        "Anterior cingulate cortex",
        "Caudate nucleus",
        "Frontal cortex",
        "DLPFC",
        "Hippocampus",
        "Hypothalamus",
        "Nucleus accumbens",
        "Putamen",
        "Spinal cord",
        "Substantia nigra",
        "DLPFC",
        "Posterior cingulate cortex",
        "Caudate nucleus",
        "Frontal cortex",
        "Temporal cortex",
        "Superior temporal gyrus",
        "Parahippocampal gyrus",
        "Frontal pole",
        "Inferior frontal gyrus",
        "Prefrontal cortex",
        "Cortex mixed",
        "Cortex mixed",
        "DLPFC",
        "Hippocampus",
        "Frontal cortex"
)

tissue_dict$macro_region <- c(
        "subcortical",
        "frontal",
        "subcortical",
        "frontal",
        "frontal",
        "temporal",
        "subcortical",
        "subcortical",
        "subcortical",
        "spinal cord",
        "substantia nigra",
        "frontal",
        "PCC",
        "subcortical",
        "frontal",
        "temporal",
        "temporal",
        "temporal",
        "frontal",
        "frontal",
        "frontal",
        "cortex_mixed",
        "cortex_mixed",
        "temporal",
        "frontal",
        "frontal"
)

table(res_pred$tissue)

res_pred$tissue_simplified <- tissue_dict$tissue_simplified[
        match(res_pred$tissue, tissue_dict$tissue)
]

res_pred$tissue_corrected <- tissue_dict$tissue_corrected[
        match(res_pred$tissue, tissue_dict$tissue)
]

res_pred$macro_region <- tissue_dict$macro_region[
        match(res_pred$tissue, tissue_dict$tissue)
]

summary(lmerTest::lmer(
        acceleration_nonlin ~ Braak + macro_region + (1 | individualID),
        data = res_pred
))

ggplot(res_pred, mapping = aes(x = tissue, y = acceleration_nonlin)) +
        geom_boxplot() +
        theme_bw() +
        theme(
                axis.text.x = element_text(
                        angle = 90,
                        vjust = 0.5,
                        hjust = 1,
                        size = 8
                )
        ) +
        xlab("Tissue") +
        ylab("Acceleration (nonlinear)")

tissue_profile <- res_pred %>%
        group_by(tissue_simplified) %>%
        summarise(
                mean_acceleration = mean(acceleration_nonlin, na.rm = TRUE),
                sd_acceleration = sd(acceleration_nonlin, na.rm = TRUE),
                n = n()
        )

tissue_disease_profile <- res_pred %>%
        group_by(tissue_simplified, diagnosis) %>%
        summarise(
                mean_acc = mean(acceleration_nonlin, na.rm = TRUE),
                .groups = "drop"
        )

mat <- tissue_disease_profile %>%
        tidyr::pivot_wider(names_from = diagnosis, values_from = mean_acc)

mat <- as.data.frame(mat)

rownames(mat) <- mat$tissue
mat <- mat[, -1]

#mat <- t(
#        apply(
#                mat, 
#                1,
#                function(x) {
#                        x[is.na(x)] <- mean(x, na.rm = TRUE)
#                        x
#                }
#        )
#)

library(cluster)
#dist_mat <- dist(scale(mat), method = "euclidean", na.rm = TRUE)
d <- daisy(mat, metric = "gower")
d <- as.dist(1 - cor(t(mat), use = "pairwise.complete.obs"))
which(is.na(as.matrix(d)), arr.ind = TRUE)

table(res_pred$substudy[res_pred$tissue == "Cortex"])
table(res_pred$substudy[res_pred$tissue == "frontal pole"])
table(res_pred$substudy[res_pred$tissue == "inferior frontal gyrus"])
table(res_pred$substudy[res_pred$tissue == "parahippocampal gyrus"])
table(res_pred$substudy[res_pred$tissue == "prefrontal cortex"])
table(res_pred$substudy[res_pred$tissue == "superior temporal gyrus"])
table(res_pred$substudy[res_pred$tissue == "Caudate (basal ganglia)"])

d_mat <- as.matrix(d)
d_mat[is.na(d_mat)] <- 1
hc <- hclust(d)

plot(hc)

na_prop_tissue <- apply(mat, 1, function(x) sum(is.na(x))/length(x))

groups <- cutree(hc, k = 4)

na_prop_tissue[names(groups)[groups == 1]]
na_prop_tissue[names(groups)[groups == 2]]
na_prop_tissue[names(groups)[groups == 3]]
na_prop_tissue[names(groups)[groups == 4]]



table(res_pred$diagnosis, res_pred$tissue_simplified)

if (!require("effsize", quietly = TRUE)){
        install.packages("effsize")
}

library(effsize)

library(lme4)

disease_accel_stats <- list()
for(d in diseases){
        res_pred_d <- res_pred[res_pred$diagnosis %in% c("control", d), ]
        res_pred_d$diagnosis <- factor(res_pred_d$diagnosis, levels = c("control", d))
        
        fit <- lmerTest::lmer(
                acceleration_nonlin ~ diagnosis + macro_region + (1 | individualID),
                data = res_pred_d
        )
        res_pred_d$age_c <- scale(
                res_pred_d$age_numeric,
                center = TRUE,
                scale = FALSE
        )
        
        fit_2 <- lmerTest::lmer(
                acceleration_nonlin ~ diagnosis + age_c + age_c * diagnosis + macro_region + (1 | individualID),
                data = res_pred_d
        )
        res_fit <- coef(summary(fit))
        res_fit_2 <- coef(summary(fit_2))
        
        tobind = as.data.frame(matrix(res_fit[2, ], nrow = 1, dimnames = list(d, names(res_fit[2, ]))))
        tobind$disease <- d
        
        tobind <- tobind[, c(ncol(tobind), 1:(ncol(tobind) - 1))]
        colnames(tobind) <- c("disease", "effect_size", "std_err", "df", "t_value", "p_value")
        
        disease_accel_stats[[d]] <- tobind
}

disease_accel_stats <- do.call(rbind.data.frame, disease_accel_stats)

disease_accel_stats$p_adj <- p.adjust(disease_accel_stats$p_value)

disease_accel_stats

# Do Forest plot of the acceleration
library(ggplot2)

# From your per-disease lmer results
disease_accel_stats$lower <- disease_accel_stats$effect_size - 
        1.96 * disease_accel_stats$std_err
disease_accel_stats$upper <- disease_accel_stats$effect_size + 
        1.96 * disease_accel_stats$std_err
disease_accel_stats$sig <- disease_accel_stats$p_adj < 0.05

# Order by effect size
disease_accel_stats <- disease_accel_stats[
        order(disease_accel_stats$effect_size),
]
disease_accel_stats$disease <- factor(
        disease_accel_stats$disease, 
        levels = disease_accel_stats$disease
)

forest_plot_disease_accel <- ggplot(disease_accel_stats, 
       aes(x = effect_size, y = disease, color = sig)) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3) +
        geom_point(size = 3) +
        scale_color_manual(values = c("grey60", "firebrick"),
                           labels = c("n.s.", "FDR < 0.05")) +
        labs(x = "Age acceleration (years)", y = NULL,
             color = NULL,
             title = "Brain transcriptomic age acceleration by disease") +
        theme_bw() +
        theme(legend.position = "bottom")

ggsave(filename = file.path(out_dir, "forest_plot_disease_accel.pdf"),
       plot = forest_plot_disease_accel, width = 4, height = 4)

res_pred_average <- res_pred %>%
        group_by(individualID, diagnosis) %>%
        summarise(
                mean_acceleration = mean(acceleration, na.rm = TRUE),
                mean_acceleration_nonlin = mean(acceleration_nonlin, na.rm = TRUE),
                .groups = "drop"
        )

disease_accel_stats_average <- list()
for(d in diseases){
        res_pred_average_d <- res_pred_average[
                res_pred_average$diagnosis %in% c("control", d),
        ]
        res_pred_average_d$diagnosis <- factor(
                res_pred_average_d$diagnosis, levels = c("control", d)
        )
        
        ttst <- t.test(
                res_pred_average_d$mean_acceleration_nonlin[res_pred_average_d$diagnosis != "control"],
                res_pred_average_d$mean_acceleration_nonlin[res_pred_average_d$diagnosis == "control"]
        )
        eff_size <- ttst$estimate[1] - ttst$estimate[2]
        pval <- ttst$p.value
        cohen_D <- cohen.d(
                res_pred_average_d$mean_acceleration_nonlin[res_pred_average_d$diagnosis != "control"],
                res_pred_average_d$mean_acceleration_nonlin[res_pred_average_d$diagnosis == "control"]
        )
        tobind <- data.frame(
                disease = d,
                effect_size = eff_size,
                cohen_d = cohen_D$estimate,
                cohen_d_95CI = sprintf(
                        "[%s, %s]",
                        round(cohen_D$conf.int[1], digits = 3),
                        round(cohen_D$conf.int[2], digits = 3)
                ),
                p_value = pval
        )
        disease_accel_stats_average[[d]] <- tobind
}

#sum(res_pred_average$individualID %in% ind_mdat$individualID[ind_mdat$has_RNAseq & ind_mdat$has_WGS & !ind_mdat$used_in_clock_training])

disease_accel_stats_average <- do.call(rbind.data.frame, disease_accel_stats_average)

disease_accel_stats_average$p_adj <- p.adjust(disease_accel_stats_average$p_value)

disease_accel_stats_average

# Compute acceleration per individual after accounting for tissue heterogeneity
################################################################################

m_lin <- lmer(acceleration ~ macro_region + (1 | individualID), data = res_pred)
m_nonlin <- lmer(acceleration_nonlin ~ macro_region + (1 | individualID), data = res_pred)

indiv_acceleration_lin <- ranef(m_lin)$individualID
indiv_acceleration_nonlin <- ranef(m_nonlin)$individualID
dim(indiv_acceleration_lin)

all(rownames(indiv_acceleration_nonlin) == res_pred_average$individualID)

cor(indiv_acceleration_nonlin[, 1], res_pred_average$mean_acceleration_nonlin)

res_pred_average$lmer_acceleration <- indiv_acceleration_lin[, 1]
res_pred_average$lmer_acceleration_nonlin <- indiv_acceleration_nonlin[, 1]

table(res_pred_average$diagnosis)

disease_accel_stats_lmer_perind <- list()
for(d in diseases){
        res_pred_average_d <- res_pred_average[
                res_pred_average$diagnosis %in% c("control", d),
        ]
        res_pred_average_d$diagnosis <- factor(
                res_pred_average_d$diagnosis, levels = c("control", d)
        )
        
        ttst <- t.test(
                res_pred_average_d$lmer_acceleration_nonlin[res_pred_average_d$diagnosis != "control"],
                res_pred_average_d$lmer_acceleration_nonlin[res_pred_average_d$diagnosis == "control"]
        )
        eff_size <- ttst$estimate[1] - ttst$estimate[2]
        pval <- ttst$p.value
        cohen_D <- cohen.d(
                res_pred_average_d$lmer_acceleration_nonlin[res_pred_average_d$diagnosis != "control"],
                res_pred_average_d$lmer_acceleration_nonlin[res_pred_average_d$diagnosis == "control"]
        )
        tobind <- data.frame(
                disease = d,
                effect_size = eff_size,
                cohen_d = cohen_D$estimate,
                cohen_d_95CI = sprintf(
                        "[%s, %s]",
                        round(cohen_D$conf.int[1], digits = 3),
                        round(cohen_D$conf.int[2], digits = 3)
                ),
                p_value = pval
        )
        disease_accel_stats_lmer_perind[[d]] <- tobind
}


disease_accel_stats_lmer_perind <- do.call(rbind.data.frame, disease_accel_stats_lmer_perind)

disease_accel_stats_lmer_perind$p_adj <- p.adjust(disease_accel_stats_lmer_perind$p_value)

disease_accel_stats_lmer_perind

res_pred_average <- as.data.frame(res_pred_average)

write.csv(res_pred_average, file = file.path(out_dir, "age_acceleration_per_individual.csv"))

res_pred_average$diagnosis <- factor(res_pred_average$diagnosis, levels = c("control", levels(res_pred_average$diagnosis)[levels(res_pred_average$diagnosis) != "control"]))
summary(lm(lmer_acceleration_nonlin ~ diagnosis, data = res_pred_average))


# Fit a global LME model
################################################################################

res_pred$diagnosis <- factor(
        res_pred$diagnosis,
        levels = c(
                "control",
                levels(res_pred$diagnosis)[
                        levels(res_pred$diagnosis) != "control"
                ]
        )
)

lme_all_samps <- lmerTest::lmer(
        acceleration_nonlin ~ diagnosis + macro_region + (1 | individualID),
        data = res_pred
)

coef(summary(lme_all_samps))

if (!require("emmeans", quietly = TRUE)){
        install.packages("emmeans")
}
if (!require("tidydr", quietly = TRUE)){
        install.packages("tidydr")
}
library(tidydr)

emmeans_res <- emmeans(lme_all_samps, "diagnosis")

emmeans_res

#counts_basonorm_pca <- big_PCA(t(baso$norm_counts), nComps = 20, stand = TRUE)


lme_all_samps_inter <- lmerTest::lmer(
        acceleration_nonlin ~ diagnosis * macro_region + (1 | individualID),
        data = res_pred
)

coef(summary(lme_all_samps_inter))


anova(lme_all_samps, lme_all_samps_inter)



# Assess whether MCI is an intemrediate state of AD in ROSMAP

res_pred %>%
        filter(
                substudy == "ROSMAP",
                diagnosis %in% c("control", "mild cognitive impairment", "Alzheimer Disease")
        ) %>%
        ggplot(aes(x = age_numeric, y = acceleration_nonlin, col = diagnosis)) +
        geom_point(alpha = 0.4) +
        geom_smooth(method = "lm", se = FALSE) +
        facet_wrap(~ macro_region)

res_pred %>%
        filter(substudy == "ROSMAP",
               diagnosis %in% c("control", "mild cognitive impairment", "Alzheimer Disease")) %>%
        ggplot(aes(x = acceleration_nonlin, fill = diagnosis)) +
        geom_density(alpha = 0.4) +
        facet_wrap(~ macro_region)




#pca_substudy <- pcBiplot(counts_basonorm_pca, colBy = "substudy", biPlot = F,
#                         colROSMAPBatch = F)


# Check after filtering ambiguous ages and bad RIN how many individuals
# with WGS not used in training do we end up with


sum(unique(res_pred$individualID) %in% ind_mdat$individualID[ind_mdat$has_RNAseq & ind_mdat$has_WGS & !ind_mdat$used_in_clock_training])


traintest_clock <- read_table_fast("/Users/guillem.santamaria/Downloads/train_test.csv")
clock_mdat <- read_table_fast("/Users/guillem.santamaria/Downloads/merged_metdat.csv")
table(traintest_clock$substudy)

hist(traintest_clock$ageDeath[traintest_clock$train_test_split == "train"])

hist(clock_mdat$ageDeath[clock_mdat$substudy == "brainSeq_pI"])

# Stratify per region

regs <- unique(res_pred$macro_region)

disease_accel_stats_strat_reg <- list()
for(i in seq_along(regs)){
        reg <- regs[i]
        res_pred_r <- res_pred[res_pred$macro_region == reg, ]
        diseases <- as.character(unique(res_pred_r$diagnosis))
        res_pred_r$diagnosis <- factor(
                res_pred_r$diagnosis,
                levels = c("control", diseases[diseases != "control"])
        )
        samples_per_id <- table(res_pred_r$individualID)
        can_fit_mixed <- any(samples_per_id >= 2) && length(unique(res_pred_r$individualID)) > 1
        if (can_fit_mixed){
                fit <- lmerTest::lmer(
                        acceleration_nonlin ~ diagnosis + (1 | individualID),
                        data = res_pred_r
                )
                res_fit <- as.data.frame(coef(summary(fit)))
                res_fit$model_type <- "lmer"
        }else{
                fit <- lm(
                        acceleration_nonlin ~ diagnosis,
                        data = res_pred_r
                )
                res_fit <- as.data.frame(coef(summary(fit)))
                res_fit$model_type <- "lm"
        }
        
        res_fit <- res_fit[2:nrow(res_fit), ]
        
        res_fit$disease <- gsub("diagnosis", "", rownames(res_fit))
        res_fit$tissue <- reg
        
        res_fit <- res_fit[, c(ncol(res_fit), ncol(res_fit) - 1, 1:(ncol(res_fit) - 2))]
        res_fit <- res_fit[, colnames(res_fit) != "df"]
        
        colnames(res_fit) <- c(
                "tissue",
                "disease",
                "effect_size",
                "std_err",
                "t_value",
                "p_value",
                "model_type"
        )
        
        disease_accel_stats_strat_reg[[reg]] <- res_fit
}

disease_accel_stats_strat_reg <- do.call(rbind.data.frame, disease_accel_stats_strat_reg)

disease_accel_stats_strat_reg[order(disease_accel_stats_strat_reg$disease), ]


library(tidyr)

# Prepare matrix from stratified results
heat_dat <- disease_accel_stats_strat_reg %>%
        filter(disease %in% c(
                "Alzheimer Disease", "Parkinson's disease",
                "progressive supranuclear palsy",
                "Cerebrovascular disease"#,
                #"mild cognitive impairment"
        )) %>%
        mutate(
                sig_label = case_when(
                        p_value < 0.001 ~ "***",
                        p_value < 0.01  ~ "**",
                        p_value < 0.05  ~ "*",
                        TRUE            ~ ""
                ),
                # Cap effect sizes for color scale
                effect_capped = pmax(pmin(effect_size, 6), -6)
        )

heatmap_age_accel_disease_brain_reg <- ggplot(heat_dat, 
       aes(x = tissue, y = disease, fill = effect_capped)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sig_label), size = 4, vjust = 0.75) +
        scale_fill_gradient2(
                low = "steelblue", mid = "white", high = "firebrick",
                midpoint = 0, name = "Effect size\n(years)",
                limits = c(-6, 6)
        ) +
        labs(x = "Brain region", y = NULL#,
             #title = "Age acceleration by disease and brain region"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(out_dir, "disease_reg_accel_hm.pdf"),
       plot = heatmap_age_accel_disease_brain_reg, width = 6, height = 2.5)



res_pred %>%
        filter(diagnosis %in% c(
                "control", "Alzheimer Disease", 
                "Parkinson's disease"
        )) %>%
        mutate(diagnosis = factor(diagnosis, 
                                  levels = c("control", "Alzheimer Disease", 
                                             "Parkinson's disease"))) %>%
        ggplot(aes(x = age_numeric, y = acceleration_nonlin,
                   color = diagnosis)) +
        geom_point(alpha = 0.2, size = 0.8) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("grey50", "firebrick", 
                                      "steelblue")) +
        labs(x = "Chronological age", 
             y = "Age acceleration (years)",
             color = NULL) +
        theme_bw() +
        theme(legend.position = "bottom")

# Check if subcortical is significantly higher in age in ROSMAP or in 
# all the cohort
res_pred %>%
        filter(substudy == "ROSMAP",
               macro_region == "subcortical",
               diagnosis %in% c("control", "mild cognitive impairment")) %>%
        group_by(diagnosis) %>%
        summarise(
                mean_age = mean(age_numeric),
                mean_RIN = mean(RIN, na.rm=TRUE),
                n = n()
        )

# Check if subcortical is higher in age in MCI than in controls
res_pred %>%
        filter(macro_region == "subcortical",
               diagnosis %in% c("control", "mild cognitive impairment")) %>%
        group_by(diagnosis) %>%
        summarise(
                mean_age = mean(age_numeric),
                mean_RIN = mean(RIN, na.rm=TRUE),
                n = n()
        )


# Do subcortical analysis only with rosmap to see if results prevail
res_pred_r_rosmap <- res_pred[res_pred$macro_region == "subcortical" & res_pred$substudy == "ROSMAP", ]
diseases <- as.character(unique(res_pred_r_rosmap$diagnosis))
res_pred_r_rosmap$diagnosis <- factor(
        res_pred_r_rosmap$diagnosis,
        levels = c("control", diseases[diseases != "control"])
)

fit <- lmerTest::lmer(
        acceleration_nonlin ~ diagnosis + (1 | individualID),
        data = res_pred_r_rosmap
)
res_fit <- as.data.frame(coef(summary(fit)))
res_fit

# When focusing only on ROSMAP MCI is not significant anymore in subcortical
# Might it be because of the controls being significantly younger.

table(res_pred$diagnosis[res_pred$substudy == "ROSMAP"])

mean(res_pred$age_numeric[res_pred$substudy == "ROSMAP" & res_pred$diagnosis == "control"])

table(res_pred$diagnosis[res_pred$substudy == "GTEx"])

mean(res_pred$age_numeric[res_pred$substudy == "GTEx" & res_pred$diagnosis == "control"])

table(res_pred$diagnosis[res_pred$substudy == "MSBB"])

mean(res_pred$age_numeric[res_pred$substudy == "MSBB" & res_pred$diagnosis == "control"])

table(res_pred$diagnosis[res_pred$substudy == "Mayo"])

mean(res_pred$age_numeric[res_pred$substudy == "Mayo" & res_pred$diagnosis == "control"])

table(res_pred$diagnosis[res_pred$substudy == "LBP"])

mean(res_pred$age_numeric[res_pred$substudy == "LBP" & res_pred$diagnosis == "control"])

table(res_pred$macro_region[res_pred$diagnosis == "control"],
      res_pred$substudy[res_pred$diagnosis == "control"])

table(res_pred$macro_region[res_pred$diagnosis == "mild cognitive impairment"],
      res_pred$substudy[res_pred$diagnosis == "mild cognitive impairment"])



library(ggplot2)
library(dplyr)

# Compute mean age and n per disease/region/cohort
age_summary <- res_pred %>%
        group_by(macro_region, diagnosis, substudy) %>%
        summarise(
                mean_age = mean(age_numeric, na.rm=TRUE),
                n = n(),
                .groups = "drop"
        )

# Flag combinations with few samples
age_summary$label <- ifelse(age_summary$n < 20, 
                            as.character(age_summary$n), "")

# Heatmap of mean age per disease x region
ggplot(age_summary, 
       aes(x = macro_region, y = diagnosis, fill = mean_age)) +
        geom_tile(color = "white") +
        geom_text(aes(label = n), size = 2.5) +
        scale_fill_gradient2(
                low  = "steelblue", 
                mid  = "white", 
                high = "firebrick",
                midpoint = median(res_pred$age_numeric, na.rm=TRUE),
                name = "Mean age"
        ) +
        facet_wrap(~ substudy) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(title = "Mean age by disease x region x cohort",
             subtitle = "Numbers = sample size, color = mean age",
             x = "Region", y = "Diagnosis")

# For each disease, plot age distributions vs control
res_pred %>%
        filter(diagnosis %in% c("control", diseases)) %>%
        mutate(is_control = diagnosis == "control") %>%
        ggplot(aes(x = age_numeric, fill = is_control)) +
        geom_density(alpha = 0.4) +
        scale_fill_manual(values = c("firebrick", "steelblue"),
                          labels = c("Disease", "Control")) +
        facet_wrap(~ diagnosis, scales = "free_y") +
        theme_minimal() +
        labs(title = "Age distribution overlap: cases vs controls",
             x = "Age", fill = "")

table(res_pred$macro_region,
      res_pred$diagnosis)

print(age_summary, n = 70)

################################################################################
# Try per-cohort residualization
################################################################################

cohorts <- unique(res_pred$substudy)

res_pred$acceleration_nonlin_v3 <- NA

for(coh in cohorts){
        # Controls in this cohort
        ctrl_idx <- res_pred$substudy == coh & res_pred$diagnosis == "control"
        all_idx  <- res_pred$substudy == coh
        
        # Need enough controls to fit spline
        n_ctrl <- sum(ctrl_idx)
        df_use <- ifelse(n_ctrl >= 100, 4, ifelse(n_ctrl >= 50, 3, 2))
        
        if(n_ctrl < 20){
                # Too few controls — use linear
                mod <- lm(accel_raw ~ age_numeric, data = res_pred[ctrl_idx,])
        } else {
                if (length(unique(res_pred$substudy_batch[ctrl_idx])) == 1){
                        mod <- lm(
                                accel_raw ~ splines::ns(age_numeric, df=df_use),
                                data = res_pred[ctrl_idx,]
                        )
                }else{
                        mod <- lm(
                                accel_raw ~ splines::ns(age_numeric, df=df_use) + substudy_batch,
                                data = res_pred[ctrl_idx,]
                        )
                }
        }
        
        res_pred$acceleration_nonlin_v3[all_idx] <- 
                res_pred$accel_raw[all_idx] - predict(mod, newdata=res_pred[all_idx,])
}

# Verify
res_pred %>%
        group_by(substudy, diagnosis) %>%
        summarise(mean_accel = mean(acceleration_nonlin_v3, na.rm=TRUE),
                  n = n(), .groups="drop") %>%
        filter(diagnosis == "control")

# Visual check per cohort
res_pred %>%
        filter(diagnosis == "control") %>%
        ggplot(aes(x = age_numeric, y = accel_raw)) +
        geom_point(alpha=0.3) +
        geom_smooth(method="lm", 
                    formula = y ~ splines::ns(x, df=4),
                    se=TRUE) +
        facet_wrap(~ substudy, scales="free") +
        labs(title="Per-cohort aging trajectories in controls")

lme_all_samps_accel_by_cohort <- lmerTest::lmer(
        acceleration_nonlin_v3 ~ diagnosis + macro_region + (1 | individualID),
        data = res_pred
)

coef(summary(lme_all_samps_accel_by_cohort))


# Within LBP: PD vs controls
res_pred %>%
        filter(substudy == "LBP",
               diagnosis %in% c("control", "Parkinson's disease")) %>%
        group_by(diagnosis) %>%
        summarise(
                mean_age = mean(age_numeric),
                mean_accel = mean(acceleration_nonlin_v3, na.rm=TRUE),
                n = n()
        )

# Within GTEx: PD vs controls  
res_pred %>%
        filter(substudy == "GTEx",
               diagnosis %in% c("control", "Parkinson's disease")) %>%
        group_by(diagnosis) %>%
        summarise(
                mean_age = mean(age_numeric),
                mean_accel = mean(acceleration_nonlin_v3, na.rm=TRUE),
                n = n()
        )


regs <- unique(res_pred$macro_region)

disease_accel_stats_strat_reg_per_cohort <- list()
for(i in seq_along(regs)){
        reg <- regs[i]
        res_pred_r <- res_pred[res_pred$macro_region == reg, ]
        diseases <- as.character(unique(res_pred_r$diagnosis))
        res_pred_r$diagnosis <- factor(
                res_pred_r$diagnosis,
                levels = c("control", diseases[diseases != "control"])
        )
        samples_per_id <- table(res_pred_r$individualID)
        can_fit_mixed <- any(samples_per_id >= 2) && length(unique(res_pred_r$individualID)) > 1
        if (can_fit_mixed){
                fit <- lmerTest::lmer(
                        acceleration_nonlin_v3 ~ diagnosis + (1 | individualID),
                        data = res_pred_r
                )
                res_fit <- as.data.frame(coef(summary(fit)))
                res_fit$model_type <- "lmer"
        }else{
                fit <- lm(
                        acceleration_nonlin_v3 ~ diagnosis,
                        data = res_pred_r
                )
                res_fit <- as.data.frame(coef(summary(fit)))
                res_fit$model_type <- "lm"
        }
        
        res_fit <- res_fit[2:nrow(res_fit), ]
        
        res_fit$disease <- gsub("diagnosis", "", rownames(res_fit))
        res_fit$tissue <- reg
        
        res_fit <- res_fit[, c(ncol(res_fit), ncol(res_fit) - 1, 1:(ncol(res_fit) - 2))]
        res_fit <- res_fit[, colnames(res_fit) != "df"]
        
        colnames(res_fit) <- c(
                "tissue",
                "disease",
                "effect_size",
                "std_err",
                "t_value",
                "p_value",
                "model_type"
        )
        
        disease_accel_stats_strat_reg_per_cohort[[reg]] <- res_fit
}

plot(res_pred$age_numeric[res_pred$substudy == "GSE216281" & res_pred$diagnosis == "control"],
     res_pred$predicted_age[res_pred$substudy == "GSE216281" & res_pred$diagnosis == "control"])

cor(res_pred$age_numeric[res_pred$substudy == "GSE216281" & res_pred$diagnosis == "control"],
res_pred$predicted_age[res_pred$substudy == "GSE216281" & res_pred$diagnosis == "control"])

table(res_pred$substudy)

disease_accel_stats_strat_reg_per_cohort <- do.call(rbind.data.frame, disease_accel_stats_strat_reg_per_cohort)

disease_accel_stats_strat_reg_per_cohort[order(disease_accel_stats_strat_reg_per_cohort$disease), ]

# Add the used_in_training column
traintest_clock <- read.csv("/Users/guillem.santamaria/Documents/postdoc/comput/clock_gwas/data/metadata/brain_clock/train_test.csv")

res_pred$used_in_clock_training <- res_pred$specimenID %in% traintest_clock$specimenID[traintest_clock$train_test_split == "train"]

write.csv(res_pred, file = file.path(out_dir, "age_prediction_results.csv"))
