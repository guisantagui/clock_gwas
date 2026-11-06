################################################################################
# Brain clock: build predictive model of surrogate variables from integrated.  #
# database, with the aim of using it with new samples user might give, to      #
# estimate surrogate variables in new data and be able to adjust for new batch #
# effects.                                                                     #
################################################################################
library(dplyr)

svObj <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/wTBI_lincsLndmrk/combined_TPMs_wTBI_lincsLndmrk_preproc_log10_noCerebell_onlyAge_svobj.rds"
preprocDat <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/wTBI_lincsLndmrk/combined_TPMs_wTBI_lincsLndmrk_preproc_log10_noCerebell.csv"
rawDat <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_counts_wTBI_lincsLndmrk.csv"


svObj <- readRDS(svObj)
svObj$sv

readCsvFst <- function(path){
        df <- data.frame(data.table::fread(path))
        rownames(df) <- df$V1
        df <- df[, colnames(df) != "V1"]
        return(df)
}

dat <- readCsvFst(preprocDat)

sv_mat <- svObj$sv

colnames(sv_mat) <- paste0("SV", as.character(1:ncol(sv_mat)))

#modMat <- 

dim(sv_mat)
dim(dat)
       
surrogate_model <- lm(sv_mat ~ ., data = as.data.frame(dat))
summary(surrogate_model)


# Given a counts matrix with genes in rows and samples in columns, and a 
# geneInfo dataframe containing gene sizes, returns a matrix of TPMs.
counts2TPMs <- function(countDF,
                        geneInfo,
                        IDsIn = "ensembl_gene_id_version_orig",
                        genesInCols = F){
        if(genesInCols){
                inDF <- t(countDF)
        }else{
                inDF <- countDF
        }
        # Filter out genes that are not in geneInfo DF (don't have size)
        inDF <- inDF[rownames(inDF) %in% geneInfo[, IDsIn], ]
        # Get a vector of gene sizes
        sizeVec <- geneInfo$size[match(rownames(inDF), geneInfo[, IDsIn])]
        # Obtain TPMs
        preTPMs <- inDF/sizeVec
        TPMs <- t(t(preTPMs) * 1e6 / colSums(preTPMs))
        if(genesInCols){
                TPMs <- t(TPMs)
        }
        return(TPMs)
}


# Predict the age in the test dataset
predictAge <- function(model,
                       predExpr
){
        predExpr <- predExpr[model@parameters$x,]
        predAges <- c()
        maxiter <- ceiling(ncol(predExpr)/1000)
        for(i in 1:maxiter){
                #i <- 1
                from <- ((i-1)*1000+1)
                to <- min(i*1000,ncol(predExpr))
                
                inp_h2o <- as.h2o(t(predExpr[,from:to]))
                
                pred <- h2o.predict(model,inp_h2o)
                pred <- as.data.frame(pred)$predict
                names(pred) <- colnames(predExpr)[from:to]
                
                predAges <- c(predAges,pred)
        }
        return(predAges)
}

# Preprocess lincs samples using surrogate model
dat_raw <- readCsvFst(rawDat)
dat_svaAdjFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/preprocessing/wTBI_lincsLndmrk/combined_TPMs_wTBI_lincsLndmrk_preproc_log10_noCerebell_onlyAge_svaAdj.csv"
dat_svaAdj <- readCsvFst(dat_svaAdjFile)
geneInfoFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_geneInfo_wTBI_lincsLndmrk.csv"
geneInfo <- read.csv(geneInfoFile, row.names = 1)


lincs_files1 <- list.files("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/lincs/divided_by_plate/exp_by_plate_1", full.names = T)
lincs_ctrls <- list.files("/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/data/perturbation/lincs/divided_by_plate/all_lincs_control_by_plate", full.names = T)

xpr <- readCsvFst(lincs_files1[1])

# Keep only features in training dataset, and transpose
xpr <- t(xpr[rownames(xpr) %in% colnames(dat_raw), ])

# Impute variables that might not be in the input data
library(impute)
genesToImpute <- colnames(dat_raw)[!colnames(dat_raw) %in% colnames(xpr)]

print(sprintf("There are %s genes that need to be imputed. This represents the %s %% of the genes in the training data.",
              length(genesToImpute),
              round(length(genesToImpute)/ncol(dat_raw) * 100, digits = 2)))

print("These are the genes that need to be imputed:")
print(genesToImpute)

imputeTo <- "knn"
if(imputeTo == "median"){
        impVals <- apply(dat_raw[, genesToImpute], 2, median)
        
        impCols <- matrix(rep(impVals, each = nrow(xpr)),
                          ncol = length(genesToImpute),
                          nrow = nrow(xpr),
                          dimnames = list(rownames(xpr), genesToImpute))
        
        # Merge and reorder to match order of variables of dat
        xpr_imp <- cbind(xpr, impCols)
        xpr_imp <- xpr_imp[, match(colnames(dat_raw), colnames(xpr_imp))]
        xpr_imp <- as.data.frame(xpr_imp)
}else if(imputeTo == "median"){
        impVals <- apply(dat_raw[, genesToImpute], 2, mean)
        
        impCols <- matrix(rep(impVals, each = nrow(xpr)),
                          ncol = length(genesToImpute),
                          nrow = nrow(xpr),
                          dimnames = list(rownames(xpr), genesToImpute))
        
        # Merge and reorder to match order of variables of dat
        xpr_imp <- cbind(xpr, impCols)
        xpr_imp <- xpr_imp[, match(colnames(dat_raw), colnames(xpr_imp))]
        xpr_imp <- as.data.frame(xpr_imp)
}else if(imputeTo == "zero"){
        impVals <- rep(0, length(genesToImpute))
        names(impVals) <- genesToImpute
        
        impCols <- matrix(rep(impVals, each = nrow(xpr)),
                          ncol = length(genesToImpute),
                          nrow = nrow(xpr),
                          dimnames = list(rownames(xpr), genesToImpute))
        
        # Merge and reorder to match order of variables of dat
        xpr_imp <- cbind(xpr, impCols)
        xpr_imp <- xpr_imp[, match(colnames(dat_raw), colnames(xpr_imp))]
        xpr_imp <- as.data.frame(xpr_imp)
}else if(imputeTo == "knn"){
        impCols <- matrix(rep(NA, nrow(xpr) * length(genesToImpute)),
                          ncol = length(genesToImpute),
                          nrow = nrow(xpr),
                          dimnames = list(rownames(xpr), genesToImpute))
        xpr_toImp <- cbind(xpr, impCols)
        xpr_toImp <- xpr_toImp[, match(colnames(dat_raw),
                                       colnames(xpr_toImp))]
        xpr_toImp <- data.frame(xpr_toImp)
        xprToImpWRaw <- rbind.data.frame(dat_raw, xpr_toImp)
        xprToImpWRaw_imputed <- impute.knn(t(xprToImpWRaw))
        xpr_imp <- t(xprToImpWRaw_imputed$data)
        xpr_imp <- xpr_imp[rownames(xpr_imp) %in% rownames(xpr_toImp), ]
        xpr_imp <- as.data.frame(xpr_imp)
}




# Convert counts to TPMs and log10 transform
xpr_imp_tpms <- counts2TPMs(xpr_imp,
                            geneInfo = geneInfo,
                            IDsIn = "ensembl_gene_id",
                            genesInCols = T)
xpr_imp_tpms_log10 <- data.frame(log10(xpr_imp_tpms + 0.00001))

# Predict surrogate variables and regress them out
sv_pred <- predict(surrogate_model, xpr_imp_tpms_log10)
xpr_imp_tpms_log10_svaAdj <- t(removeBatchEffect(t(xpr_imp_tpms_log10), covariates = sv_pred))

plot(density(as.matrix(xpr_imp_tpms_log10_svaAdj)))
plot(density(as.matrix(dat)))
plot(density(as.matrix(dat_svaAdj)))

mean(as.matrix(xpr_imp_tpms_log10_svaAdj))
sd(as.matrix(xpr_imp_tpms_log10_svaAdj))
mean(as.matrix(dat))
sd(as.matrix(dat))
mean(as.matrix(dat_svaAdj))
sd(as.matrix(dat_svaAdj))

pred_valid <- predictAge(mod, t(xpr_imp_tpms_log10_svaAdj))
max(pred_valid)
min(pred_valid)


# Now let's check brain clock predictions to see if we get predictions that
# fall within 0-1 range.
modFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/models/modAllGenes_lincsLndmrk_wTBI/modFuncsAlpha0.5/GLM_model_R_1720454074420_1"
library(h2o)

conn <- h2o.init(max_mem_size="16G")

mod <- h2o.loadModel(modFile)

pred_valid <- predictAge(mod, t(xpr_imp_tpms_log10_svaAdj))
max(pred_valid)
min(pred_valid)



# Try to apply SVA on new data coupled with database to see if results are 
# better.
library(sva)
xpr




xpr[1:10, 1:10]

sum(colnames(dat) %in% rownames(xpr))


ctr <- readCsvFst(lincs_ctrls[40])
head(colnames(ctr))
head(colnames(xpr))
table(sapply(colnames(xpr), function(x) paste(strsplit(x, split = "_", fixed = T)[[1]][6],
                                        strsplit(x, split = "_", fixed = T)[[1]][7],
                                        sep = "_")))

drug <- names(table(sapply(rownames(xpr_imp_tpms_log10), function(x) strsplit(x, split = "_")[[1]][8])))[33]

xpr_imp_tpms_log10_drug <- xpr_imp_tpms_log10[grep(drug, rownames(xpr_imp_tpms_log10)), ]
dat_xpr_imp_tpms_log10_drug <- rbind.data.frame(dat, xpr_imp_tpms_log10_drug)

metDatFile <- "/Users/guillem.santamaria/Documents/postdoc/comput/brain_clock/results/parsed/combined_metDat_wTBI.csv"
metDat <- read.csv(metDatFile, row.names = 1)

input <- dat_xpr_imp_tpms_log10_drug
pheno <- data.frame(age_death = metDat$ageDeath[match(rownames(input),
                                                      make.names(metDat$specimenID))],
                    pmi = metDat$pmi[match(rownames(input),
                                           make.names(metDat$specimenID))],
                    sex = metDat$sex[match(rownames(input),
                                           make.names(metDat$specimenID))],
                    substudy = metDat$substudy[match(rownames(input),
                                                     make.names(metDat$specimenID))])

pheno$age_death[is.na(pheno$age_death)] <- median(pheno$age_death[!is.na(pheno$age_death)])

mod_onlyAge = model.matrix(~age_death, data=pheno)
mod0_onlyAge = model.matrix(~1, data=pheno)

edata <- t(input)

n.sv = num.sv(edata, mod_onlyAge, method="leek")

svobj = sva(edata,
            mod_onlyAge,
            mod0_onlyAge,
            n.sv=n.sv)

mod4adj <- model.matrix(~age_death,
                        data = data.frame(age_death = mod_onlyAge[, "age_death"]))

edata_adj <- removeBatchEffect(edata, covariates = svobj$sv, design = mod4adj)

edata_adj <- t(edata_adj)


xpr_imp_tpms_log10_drug_adj <- edata_adj[rownames(edata_adj) %in% rownames(xpr_imp_tpms_log10_drug), ]

pred_xpr_imp_tpms_log10_drug_adj <- predictAge(mod, t(xpr_imp_tpms_log10_drug_adj))

# It seems that this approach works, at least with the first dataset
# and with these genes.

plot(density(pred_xpr_imp_tpms_log10_drug_adj))

min(pred_xpr_imp_tpms_log10_drug_adj)
max(pred_xpr_imp_tpms_log10_drug_adj)
dat_xpr_imp_tpms_log10_drug

h2o.shutdown(prompt = F)
