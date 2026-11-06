################################################################################
# Brain clock cross-fitting, step 2: preprocess, train and predict one fold.   #
#                                                                              #
# For outer fold k of repeat r:                                                #
#   train = eligible controls (always_train + cv individuals not in fold k)    #
#   test  = cv individuals in fold k (cases and controls) + external samples   #
# --fold 0 trains the final clock on all eligible controls and predicts every #
# sample (training samples are then in-sample and flagged as such).           #
#                                                                              #
# Everything is fit on the training samples only and then applied to the     #
# held-out ones sample by sample, without their ages:                         #
#   1. gene filter (zeros / low expression per age bin)                        #
#   2. log2(counts + 1) + quantile normalisation to the training reference    #
#   3. SVA on training (mod = ~age + region, mod0 = ~region) and per-sample    #
#      projection of SVs (see sva_fit() in cv_utils.R)                         #
#   4. elastic net, alpha and lambda tuned by inner CV grouped by individual   #
#      (inner folds reuse this outer fold's preprocessing).                    #
################################################################################

suppressPackageStartupMessages({
        library(argparser)
        library(glmnet)
})
source(file.path((function() {
        a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
        if (length(f)) dirname(normalizePath(f[1])) else getwd()
})(), "cv_utils.R"))

p <- arg_parser("Train and predict one cross-fitting fold.")
p <- add_argument(p, "--fold", help = "Outer fold (0 = final model on all controls).", type = "integer")
p <- add_argument(p, "--rep", help = "Repeat (column fold_r<rep> of folds.csv).", default = 1L)
p <- add_argument(p, "--folds", help = "folds.csv from 01_make_folds.R.", default = "results/clock_cv/folds.csv")
p <- add_argument(p, "--counts", help = "Merged counts (genes x samples).",
                  default = "results/counts_parsed/counts_merged.csv")
p <- add_argument(p, "--samp_mdat", help = "Sample metadata.",
                  default = "results/counts_parsed/sample_metadata_merged.csv")
p <- add_argument(p, "--ind_mdat", help = "Individual metadata.",
                  default = "results/counts_parsed/individual_metadata_merged_rnaseq_filt.csv")
p <- add_argument(p, "--region_map", help = "tissue -> region table.",
                  default = file.path(script_dir(), "region_map.csv"))
p <- add_argument(p, "--gene_universe", help = "'lincs', 'all' or a file with one Ensembl ID per line.",
                  default = "lincs")
p <- add_argument(p, "--lincs_genes", help = "LINCS gene list.",
                  default = "../brain_clock/data/LINCS/lincs_genes.txt")
p <- add_argument(p, "--exclude_regions", help = "Regions left out (must match 01_make_folds.R).",
                  default = "cerebellum")
p <- add_argument(p, "--external", help = "Substudies never used for training.", default = "brainSeq_pI")
p <- add_argument(p, "--min_train_age", help = "Minimum age of training controls.", default = 20)
p <- add_argument(p, "--min_rin", help = "Minimum RIN of training samples (NA RIN allowed).", default = 6)
p <- add_argument(p, "--sva_mod", help = "Full model for SVA.", default = "~ age + region")
p <- add_argument(p, "--sva_mod0", help = "Null model for SVA.", default = "~ region")
p <- add_argument(p, "--n_sv", help = "Number of SVs, or 'leek' (num.sv) or 'rmt' (EstDimRMT + 1).", default = "leek")
p <- add_argument(p, "--sva_method", help = "'smartsva' (SmartSVA::smartsva.cpp, as in the paper) or 'sva' (sva::sva, slower).", default = "smartsva")
p <- add_argument(p, "--alphas", help = "Comma-separated elastic net alphas.", default = "0.05,0.1,0.25,0.5,1")
p <- add_argument(p, "--n_inner", help = "Inner CV folds.", default = 10L)
p <- add_argument(p, "--lambda", help = "'lambda.min' or 'lambda.1se'.", default = "lambda.min")
p <- add_argument(p, "--region_covariates", help = "Add unpenalised region indicators to the model.",
                  flag = TRUE)
p <- add_argument(p, "--cores", help = "Cores for cv.glmnet.", default = 1L)
p <- add_argument(p, "--debug_n_genes", help = "Debug: random subset of genes.", type = "integer", default = NA)
p <- add_argument(p, "--debug_n_train", help = "Debug: random subset of training individuals.", type = "integer", default = NA)
p <- add_argument(p, "--outdir", help = "Output root.", default = "results/clock_cv")
a <- parse_args(p)

t0 <- Sys.time()
k <- a$fold; r <- a$rep
seed <- 1000L * r + k
out <- file.path(a$outdir, sprintf("rep%d", r), sprintf("fold%d", k))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_f <- file.path(out, "log.txt")
logmsg <- function(...) {
        m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
        message(m); cat(m, "\n", file = log_f, append = TRUE)
}
cat("", file = log_f)
arg_show <- a[setdiff(names(a), c("", "help", "opts"))]
logmsg("rep ", r, " fold ", k, " | args: ",
       paste(names(arg_show), sapply(arg_show, paste, collapse = ","), sep = "=", collapse = " "))

# Samples and roles
################################################################################
exclude_regions <- strsplit(a$exclude_regions, ",")[[1]]
external <- strsplit(a$external, ",")[[1]]

s <- load_metadata(a$samp_mdat, a$ind_mdat, a$region_map)
s <- s[!s$region %in% exclude_regions, ]
folds <- read.csv(a$folds, stringsAsFactors = FALSE)
fold_col <- paste0("fold_r", r)
if (!fold_col %in% colnames(folds)) stop(fold_col, " not in ", a$folds)
m <- match(s$individualID, folds$individualID)
if (anyNA(m)) stop("Samples whose individual is not in folds.csv; rerun 01_make_folds.R with the same settings.")
s$role <- folds$role[m]
s$fold <- folds[[fold_col]][m]
s$group <- ifelse(is.na(folds$group[m]), s$individualID, folds$group[m])
s$train_eligible <- flag_train_eligible(s, a$min_train_age, a$min_rin, external)

if (k == 0) {
        is_train <- s$train_eligible & s$role != "external"
        is_test <- rep(TRUE, nrow(s))
} else {
        if (!k %in% s$fold) stop("Fold ", k, " not found in ", fold_col)
        is_train <- s$train_eligible &
                (s$role == "always_train" | (s$role == "cv" & s$fold != k))
        is_test <- (s$role == "cv" & s$fold == k) | s$role == "external"
        if (length(intersect(s$group[is_train], s$group[is_test])) > 0) {
                stop("Leakage: a group has samples in both train and test.")
        }
}

if (!is.na(a$debug_n_train)) {
        set.seed(seed)
        keep_ind <- sample(unique(s$individualID[is_train]),
                           min(a$debug_n_train, length(unique(s$individualID[is_train]))))
        is_train <- is_train & s$individualID %in% keep_ind
}

s_train <- s[is_train, ]
s_test <- s[is_test, ]
logmsg(sprintf("train: %d samples / %d individuals | test: %d samples / %d individuals",
               nrow(s_train), length(unique(s_train$individualID)),
               nrow(s_test), length(unique(s_test$individualID))))

# Counts
################################################################################
genes_univ <- read_gene_universe(a$gene_universe, a$lincs_genes)
cache_f <- file.path(a$outdir, "cache",
                     sprintf("counts_%s.rds", gsub("[^A-Za-z0-9]", "_", basename(a$gene_universe))))
counts <- load_counts(a$counts, genes_univ, cache_f)
need <- unique(c(s_train$sample_id, s_test$sample_id))
miss <- setdiff(need, colnames(counts))
if (length(miss) > 0) stop(length(miss), " samples missing from counts, e.g. ", miss[1])

if (!is.na(a$debug_n_genes)) {
        set.seed(seed)
        counts <- counts[sample(rownames(counts), min(a$debug_n_genes, nrow(counts))), ]
}

# 1. Gene filter on training samples
genes <- select_genes(counts[, s_train$sample_id], s_train$age)
logmsg(length(genes), " of ", nrow(counts), " genes kept after expression filter")

x_train <- log2(counts[genes, s_train$sample_id] + 1)
x_test <- log2(counts[genes, s_test$sample_id, drop = FALSE] + 1)
rm(counts); invisible(gc())

# 2. Quantile normalisation to the training reference
qn_target <- qn_fit(x_train)
x_train <- qn_apply(x_train, qn_target)
x_test <- qn_apply(x_test, qn_target)
logmsg("quantile normalisation done")

# 3. SVA
################################################################################
s_train$region_f <- droplevels(factor(s_train$region))
sva_mod <- a$sva_mod; sva_mod0 <- a$sva_mod0
if (nlevels(s_train$region_f) < 2) {
        sva_mod <- gsub("\\+\\s*region", "", sva_mod)
        sva_mod0 <- "~ 1"
}
mdat_sva <- data.frame(age = s_train$age, region = s_train$region_f)
mod <- model.matrix(as.formula(sva_mod), data = mdat_sva)
mod0 <- model.matrix(as.formula(sva_mod0), data = mdat_sva)
n_sv <- if (a$n_sv %in% c("leek", "rmt", "auto")) a$n_sv else as.integer(a$n_sv)
sv <- sva_fit(x_train, mod, mod0, n_sv = n_sv, method = a$sva_method)
logmsg("SVA (", a$sva_method, "): ", sv$n_sv, " SVs; |cor| projection vs SVA SVs: ",
       if (sv$n_sv > 0) sprintf("%.3f", min(abs(sv$cor_with_sva))) else "NA")
x_train <- sva_apply(sv, x_train)
x_test <- sva_apply(sv, x_test)

# 4. Elastic net with inner CV grouped by individual
################################################################################
x_tr <- t(x_train); x_te <- t(x_test)
y_tr <- s_train$age
pen <- rep(1, ncol(x_tr))
if (a$region_covariates) {
        lv <- levels(s_train$region_f)
        dummies <- function(reg) {
                d <- sapply(lv[-1], function(l) as.numeric(reg == l))
                matrix(d, nrow = length(reg), dimnames = list(NULL, paste0("region_", lv[-1])))
        }
        x_tr <- cbind(x_tr, dummies(s_train$region))
        x_te <- cbind(x_te, dummies(s_test$region))
        pen <- c(pen, rep(0, length(lv) - 1))
}

inner_units <- unique(data.frame(unit = s_train$group, stratum = s_train$substudy))
inner_units$n_samples <- as.numeric(table(s_train$group)[inner_units$unit])
inner_units <- inner_units[!duplicated(inner_units$unit), ]
inner_fold <- assign_folds(inner_units, a$n_inner, seed = seed)
foldid <- as.integer(inner_fold[s_train$group])

if (a$cores > 1) {
        doParallel::registerDoParallel(cores = a$cores)
}
alphas <- as.numeric(strsplit(a$alphas, ",")[[1]])
cv_fits <- list(); tuning <- data.frame()
for (al in alphas) {
        set.seed(seed)
        cvf <- cv.glmnet(x_tr, y_tr, alpha = al, foldid = foldid, standardize = TRUE,
                         penalty.factor = pen, parallel = a$cores > 1)
        cv_fits[[as.character(al)]] <- cvf
        i_l <- which(cvf$lambda == cvf[[a$lambda]])
        tuning <- rbind(tuning, data.frame(alpha = al, lambda = cvf[[a$lambda]],
                                           cv_mse = cvf$cvm[i_l], cv_se = cvf$cvsd[i_l],
                                           n_nonzero = cvf$nzero[i_l]))
        logmsg(sprintf("alpha %.2f: lambda %.4g, inner CV RMSE %.2f, %d genes",
                       al, cvf[[a$lambda]], sqrt(cvf$cvm[i_l]), cvf$nzero[i_l]))
}
best <- which.min(tuning$cv_mse)
tuning$selected <- seq_len(nrow(tuning)) == best
fit <- cv_fits[[best]]
lam <- fit[[a$lambda]]
logmsg(sprintf("selected alpha %.2f", tuning$alpha[best]))

pred_test <- as.numeric(predict(fit, newx = x_te, s = lam))
pred_train <- as.numeric(predict(fit, newx = x_tr, s = lam))

# Outputs
################################################################################
pred <- data.frame(
        sample_id = s_test$sample_id,
        specimenID = s_test$specimenID,
        individualID = s_test$individualID,
        rep = r, fold = k,
        predicted_age = pred_test,
        in_training = s_test$sample_id %in% s_train$sample_id
)
write.csv(pred, file.path(out, "predictions.csv"), row.names = FALSE)

cf <- as.matrix(coef(fit, s = lam))
cf <- data.frame(feature = rownames(cf), coefficient = cf[, 1])
write.csv(cf[cf$coefficient != 0, ], file.path(out, "coefficients.csv"), row.names = FALSE)
write.csv(tuning, file.path(out, "tuning.csv"), row.names = FALSE)

# Everything needed to apply this fold's clock to new samples.
sv$sv_train <- NULL
saveRDS(list(genes = genes, qn_target = qn_target, sva = sv, glmnet = fit$glmnet.fit,
             lambda = lam, alpha = tuning$alpha[best],
             region_levels = levels(s_train$region_f),
             region_covariates = a$region_covariates,
             train_samples = s_train$sample_id, args = a),
        file.path(out, "clock.rds"))

trm <- metrics(y_tr, pred_train)
tem <- metrics(s_test$age[s_test$is_control], pred_test[s_test$is_control])
logmsg(sprintf("training (in-sample) MAE %.2f r %.3f | held-out controls MAE %.2f r %.3f (n=%d)",
               trm$mae, trm$r, tem$mae, tem$r, tem$n))
logmsg("done in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
