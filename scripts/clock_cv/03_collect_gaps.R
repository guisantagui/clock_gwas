################################################################################
# Brain clock cross-fitting, step 3: collect out-of-fold predictions and       #
# compute age gaps.                                                            #
#                                                                              #
# * Out-of-fold (OOF) predictions of cv samples are averaged over repeats;     #
#   external samples are averaged over all fold models.                        #
# * Raw gap = predicted - chronological age. Top-coded ages ("90+", "89+") get #
#   the sex-specific imputed age given in --topcode_age_female/male; if not    #
#   given, their gaps are left NA (age_topcoded flags them either way).        #
# * Residualised gap: lm(--resid_formula) fit on the OOF predictions of        #
#   controls with known, non-top-coded ages >= --min_train_age, applied to     #
#   everyone. A version without cohort/region terms is also given, so that     #
#   samples from cohorts/regions without OOF controls (e.g. TBI cases) still   #
#   get a gap.                                                                 #
# * Individual level: mean residualised gap (+ random-effect BLUP if lme4 is   #
#   installed); regional deviation = sample gap - individual's mean gap.       #
################################################################################

suppressPackageStartupMessages({
        library(argparser)
        library(splines)
})
source(file.path((function() {
        a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
        if (length(f)) dirname(normalizePath(f[1])) else getwd()
})(), "cv_utils.R"))

p <- arg_parser("Collect OOF predictions and compute age gaps.")
p <- add_argument(p, "--outdir", help = "Output root used by 02_run_fold.R.", default = "results/clock_cv")
p <- add_argument(p, "--folds", help = "folds.csv.", default = "results/clock_cv/folds.csv")
p <- add_argument(p, "--samp_mdat", help = "Sample metadata.",
                  default = "results/counts_parsed/sample_metadata_merged.csv")
p <- add_argument(p, "--ind_mdat", help = "Individual metadata.",
                  default = "results/counts_parsed/individual_metadata_merged_rnaseq_filt.csv")
p <- add_argument(p, "--region_map", help = "tissue -> region table.",
                  default = file.path(script_dir(), "region_map.csv"))
p <- add_argument(p, "--exclude_regions", help = "Regions left out.", default = "cerebellum")
p <- add_argument(p, "--min_train_age", help = "Minimum age of controls in the residualisation fit.", default = 20)
p <- add_argument(p, "--resid_formula", help = "Residualisation model (fit on OOF controls).",
                  default = "predicted_age ~ ns(age_used, df = 4) + substudy_batch + region")
p <- add_argument(p, "--topcode_age_female", help = "Imputed age for top-coded females (e.g. life-table E[age | >= 90]).",
                  type = "numeric", default = NA)
p <- add_argument(p, "--topcode_age_male", help = "Imputed age for top-coded males.",
                  type = "numeric", default = NA)
a <- parse_args(p)

# Predictions
################################################################################
pred_files <- list.files(a$outdir, pattern = "^predictions\\.csv$", recursive = TRUE, full.names = TRUE)
pred_files <- pred_files[grepl("rep[0-9]+/fold[1-9][0-9]*/predictions\\.csv$", pred_files)]
if (length(pred_files) == 0) stop("No fold predictions under ", a$outdir)
pr <- do.call(rbind, lapply(pred_files, read.csv, stringsAsFactors = FALSE))
pr <- pr[!pr$in_training, ]
message(length(pred_files), " fold prediction files, repeats: ",
        paste(sort(unique(pr$rep)), collapse = ","))

s <- load_metadata(a$samp_mdat, a$ind_mdat, a$region_map)
s <- s[!s$region %in% strsplit(a$exclude_regions, ",")[[1]], ]
folds <- read.csv(a$folds, stringsAsFactors = FALSE)
s$role <- folds$role[match(s$individualID, folds$individualID)]

# Each cv sample must have exactly one OOF prediction per repeat.
per_rep <- table(pr$sample_id[pr$sample_id %in% s$sample_id[s$role == "cv"]],
                 pr$rep[pr$sample_id %in% s$sample_id[s$role == "cv"]])
if (any(per_rep > 1)) stop("Some cv samples have more than one OOF prediction in a repeat.")
n_missing <- sum(!s$sample_id[s$role == "cv"] %in% pr$sample_id)
if (n_missing > 0) warning(n_missing, " cv samples have no OOF prediction (folds not run yet?).")

agg <- aggregate(predicted_age ~ sample_id, data = pr, FUN = mean)
agg$predicted_age_sd <- aggregate(predicted_age ~ sample_id, data = pr, FUN = sd)$predicted_age
agg$n_predictions <- aggregate(predicted_age ~ sample_id, data = pr, FUN = length)$predicted_age

d <- merge(s, agg, by = "sample_id")

# Ages used for gaps
################################################################################
d$age_used <- d$age
tc <- d$age_topcoded
d$age_used[tc & d$sex == "female"] <- a$topcode_age_female
d$age_used[tc & d$sex == "male"] <- a$topcode_age_male
d$gap_raw <- d$predicted_age - d$age_used
d$age_below_train_range <- !is.na(d$age) & d$age < a$min_train_age

# Residualisation on OOF controls
################################################################################
d$in_resid_fit <- d$is_control & !is.na(d$age) & !d$age_topcoded & d$age >= a$min_train_age
fit_d <- d[d$in_resid_fit, ]
f_full <- as.formula(a$resid_formula)
f_age <- update(f_full, . ~ ns(age_used, df = 4))

resid_with <- function(f) {
        fit <- lm(f, data = fit_d)
        vars <- intersect(all.vars(f)[-1], c("substudy", "substudy_batch", "region", "sex"))
        ok <- !is.na(d$age_used)
        for (v in vars) ok <- ok & d[[v]] %in% unique(fit_d[[v]])
        g <- rep(NA_real_, nrow(d))
        g[ok] <- d$predicted_age[ok] - predict(fit, newdata = d[ok, ])
        list(gap = g, fit = fit, unfit_levels = unique(d[!ok & !is.na(d$age_used), vars, drop = FALSE]))
}
r_full <- resid_with(f_full)
r_age <- resid_with(f_age)
d$gap_resid <- r_full$gap
d$gap_resid_ageonly <- r_age$gap
if (nrow(r_full$unfit_levels) > 0) {
        message("No OOF controls for these levels; gap_resid is NA there (use gap_resid_ageonly):")
        print(r_full$unfit_levels)
}

# Individual-level phenotypes and regional deviations
################################################################################
ind_mean <- tapply(d$gap_resid, d$individualID, mean, na.rm = TRUE)
ind_mean[is.nan(ind_mean)] <- NA
d$gap_ind_mean <- as.numeric(ind_mean[d$individualID])
n_reg <- tapply(d$region[!is.na(d$gap_resid)], d$individualID[!is.na(d$gap_resid)],
                function(x) length(unique(x)))
d$gap_regional_dev <- ifelse(d$individualID %in% names(n_reg)[n_reg >= 2],
                             d$gap_resid - d$gap_ind_mean, NA)

keep_cols <- c("sample_id", "specimenID", "individualID", "substudy", "substudy_batch",
               "tissue", "region", "sex", "race", "diagn_4BrainClck", "diagnosis",
               "has_genomics", "role", "RIN", "pmi", "age_label", "age", "age_topcoded",
               "age_used", "age_below_train_range", "predicted_age", "predicted_age_sd",
               "n_predictions", "gap_raw", "gap_resid", "gap_resid_ageonly",
               "gap_ind_mean", "gap_regional_dev", "in_resid_fit")
d <- d[order(d$substudy, d$individualID, d$region), keep_cols]
write.csv(d, file.path(a$outdir, "sample_gaps.csv"), row.names = FALSE)

ind <- as.data.frame(data.table::as.data.table(d)[, .(
        substudy = substudy[1], diagn = diagn_4BrainClck[1], has_genomics = has_genomics[1],
        sex = sex[1], age_label = age_label[1], age_topcoded = age_topcoded[1],
        n_samples = sum(!is.na(gap_resid)), n_regions = data.table::uniqueN(region[!is.na(gap_resid)]),
        regions = paste(sort(unique(region)), collapse = ";"),
        gap_mean = mean(gap_resid, na.rm = TRUE),
        gap_ageonly_mean = mean(gap_resid_ageonly, na.rm = TRUE)
), by = individualID])
ind$gap_mean[is.nan(ind$gap_mean)] <- NA
ind$gap_ageonly_mean[is.nan(ind$gap_ageonly_mean)] <- NA

if (requireNamespace("lme4", quietly = TRUE)) {
        dd <- d[!is.na(d$gap_resid), ]
        mm <- lme4::lmer(gap_resid ~ region + (1 | individualID), data = dd)
        re <- lme4::ranef(mm)$individualID
        ind$gap_blup <- re[match(ind$individualID, rownames(re)), 1]
        message(sprintf("Mixed model: individual SD %.2f, residual SD %.2f",
                        attr(lme4::VarCorr(mm)$individualID, "stddev"), sigma(mm)))
}
write.csv(ind, file.path(a$outdir, "individual_gaps.csv"), row.names = FALSE)

# Performance and sanity summaries
################################################################################
ctrl <- d[d$in_resid_fit, ]
by_grp <- function(v) {
        do.call(rbind, lapply(split(ctrl, ctrl[[v]]), function(x)
                cbind(group = x[[v]][1], metrics(x$age, x$predicted_age))))
}
perf <- rbind(cbind(level = "all", group = "all", metrics(ctrl$age, ctrl$predicted_age)),
              cbind(level = "substudy", by_grp("substudy")),
              cbind(level = "region", by_grp("region")))
write.csv(perf, file.path(a$outdir, "oof_performance_controls.csv"), row.names = FALSE)

cc <- aggregate(gap_resid ~ substudy + diagn_4BrainClck, data = d,
                FUN = function(x) c(n = length(x), mean = mean(x), sd = sd(x)))
cc <- do.call(data.frame, cc)
write.csv(cc, file.path(a$outdir, "gap_by_substudy_diagnosis.csv"), row.names = FALSE)

cat("\nOOF performance on controls (residualisation set):\n")
print(perf, digits = 3, row.names = FALSE)
cat("\nResidualised gap by substudy x diagnosis:\n")
print(cc, digits = 3, row.names = FALSE)
cat("\nGenotyped individuals with a gap:", sum(ind$has_genomics & !is.na(ind$gap_mean)), "\n")
cat("Top-coded genotyped individuals without a gap (no imputed age given):",
    sum(ind$has_genomics & ind$age_topcoded & is.na(ind$gap_mean)), "\n")

# Fold-to-fold consistency of the clocks
################################################################################
# Each held-out sample is predicted by a different fold model. These checks
# quantify how much the models differ:
#   * per-fold performance on held-out controls (before averaging repeats);
#   * agreement of the fold models on samples that ALL of them predict out of
#     sample (external samples, and cases vs the final model, which never
#     trains on cases);
#   * correlation of coefficient vectors between fold models.
ctrl_ids <- d$sample_id[d$in_resid_fit]
pf <- merge(pr, d[, c("sample_id", "age")], by = "sample_id")
pf <- pf[pf$sample_id %in% ctrl_ids & !pf$sample_id %in% d$sample_id[d$role == "external"], ]
per_fold <- do.call(rbind, lapply(split(pf, list(pf$rep, pf$fold), drop = TRUE), function(x)
        cbind(rep = x$rep[1], fold = x$fold[1], metrics(x$age, x$predicted_age))))
write.csv(per_fold, file.path(a$outdir, "oof_performance_per_fold.csv"), row.names = FALSE)
cat("\nPer-fold performance on held-out controls:\n")
print(per_fold, digits = 3, row.names = FALSE)

ext <- pr[pr$sample_id %in% d$sample_id[d$role == "external"], ]
if (nrow(ext) > 0) {
        ext$model <- paste0("r", ext$rep, "f", ext$fold)
        w <- reshape(ext[, c("sample_id", "model", "predicted_age")], idvar = "sample_id",
                     timevar = "model", direction = "wide")
        wm <- as.matrix(w[, -1])
        cc_m <- cor(wm, use = "pairwise.complete.obs")
        cat(sprintf("\nExternal samples (n=%d) predicted by %d fold models: median SD across models %.2f y; pairwise r min %.3f, median %.3f\n",
                    nrow(wm), ncol(wm), median(apply(wm, 1, sd)),
                    min(cc_m[upper.tri(cc_m)]), median(cc_m[upper.tri(cc_m)])))
}

f0 <- file.path(a$outdir, "rep1", "fold0", "predictions.csv")
if (file.exists(f0)) {
        p0 <- read.csv(f0, stringsAsFactors = FALSE)
        cases <- d[d$diagn_4BrainClck != "Control" & d$role == "cv", c("sample_id", "predicted_age")]
        cmp <- merge(cases, p0[, c("sample_id", "predicted_age")], by = "sample_id",
                     suffixes = c("_oof", "_final"))
        dif <- cmp$predicted_age_oof - cmp$predicted_age_final
        cat(sprintf("Cases (n=%d), OOF (mean of repeats) vs final all-controls clock: r %.3f, mean diff %.2f, SD diff %.2f y\n",
                    nrow(cmp), cor(cmp$predicted_age_oof, cmp$predicted_age_final),
                    mean(dif), sd(dif)))
        write.csv(cmp, file.path(a$outdir, "cases_oof_vs_final.csv"), row.names = FALSE)
}

cf_files <- list.files(a$outdir, pattern = "^coefficients\\.csv$", recursive = TRUE, full.names = TRUE)
cf_files <- cf_files[grepl("rep[0-9]+/fold[0-9]+/", cf_files)]
if (length(cf_files) > 1) {
        cfs <- lapply(cf_files, read.csv, stringsAsFactors = FALSE)
        names(cfs) <- sub(".*(rep[0-9]+)/(fold[0-9]+)/.*", "\\1_\\2", cf_files)
        feats <- setdiff(unique(unlist(lapply(cfs, `[[`, "feature"))), "(Intercept)")
        cm <- sapply(cfs, function(x) { v <- setNames(rep(0, length(feats)), feats)
                     v[intersect(x$feature, feats)] <- x$coefficient[match(intersect(x$feature, feats), x$feature)]; v })
        ccf <- cor(cm)
        write.csv(ccf, file.path(a$outdir, "coefficient_correlation_between_folds.csv"))
        sel <- colSums(cm != 0)
        cat(sprintf("Coefficient vectors across %d models: pairwise r min %.3f, median %.3f; %d genes selected by every model, %d by any\n",
                    ncol(cm), min(ccf[upper.tri(ccf)]), median(ccf[upper.tri(ccf)]),
                    sum(rowSums(cm != 0) == ncol(cm)), sum(rowSums(cm != 0) > 0)))
}
