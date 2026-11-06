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
                  default = "predicted_age ~ ns(age_used, df = 4) + substudy + region")
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
