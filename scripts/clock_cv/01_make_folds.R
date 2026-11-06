################################################################################
# Brain clock cross-fitting, step 1: assign individuals to CV folds.           #
#                                                                              #
# Roles:                                                                       #
#   always_train : controls without genomics -> in the training set of every   #
#                  fold (they never need an out-of-fold gap for the GWAS).     #
#   external     : substudies kept out of training entirely (default           #
#                  brainSeq_pI, as in the published clock); predicted by every #
#                  fold model and by the final model.                          #
#   cv           : everyone else (genotyped individuals, cases and controls,   #
#                  plus cases without genomics). Assigned to K folds, grouped  #
#                  by individual (all regions of a donor together; related or  #
#                  duplicated donors together if --related is given) and       #
#                  stratified by substudy x diagnosis x age bin.               #
#                                                                              #
# Output: <outdir>/folds.csv with one row per individual and one fold column   #
# per repeat (fold_r1, fold_r2, ...). fold = 0 for always_train / external.    #
################################################################################

suppressPackageStartupMessages(library(argparser))
source(file.path((function() {
        a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
        if (length(f)) dirname(normalizePath(f[1])) else getwd()
})(), "cv_utils.R"))

p <- arg_parser("Assign individuals to cross-fitting folds.")
p <- add_argument(p, "--samp_mdat", help = "Sample metadata (aligned with counts).",
                  default = "results/counts_parsed/sample_metadata_merged.csv")
p <- add_argument(p, "--ind_mdat", help = "Individual metadata with has_WGS/has_snpArray.",
                  default = "results/counts_parsed/individual_metadata_merged_rnaseq_filt.csv")
p <- add_argument(p, "--region_map", help = "tissue -> region table.",
                  default = file.path(script_dir(), "region_map.csv"))
p <- add_argument(p, "--exclude_regions", help = "Comma-separated regions left out of the analysis.",
                  default = "cerebellum")
p <- add_argument(p, "--external", help = "Comma-separated substudies held out from training.",
                  default = "brainSeq_pI")
p <- add_argument(p, "--related", help = "Optional two-column file (ID1, ID2) of duplicated/related individuals (e.g. from KING) to force into the same fold.",
                  default = NA)
p <- add_argument(p, "--K", help = "Number of folds.", default = 5L)
p <- add_argument(p, "--n_rep", help = "Number of repeats with different seeds.", default = 1L)
p <- add_argument(p, "--seed", help = "Base seed.", default = 2026L)
p <- add_argument(p, "--outdir", help = "Output directory.", default = "results/clock_cv")
a <- parse_args(p)

exclude_regions <- strsplit(a$exclude_regions, ",")[[1]]
external <- strsplit(a$external, ",")[[1]]
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

s <- load_metadata(a$samp_mdat, a$ind_mdat, a$region_map)
s <- s[!s$region %in% exclude_regions, ]
message(sprintf("%d samples / %d individuals after excluding regions: %s",
                nrow(s), length(unique(s$individualID)), paste(exclude_regions, collapse = ", ")))

# Individual table
ind <- as.data.table(s)[, .(
        substudy = substudy[1],
        diagn = diagn_ind[1],
        has_genomics = has_genomics[1],
        age = age[1],
        age_topcoded = age_topcoded[1],
        n_samples = .N,
        n_regions = uniqueN(region),
        regions = paste(sort(unique(region)), collapse = ";")
), by = individualID]
ind <- as.data.frame(ind)

ind$role <- ifelse(ind$substudy %in% external, "external",
                   ifelse(!ind$has_genomics & ind$diagn == "Control", "always_train", "cv"))
ind$age_bin <- age_bin(ind$age, ind$age_topcoded)
ind$stratum <- paste(ind$substudy, ind$diagn, ind$age_bin, sep = "|")

# Related / duplicated individuals
pairs <- NULL
if (!is.na(a$related)) {
        pairs <- read.table(a$related, header = TRUE, stringsAsFactors = FALSE)[, 1:2]
        message(nrow(pairs), " related pairs read from ", a$related)
}
cv <- ind$role == "cv"
ind$group <- NA_character_
ind$group[cv] <- make_groups(ind$individualID[cv], pairs)
if (!is.null(pairs)) {
        cross_role <- pairs[(pairs[[1]] %in% ind$individualID[cv]) !=
                                    (pairs[[2]] %in% ind$individualID[cv]) &
                                    pairs[[1]] %in% ind$individualID &
                                    pairs[[2]] %in% ind$individualID, ]
        if (nrow(cross_role) > 0) {
                warning(nrow(cross_role), " related pairs link a CV individual with an ",
                        "always_train/external one; that relative will be in training ",
                        "when the other is held out. Review these:")
                print(cross_role)
        }
}

# A group takes the stratum of its member with most samples.
units <- as.data.table(ind[cv, ])[order(-n_samples), .(
        stratum = stratum[1], n_samples = sum(n_samples)), by = group]
setnames(units, "group", "unit")
units <- as.data.frame(units)

for (r in seq_len(a$n_rep)) {
        f <- assign_folds(units, a$K, seed = a$seed + r)
        col <- paste0("fold_r", r)
        ind[[col]] <- 0L
        ind[[col]][cv] <- f[ind$group[cv]]
}

write.csv(ind, file.path(a$outdir, "folds.csv"), row.names = FALSE)

# Summaries
################################################################################
sink(file.path(a$outdir, "folds_summary.txt"), split = TRUE)
cat("Individuals by role and substudy:\n")
print(table(ind$substudy, ind$role))
cat("\nSamples by role and substudy:\n")
print(xtabs(n_samples ~ substudy + role, data = ind))
cat("\nCV individuals by fold (repeat 1) and substudy x diagnosis:\n")
print(table(paste(ind$substudy, ind$diagn)[cv], ind$fold_r1[cv]))
cat("\nCV samples by fold (repeat 1):\n")
print(tapply(ind$n_samples[cv], ind$fold_r1[cv], sum))
cat("\nCV individuals by fold (repeat 1) and age bin:\n")
print(table(ind$age_bin[cv], ind$fold_r1[cv]))
cat("\nTop-coded ages among genotyped individuals (substudy x diagnosis):\n")
g <- ind$has_genomics & ind$age_topcoded
print(table(ind$substudy[g], ind$diagn[g]))
cat("\nGenotyped individuals per region (for per-region GWAS power):\n")
reg <- s[s$has_genomics, c("individualID", "region", "is_control")]
reg <- unique(reg)
print(table(reg$region, ifelse(reg$is_control, "Control", "Rest")))
sink()
message("Saved ", file.path(a$outdir, "folds.csv"))
