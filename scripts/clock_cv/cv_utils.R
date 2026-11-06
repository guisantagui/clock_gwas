################################################################################
# Brain clock cross-fitting: shared helper functions.                          #
#                                                                              #
# Sourced by 01_make_folds.R, 02_run_fold.R and 03_collect_gaps.R.             #
################################################################################

suppressPackageStartupMessages({
        library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# Directory of the script being run with Rscript (used to locate region_map.csv
# and this file regardless of the working directory).
script_dir <- function() {
        args <- commandArgs(trailingOnly = FALSE)
        f <- sub("^--file=", "", args[grepl("^--file=", args)])
        if (length(f) == 0) return(getwd())
        dirname(normalizePath(f[1]))
}

# Ages
################################################################################

# Parses the character ageDeath column. Numeric strings are kept as numbers.
# Top-coded labels ("90+", "89+") are flagged and their numeric value is set to
# NA; the lower bound is kept in age_floor. Empty strings / NA are missing.
parse_age <- function(x) {
        x <- trimws(as.character(x))
        x[x == ""] <- NA
        topcoded <- !is.na(x) & grepl("^[0-9.]+\\+$", x)
        num <- suppressWarnings(as.numeric(x))
        num[topcoded] <- NA
        floor_ <- rep(NA_real_, length(x))
        floor_[topcoded] <- as.numeric(sub("\\+$", "", x[topcoded]))
        unknown <- !is.na(x) & !topcoded & is.na(num)
        if (any(unknown)) {
                stop("Unrecognised ageDeath values: ",
                     paste(unique(x[unknown]), collapse = ", "))
        }
        data.frame(age = num, age_topcoded = topcoded, age_floor = floor_,
                   age_label = x, stringsAsFactors = FALSE)
}

# Metadata
################################################################################

# Loads sample and individual metadata, adds parsed ages, harmonised region and
# genomics availability. sample_id is make.names(specimenID), which is how the
# samples are named in the columns of counts_merged.csv.
load_metadata <- function(samp_f, ind_f, region_map_f) {
        s <- read.csv(samp_f, row.names = 1, stringsAsFactors = FALSE)
        i <- read.csv(ind_f, row.names = 1, stringsAsFactors = FALSE)
        if (any(duplicated(i$individualID))) stop("Duplicated individualID in ", ind_f)
        if (!all(s$individualID %in% i$individualID)) {
                stop("Samples whose individualID is not in ", ind_f)
        }
        s$sample_id <- make.names(s$specimenID)
        if (any(duplicated(s$sample_id))) stop("Duplicated sample IDs after make.names()")
        s <- cbind(s, parse_age(s$ageDeath))

        rmap <- read.csv(region_map_f, stringsAsFactors = FALSE)
        s$region <- rmap$region[match(s$tissue, rmap$tissue)]
        if (any(is.na(s$region))) {
                stop("Tissues missing from ", region_map_f, ": ",
                     paste(unique(s$tissue[is.na(s$region)]), collapse = ", "))
        }

        m <- match(s$individualID, i$individualID)
        s$has_genomics <- (i$has_WGS | i$has_snpArray)[m] %in% TRUE
        s$diagn_ind <- i$diagn_4BrainClck[m]
        s$is_control <- s$diagn_4BrainClck == "Control"
        s
}

# Flags samples whose labels can be used for training (independently of folds).
flag_train_eligible <- function(s, min_train_age, min_rin, external_substudies) {
        s$is_control &
                !is.na(s$age) &
                s$age >= min_train_age &
                (is.na(s$RIN) | s$RIN >= min_rin) &
                !s$substudy %in% external_substudies
}

# Counts
################################################################################

# Returns the gene universe: "all", "lincs" (LINCS L1000 genes) or the path to a
# text file with one Ensembl gene ID per line.
read_gene_universe <- function(universe, lincs_f) {
        if (universe == "all") return(NULL)
        f <- if (universe == "lincs") lincs_f else universe
        g <- trimws(readLines(f))
        g[g != ""]
}

# Reads counts_merged.csv (genes in rows, samples in columns) restricted to the
# gene universe, as a numeric matrix. The result is cached as RDS because the
# CSV takes a while to parse.
load_counts <- function(counts_f, genes = NULL, cache_f = NULL) {
        if (!is.null(cache_f) && file.exists(cache_f) &&
            file.mtime(cache_f) < file.mtime(counts_f)) {
                message("Counts file is newer than the cache; rebuilding ", cache_f)
                file.remove(cache_f)
        }
        if (!is.null(cache_f) && file.exists(cache_f)) {
                message("Loading cached counts: ", cache_f)
                return(readRDS(cache_f))
        }
        message("Reading ", counts_f, " (this takes a few minutes)...")
        dt <- fread(counts_f, header = TRUE)
        gene_ids <- dt[[1]]
        dt[, 1 := NULL]
        keep <- if (is.null(genes)) rep(TRUE, length(gene_ids)) else gene_ids %in% genes
        m <- as.matrix(dt[keep])
        rownames(m) <- gene_ids[keep]
        storage.mode(m) <- "double"
        rm(dt); gc()
        if (!is.null(cache_f)) {
                dir.create(dirname(cache_f), recursive = TRUE, showWarnings = FALSE)
                saveRDS(m, cache_f)
        }
        m
}

# Gene filtering (fit on training samples only)
################################################################################

# Same criteria as filt_dataset.R in the published pipeline, but computed on the
# training samples only: drop genes with >= prop_zero_max zeros, and genes whose
# mean count is < 1 in every 5-year age bin.
select_genes <- function(counts_train, ages_train, prop_zero_max = 0.8, bin_width = 5) {
        prop_zero <- rowMeans(counts_train == 0)
        keep <- prop_zero < prop_zero_max
        bins <- floor(ages_train / bin_width)
        low_all <- rep(TRUE, nrow(counts_train))
        for (b in unique(bins)) {
                idx <- which(bins == b)
                if (length(idx) < 2) next
                low_all <- low_all & (rowMeans(counts_train[, idx, drop = FALSE]) < 1)
        }
        rownames(counts_train)[keep & !low_all]
}

# Quantile normalisation (reference from training samples)
################################################################################

qn_fit <- function(log2_train) {
        preprocessCore::normalize.quantiles.determine.target(log2_train)
}

# Each sample is mapped onto the reference independently of the others.
qn_apply <- function(log2_x, target) {
        out <- preprocessCore::normalize.quantiles.use.target(log2_x, target)
        dimnames(out) <- dimnames(log2_x)
        out
}

# Surrogate variable analysis with per-sample projection
################################################################################

# Estimates SVs on the training samples with sva (or SmartSVA), and builds a
# projection that computes the SVs of any sample from its expression alone,
# without its age. It follows the "fast" branch of sva::fsva, with two
# changes so that a sample's correction depends only on itself and the training
# set, not on the other held-out samples:
#   * new samples are centred with the training gene means;
#   * the new SVs are NOT re-centred/re-scaled within the held-out batch.
# The same projection is applied to the training samples (for sva it
# reproduces sva's SVs exactly), so training and held-out data are corrected
# with the same map.
sva_fit <- function(dat, mod, mod0, n_sv = "auto", method = c("sva", "smartsva"), B = 5) {
        method <- match.arg(method)
        if (identical(n_sv, "auto")) {
                n_sv <- sva::num.sv(dat, mod, method = "leek")
                message("num.sv (leek): ", n_sv)
        }
        n_sv <- as.integer(n_sv)
        if (n_sv == 0) return(list(n_sv = 0L, genes = rownames(dat)))
        if (method == "sva") {
                fit <- sva::sva(dat, mod, mod0, n.sv = n_sv, B = B)
        } else {
                fit <- SmartSVA::smartsva.cpp(dat, mod, mod0, n.sv = n_sv, B = B)
        }
        wts <- fit$pprob.gam * (1 - fit$pprob.b)
        mu <- rowMeans(dat)
        xw <- wts * (dat - mu)
        sv_dec <- irlba::irlba(xw, nv = n_sv, nu = n_sv, tol = 1e-8, maxit = 2000)
        U <- sv_dec$u
        D <- sv_dec$d
        V <- sv_dec$v
        # Align signs with the SVs returned by sva (cosmetic).
        sgn <- sign(colSums(V * fit$sv[, seq_len(n_sv), drop = FALSE]))
        sgn[sgn == 0] <- 1
        U <- sweep(U, 2, sgn, "*")
        V <- sweep(V, 2, sgn, "*")
        mod_full <- cbind(mod, V)
        coefs <- dat %*% mod_full %*% solve(crossprod(mod_full))
        gammahat <- coefs[, (ncol(mod) + 1):ncol(mod_full), drop = FALSE]
        list(n_sv = n_sv, method = method, genes = rownames(dat), wts = wts,
             mu = mu, U = U, D = D, gammahat = gammahat,
             sv_train = V, cor_with_sva = diag(cor(V, fit$sv[, seq_len(n_sv), drop = FALSE])))
}

sva_project <- function(fit, x) {
        x <- x[fit$genes, , drop = FALSE]
        t(crossprod(fit$U, fit$wts * (x - fit$mu)) / fit$D)
}

sva_apply <- function(fit, x) {
        x <- x[fit$genes, , drop = FALSE]
        if (fit$n_sv == 0) return(x)
        x - fit$gammahat %*% t(sva_project(fit, x))
}

# Fold assignment
################################################################################

# Groups individuals that must share a fold (duplicates/relatives). pairs is a
# two-column data.frame of individualIDs. Returns a group label per individual.
make_groups <- function(ids, pairs = NULL) {
        grp <- setNames(seq_along(ids), ids)
        if (is.null(pairs) || nrow(pairs) == 0) return(paste0("g", grp))
        pairs <- pairs[pairs[[1]] %in% ids & pairs[[2]] %in% ids, , drop = FALSE]
        repeat {
                changed <- FALSE
                for (r in seq_len(nrow(pairs))) {
                        a <- pairs[[1]][r]; b <- pairs[[2]][r]
                        m <- min(grp[a], grp[b])
                        if (grp[a] != m || grp[b] != m) {
                                grp[grp == grp[a] | grp == grp[b]] <- m
                                changed <- TRUE
                        }
                }
                if (!changed) break
        }
        paste0("g", grp)
}

# Greedy stratified assignment of units (groups of individuals) to K folds.
# Units are processed stratum by stratum, largest first; each goes to the fold
# with fewest samples of its stratum, ties broken by total samples, then
# randomly. units: data.frame(unit, stratum, n_samples).
assign_folds <- function(units, K, seed) {
        set.seed(seed)
        units <- units[sample(nrow(units)), , drop = FALSE]
        units <- units[order(units$stratum, -units$n_samples), , drop = FALSE]
        strata <- unique(units$stratum)
        by_stratum <- matrix(0, nrow = length(strata), ncol = K, dimnames = list(strata, NULL))
        total <- numeric(K)
        fold <- integer(nrow(units))
        for (u in seq_len(nrow(units))) {
                st <- units$stratum[u]
                cand <- which(by_stratum[st, ] == min(by_stratum[st, ]))
                cand <- cand[total[cand] == min(total[cand])]
                f <- if (length(cand) == 1) cand else sample(cand, 1)
                fold[u] <- f
                by_stratum[st, f] <- by_stratum[st, f] + units$n_samples[u]
                total[f] <- total[f] + units$n_samples[u]
        }
        setNames(fold, units$unit)
}

age_bin <- function(age, topcoded) {
        b <- as.character(cut(age, c(-Inf, 40, 60, 70, 80, 90, Inf), right = FALSE,
                              labels = c("<40", "40-60", "60-70", "70-80", "80-90", "90+")))
        b[topcoded] <- "90+"
        b[is.na(b)] <- "NA"
        b
}

# Metrics
################################################################################

metrics <- function(obs, pred) {
        ok <- !is.na(obs) & !is.na(pred)
        obs <- obs[ok]; pred <- pred[ok]
        data.frame(n = length(obs),
                   mae = mean(abs(obs - pred)),
                   rmse = sqrt(mean((obs - pred)^2)),
                   r = if (length(obs) > 2) cor(obs, pred) else NA,
                   r2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
                   bias = mean(pred - obs))
}
