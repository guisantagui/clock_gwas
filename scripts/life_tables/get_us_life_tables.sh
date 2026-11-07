#!/bin/bash
################################################################################
# Download the official US life tables by sex (CDC/NCHS, National Vital        #
# Statistics Reports) and compute E[age at death | death at age >= 89 / 90]    #
# for each year and sex, used to impute the age of top-coded donors            #
# ("89+", "90+") in scripts/clock_cv/03_collect_gaps.R.                        #
#                                                                              #
# E[age at death | alive at x] = x + e(x), where e(x) is the life expectancy   #
# at exact age x in the table, so E[age | >= 90] = 90 + e(90).                 #
#                                                                              #
# Source: https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/    #
#   Table02 = males, Table03 = females (one NVSR report per data year).        #
#   The SSA actuarial tables (ssa.gov/oact/STATS/table4c6.html) are equivalent #
#   but ssa.gov refuses scripted downloads.                                    #
#                                                                              #
# Usage (from the repository root):                                           #
#   bash scripts/life_tables/get_us_life_tables.sh [outdir]                    #
# Output: <outdir>/raw/*.xls(x) and <outdir>/us_life_tables_e89_e90.csv        #
################################################################################
set -euo pipefail

OUTDIR=${1:-data/life_tables}
BASE=https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR
mkdir -p "$OUTDIR/raw"

# Every NVSR report directory on the server.
dirs=$(curl -sS "$BASE/" | grep -o 'NVSR/[0-9][0-9]*[-_][0-9]*/' | sed 's#NVSR/##; s#/##' | sort -u)

for d in $dirs; do
        listing=$(curl -sS "$BASE/$d/" | grep -o 'HREF="[^"]*"' | sed 's#.*/##; s#"##')
        for t in 02 03; do
                # Prefer the regular table; skip intercensal re-estimates.
                f=$(echo "$listing" | grep -i "^table$t\.xlsx\?$" | head -1 || true)
                [[ -z "$f" ]] && continue
                out="$OUTDIR/raw/${d}_$f"
                [[ -s "$out" ]] || curl -sS -o "$out" "$BASE/$d/$f"
                echo "downloaded $d/$f"
        done
done

# Parse: keep only US life tables for males/females (other NVSR reports, e.g.
# state tables, have different titles), extract e(89) and e(90).
Rscript - "$OUTDIR" <<'RS'
outdir <- commandArgs(TRUE)[1]
files <- list.files(file.path(outdir, "raw"), full.names = TRUE)
res <- list()
for (f in files) {
        x <- tryCatch(suppressMessages(readxl::read_excel(f, col_names = FALSE, .name_repair = "minimal")),
                      error = function(e) NULL)
        if (is.null(x)) next
        title <- as.character(x[[1]][1])
        m <- regmatches(title, regexec("Life table for (males|females): United States, ([0-9]{4})$", trimws(title)))[[1]]
        if (length(m) == 0) next
        age_col <- as.character(x[[1]])
        e_col <- suppressWarnings(as.numeric(x[[ncol(x)]]))
        get_e <- function(a) e_col[grep(sprintf("^%d\\D", a), age_col)[1]]
        res[[f]] <- data.frame(year = as.integer(m[3]), sex = sub("s$", "", m[2]),
                               e89 = get_e(89), e90 = get_e(90),
                               nvsr = sub("_.*", "", basename(f)), file = basename(f))
}
d <- do.call(rbind, res)
d$E_age_ge89 <- 89 + d$e89
d$E_age_ge90 <- 90 + d$e90
# If a year appears in several reports, keep the latest report (final revision).
d <- d[order(d$year, d$sex, d$nvsr, decreasing = TRUE), ]
d <- d[!duplicated(d[, c("year", "sex")]), ]
d <- d[order(d$year, d$sex), ]
write.csv(d, file.path(outdir, "us_life_tables_e89_e90.csv"), row.names = FALSE)
print(d[, c("year", "sex", "e89", "e90", "E_age_ge89", "E_age_ge90")], row.names = FALSE)
pre <- d[d$year <= 2019, ]
cat("\nMean over years <= 2019 (excludes COVID-19 years):\n")
print(aggregate(cbind(E_age_ge89, E_age_ge90) ~ sex, data = pre, FUN = function(v) round(mean(v), 2)))
RS
