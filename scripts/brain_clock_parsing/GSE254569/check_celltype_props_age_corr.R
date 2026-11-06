library(plotUtils)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)

mdat_f <- "/work/projects/age_sorter/data/counts/GSE254569/GSE254569/extracted_data/metadata.csv"
mdat <- read_table_fast(mdat_f)
mdat <- mdat[mdat$Classification == "Control", ]

out_dir <- "/work/projects/age_sorter/results/explor_counts/GSE254569/"
create_dir_if_not(out_dir)

cell_type_counts <- mdat %>%
  group_by(Donor, celltypes_final) %>%
  summarise(count = n(), .groups = 'drop')

  individual_totals <- mdat %>%
    group_by(Donor) %>%
    summarise(total_cells = n(), .groups = 'drop')

merged_data <- left_join(cell_type_counts, individual_totals, by = "Donor")

merged_data <- merged_data %>%
  mutate(proportion = count / total_cells)

pivoted_data <- merged_data %>%
  select(Donor, celltypes_final, proportion) %>%
  spread(key = celltypes_final, value = proportion, fill = 0)

pivoted_data <- as.data.frame(pivoted_data)
rownames(pivoted_data) <- pivoted_data$Donor
pivoted_data <- pivoted_data[, 2:ncol(pivoted_data)]
pivoted_data$age <- mdat$Age[match(rownames(pivoted_data), mdat$Donor)]
pivoted_data$RIN <- mdat$RIN[match(rownames(pivoted_data), mdat$Donor)]
pivoted_data$Brain.pH <- mdat$Brain.pH[match(rownames(pivoted_data), mdat$Donor)]
pivoted_data$PMI <- mdat$PMI[match(rownames(pivoted_data), mdat$Donor)]
pivoted_data$lib_batch <- as.factor(mdat$lib_batch[match(rownames(pivoted_data), mdat$Donor)])
pivoted_data$Sex <- mdat$Sex[match(rownames(pivoted_data), mdat$Donor)]
colnames(pivoted_data) <- make.names(colnames(pivoted_data))
c_types <- colnames(pivoted_data)
c_types <- c_types[!c_types %in% c("age", "RIN", "Brain.pH", "PMI", "lib_batch", "Sex")]
pvals <- c()
coefs <- c()
for (cell in c_types){
    l_mod <- lm(as.formula(sprintf("age ~ %s + RIN + Brain.pH + PMI + lib_batch + Sex", cell)),
                pivoted_data)
    #l_mod <- lm(as.formula(sprintf("age ~ %s", cell)),
    #            pivoted_data)
    pval <- summary(l_mod)$coefficients[cell, "Pr(>|t|)"]
    coef <- summary(l_mod)$coefficients[cell, "Estimate"]
    pvals <- c(pvals, pval)
    coefs <- c(coefs, coef)
}
pvals_df <- data.frame(cell_type = c_types,
                       p_val = pvals,
                       coef = coefs)
pvals_df$p_adj <- p.adjust(pvals_df$p_val, method = "BH")
write.csv(pvals_df, file = sprintf("%scellprop_vs_age_sign.csv", out_dir))

sign_cells <- pvals_df[pvals_df$p_adj <= 0.05, ]
for (cell in sign_cells$cell_type){
    plt <- ggplot(pivoted_data,
                   mapping = aes(y = !!sym(cell),
                                 x = age)) +
            geom_point()
    ggsave(plot = plt, filename = sprintf("%s%s.pdf", out_dir, cell))
}