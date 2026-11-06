################################################################################
# Brain clock: divide samples in train and test, ensuring that the donors are  #
# stratified, i.e. all the samples coming from a single donor are either in    #
# train set or in test set.                                                    #
################################################################################

if (!require("dplyr",quietly = T)){
    install.packages("dplyr",
                     repos = 'http://cran.us.r-project.org')
}
if (!require("devtools",quietly = T)){
    install.packages("devtools",
                     repos = 'http://cran.us.r-project.org')
}
if (!require(plotUtils, quietly = T)){
        devtools::install_github("guisantagui/plotUtils", upgrade = "never")
}
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)
library(dplyr)
library(plotUtils)

# Terminal argument parser
################################################################################
parser <- arg_parser("Divide the samples in train and test, stratifying by individualID.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--metDat",
                               "--trainProp",
                               "--exclude_substudy",
                               "--outdir"),
                       help = c("Expression dataframe that is to be split.",
                                "Metadata file",
                                "Proportion of controls to be kept in training set.",
                                "Substudy to be excluded from the train-test split, with the purpose of doing an additional validation with an entire substudy.",
                                "Output directory."),
                       flag = c(F,
                                F,
                                F,
                                F,
                                F))

parsed <- parse_args(parser)

ind_mdat_filt[!(!ind_mdat_filt$has_WGS | !ind_mdat_filt$has_snpArray) & ind_mdat_filt$diagn_4BrainClck == "Control", ]
ind_mdat_filt[(!ind_mdat_filt$has_WGS & !ind_mdat_filt$has_snpArray) & ind_mdat_filt$diagn_4BrainClck == "Control", ]


controls_wo_genomics <- ind_mdat_filt$individualID[
        (!ind_mdat_filt$has_WGS & !ind_mdat_filt$has_snpArray) & ind_mdat_filt$diagn_4BrainClck == "Control"
]

dim(ind_mdat_filt[(!ind_mdat_filt$has_WGS & !ind_mdat_filt$has_snpArray) & ind_mdat_filt$diagn_4BrainClck == "Control", ])

table(ind_mdat_filt[(!ind_mdat_filt$has_WGS & !ind_mdat_filt$has_snpArray) & ind_mdat_filt$diagn_4BrainClck == "Control", ]$substudy)

sum(make.names(samp_mdat_rnaseq$individualID) %in% controls_wo_genomics)


controls_w_genomics <- ind_mdat_filt$individualID[
        (ind_mdat_filt$has_WGS | ind_mdat_filt$has_snpArray) & ind_mdat_filt$diagn_4BrainClck == "Control"
]

dim(ind_mdat_filt[(ind_mdat_filt$has_WGS | ind_mdat_filt$has_snpArray) & ind_mdat_filt$diagn_4BrainClck == "Control", ])

table(ind_mdat_filt[(ind_mdat_filt$has_WGS | ind_mdat_filt$has_snpArray) & ind_mdat_filt$diagn_4BrainClck == "Control", ]$substudy)

sum(samp_mdat_rnaseq$individualID %in% controls_w_genomics)


indiv_w_genomics <- ind_mdat_filt$individualID[ind_mdat_filt$has_WGS | ind_mdat_filt$has_snpArray]
table(ind_mdat_filt$substudy[ind_mdat_filt$individualID %in% indiv_w_genomics])
sum(samp_mdat_rnaseq$individualID %in% indiv_w_genomics)


ind_mdat_filt[ind_mdat_filt$has_snpArray & !ind_mdat_filt$has_WGS, ]

table(ind_mdat_filt$diagnosis[!(!ind_mdat_filt$has_WGS | !ind_mdat_filt$has_snpArray)])

table(ind_mdat_filt$diagnosis[ind_mdat_filt$has_WGS | ind_mdat_filt$has_snpArray])

# Directory stuff
################################################################################
exprsn_f <- "../../results/preproc/test_no_lincs/merged_counts_log2_qnorm_noCerebell_onlyAge_svaAdj.csv"
metdat_f <- "../../results/parsed/merged/merged_metdat.csv"
trainProp <- .66
exclude_substudy <- "brainSeq_pI"
#exclude_tissue <- "cerebellum,cerebellar hemisphere"
outdir <- "../../results/parsed/merged/"

exprsn_f <- parsed$input
metdat_f <- parsed$metDat
trainProp <- as.numeric(parsed$trainProp)
exclude_substudy <- parsed$exclude_substudy
outdir <- add_slash_if_not(parsed$outdir)

create_dir_if_not(outdir)

# Load data
################################################################################
metdat <- read_table_fast(metdat_f, row.names = 1)
exprsn <- read_table_fast(exprsn_f, row.names = 1)

# Do train/test split
################################################################################

# If indicated, remove out samples from a given tissue before doing split.
#if (!is.null(exclude_tissue)){
#        exclude_tissue <- strsplit(exclude_tissue, split = ",")[[1]]
#        boolVec <- !metdat$tissue[match(rownames(exprsn),
#                                        make.names(metdat$specimenID))] %in% exclude_tissue
#        exprsn <- exprsn[boolVec, ]
#}

# Filter metadata to include only what is in the expression matrix (we filtered
# out cerebellum samples)
metdat <- metdat[match(make.names(rownames(exprsn)),
                       make.names(metdat$specimenID)), ]

# Filter metadata to keep only the controls
metdat <- metdat[metdat$diagn_4BrainClck == "Control", ]

# Keep out brainSeq_pI from the dataset, as this will be the final external
# validation
metdat <- metdat[metdat$substudy != exclude_substudy, ]

# Do train/test splitting ensuring that all the 5 year intervals are represented
# at the desired proportion.
agesBins <- 20:100
agesBins <- agesBins[20:100 %% 5 == 0]

inTrain <- c()
inTest <- c()
seed <- 222 # This seed is the good one for now, but the resulting model doesnt_work great on mice data
for(i in 1:(length(agesBins) - 1)){
        intLow <- agesBins[i]
        intHigh <- agesBins[i + 1]
        toSamp <- metdat[metdat$ageDeath >= intLow & metdat$ageDeath < intHigh, ]
        donor_counts <- toSamp %>%
                count(individualID, name = "n_samples")

        # Sort donors randomly to avoid bias
        set.seed(seed + i)
        donor_counts <- donor_counts %>% slice_sample(prop = 1)

        # Accumulate sample counts to get as close as possible to trainProp
        donor_counts <- donor_counts %>%
                mutate(cum_samples = cumsum(n_samples),
                       total_samples = sum(n_samples),
                       target_train_samples = round(total_samples * trainProp))
        split_idx <- which.min(abs(donor_counts$cum_samples - donor_counts$target_train_samples[1]))
        train_donors <- donor_counts$individualID[1:split_idx]
        test_donors  <- donor_counts$individualID[(split_idx + 1):nrow(donor_counts)]
        inTrain <- c(inTrain, train_donors)
        inTest <- c(inTest, test_donors)
}

train_test_df <- metdat[, c("specimenID",
                            "individualID",
                            "tissue",
                            "sex",
                            "race",
                            "ageDeath",
                            "diagn_4BrainClck",
                            "substudy")]
train_test_df <- train_test_df %>%
        mutate(train_test_split = case_when(
                individualID %in% inTrain ~ "train",
                individualID %in% inTest ~ "test"
        ))

write_table_fast(train_test_df, f = sprintf("%strain_test.csv", outdir))
print(sprintf("train_test.csv saved at %s", outdir))
