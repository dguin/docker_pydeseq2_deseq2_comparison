# R script to run DESeq2 on count data.
# Usage Rscript deseq2.R $input_directory $out_directory
# Input directory must contain counts.csv with sample names as columns and gene
# names as indices and metadata.csv with index as the sample names that match
# column names of counts. Columns expected = ['condition'] which annotates
# control or test per sample

# suppressMessages(
#   library("DESeq2")
# )
library("DESeq2")

args <- commandArgs(trailingOnly = TRUE)
input_directory <- args[1]
out_directory <- args[2]

# Load counts of the shape genesXsamples
countData <- read.csv(file.path(input_directory, "counts.csv"),
  header = TRUE, sep = ","
)
# Load metadata dataframe - index must be sample names that match column names
# of counts. Columnd expected = ['condition'] which annotates control or test
# per sample
metaData <- read.csv(file.path(input_directory, "metadata.csv"),
  header = TRUE, sep = ","
)

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = metaData,
  design = ~condition, tidy = TRUE
)

dds <- DESeq(dds)
size_factors <- sizeFactors(dds)
res <- results(dds, name = "condition_test_vs_control")
write.csv(res, file.path(out_directory, "r_results_pre_shrunk.csv"),
  row.names = TRUE
)
write.csv(size_factors, file.path(out_directory, "r_size_factors.csv"),
  row.names = TRUE
)
# or to shrink log fold changes association with condition:
print("shrinking")
res <- lfcShrink(dds, coef = "condition_test_vs_control")
write.csv(res, file.path(out_directory, "r_results_post_shrunk.csv"),
  row.names = TRUE
)
