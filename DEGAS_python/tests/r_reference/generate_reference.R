# Generate R reference output for Python preprocess_counts test.
#
# Run from this directory:
#   module load r/4.4.1
#   Rscript generate_reference.R
#
# Output CSVs are committed to the repo; this script only needs to be
# re-run if the fixture spec changes.

source("../../../DEGAS_R/R/DEGAS_preprocessing.R")

set.seed(42)
n_genes <- 50
n_samples <- 20
# DEGAS R convention: (genes x samples). Integer counts in [0, 9999].
X <- matrix(
  sample(0:9999, n_genes * n_samples, replace = TRUE),
  nrow = n_genes, ncol = n_samples
)

write.table(X, "fixture_input.csv",
            sep = ",", row.names = FALSE, col.names = FALSE)

# R output shape is (n_samples x n_genes) — preprocessCounts transposes.
Y <- preprocessCounts(X)
stopifnot(nrow(Y) == n_samples, ncol(Y) == n_genes)

write.table(Y, "fixture_expected_output.csv",
            sep = ",", row.names = FALSE, col.names = FALSE)

cat(sprintf(
  "wrote fixture_input.csv (%d x %d) and fixture_expected_output.csv (%d x %d)\n",
  nrow(X), ncol(X), nrow(Y), ncol(Y)
))
