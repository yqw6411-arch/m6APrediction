# m6APrediction

**m6APrediction** is an R package for predicting N6-methyladenosine (m6A) modification sites in RNA sequences using a pre-trained **Random Forest** model. The package provides functions to predict m6A status and probability for single or multiple RNA samples, with example datasets included for demonstration.

## Installation

You can install the package directly from GitHub:

```r
# If you don't have remotes installed
install.packages("remotes")

# Install m6APrediction from GitHub
remotes::install_github("yqw6411-arch/m6APrediction")

Usage
Load the package:



library(m6APrediction)
Predict multiple RNA sites


# Load example model and data from the package
model_path <- system.file("extdata", "rf_fit.rds", package="m6APrediction")
data_path  <- system.file("extdata", "m6A_input_example.csv", package="m6APrediction")

ml_fit <- readRDS(model_path)
feature_df <- read.csv(data_path, stringsAsFactors = FALSE)

# Convert factor columns to match training model levels
feature_df$RNA_type <- factor(feature_df$RNA_type, levels = ml_fit$forest$xlevels$RNA_type)
feature_df$RNA_region <- factor(feature_df$RNA_region, levels = ml_fit$forest$xlevels$RNA_region)

# Predict m6A sites
res <- prediction_multiple(ml_fit, feature_df)
head(res)
Predict a single RNA site


single_res <- prediction_single(
  ml_fit,
  gc_content = 0.55,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 10,
  distance_to_junction = 5,
  evolutionary_conservation = 0.2,
  DNA_5mer = "ATGCG"
)
single_res
Functions
dna_encoding(dna_strings) – Encode DNA 5-mer sequences into categorical nucleotide positions.

prediction_multiple(ml_fit, feature_df, positive_threshold = 0.5) – Predict m6A probabilities and status for multiple RNA samples.

prediction_single(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5) – Predict m6A probability and status for a single RNA site.

Example Data
The package includes:

rf_fit.rds – Pre-trained Random Forest model.

m6A_input_example.csv – Example feature dataset.
