#' Encode DNA sequences into categorical features
#'
#' This function converts a vector of DNA strings into a data frame of categorical
#' nucleotide positions (A, T, C, G), which can be used as input features for
#' machine learning models.
#'
#' @param dna_strings A character vector containing DNA sequences (e.g., 5-mers)
#'
#' @return A data frame where each column represents a nucleotide position
#' (e.g., nt_pos1, nt_pos2, ...) and each row corresponds to an input sequence.
#'
#' @examples
#' # 简单示例：直接编码两条 5-mer 序列
#' dna_encoding(c("ATGCG", "TTGCA"))
#'
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Make multiple m6A site predictions using a trained model
#'
#' Predict m6A modification probabilities and statuses for multiple RNA samples.
#' Automatically encodes the DNA 5-mer sequences and appends the predictions.
#'
#' @param ml_fit A trained random forest model object
#' @param feature_df A data frame containing features:
#' gc_content, RNA_type, RNA_region, exon_length,
#' distance_to_junction, evolutionary_conservation, DNA_5mer
#' @param positive_threshold Numeric (0–1) classification threshold
#'
#' @return A data frame identical to input \code{feature_df} but with two
#' additional columns: \code{predicted_m6A_prob} and \code{predicted_m6A_status}.
#'
#' @examples
#' \dontrun{
#' model_path <- system.file("extdata", "rf_fit.rds", package="m6APrediction")
#' data_path  <- system.file("extdata", "m6A_input_example.csv", package="m6APrediction")
#' ml_fit <- readRDS(model_path)
#' feature_df <- read.csv(data_path, stringsAsFactors = FALSE)
#'
#' # 转 factor 并指定 levels 与训练模型一致
#' feature_df$RNA_type <- factor(feature_df$RNA_type, levels = ml_fit$forest$xlevels$RNA_type)
#' feature_df$RNA_region <- factor(feature_df$RNA_region, levels = ml_fit$forest$xlevels$RNA_region)
#'
#' # 调用函数
#' res <- prediction_multiple(ml_fit, feature_df)
#' head(res)
#' }
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length",
                  "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
                %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = ml_fit$forest$xlevels$RNA_type)
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = ml_fit$forest$xlevels$RNA_region)
  feature_df <- cbind(feature_df[, -which(names(feature_df) == 'DNA_5mer')],
                      dna_encoding(feature_df$DNA_5mer))
  predicted_prob <- predict(ml_fit, newdata = feature_df, type = "prob")[, "Positive"]
  feature_df$predicted_m6A_prob <- predicted_prob
  feature_df$predicted_m6A_status <- ifelse(predicted_prob > positive_threshold,
                                            "Positive", "Negative")
  return(feature_df)
}

#' Make a single m6A site prediction
#'
#' Predict the m6A probability and class for a single RNA site.
#'
#' @param ml_fit A trained random forest model object
#' @param gc_content Numeric GC content
#' @param RNA_type Character: "mRNA", "lincRNA", "lncRNA", "pseudogene"
#' @param RNA_region Character: "CDS", "intron", "3'UTR", "5'UTR"
#' @param exon_length Numeric exon length
#' @param distance_to_junction Numeric distance to exon junction
#' @param evolutionary_conservation Numeric conservation score
#' @param DNA_5mer Character 5-mer sequence
#' @param positive_threshold Numeric classification threshold (0–1)
#'
#' @return Named vector: predicted_m6A_prob and predicted_m6A_status
#'
#' @examples
#' \dontrun{
#' model_path <- system.file("extdata", "rf_fit.rds", package="m6APrediction")
#' ml_fit <- readRDS(model_path)
#'
#' # 示例调用
#' res <- prediction_single(
#'   ml_fit,
#'   gc_content = 0.5,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 100,
#'   distance_to_junction = 10,
#'   evolutionary_conservation = 0.2,
#'   DNA_5mer = "ATGCG"
#' )
#' res
#' }
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5){
  provided_feature <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(RNA_type, levels = ml_fit$forest$xlevels$RNA_type),
    RNA_region = factor(RNA_region, levels = ml_fit$forest$xlevels$RNA_region),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )
  result_df <- prediction_multiple(ml_fit, provided_feature, positive_threshold)
  returned_vector <- c(
    predicted_m6A_prob = result_df$predicted_m6A_prob[1],
    predicted_m6A_status = result_df$predicted_m6A_status[1]
  )
  return(returned_vector)
}
