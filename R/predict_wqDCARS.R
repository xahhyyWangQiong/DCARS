#' Calculate DCAsignature
#'
#' This function calculates the DCAsignature based on a given gene expression matrix.
#' The expression matrix should have rows as genes and columns as samples.
#' If the matrix is missing required genes, an error will be thrown.
#'
#' @param my_expr A numeric matrix with rows as genes and columns as samples.
#' @return A data frame with sample names and their DCAsignature.
#' @export
predict_DCAsignature <- function(my_expr) {
  required_genes <- c("NRP1","ANAPC11","MYH9","FBXW9","ZNF645","COPS7A","NEDD8","GNAI1","DCUN1D3","TPM2")
  missing_genes <- required_genes[!required_genes %in% rownames(my_expr)]
  if (length(missing_genes) > 0) {
    stop("Error: The following genes are missing in the expression matrix and are required for DCAsignature calculation: ",
         paste(missing_genes, collapse = ", "), ". Please provide a complete expression matrix.")
  }
  load(system.file("data/my_coef.rdata", package = "wqDCARS"))
  DCAsignature <- cbind(colnames(my_expr), my_DCAsignature=colSums(my_expr[required_genes,]*my_coef$coef))
  return(DCAsignature)
}


