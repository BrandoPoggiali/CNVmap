#' This function check if .cnr and .cns files have the correct structure for the plots.
#'
#' @param cnr_data Dataframe containing .cnr file originated from the tool CNVkit.
#' @param cns_data Dataframe containing .cns file originated from the tool CNVkit.
#' @return Nothing
#' @export
check_cnvkit_data <- function(cnr_data, cns_data) {
  if (!is.data.frame(cnr_data)){
    stop(paste(substitute(cnr_data),"is not a dataframe!"), call. = FALSE)
  }
  
  if (length(cnr_data) != 7){
    stop(paste(substitute(cnr_data),"does not have 7 columns!"), call. = FALSE) 
  }
  
  if (!identical(colnames(cnr_data), c("chromosome", "start", "end", "gene", "depth", "log2", "weight"))){
    stop(paste(substitute(cnr_data), "must have the following column names in this order: chromosome, start, end, gene, depth, log2, weight"), call. = FALSE) 
  }
  
  if (!is.data.frame(cns_data)){graph
    stop(paste(cns_data,"is not a dataframe!"),call.=FALSE)
  }
  
  if (length(cns_data) != 10){
    stop(paste(substitute(cns_data),"does not have 10 columns!"),call. = FALSE) 
  }
  
  if (!identical(colnames(cns_data), c("chromosome", "start", "end", "gene", "log2", "depth", "probes", "weight", "ci_lo", "ci_hi" ))){
    stop(paste(substitute(cns_data), "must have the following column names in this order: chromosome, start, end, gene, depth, log2, weight"), call. = FALSE) 
  }
}