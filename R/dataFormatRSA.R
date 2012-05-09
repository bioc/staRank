#' dataFormatRSA
#' 
#' Creates a data.frame in the RSA input format from a matrix, where rows
#' indicate the genes and columns the replicate or siRNA values.
#' 
#' @param dataMatrix a matrix with one row per gene and columns are replicate
#' or siRNA values.
#' 
#' @return a data.frame with three columns Gene_ID, Well_ID, Score. The Well_ID
#' does not have a meaning if data are replicates, but is usefull to 
#' distinguish between siRNAs.
#' 
#' @example inst/example/staRank-example.R
#' @export 
#' 
dataFormatRSA <-
function(dataMatrix){
	d<-data.frame(Gene_ID=rep(rownames(dataMatrix),ncol(dataMatrix)),
			Well_ID=1:length(dataMatrix),
			Score=c(dataMatrix),stringsAsFactors = FALSE)
	return(d)
}