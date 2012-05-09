#' uniqueRSARanking
#' 
#' Performs RSA analysis and summarizes the results gene wise.
#' 
#' @param dataRSA a matrix with data in RSA format (can be created with the
#' function dataFormatRSA).
#' @param opts the options for the RSA ranking. This is a list of: 
#' LB: lower_bound (defaults = 0), UB: upper bound (defaults = 1)
#' reverse: boolean (if TRUE: reverse hit picking, higher scores are better,
#' default = FALSE).
#' 
#' @return a vector of ranked genes with their RSA LogP-values.
#' 
#' @example inst/example/staRank-example.R
#' @export 
#' 
uniqueRSARanking <-
function(dataRSA,opts){
	names(dataRSA)<-c("Gene_ID","Well_ID","Score")
	ranking <- runRSA(dataRSA,opts)
	#make RSA output ranking unique
	ranking <- ranking[!(duplicated(ranking$Gene_ID)),c('Gene_ID','LogP')]
	r<-ranking$LogP
	names(r)<-ranking$Gene_ID
	return(r)
}