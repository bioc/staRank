#' getRankmatrix
#' 
#' Extracts all three rankings from a RankSummary object.
#' 
#' @param sr an object of class RankSummary.
#' 
#' @return a ranking matrix where each row corresponds to one element
#' and the columns give the ranks from stability, base and averaged ranking.
#' 
#' @example inst/example/staRank-example.R
#' @export 
#' 
getRankmatrix <-
function(sr)
{
	ranks<-cbind(stability=sr@stabRank,
			base=match(names(sr@stabRank),names(sr@baseRank)),	
			averaged=match(names(sr@stabRank),names(sr@avrgRank)))
	rownames(ranks)<-names(sr@stabRank)
	colnames(ranks)<-paste(sr@method,colnames(ranks),sep='_')	
	return(ranks)
}