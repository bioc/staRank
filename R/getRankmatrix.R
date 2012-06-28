#' getRankmatrix
#' 
#' Extracts all three rankings from a RankSummary object.
#' 
#' @param rs an object of class RankSummary.
#' 
#' @return a ranking matrix where each row corresponds to one element
#' and the columns give the ranks from stability, base and averaged ranking.
#' 
#' @example inst/example/staRank-example.R
#' @export 
#' 
getRankmatrix <-
function(rs)
{
	ranks<-matrix(NA,nrow=length(rs@stabRank), ncol=3)
	ranks[,1]<-rs@stabRank
	rownames(ranks)<-names(rs@stabRank)
	ranks[,2]<-rank(rs@baseRank,na.last='keep')[rownames(ranks)]
	ranks[,3]<-rank(rs@avrgRank,na.last='keep')[rownames(ranks)]
	colnames(ranks)<-paste(rs@method,c('stability','base','averaged'),sep='_')
	
	return(ranks)
}

