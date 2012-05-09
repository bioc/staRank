#' aggregRank
#'
#' Performs rank aggregation on a chosen summary statistic that is applied on a 
#' matrix of (columnwise) rankings.
#' 
#' @param rankMatrix a matrix where each column is a ranking and the rows 
#' correspond to the elements.
#' 
#' @param method one of c('mean','median','max','min') that is used to calulate
#' the joint value of one element.
#' 
#' @return the aggregated ranking as a named vector.
#' 
#' @example inst/example/staRank-example.R
#' @export 
#' 
aggregRank <-
function(rankMatrix,method='mean')
{
	r<-apply(rankMatrix,1,function(x){
				switch(method,
						'mean'=mean(x,na.rm=TRUE),
						'median'=median(x,na.rm=TRUE),
						'max'=max(x,na.rm=TRUE),
						'min'=min(x,na.rm=TRUE))}
	)
	return(sort(r))
}