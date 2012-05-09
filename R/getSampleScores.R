#' @title Get Sample Scores
#' 
#' @description an S4 method to a sample dataset that is scored by a given method.
#' 
#' @param data can be a matrix with one row per element or a list of vectors
#' of different length, one for each element.
#' @param method one of the ranking methods: 'mean' (default), 'median', 
#' 'mwtest' (two sample one sided mann-whitney test), 
#' 'ttest'(two sample one sided t-test).
#' @param decreasing a boolean indicating the direction of the ranking.
#' @param bootstrap a boolean indicating whether bootstrapping or subsampling 
#' is used.
#' 
#' @return a vector containing the summary values by the given method for the 
#' sample dataset.
#' 
#' @export 
#' @docType methods
#' @rdname getSampleScores-methods
#' 
#' @example inst/example/staRank-example.R
#'
setGeneric("getSampleScores",function(data,method='mean',decreasing=FALSE,
				bootstrap=TRUE) standardGeneric("getSampleScores"))

#' @rdname getSampleScores-methods
#' @aliases getSampleScores,matrix-method
setMethod("getSampleScores", signature=signature(data='matrix'),
	function(data,method,decreasing,bootstrap) {
		c<-ifelse(bootstrap,ncol(data),ceiling(ncol(data)/2))
		s<-t(apply(data,1,function(r){sample(r,c,replace=bootstrap)}))		
		return(summaryStats(s,method,decreasing))
})

#' @rdname getSampleScores-methods
#' @aliases getSampleScores,list-method
setMethod("getSampleScores", signature=signature(data='list'),
		function(data,method,decreasing,bootstrap){
			s<-Reduce(lapply(data,function(r){ifelse(length(r)==1,r,function(r){
					c<-ifelse(bootstrap,length(r),ceiling(ncol(data)/2))
					return(sample(r,size=c,replace=bootstrap))
					})}),rbind
			)
			rownames(s)<-names(data)
			return(summaryStats(s,method,decreasing))
})

