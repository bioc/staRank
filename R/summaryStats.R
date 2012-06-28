#' @title Summary statistic per element
#' 
#' @description Calculates the summary statistic per element for the whole dataset.
#' 
#' @param data a matrix with one row per element or a list containing one
#' vector per element. 
#' @param method the score that is calculated per gene, one of 'mean' (default), 
#' 'median', 'mwtest' (two-sample one sided Mann-Whitney test), 
#' 'ttest'(two-sample one sided t-test), 'RSA' (redundant siRNA activity).
#' @param decreasing a boolean indicating the direction of the ranking.
#' 
#' @return a named vector of the scored elements.
#' 
#' @export 
#' @docType methods
#' @rdname summaryStats-methods
#' 
#' @example inst/example/staRank-example.R
#' 
setGeneric("summaryStats", function(data,method='mean',decreasing=FALSE) 
			standardGeneric('summaryStats')
)


#' @rdname summaryStats-methods
#' @aliases summaryStats,matrix-method
setMethod("summaryStats",signature=signature(data = "matrix"),
		function(data,method,decreasing){
			scores<-0
			a<-'less'
			if(decreasing){a<-'greater'}
			if(dim(data)[1]==1){method<-'sort'}		
			if(method=='median') {
				scores<-apply(data,1,function(x){median(x,na.rm=TRUE)})
			}
			if(method=='mean') {
				scores<-apply(data,1,function(x){mean(x,na.rm=TRUE)})
			}
			if(method=='mwtest'){
				d<-c(data);
				scores<-apply(data,1,function(x){
					#wilcox.test(x, d, alternative = a)$p.value}) # slow!
					mwTest2samp(x, d, alternative = a)$p.value})
			}
			if(method=='ttest'){
				d<-c(data); 
				scores<-apply(data,1,function(x){
					t.test(x, d, alternative = a,var.equal=FALSE)$p.value})
			}
			if(method=='RSA') {
				opts<-list(LB=min(c(data)),UB=max(c(data)),outputFile=NA,
						inputFile=NA,reverse=decreasing)
				scores<-uniqueRSARanking(dataFormatRSA(data),opts)
				scores<-scores[rownames(data)]
			}
			if(method=='sort'){
				if(dim(data)[1]==1){
					x<-colnames(data)
					scores<-c(data)
					names(scores)<-x
				}else{
					scores<-apply(data,1,function(x){mean(x)})
				}
			}
			# ... to be continued for further summary values
			scores<-sapply(scores,function(x){ifelse(is.na(x),NA,x)})
			return(scores)
		})


#' @rdname summaryStats-methods
#' @aliases summaryStats,list-method
setMethod("summaryStats",signature=signature(data = "list"),
		function(data,method,decreasing){
			scores<-0
			a<-'less'
			if(decreasing){a<-'greater'}			
			if(method=='median') {
				scores<-sapply(data,function(x){median(x,na.rm=TRUE)})
			}
			if(method=='mean') {
				scores<-sapply(data,function(x){mean(x,na.rm=TRUE)})
			}
			if(method=='mwtest'){
				scores<-sapply(data,function(x){
					#wilcox.test(x, c(data), alternative = a)$statistic})
					mwTest2samp(x, c(data), alternative = a)$p.value})
			}
			if(method=='ttest'){
				scores<-sapply(data,function(x){
					#t.test(x, c(data), alternative = a,var.equal=FALSE)$statistic})
					t.test(x, c(data), alternative = a,var.equal=FALSE)$p.value})
			}
			if(method=='RSA') {
				opts<-list(LB=min(c(data)),UB=max(c(data)),outputFile=NA,
						inputFile=NA,reverse=decreasing)
				scores<-uniqueRSARanking(dataFormatRSA(data),opts)
				scores<-scores[rownames(data)]
			}
			if(method=='sort'){
				scores<-sapply(data,function(x){mean(x)})
			}
			# ... to be continued for further summary values
			
			return(scores)
		})
