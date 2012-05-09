#' getStability
#' 
#' Performes stability selection on sample rankings with a given stability
#' threshold. Selection probabilities and stability ranking are calculated.
#' 
#' @param sr sample rankings in the form of a matrix where each row corresponds 
#' to one element and each column gives one ranking.
#' @param thr the threshold for the stability selection, indicating above which
#' frequency in the samples an element is considered stable.
#' 
#' @return a list containing:
#' \describe{
#' 		\item{\code{stabRank}}{the stable ranking.}
#' 		\item{\code{Pi}}{the frequency matrix with all values per gene and per 
#' cutoff.}
#' 		\item{\code{stableSetSize}}{a table with the number of stable genes per 
#' cutoff.}
#' 	}
#' @example inst/example/staRank-example.R
#' @export
#' 
getStability <-
function(sr,thr) 
{
	n<-ncol(sr)		# sample number
	p<-nrow(sr)		# parameter number (genes)
	
	Pi<-matrix(0,nrow=p,ncol=p)
	rownames(Pi)<-rownames(sr)
	setSize<-rep(0,p)
	rank<-c()
	for(k in 1:p){
		Pi[,k]<-(apply(sr<=k,1,sum))/n	# selection Pr. for each gene at cutoff k		
		stabK<-names(which(Pi[,k]>=thr))
		setSize[k]<-length(stabK)		
		rank<-c(rank,setdiff(stabK,rank))
	}
		
	stabRank<-1:p
	names(stabRank)<-rank
	return(list(stabRank=stabRank,stableSetSize=setSize,Pi=Pi))
}

