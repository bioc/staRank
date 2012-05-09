#
#	RankSummary class definition and standard methods for it
#

#' @title Class RankSummary
#' 
#' @description A class containing different rankings of a dataset further 
#' summary information like the stable set sizes, a correlation matrix of 
#' the rankings, the name of the base ranking method,
#' the function call and if needed the Pi matrix of the stabilityRanking.
#' 
#' \bold{Slots}
#'	\describe{
#'		\item{\code{baseRank}:}{a numeric vector containing the base ranking}
#'		\item{\code{stabRank}:}{a numeric vector containing the stability ranking}
#'		\item{\code{avrgRank}:}{a numeric vector containing the average ranking}
#'		\item{\code{stableSetSize}:}{the sizes of the stable sets for each cutoff}
#' 		\item{\code{rankCor}}{Spearman's rank correlation coefficient for the 
#' 			three rankings}
#' 		\item{\code{method}}{the used method}
#' 		\item{\code{Pi}:}{a matrix containing}
#' }
#' 
#' \bold{Accessors}
#'	\itemize{
#'		\item \code{\link{baseRank}}
#'		\item \code{\link{stabRank}}
#'		\item \code{\link{avrgRank}}
#'		\item \code{\link{stableSetSize}}
#' 		\item \code{\link{rankCor}}
#' 		\item \code{\link{method}}
#' 		\item \code{\link{Pi}}
#' }
#' 
#' \bold{Methods}
#'	\itemize{
#'		\item \code{\link{show}} 
#' 		\item \code{\link{summary}}
#' }
#' 
#' @export
#' @docType class
#' @name RankSummary-class
#' @rdname RankSummary-class
#' @aliases RankSummary-class
#' 
#' @examples 
#' # overview of RankSummary class
#' showClass("RankSummary")
#' # generate dataset
#' d<-replicate(4,sample(1:10,10,replace=FALSE))
#' rownames(d)<-letters[1:10]
#' # run default stability ranking
#' s<-stabilityRanking(d)
#' # using an accessor functions on the RankSummary object
#' stabRank(s)
#' rankCor(s)
#' 
#' @keywords classes
#' 
setClass("RankSummary", 
		representation(baseRank="numeric",stabRank="numeric",
				avrgRank="numeric",stableSetSize="numeric", rankCor="matrix",
				method="character",Pi="matrix")
)


#' \bold{Constructor for class RankSummary}
#' 
#' RankSummary objects can be created using the constructor \code{RankSummary}, 
#' however this is most likely not needed.
#' 
#' @param baseRank a numeric vector containing the base ranking
#' @param stabRank a numeric vector containing the stability ranking
#' @param avrgRank a numeric vector containing the average ranking 
#' @param stableSetSize the sizes of the stable sets for each cutoff
#' @param rankCor Spearman's rank correlation coefficient for the 
#' 			three rankings
#' @param method the used method 
#' @param Pi a matrix containing
#' 
#' @export
#' @rdname RankSummary-class
#' 
#' @examples 
#' # create a new empty object of class RankSummary using the constructor
#' RankSummary()
RankSummary<-function(baseRank=0,stabRank=0,avrgRank=0,
		stableSetSize=0,rankCor=matrix(NA,nrow=3,ncol=3),method="",
		Pi=matrix(NA,nrow=1,ncol=1))
{new("RankSummary", baseRank=baseRank,stabRank=stabRank,avrgRank=avrgRank,
			stableSetSize=stableSetSize, rankCor=rankCor,method=method,Pi=Pi)}


#' @title Get base ranking
#' 
#' @description Accessor function to obtain the base ranking from 
#' 	RankSummary objects.
#'
#' @param x a RankSummary object
#' @return a named vector with the base ranking
#' 
#' @export
#' @docType methods
#' @rdname baseRank-methods
setGeneric("baseRank",function(x) standardGeneric("baseRank"))

#' @rdname baseRank-methods
#' @aliases baseRank,RankSummary-method
setMethod("baseRank",signature=signature(x="RankSummary"),function(x) x@baseRank)


#' @title Get stability ranking
#' 
#' @description Accessor function to obtain the stability ranking from 
#' 	RankSummary objects.
#' 
#' @param x a RankSummary object
#' @return a named vector with the stability ranking
#' 
#' @export
#' @docType methods
#' @rdname stabRank-methods
setGeneric('stabRank',function(x) standardGeneric('stabRank'))

#' @rdname stabRank-methods
#' @aliases stabRank,RankSummary-method
setMethod('stabRank','RankSummary',function(x) x@stabRank)


#' @title Get average ranking
#' 
#' @description Accessor function to obtain the average ranking from 
#' 	RankSummary objects.
#' 
#' @param x a RankSummary object
#' @return a named vector with the average ranking
#' 
#' @export 
#' @rdname avrgRank-methods
#' @aliases avrgRank
setGeneric('avrgRank',function(x) standardGeneric('avrgRank'))

#' @rdname avrgRank-methods
#' @aliases avrgRank,RankSummary-method
setMethod('avrgRank','RankSummary',function(x) x@avrgRank)


#' @title Get stable set sizes
#' 
#' @description Accessor function to obtain the sizes of stable sets per cutoff 
#' 	from RankSummary objects.
#' 
#' @param x a RankSummary object
#' @return a vector with the per cutoff stable set size. 
#' 
#' @export
#' @rdname stableSetSize-methods
#' @aliases stableSetSize
setGeneric('stableSetSize',function(x) standardGeneric('stableSetSize'))

#' @rdname stableSetSize-methods
#' @aliases stableSetSize,RankSummary-method
setMethod('stableSetSize','RankSummary',function(x) x@stableSetSize)


#' @title Get rank correlations
#' 
#' @description Accessor function to obtain the rank correlation matrix from 
#' RankSummary objects.
#' 
#' @param x a RankSummary object
#' @return A matrix with Spearman's rank correlation coefficients of the three
#' rankings contained in a RankSummary object.
#' 
#' @export
#' @rdname rankCor-methods
#' @aliases rankCor
setGeneric('rankCor',function(x) standardGeneric('rankCor'))

#' @rdname rankCor-methods
#' @aliases rankCor,RankSummary-method
setMethod('rankCor','RankSummary',function(x) x@rankCor)


#' @title Get method
#' 
#' @description Accessor function to obtain the name of the base ranking method 
#' from RankSummary objects.
#' 
#' @param x a RankSummary object
#' @return a string with the name of the base ranking method
#' 
#' @export
#' @rdname method-methods
#' @aliases method
setGeneric('method',function(x) standardGeneric('method'))

#' @rdname method-methods
#' @aliases method,RankSummary-method
setMethod('method','RankSummary',function(x) x@method)


#' @title Get Pi
#' 
#' @description Accessor function to obtain the Pi matrix from 
#' 	RankSummary objects.
#' 
#' @param x a RankSummary object
#' @return a matrix with the per cutoff stability values for each gene.
#' 
#' @export
#' @rdname Pi-methods
#' @aliases Pi
setGeneric('Pi',function(x) standardGeneric('Pi'))

#' @rdname Pi-methods
#' @aliases Pi,RankSummary-method
setMethod('Pi','RankSummary',function(x) x@Pi)


#' @title Show function for a RankSummary object.
#' 
#' @description This prints the number of ranked elements, the base ranking 
#' method, the top elements of the three rankings and the top stable set sizes 
#' from a RankSummary object.
#' @param object an object of class RankSummary.  
#' 
#' @rdname show-methods
#' @aliases show,RankSummary-method 
setMethod("show",
		signature = signature(object="RankSummary"), 
		function(object){
			cat("\nStability Ranking of:",length(object@baseRank),'elements')
			cat("\nbase ranking method:",object@method,'\n')
			cat("\nTop base ranking:\n")
			print(head(object@baseRank))
			cat("\nTop stability ranking:\n")
			print(head(object@stabRank))
			cat("\nTop averaged ranking:\n")
			print(head(object@avrgRank))
			cat("\nTop stability:\n")
			print(head(object@stableSetSize))
		}
)


#' @title Summary method for RankSummary object
#' 
#' @description This function summarizes the top 10 elements for each of the 
#' three ranking types 'stability', 'base', 'averaged'. It also contains their 
#' (spearman) rank correlation coefficients and stability information about 
#' the top 1,10 and 50% of the ranking.
#' 
#' @param object a RankSummary object.  
#' @return a summary of a 'RankSummary' object.
#' 
#' @example inst/example/staRank-example.R
#' 
#' @rdname summary-methods
#' @aliases summary,RankSummary-method
setMethod("summary",
		signature = signature(object="RankSummary"),
		function(object){
			n<-length(baseRank(object))
			t<-min(n,10)
			top<-cbind(names(baseRank(object)[1:t]),names(stabRank(object)[1:t]),
					names(avrgRank(object)[1:t]))
			colnames(top)<-c('base','stability','averaged')
			perc<-c(1,10,50)
			cutoff<-ceiling(perc/100*length(stableSetSize(object)))
			res <- list(method=method(object),
					topElements=top,
					stabilityInfo=cbind(percent=perc,cutoff=cutoff,
							setSize=stableSetSize(object)[cutoff]),
					elementNumber=n,
					correlation=rankCor(object))
			return(res)
		}
)
