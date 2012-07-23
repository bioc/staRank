#' @title Wrapper to perform stabilityRanking
#' 
#' @description An S4 function to perform stability ranking on a dataset or 
#' directly on given sample rankings.
#' 
#' @param x the data that is to be ranked. It can be of different types:
#' a data matrix with one row per element to be ranked, or a ranking from an 
#' external method, or an object of class cellHTS. Depending on this the other
#' parameters can vary.
#' @param samps a matrix of scored sample data, each row corresponds to an 
#' element, the columns to a scoring. This is only needed when an external 
#' ranking method is used.
#' @param channel a string with the name of the feature (channel) to be ranked.
#' This is only needed for cellHTS objects.
#' @param replicates names or indices of the replicates (samples) to be used 
#' for the rankings (default: all samples are used).
#' @param method one of the ranking methods: 'mean' (default), 'median', 
#' 'mwtest' (two-sample one sided Mann-Whitney test), 
#' 'ttest'(two-sample one sided t-test) or 'RSA' (redundant siRNA analysis).
#' If an external ranking is used, you can specify the name of that ranking
#' method in the method argument.
#' @param decreasing a boolean indicating the direction of the ranking.
#' @param bootstrap a boolean indicating if bootstrapping or subsampling is used.
#' @param thr threshold for stability (default = 0.9).
#' @param nSamp the number of samples to generate (default = 100).
#' @param Pi boolean indicating if the Pi matrix should be returned
#' (can be very large, default=FALSE).
#' @param verbose boolean indicating if status update should be printed
#' @param ... further parameter for the stability ranking.
#' 
#' @return an object of class \code{\link{RankSummary}}.
#' 
#' @export 
#' @docType methods
#' @rdname stabilityRanking-methods
#' 
#' @example inst/example/staRank-example.R 

setGeneric('stabilityRanking', function(x,...) 
			standardGeneric('stabilityRanking'))


#' @rdname stabilityRanking-methods
#' @aliases stabilityRanking,numeric-method
setMethod("stabilityRanking", signature=signature(x = "numeric"),
function(x, samps , method = 'mean', decreasing = FALSE, bootstrap = TRUE,
		Pi = FALSE, thr = 0.9, verbose = FALSE)
{
	# replace values by ranks
	if(decreasing){
		sr<-apply(samps,2,function(s){
					rank(-s,ties.method='random',na.last='keep')})#NA at end
	}else{
		sr<-apply(samps,2,function(s){
					rank(s,ties.method='random',na.last='keep')})
	}
	
	gs<-getStability(sr,thr,Pi=Pi,verbose=verbose)
	s<-RankSummary(
			baseRank=x,
			stabRank=gs$stabRank,
			avrgRank=aggregRank(sr,'mean'), 
			stableSetSize=gs$stableSetSize,
			rankCor=as.matrix(NA),
			method=method
	)
	
	if(decreasing){avrgRank(s)<-rev(avrgRank(s))}
	s@rankCor<-cor(cbind(stability=stabRank(s),base=match(names(stabRank(s)),
							names(baseRank(s))),average=match(names(stabRank(s)),
							names(avrgRank(s)))),method='spearman')
	# reset Pi if not needed
	if(Pi){
		s@Pi<-gs$Pi
	}#else{s@Pi<-matrix(NA)}
	rm(gs)
	return(s)
})


#' @rdname stabilityRanking-methods
#' @aliases stabilityRanking,matrix-method
setMethod("stabilityRanking", signature=signature(x= "matrix"),
function(x, method = 'mean',decreasing = FALSE, bootstrap = TRUE,
		Pi = FALSE, thr = 0.9, nSamp = 100, verbose = FALSE)
{
	# rank data by method
	dataRank<-sort(summaryStats(x,method,decreasing),decreasing=decreasing,
			na.last = TRUE)
	# create samples
	if(verbose){cat('generating bs-samples ...')}
	samps<-replicate(nSamp,getSampleScores(x,method,decreasing,bootstrap))
	if(verbose){cat('done\n')}
	if(verbose){cat('running stability selection...\n')}
	s<-stabilityRanking(dataRank,samps,method=method,decreasing=decreasing,Pi=Pi,
			thr=thr,verbose=verbose)
	if(verbose){cat('done\n')}
	rm(samps)
	return(s)
})


#' @rdname stabilityRanking-methods
#' @aliases stabilityRanking,cellHTS-method
#'
 
# examples
# the method buildCellHTS2() seems not to exist anymore
# wells = sprintf("%s%02d", rep(LETTERS[1:8], each=12), 1:12)
# xd = expand.grid(well=wells, plate=1:3, replicate=1:4)
# xd$cell.number = rnorm(nrow(xd))
# xd$cell.size = rnorm(nrow(xd))
# x = buildCellHTS2(xd)
# GeneID = paste('gene_',1:(nrow(xd)/4),sep='')
# featureNames(x) = GeneID
# s<-stabilityRanking(x,channel='cell.number',replicates=1:4)
 
setMethod("stabilityRanking", signature=signature(x= "cellHTS"),
function(x, channel, replicates=NULL, method = 'mean',decreasing = FALSE, 
		bootstrap = TRUE, Pi = FALSE, thr = 0.9, nSamp = 100, verbose = FALSE)
{
	# check that selected channel is available
	if(!channel%in%channelNames(x)){
		stop('Please select valid channel')
	}
	
	# get the replicates or use all
	if(is.null(replicates)){samples<-sampleNames(x)
	}else{samples<-as.character(replicates)}
	
	# check that selected replicates are available
	if(!all(samples%in%sampleNames(x))){
		stop('Please select valid replicates')
	}
	
	# check that at least 2 to replicates are chosen
	if(!all(samples%in%sampleNames(x))){
		stop('At least two replicate measurements are needed for a stability
						analysis')
	}
	# extract assay data from 
	data<-Data(x)
		
	# get the selected channel
	d<-data[,,channel]
	# attach gene names
	rownames(d)<-featureNames(x)
	# run stabilityRanking on the matrix
	s<-stabilityRanking(d,method=method, decreasing=decreasing, 
			bootstrap=bootstrap,Pi=Pi, thr=thr, nSamp=nSamp, verbose=verbose)
	return(s)
})
