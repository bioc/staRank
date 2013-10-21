#' getStability
#' 
#' Performes stability selection on sample rankings with a given stability
#' threshold. Selection probabilities and stability ranking are calculated.
#' 
#' @param sr sample rankings in the form of a matrix where each row corresponds 
#' to one element and each column gives one ranking.
#' @param thr the threshold for the stability selection, indicating above which
#' frequency in the samples an element is considered stable.
#' @param Pi boolean indicating if the Pi matrix should be returned
#' (can be very large, default=FALSE).
#' @param verbose boolean indicating whether status updates should be printed
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
function(sr,thr,Pi=FALSE,verbose=FALSE) 
{
	n<-ncol(sr)		# sample number
	p<-nrow(sr)		# parameter number (genes)
	
	if(Pi){
		setSize<-rep(0,p)
		Pimat<-matrix(0,nrow=p,ncol=p)
		rownames(Pimat)<-rownames(sr)
		# new version using less memory
		stab<-c()
		for(g in 1:p){
			Pimat[g,1]<-(sum(sr[g,]<=1,na.rm=TRUE)/n) # selection Pr. for each gene at cutoff 1
			if(Pimat[g,1]>=thr){stab<-c(stab,g)}
		}
		setSize[1]<-length(stab)	
		if(setSize[1]==0){rank<-c()}else{rank<-stab}
		
		if(verbose){
			pks<-1:10*p/10
			for(k in 2:p){
				if(k %in% pks) cat('cutoff',k,'of',p,'\n')
				stabK<-c()
				for(g in 1:p){
				Pimat[g,k]<-sum(sr[g,]<=k,na.rm=TRUE)/n # selection Pr. for each gene at cutoff k			
				if(Pimat[g,k]>=thr & Pimat[g,k-1]<thr){stabK<-c(stabK,g)}
				}
				setSize[k]<-setSize[k-1]+length(stabK)		
				rank<-c(rank,stabK)
			}
		}else{
			for(k in 2:p){
				stabK<-c()
				for(g in 1:p){
					Pimat[g,k]<-sum(sr[g,]<=k,na.rm=TRUE)/n # selection Pr. for each gene at cutoff k			
					if(Pimat[g,k]>=thr & Pimat[g,k-1]<thr){stabK<-c(stabK,g)}
				}
				setSize[k]<-setSize[k-1]+length(stabK)		
				rank<-c(rank,stabK)
			}
		}

## old version using much memory
#	for(k in 1:p){
#		if(verbose){
#			if(k %in% pks) cat('cutoff',k,'of',p,'\n')
#		}
#		Pimat[,k]<-(apply(sr<=k,1,sum))/n # selection Pr. for each gene at cutoff k		
#		stabK<-names(which(Pimat[,k]>=thr))
#		setSize[k]<-length(stabK)		
#		rank<-c(rank,setdiff(stabK,rank))
#	}
		nas<-rownames(sr)[which(is.na(Pimat[,p]))]
		stabRank<-c(1:(length(rank)),rep(NA,length(nas)))
		names(stabRank)<-c(rownames(sr)[rank],nas)		
		
		return(list(stabRank=stabRank,stableSetSize=setSize,Pi=Pimat))
		
	}else{
		# genewise
		stabAtK<-rep(NA,p)
		setSize<-rep(0,p)
		if(verbose){
			ng<-c(1,floor(1:10*p/10))
			for(g in 1:p){	
				if(g %in% ng) cat('estimating gene',g,'of',p,'\n')
				k<-1		  	# current cutoff
				stable<-FALSE 	#is this gene stable yet
				rgs<-sr[g,] 	# all ranks of one gene
				rsK<-which(rgs<=1) # get all ranks <=1
				PiPrev<-length(rsK)/n	# get percentage
				if(PiPrev>=thr){
					stabAtK[g]<-k
					setSize[k:p]<-setSize[k:p]+1
					stable=TRUE
				}
				if(length(rsK)>0){rgs<-rgs[-rsK]}	#remove already counted ranks
				while(stable==FALSE & k<p){
					k<-k+1			# increase cutoff
					rsK<-which(rgs<=k) # get all ranks <=k
					PiK<-PiPrev + length(rsK)/n #sum percentages
					if(length(rsK)>0){rgs<-rgs[-rsK]}
					if(PiK>=thr){
						stabAtK[g]<-k
						setSize[k:p]<-setSize[k:p]+1
						stable=TRUE
					}
					PiPrev<-PiK
				}
			}
		}else{
			for(g in 1:p){	
				k<-1		  	# current cutoff
				stable<-FALSE 	#is this gene stable yet
				rgs<-sr[g,] 	# all ranks of one gene
				rsK<-which(rgs<=1) # get all ranks <=1
				PiPrev<-length(rsK)/n	# get percentage
				if(PiPrev>=thr){
					stabAtK[g]<-k
					setSize[k:p]<-setSize[k:p]+1
					stable=TRUE
				}
				if(length(rsK)>0){rgs<-rgs[-rsK]}	#remove already counted ranks
				while(stable==FALSE & k<p){
					k<-k+1			# increase cutoff
					rsK<-which(rgs<=k) # get all ranks <=k
					PiK<-PiPrev + length(rsK)/n #sum percentages
					if(length(rsK)>0){rgs<-rgs[-rsK]}
					if(PiK>=thr){
						stabAtK[g]<-k
						setSize[k:p]<-setSize[k:p]+1
						stable=TRUE
					}
					PiPrev<-PiK
				}
			}
		}
		
		names(stabAtK)<-rownames(sr)
		stabAtK<-sort(stabAtK,na.last=TRUE)
		stabRank<-rank(stabAtK,na.last='keep',ties.method='average')
		
		return(list(stabRank=stabRank,stableSetSize=setSize,Pi=matrix(NA)))
	}
}

