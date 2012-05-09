#R version: Bin Zhou, bzhou@gnf.org, May 3, 2007
#		
# Usage of RSA.R:
# ===============
# --l: lower_bound, defaults to 0
# --u: upper bound, defaults to 1
# --r: reverse hit picking, the higher the score the better
# if -r flag is off, the lower the score the better
# --i: input file name
# --o: output file name, STDOUT if not specified
#
# Examples 
# R CMD BATCH --vanilla --slave --args --l=0.2 --u=0.8 --i=input.csv --o=output.csv RSA.R
# R CMD BATCH --vanilla --slave --args --l=1.2 --u=2.0 --r --i=input.csv -o=output.csv RSA.R
#
# Input and Output Format:
# ========================
# Gene_ID,Well_ID,Score: columns from input spreadsheet
# LogP: OPI p-value in log10, i.e., -2 means 0.01;
# OPI_Hit: whether the well is a hit, 1 means yes, 0 means no;
# #hitWell: number of hit wells for the gene
# #totalWell: total number of wells for the gene
# if gene A has three wells w1, w2 and w3, and w1 and w2 are hits,
# #totalWell should be 3, #hitWell should be 2, w1 and w2 should have OPI_Hit 
# set as 1 and w3 should have OPI_Hit set as 0.
# OPI_Rank: ranking column to sort all wells for hit picking
# Cutoff_Rank: ranking column to sort all wells based on Score in the simple 
# activity-based method
#
# Note: a rank value of 999999 means the well is not a hit. We put a large rank 
# number here for the convenient of spreadsheet sorting.
#
# Examples A
# -------------------------
# 1221200,7_O20,0.0541,-6.810,1,3,3,1,33
# 1221200,18_A21,0.0626,-6.810,1,3,3,2,43
# 1221200,41_A21,0.0765,-6.810,1,3,3,3,72
#
# Gene ID 1221200 has three wells, 7_O20, 18_A21 and 41_A21. All show good scores.
# Therefore 3 out of 3 wells are hits (#totalWell=3, #hitWell=3, OPI_Hit=1 for 
# all three wells)
# LogP is -6.810. These three wells are ranked as the best three wells by OPI.
# However, they are ranked as the 33th, 43th and 73th well by the traditional 
# cutoff method.
#
# Examples B
# -------------------------
# 3620,21_I17,0.0537,-2.344,1,1,2,162,31
# 3620,44_I17,0.7335,-2.344,0,1,2,999999,4113
#
# Gene ID 3620 has two wells, 21_I17 is active, while 44_I17 is relative inactive.
# OPI decides that only 1 out of the 2 wells is a hit. Therefore one well has 
# OPI_Hit set as 1, and the other 0. #totalWell=2, but #hitWell=1.
# The first well is the 162th hit by OPI, 31th by cutoff method.
# The second well is not a hit by OPI, 4113th by cutoff method.

#' checkOptions
#' 
#' Checks the input options for RSA.
#' 
#' @param args the input options. This is a list of: 
#' LB: lower_bound (defaults = 0), UB: upper bound (defaults = 1)
#' reverse: boolean (if TRUE: reverse hit picking, higher scores are better).
#' 
#' @return a list of the input options.
#' 
#' @keywords internal
#' @noRd
checkOptions <- function(args)
{
	i = grep("--l",args); 
	LB <- ifelse( !any(i), 0, as.numeric(strsplit(args[i],"=")[[1]][2]))
	i = grep("--u",args); 
	UB = ifelse( !any(i), 1, as.numeric(strsplit(args[i],"=")[[1]][2]))
	i = grep("--o",args); 
	outputFile = ifelse( !any(i), NA, strsplit(args[i],"=")[[1]][2])
	i = grep("--i",args); 
	inputFile = ifelse( !any(i), NA, strsplit(args[i],"=")[[1]][2])
	i = grep("--r$",args); 
	reverse  = ifelse( !any(i), FALSE,TRUE)
	bb = list(LB=LB,UB=UB,outputFile=outputFile,inputFile=inputFile,reverse=reverse);
	return (bb);
}

#' handleOneGroup
#' 
#' Performs calculation of the OPI score on one group of siRNAs.
#' 
#' @param i one Gene group.
#' @param dataset the input dataset.
#' @param optsb the options.
#' 
#' @return a list of the input options.
#' 
#' @keywords internal
#' @noRd
handleOneGroup <- function(i,dataset, optsb) 
{
	if(optsb$reverse)
	{
  		i_max = sum(dataset$Score[i]>=optsb$LB)
  		i_min = max(1,sum(dataset$Score[i]>=optsb$UB))
	}else
	{
  		i_max = sum(dataset$Score[i]<=optsb$UB)
  		i_min = max(1,sum(dataset$Score[i]<=optsb$LB))
	}
	r = OPIScore(i,nrow(dataset),i_min,i_max)
	return ( cbind(
  				LogP = r["logp"]
				,OPI_Hit=as.numeric(seq(length(i))<=r["cutoff"])
  				,"#hitWell"=i_max
				,"#totalWell"=length(i)
				,rank = i))
}

#' OPIScore
#' 
#' Performes the OPI score calculation.
#' 
#' @param I_rank one group rank.
#' @param N row number of dataset.
#' @param i_min=1 the minimal value of group I.
#' @param i_max=-1args the maximal value of group I.
#' 
#' @return a list of OPI p-values and cutoffs.
#' 
#' @keywords internal
#' @noRd
OPIScore <- function(I_rank, N, i_min=1, i_max=-1)
{
	n_drawn = length(I_rank) # number of black
	if(i_max == -1)
	{
    	i_max=n_drawn
	}
	r1 = c(logp=1.0,cutoff=0)
	if( i_max < i_min) return (r1)
	# phyper(x, lower.tail = F), x = x-1, when lower.tail = F
	logp =  apply(cbind(seq(i_min,i_max),I_rank[i_min:i_max]),1,function(x) {
				phyper(x[1]-1,x[2] ,N-x[2], n_drawn,lower.tail = F,log.p=T)})
        
	logp = logp/log(10)
	logp[logp<(-100)] = -100
	if(all(is.na(logp))) {
		return  (r1)
	}else
		return   ( c(logp=min(logp),cutoff = i_min-1+which.min(logp)))
}

#' OPI
#' 
#' Creates the final table with OPI scores and hit lists.
#' 
#' @param Groups the input gene identifiers.
#' @param Scores the input scores.
#' @param opts the input options.
#' @param Data the input datamatrix.
#' 
#' @return a table with all ranking information.
#' 
#' @keywords internal
#' @noRd
OPI<-function(Groups,Scores,opts,Data=NULL)
{
	t = data.frame(cbind(Gene_ID = Groups, Score= Scores))
	t$Score<-as.numeric(Scores) #this line is removed in newest version -> error
	Sorted_Order = order(t$Score,decreasing=opts$reverse);
	Data = Data[Sorted_Order,]
	t = t[Sorted_Order,]
	t = do.call("rbind", tapply(seq(nrow(t)), list(t$Gene_ID), handleOneGroup, 
					dataset = t, opts))
	t = cbind(Data, t[order(t[,"rank"]),])

	# add OPI_Rank
	t = t[order(t$LogP,t$Score*ifelse(opts$reverse,-1,1)),]
	t$OPI_Rank = cumsum(t$OPI_Hit)
	t$OPI_Rank[t$OPI_Hit == 0] = 999999
	
	# add Cutoff_Rank
	t = t[order(t$Score*(ifelse(opts$reverse,-1,1)),t$LogP),]
	
	if(opts$reverse){tmp = t$Score>=opts$LB} else {tmp = t$Score<=opts$UB}
	t$Cutoff_Rank = cumsum(tmp)
	t$Cutoff_Rank[!tmp] = 999999
	
	# add EXP_Rank
	t$EXP_Rank = pmin(t$OPI_Rank,t$Cutoff_Rank)
        t$EXP_Rank = pmin(t$OPI_Rank,t$Cutoff_Rank)
        if(opts$reverse) {
                return(t[order(t$OPI_Rank, -t$Score),])
        } else {
                return(t[order(t$OPI_Rank, t$Score),])
        }
}	

#' runRSA
#' 
#' Performs RSA analysis from within R. It does exactly the same as the version
#' from BinZhou 2007 but data are parsed from within R.
#' 
#' @param t a data.frame in RSA format (can be created with the
#' function dataFormatRSA).
#' @param opts the options for the RSA ranking. This is a list of: 
#' LB: lower_bound (defaults = 0), UB: upper bound (defaults = 1)
#' reverse: boolean (if TRUE: reverse hit picking, higher scores are better,
#' default = FALSE).
#' 
#' @return a matrix containing the RSA analysis results. The gene wise 
#' LogP-value is considered for further ranking analysis.
#' 
#' @example inst/example/staRank-example.R
#' 
#' @export 
#' 
runRSA<-function(t,opts){
	options(echo=F);		#do not display input commands
#opts = checkOptions(commandArgs())
#t = read.csv(opts$inputFile);
	colNames = dimnames(t)[[2]]			# make sure the column names are right
	if(!( ("Gene_ID" %in% colNames) & ("Well_ID" %in% colNames) 
				&("Score" %in% colNames)))
	{
		warnings(" not colnmae")
		return(NA)
	}
#filter out bad record
	#t = subset(t, !is.na(Gene_ID) & Gene_ID != "" & !is.na(Score))
	t = subset(t, !is.na(t$Gene_ID) & t$Gene_ID != "" & !is.na(t$Score))
	r = OPI(t$Gene_ID,t$Score,opts,t) #the final results

#output result
#	if(is.na(opts$outputFile))
#	{
#		#print(r)
#	}else
#	{
#		write.csv(r,opts$outputFile,row.names = FALSE,quote = FALSE)
#	}
#  print(paste("Total #Genes = ", length(unique(r$Gene_ID)), collapse=""));
#  print(paste("Total #Wells = ", nrow(r), collapse=""));
#  print(paste("Total #Hit Genes = ", length(
#		unique(r$Gene_ID[r$OPI_Rank<999999])), collapse=""));
#  print(paste("Total #Hit Wells = ", sum(r$OPI_Rank<999999), collapse=""));
	return(r)
}

