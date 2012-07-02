\name{uniqueRSARanking}
\alias{uniqueRSARanking}
\title{uniqueRSARanking}
\usage{
  uniqueRSARanking(dataRSA, opts)
}
\arguments{
  \item{dataRSA}{a matrix with data in RSA format (can be
  created with the function dataFormatRSA).}

  \item{opts}{the options for the RSA ranking. This is a
  list of: LB: lower_bound (defaults = 0), UB: upper bound
  (defaults = 1) reverse: boolean (if TRUE: reverse hit
  picking, higher scores are better, default = FALSE).}
}
\value{
  a vector of ranked genes with their RSA LogP-values.
}
\description{
  Performs RSA analysis and summarizes the results gene
  wise.
}
\examples{
# generate dataset
d<-replicate(4,sample(1:10,10,replace=FALSE))
rownames(d)<-letters[1:10]

# rank aggregation on the dataset using two base methods
aggregRank(d, method='mean')
aggregRank(d, method='median')

# calculate summary statistic from the data
summaryStats(d, method='mean')
summaryStats(d, method='RSA')

# calculating replicate scores from different summary statistics
scores<-getSampleScores(d,'mean',decreasing=FALSE,bootstrap=TRUE)
scores<-getSampleScores(d,'mwtest',decreasing=FALSE,bootstrap=TRUE)

# perform RSA analysis

# get RSA format of data
rsaData<-dataFormatRSA(d)
# set RSA options
opts<-list(LB=min(d),UB=max(d),reverse=FALSE)
# run the RSA analysis
r<-runRSA(rsaData,opts)
# directly obtain the per gene RSA ranking from the data
r<-uniqueRSARanking(rsaData,opts)

# get stable Ranking, stable setsizes and the Pi matrix for default settings
# and stability threshold of 0.9
s<-getStability(d,0.9)

# run default stability ranking
s<-stabilityRanking(d)

# using an accessor function on the RankSummary object
stabRank(s)

# summarize a RankSummary object
summary(s)

# generate a rank matrix from a RankSummary object
getRankmatrix(s)
}

