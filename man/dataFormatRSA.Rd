\name{dataFormatRSA}
\alias{dataFormatRSA}
\title{dataFormatRSA}
\usage{
  dataFormatRSA(dataMatrix)
}
\arguments{
  \item{dataMatrix}{a matrix with one row per gene and
  columns are replicate or siRNA values.}
}
\value{
  a data.frame with three columns Gene_ID, Well_ID, Score.
  The Well_ID does not have a meaning if data are
  replicates, but is usefull to distinguish between siRNAs.
}
\description{
  Creates a data.frame in the RSA input format from a
  matrix, where rows indicate the genes and columns the
  replicate or siRNA values.
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

