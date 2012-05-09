\docType{methods}
\name{summaryStats}
\alias{summaryStats}
\alias{summaryStats,list-method}
\alias{summaryStats,matrix-method}
\title{Summary statistic per element}
\arguments{
  \item{data}{a matrix with one row per element or a list
  containing one vector per element.}

  \item{method}{the score that is calculated per gene, one
  of 'mean' (default), 'median', 'mwtest' (two-sample one
  sided Mann-Whitney test), 'ttest'(two-sample one sided
  t-test), 'RSA' (redundant siRNA activity).}

  \item{decreasing}{a boolean indicating the direction of
  the ranking.}
}
\value{
  a named vector of the scored elements.
}
\description{
  Calculates the summary statistic per element for the
  whole dataset.
}
\examples{
# generate dataset
d<-replicate(4,sample(1:10,10,replace=FALSE))
rownames(d)<-letters[1:10]

# rank aggregation on the dataset using two base methods
aggregRank(d, method='mean')
aggregRank(d, method='median')

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

