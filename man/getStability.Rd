\name{getStability}
\alias{getStability}
\title{getStability}
\usage{
  getStability(sr, thr)
}
\arguments{
  \item{sr}{sample rankings in the form of a matrix where
  each row corresponds to one element and each column gives
  one ranking.}

  \item{thr}{the threshold for the stability selection,
  indicating above which frequency in the samples an
  element is considered stable.}
}
\value{
  a list containing: \describe{ \item{\code{stabRank}}{the
  stable ranking.} \item{\code{Pi}}{the frequency matrix
  with all values per gene and per cutoff.}
  \item{\code{stableSetSize}}{a table with the number of
  stable genes per cutoff.} }
}
\description{
  Performes stability selection on sample rankings with a
  given stability threshold. Selection probabilities and
  stability ranking are calculated.
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

