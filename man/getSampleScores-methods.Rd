\docType{methods}
\name{getSampleScores}
\alias{getSampleScores}
\alias{getSampleScores,list-method}
\alias{getSampleScores,matrix-method}
\title{Get Sample Scores}
\arguments{
  \item{data}{can be a matrix with one row per element or a
  list of vectors of different length, one for each
  element.}

  \item{method}{one of the ranking methods: 'mean'
  (default), 'median', 'mwtest' (two sample one sided
  mann-whitney test), 'ttest'(two sample one sided
  t-test).}

  \item{decreasing}{a boolean indicating the direction of
  the ranking.}

  \item{bootstrap}{a boolean indicating whether
  bootstrapping or subsampling is used.}
}
\value{
  a vector containing the summary values by the given
  method for the sample dataset.
}
\description{
  an S4 method to a sample dataset that is scored by a
  given method.
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

