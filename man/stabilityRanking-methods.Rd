\docType{methods}
\name{stabilityRanking}
\alias{stabilityRanking}
\alias{stabilityRanking,cellHTS-method}
\alias{stabilityRanking,matrix-method}
\alias{stabilityRanking,numeric-method}
\title{Wrapper to perform stabilityRanking}
\arguments{
  \item{x}{the data that is to be ranked. It can be of
  different types: a data matrix with one row per element
  to be ranked, or a ranking from an external method, or an
  object of class cellHTS. Depending on this the other
  parameters can vary.}

  \item{samps}{a matrix of scored sample data, each row
  corresponds to an element, the columns to a scoring. This
  is only needed when an external ranking method is used.}

  \item{channel}{a string with the name of the feature
  (channel) to be ranked. This is only needed for cellHTS
  objects.}

  \item{replicates}{names or indices of the replicates
  (samples) to be used for the rankings (default: all
  samples are used).}

  \item{method}{one of the ranking methods: 'mean'
  (default), 'median', 'mwtest' (two-sample one sided
  Mann-Whitney test), 'ttest'(two-sample one sided t-test)
  or 'RSA' (redundant siRNA analysis). If an external
  ranking is used, you can specify the name of that ranking
  method in the method argument.}

  \item{decreasing}{a boolean indicating the direction of
  the ranking.}

  \item{bootstrap}{a boolean indicating if bootstrapping or
  subsampling is used.}

  \item{thr}{threshold for stability (default = 0.9).}

  \item{nSamp}{the number of samples to generate (default =
  100).}

  \item{Pi}{boolean indicating if the Pi matrix should be
  returned (can be very large, default=FALSE).}

  \item{verbose}{boolean indicating if status update should
  be printed}

  \item{...}{further parameter for the stability ranking.}
}
\value{
  an object of class \code{\link{RankSummary}}.
}
\description{
  An S4 function to perform stability ranking on a dataset
  or directly on given sample rankings.
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

