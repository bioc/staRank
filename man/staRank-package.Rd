\docType{package}
\name{staRank-package}
\alias{staRank}
\alias{staRank-package}
\title{Stability ranking}
\description{
  Ranking variables based on their stability
}
\details{
  Detecting all relevant variables from a data set is
  challenging, especially when only few samples are
  available and data is noisy. Stability ranking provides
  improved variable rankings of increased robustness using
  resampling or subsampling.
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
\author{
  \email{juliane.siebourg@bsse.ethz.ch}
}
\keyword{package}

