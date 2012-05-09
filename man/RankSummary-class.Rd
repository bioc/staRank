\docType{class}
\name{RankSummary-class}
\alias{RankSummary}
\alias{RankSummary-class}
\title{Class RankSummary}
\usage{
  RankSummary(baseRank = 0, stabRank = 0, avrgRank = 0,
    stableSetSize = 0,
    rankCor = matrix(NA, nrow = 3, ncol = 3), method = "",
    Pi = matrix(NA, nrow = 1, ncol = 1))
}
\arguments{
  \item{baseRank}{a numeric vector containing the base
  ranking}

  \item{stabRank}{a numeric vector containing the stability
  ranking}

  \item{avrgRank}{a numeric vector containing the average
  ranking}

  \item{stableSetSize}{the sizes of the stable sets for
  each cutoff}

  \item{rankCor}{Spearman's rank correlation coefficient
  for the three rankings}

  \item{method}{the used method}

  \item{Pi}{a matrix containing}
}
\description{
  A class containing different rankings of a dataset
  further summary information like the stable set sizes, a
  correlation matrix of the rankings, the name of the base
  ranking method, the function call and if needed the Pi
  matrix of the stabilityRanking.

  \bold{Slots} \describe{ \item{\code{baseRank}:}{a numeric
  vector containing the base ranking}
  \item{\code{stabRank}:}{a numeric vector containing the
  stability ranking} \item{\code{avrgRank}:}{a numeric
  vector containing the average ranking}
  \item{\code{stableSetSize}:}{the sizes of the stable sets
  for each cutoff} \item{\code{rankCor}}{Spearman's rank
  correlation coefficient for the three rankings}
  \item{\code{method}}{the used method} \item{\code{Pi}:}{a
  matrix containing} }

  \bold{Accessors} \itemize{ \item \code{\link{baseRank}}
  \item \code{\link{stabRank}} \item \code{\link{avrgRank}}
  \item \code{\link{stableSetSize}} \item
  \code{\link{rankCor}} \item \code{\link{method}} \item
  \code{\link{Pi}} }

  \bold{Methods} \itemize{ \item \code{\link{show}} \item
  \code{\link{summary}} }

  RankSummary objects can be created using the constructor
  \code{RankSummary}, however this is most likely not
  needed.
}
\examples{
# overview of RankSummary class
showClass("RankSummary")
# generate dataset
d<-replicate(4,sample(1:10,10,replace=FALSE))
rownames(d)<-letters[1:10]
# run default stability ranking
s<-stabilityRanking(d)
# using an accessor functions on the RankSummary object
stabRank(s)
rankCor(s)
# create a new empty object of class RankSummary using the constructor
RankSummary()
}
\keyword{classes}

