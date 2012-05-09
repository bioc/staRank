\name{mwTest2samp}
\alias{mwTest2samp}
\title{mwTest2samp}
\usage{
  mwTest2samp(x, y,
    alternative = c("two.sided", "less", "greater"),
    correct = TRUE)
}
\arguments{
  \item{x}{a numeric vector.}

  \item{y}{a numeric vector.}

  \item{alternative}{one of c("two.sided", "less",
  "greater").}

  \item{correct}{boolean.}
}
\value{
  the test statistic and the p-value of the test.
}
\description{
  Modified version of wilcox.test (taken from
  src/library/stats/R/wilcox.test.R) which performes a
  two-sample Mann-Withney test only and therefore is faster
  than the original version.
}
\examples{
x<-rnorm(100)
y<-rnorm(100,mean=1)
mwTest2samp(x,y,alternative='two.sided')
}

