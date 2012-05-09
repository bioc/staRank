\docType{data}
\name{salmonellaInfection}
\alias{salmonellaInfection}
\title{RNAi screen of human cells during Salmonella infection}
\format{a numeric matrix of 6860 genes (rows) with 3 siRNA values each
(columns).}
\description{
  The package contains an RNAi screen on human cells under
  Salmonella infection. Briefly, in HeLa cells all genes
  were knocked-down individually by RNA interference.
  Subsequently cells were infected with Salmonella. The
  cells were imaged with a microscope and from the images
  several features were extracted (for details see original
  paper (Misselwitz2011)). This dataset contains normalized
  infection ratios. These are computed as the logarithm of
  the fraction of infected cells per knock-down.
  Subsequently they are z-scored per plate. Duplicated
  genes were removed and also genes containing NA values.
}
\references{
  Misselwitz B. et al. (2011) \emph{RNAi screen of
  Salmonella invasion shows role of COPI in membrane
  targeting of cholesterol and Cdc42.} Molecular Systems
  Biology.
}
\keyword{datasets}

