#' Ranking variables based on their stability
#'
#' Detecting all relevant variables from a data set is challenging, especially 
#' when only few samples are available and data is noisy. Stability ranking 
#' provides improved variable rankings of increased robustness using resampling 
#' or subsampling.
#' 
#' @name staRank-package
#' @aliases staRank
#' @docType package
#' @title Stability ranking
#' @author \email{juliane.siebourg@@bsse.ethz.ch}
#' 
#' @example inst/example/staRank-example.R
#' @import methods cellHTS2
#' @keywords package
#'
NULL

#' The package contains an RNAi screen on human cells under 
#' Salmonella infection. Briefly, in HeLa cells all genes were knocked-down 
#' individually by RNA interference. Subsequently cells were infected with 
#' Salmonella. The cells were  imaged with a microscope and from the images 
#' several features were extracted 
#' (for details see original paper (Misselwitz2011)).
#' This dataset contains normalized infection ratios. These are computed as
#' the logarithm of the fraction of infected cells per knock-down. Subsequently
#' they are z-scored per plate.
#' Duplicated genes were removed and also genes containing NA values.
#' 
#' @name salmonellaInfection
#' @docType data
#' @title RNAi screen of human cells during Salmonella infection
#' 
#' @usage salmonellaInfection
#' @format a numeric matrix of 6860 genes (rows) with 3 siRNA values each 
#' (columns).
#' @references Misselwitz B. et al. (2011) \emph{RNAi screen of Salmonella 
#' invasion shows role of COPI in membrane targeting of cholesterol and Cdc42.}
#' Molecular Systems Biology.
#' 
#' @keywords datasets
#' 
NULL

#' The package contains the results of an example stability analysis for the
#' Salmonella infection dataset.
#' RankSummary objects were created for five base ranking methods on a subset 
#' of the genes.
#' 
#' @name salmonellaStability
#' @docType data
#' @title Example stability analysis on a subset of the Salmonella infection 
#' dataset
#' 
#' @usage salmonellaStability
#' @format a list of RankSummary objects that use different base ranking methods
#' 
#' @keywords datasets
#' 
NULL