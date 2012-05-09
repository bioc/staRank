#' mwTest2samp
#' 
#' Modified version of wilcox.test (taken from src/library/stats/R/wilcox.test.R)
#' which performes a two-sample Mann-Withney test only and therefore is faster
#' than the original version.
#' 
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param alternative one of c("two.sided", "less", "greater").
#' @param correct boolean.
#' 
#' @return the test statistic and the p-value of the test.
#' 
#' @examples 
#' x<-rnorm(100)
#' y<-rnorm(100,mean=1)
#' mwTest2samp(x,y,alternative='two.sided')
#' 
#' @export 
#' 
mwTest2samp <-
function(x, y, alternative = c("two.sided", "less", "greater"),
		correct = TRUE)
{
	x <- x[is.finite(x)]
	y <- y[is.finite(y)]
	CORRECTION <- 0
	r <- rank(c(x, y))
	n.x <- as.double(length(x))
	n.y <- as.double(length(y))
	exact <- (n.x < 50) && (n.y < 50)
	STATISTIC <- sum(r[seq_along(x)]) - n.x * (n.x + 1) / 2
	names(STATISTIC) <- "W"
	ur<-unique(r)
	TIES <- (length(r) != length(ur))
	if(exact && !TIES) {
		PVAL <- switch(alternative,
				"two.sided" = {p <- if(STATISTIC > (n.x * n.y / 2))
								pwilcox(STATISTIC - 1, n.x, n.y, lower.tail = FALSE)
							else
								pwilcox(STATISTIC, n.x, n.y)
					min(2 * p, 1)},
				"greater" = {pwilcox(STATISTIC - 1, n.x, n.y, lower.tail = FALSE)},
				"less" = pwilcox(STATISTIC, n.x, n.y))			
	}
	else {
		NTIES <- tabulate(r*10)		
		z <- STATISTIC - n.x * n.y / 2
		SIGMA <- sqrt((n.x * n.y / 12) * ((n.x + n.y + 1)
							- sum(NTIES^3 - NTIES) / ((n.x + n.y) * (n.x + n.y - 1))))
		if(correct) {
			CORRECTION <- switch(alternative,
					"two.sided" = sign(z) * 0.5,
					"greater" = 0.5,
					"less" = -0.5)
		}
		z <- (z - CORRECTION) / SIGMA
		PVAL <- switch(alternative,
				"less" = pnorm(z),
				"greater" = pnorm(z, lower.tail=FALSE),
				"two.sided" = 2 * min(pnorm(z),pnorm(z, lower.tail=FALSE)))
	}	
	RVAL <- list(statistic = STATISTIC,p.value = as.numeric(PVAL))
	return(RVAL)
}

