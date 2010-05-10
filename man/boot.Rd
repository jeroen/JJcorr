\name{boot}
\Rdversion{1.1}
\alias{boot}

\title{
Non parametric bootstrap for a fitted correlation.
}
\description{
	This function performs a non parametric bootstrap for a fitted correlation. The argument should be an object returned by PC, empPC or spearmanR.
}
\usage{
	boot(model, size)
}

\arguments{
	  \item{model}{
		a correlation object as returned by functions PC, empPC or spearmanR.
	}
	  \item{size}{
		the number of resamples
	}
}

\value{
 \item{model}{copy of the object argument model}
 \item{r}{vector with correlation (r) values for every sample}
 \item{theta}{only for PC models: vector of fitted theta values for the Copula distribution}
 \item{loss}{only for PC models: vector of loss values for the fitted Copula distribution}
}


\examples{

#load contingency table
data(table4);

#fit a copula correlation
myOC <- PC(table4,"gauss","MH");
bs1 <- boot(myOC, size=10);
print(bs1);

#fit an empirical correlation
myEC <- empPC(table4);
bs2 <- boot(myEC);

#fit a spearman correlation
mySR <- spearmanR(table4);
bs3 <- boot(mySR);
}
