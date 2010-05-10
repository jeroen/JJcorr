\name{PC}
\Rdversion{1.1}
\alias{PC}
\title{
Fit a copula distribution for a contingency table.
}
\description{
PC fits the specified bivariate copula distribution for a 2 dimensional contingency table, and returns the fitted correlation (r) value. 
}
\usage{
PC(P, cop, loss, domain = NULL, subdomains = 1)
}

\arguments{
  \item{P}{
	a 2 dimensional contingency table. 
}
  \item{cop}{
	A string specifying the copula family to fit. Should be one of "gauss", "clayton", "frank", "genest" or "nelsen2".
}
  \item{loss}{
	A string specifying which loss function is to be used for fitting the distribution. Currently "MH" (likelihood) and "SS" (sum of squares) are implemented.
}
  \item{domain}{
	optional: vector of lenght 2, specifying the domain for theta over which to optimize.
}
  \item{subdomains}{
	optional: an integer specifying the number of intervals to split the domain in for the optimization.
}
}
\details{
	PC uses numerical optimization to find an optimal value for the association parameter theta for the given distribution. The loss is calculated by comparing the observed contingency table with the fitted one as derived from the fitted distribution. 
}
\value{
	The value for the correlation (r). Has attributes with additional information. See examples.
}

\note{
	If you suspect numerical optimization to find local minima, try a smaller domain for theta, or specify a higher number of subdomains to perform numerical optimization over.
}

\examples{
#load contingency table
data(table4);

#fit a copula correlation
myOC <- PC(table4,"gauss","MH");
print(myOC);
attributes(myOC);
(bs1 <- boot(myOC));
}