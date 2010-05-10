\name{multiPC}
\Rdversion{1.1}
\alias{multiPC}
\title{
Calculate generalized correlation coefficients for multiple families and loss functions.
}
\description{
	A wrapper for PC and GF functions, that calculates correlation coefficients and goodness of fit for multiple copula families and loss functions at once. Usefull to explore which copula family might be appropriate.
}
\usage{
multiPC(P, cop = c("gauss", "frank", "clayton", "nelsen2", "genest"), 
        loss = c("MH", "SS"), fitLoss = "modelLoss", subdomains = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
	A 2 dimensional contingency table.
}
  \item{cop}{
	A vector of strings, specifying which families to try. Defaults to all implemented copula families.
}
  \item{loss}{
	A vector of strings, specifying which loss functions to try. Defaults to both "MH" (likelihood) and "SS" (sum of squares).
}
  \item{fitLoss}{
	A character string, specifying which loss function to use for goodness of fit testing. Defaults to "modelLoss", which means the same loss function as was used to fit the model.
}
  \item{subdomains}{
	The number of subdomains to perform numerical integration over. 
}
}

\value{
	\item{P}{a copy of the argument 'P' (the contingency table)}
	\item{models}{A list with all the fitted models.}
	\item{cop}{A copy of the argument 'cop'}
	\item{loss}{A copy of the argument 'loss'}
}

\examples{
data(table4);
myPC <- multiPC(table4,loss=c("MH","SS"),fitLoss="SS",subdomains=2);
}

