\name{corrMatrix}
\Rdversion{1.1}
\alias{corrMatrix}
\title{
Calculate a correlation matrix for a dataframe.
}
\description{
	This function calculates a full correlation matrix for a given dataframe. Note that all numeric variables will be converted to ordered factors.
}
\usage{
corrMatrix(dataframe, cop = "gauss", loss = "MH", domain = NULL, 
          subdomains = 1, method = "copula")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{
	dataframe with variables.
}
  \item{cop}{
	Only if method="copula": character string specifying which copula family to use. 
}
  \item{loss}{
	Only if method="copula": character string specifying which loss function to use. 
}
  \item{domain}{
	Only if method="copula": vector of length two specifying the allowed domain for theta. Defaults to the maximum domain. 
}
  \item{subdomains}{
	Only if method="copula": Integer specifying the number subdomains to perform numerical integration over. Used to avoid local optima.
}
  \item{method}{
	Correlation method. Should be one of "copula", "empirical", or "spearman".
}
}

\value{
		k by k correlation matrix. 
}

\examples{

#Correlation matrix for multiple variables in a dataframe
mydf <- data.frame(W=rbinom(100,1,.5), X=rbinom(100,2,.5), Y=rbinom(100,3,.5), 
        Z=rbinom(100,4,.5));
corrMatrix(mydf,cop="gauss",loss="MH");
corrMatrix(mydf,method="spearman");
#should be equal to:
cor(mydf,method="spearman");
}