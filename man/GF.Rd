\name{GF}
\Rdversion{1.1}
\alias{GF}
\title{
Perform a Goodness of Fit test for a fitted copula model.
}
\description{
Resamples from the fitted model (parametric bootstrap) and calculates the loss with respect to this model. Compares the loss of the resamples with the observed data. 
}
\usage{
GF(model, size = 10000, fitLoss = "modelLoss")
}

\arguments{
	  \item{model}{
		a correlation object as returned by the function PC.
	}
	  \item{size}{
		the number of resamples from the fitted model.
	}
	  \item{fitLoss}{
		the loss function used to calculate the loss of the observed and resampled data with respect to the model. Should be one of "modelLoss", "MH" or "SS". Defaults to "modelLoss", which means the same loss function is used as was used to fit the model.
	}
}

\value{
 \item{model}{
	A copy of the model argument
 }
 \item{fitLossValue}{
	The loss value for the observed data with respect to the model
 }
 \item{fitLossDist}{
	A vector with fitLoss values for the resampled data with respect to the model
 }
 \item{p}{
	The p value for the model fit; i.e. the proportion of resamples with a higher loss than the observed data
 }
}

\examples{
#load contingency table
 data(table4);

#fit a copula correlation
 myOC <- PC(table4,"gauss","MH");
 
#perform goodness of fit:
 GF(myOC,fitLoss="SS");
}