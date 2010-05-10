\name{empPC}
\Rdversion{1.1}
\alias{empPC}
\title{
calculate the empirical correlation for a 2 dimensional contingency table.
}
\description{
	calculate the empirical correlation for a 2 dimensional contingency table.
}
\usage{
empPC(P)
}

\arguments{
  \item{P}{
	a 2 dimensional contingency table.
}
}

\value{
  the empirical correlation value (r) for the observed data.
}

\examples{
data(table4);

#fit an empirical correlation
myEC <- empPC(table4);
bs2 <- boot(myEC);
}

