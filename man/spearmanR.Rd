\name{spearmanR}
\Rdversion{1.1}
\alias{spearmanR}
\title{
calculate the Spearman correlation for a 2 dimensional contingency table.
}
\description{
	calculate the Spearman correlation for a 2 dimensional contingency table.
}
\usage{
spearmanR(P)
}

\arguments{
  \item{P}{
	a 2 dimensional contingency table.
}
}

\value{
  the Spearman correlation value (r) for the observed data.
}

\examples{
data(table4);
#fit an empirical correlation
mySR <- spearmanR(table4);
(boot(mySR));

#example with a dataframe:
mydf <- data.frame(W=rbinom(100,1,.5), X=rbinom(100,2,.5), Y=rbinom(100,3,.5), 
        Z=rbinom(100,4,.5));
myMC <- corrMatrix(mydf,method="spearman");

#results should be exactly identical to:
cor(mydf,method="spearman")	
}

