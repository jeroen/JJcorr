#require:
library(JJcor);
data(table4);
data(table6);

#correlations and bootstraps
myOC <- PC(table4,"gauss","MH");
myEC <- empPC(table4);
mySR <- spearmanR(table4);

#Goodness of Fit
GF(myOC);

#Bootstrapping
bs1 <- boot(myOC, size=10);
bs2 <- boot(myEC);
bs3 <- boot(mySR);

#polychoric correlations for multiple Copula families and loss functions
myPC <- multiPC(table4,loss=c("MH","SS"),fitLoss="SS",subdomains=2);
myPC <- multiPC(table6,loss=c("MH","SS"),fitLoss="SS",subdomains=2);

#Correlation matrix for multiple variables in a dataframe
mydf <- data.frame(W=rbinom(100,1,.5), X=rbinom(100,2,.5), Y=rbinom(100,3,.5), Z=rbinom(100,4,.5));
myMC <- corrMatrix(mydf,cop="gauss",loss="MH");
myMC <- corrMatrix(mydf,method="spearman");
cor(mydf,method="spearman");

