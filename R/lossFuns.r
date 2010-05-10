################## Loss functions #########################
#                                                         #
# These functions compare a fitted contingency table with #
# an observed contingency table, and return a loss value. #
# Optimization functions minimize this value so a bad fit #
# should return a higher value.                           #
#                                                         #
###########################################################

lossMH <- function(Fit,Obs) {
    Fit <- Fit/sum(Fit);
    Obs <- Obs/sum(Obs);
    return(-1*(prod(Fit^Obs)));
}

lossSS <- function(Fit,Obs){
    Fit <- Fit/sum(Fit);
    Obs <- Obs/sum(Obs);
    return(sum((Fit-Obs)^2));
}

lossLL <- function(Fit,Obs){

    Fit <- round(Fit,15);
    return(-2*dmultinom(x=Obs,prob=Fit,log=T));
}

lossChiSq <- function(Fit,Obs){

    #ignore cells in which both Fit and Obs are zero
    O <- Obs[(Fit!=0) && (Obs!=0)];
    E <- Fit[(Fit!=0) && (Obs!=0)];
    
    chiSq <- sum(((O-E)^2)/((O+E)/2));
    return(chiSq);
}

lossPnorm <- function(Fit,Obs,p){
    #note: an additional parameter p is needed. If p=2, this equals lossSS.
    
    Fit <- Fit/sum(Fit);
    Obs <- Obs/sum(Obs);
    return(sum(abs(Fit-Obs)^p));
}

