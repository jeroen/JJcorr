rhos <- function(cop,thet){

    copulaFunction <- get(paste(cop,"c",sep=""));
    
    integral <- adapt(2,c(0,0),c(1,1),functn=copulaFunction,thet=thet);
    
    rho <- 12*integral$value-3;
    return(rho);    
}