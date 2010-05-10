########################## copLoss ###########################
#                                                            #
# This function calculates volumes for a given copula family #
# and corresponding theta parameter using H_Vol. Then it     #
# calculates the loss of this contingency table with respect #
# to the observed contingency table P.                       #
#                                                            #
##############################################################

copLoss <- function(thet,cop,P,loss){

    copulaFunction <- get(paste(cop,"c",sep=""));
    lossFunction <- get(paste("loss",loss,sep=""));

    pNormalized <- P/sum(P);

    dim_x <- dim(pNormalized)[1];
    dim_y <- dim(pNormalized)[2];
    
    marg_x <- apply(pNormalized,1,sum);
    marg_y <- apply(pNormalized,2,sum);
    
    cum_x <- c(0,cumsum(marg_x));
    cum_y <- c(0,cumsum(marg_y));
    
    cum_x[cum_x > 1] <- 1;
    cum_y[cum_y > 1] <- 1;

    V <- matrix(NA,dim_x,dim_y);
    
    for(i in 1:dim_x){
        for(j in 1:dim_y){
            V[i,j] = H_Vol(cop,c(cum_x[i],cum_x[i+1]),c(cum_y[j],cum_y[j+1]),thet);
        }
    }
    
    out <- lossFunction(V,P);

    ##Verbose Debug info:##
    #cat("Optimizing... theta:",thet,", loss:",out,"\n");
    ###############
    
    attr(out,"fitted") <- V;
    return(out);
}

