corrMatrix <- function(dataframe,cop="gauss",loss="MH",domain=NULL,subdomains=1,method="copula"){

    k = length(dataframe);
    output <- matrix(NA,k,k);
    diag(output) <- 1;
    
    for(ki in 1:k){
    
        if(!is.factor(dataframe[[ki]])){
            message(paste("Converted variable '",names(dataframe[ki]),"' to a factor!",sep=""));
            dataframe[[ki]] <- factor(dataframe[[ki]],ordered=T);
        }
    }
    
    for(thisVar1 in 1:(k-1)){
    
        var1 <- dataframe[[thisVar1]];

        for(thisVar2 in (thisVar1+1):k){

            var2 <- dataframe[[thisVar2]];

            P <- table(var1,var2);

            try({
                if(method=="copula"){
                    myOC <- PC(P,cop,loss,domain,subdomains);
                    output[thisVar2,thisVar1] <- output[thisVar1,thisVar2] <- attr(myOC,"r");
                } else if(method=="empirical"){
                    output[thisVar2,thisVar1] <- output[thisVar1,thisVar2] <- empPC(P);
                } else if(method=="spearman"){
                    output[thisVar2,thisVar1] <- output[thisVar1,thisVar2] <- spearmanR(P);
                }
            });
        }
    }
    colnames(output) <- names(dataframe);
    rownames(output) <- names(dataframe);
    return(output);
}