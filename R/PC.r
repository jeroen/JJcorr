PC <- function(P,cop,loss,domain=NULL,subdomains=1){

    if(is.null(domain)){
        copulaFunction <- get(paste(cop,"c",sep=""));
        domain <- copulaFunction()$domain;
    }
    
    ## add 1 because 2 points make only 1 interval
    subdomains <- subdomains+1;

    intervals <- seq(domain[1],domain[2],length=subdomains);

    lossValue <- Inf;
    minimum <- NULL;
    
    for(i in 2:subdomains){

        tempMinimum <- optimize(f=copLoss,interval=c(intervals[i-1],intervals[i]),cop=cop,P=P,loss=loss);
        if(tempMinimum$objective < lossValue){
        
            lossValue <- tempMinimum$objective;
            minimum <- tempMinimum;
        }
    }

    b <- minimum$minimum;

    if(cop=="gauss"){
        out = b;
    } else {
        out = 2*sin(rhos(cop,b)*pi/6);
    }

    if (((b-domain[1]) < ((domain[2]-domain[1])/1000)) || ((domain[2]-b) < ((domain[2]-domain[1])/1000))) {
        warning("The minimum is very close the the boundary of the parameter domain.\n\nIt may be a numerical optimization problem. Try altering the domain.");
    }
    
    attr(out,"theta") <- minimum$minimum;
    result <- minimum$objective;
    
    attr(out,"fitted") <- attr(result,"fitted");
    attr(result,"fitted") <- NULL;
    attr(out,"lossValue") <- result;
    if(cop=="gauss"){
        attr(out,"r") <- attr(out,"theta");
        #for gaussian empPCs, r equals theta
    } else {
        attr(out,"r") <- 2*sin(rhos(cop,minimum$minimum)*pi/6);
    }
    
    attr(out,"P") <- P;
    attr(out,"cop") <- cop;
    attr(out,"loss") <- loss;
    attr(out,"subdomains") <- subdomains;
    attr(out,"domain") <- domain;
    
    class(out) <- "PC";
    return(out);
}

print.PC <- function(x, ...){

    cat("\n--------- PC Report ------------\n");
    cat("empPC Family:", attr(x,"cop"),"\n");
    cat("Loss Function Used:", attr(x,"loss"),"\n");
    cat("Optimal Theta", attr(x,"theta"),"\n");
    cat("Fitted Contingency Table:\n\n");
    print(attr(x,"fitted"));
    cat("\nLoss", attr(x,"lossValue"),"\n");
    cat("Correlation (r)", attr(x,"r"),"\n");
    cat("---------------------------------------\n");
    

}