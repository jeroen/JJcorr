boot.empPC <- boot.spearmanR <- function(model,size=1000){

    if((class(model)!="empPC") && (class(model)!="spearmanR")){
        stop("model needs to be empirical or spearman correlation");
    }

    P <- attr(model,"P");
    sampleSize <- sum(P);
    Pnormalized <- P/sampleSize;
    
    rVector <- rep(NA,size);
    
    for(i in 1:size){

        thisDraw <- rmultinom(1,sampleSize,Pnormalized);
        dim(thisDraw) <- dim(P);
        
        if(class(model)=="spearmanR"){
            rVector[i] <- spearmanR(thisDraw);
        } else if(class(model)=="empPC") {
            rVector[i] <- empPC(thisDraw);
        }
    }
    output <- list(model=model,r=rVector);
    class(output) <- "empBoot";
    return(output);
}

print.empBoot <- function(x, ...){

    model <- x$model;
    rquantiles <- round(quantile(x$r,c(0.01,0.025,0.05,0.10,0.50,0.90,0.97,0.975,0.99)),3);

    cat("\n--------- empBoot Report ------------\n");
    cat("Correlation method:", attr(model,"method"),"\n");
    cat("Estimated Correlation (r):", model,"\n");
    cat("Standard Deviation for r:", sd(x$r),"\n");
    cat("\n------ Bootstrapped empPC ------\n");
    cat("Quantile\t r\n");
    cat("0.01\t\t",rquantiles[1],"\n");
    cat("0.025\t\t",rquantiles[2],"\n");
    cat("0.05\t\t",rquantiles[3],"\n");
    cat("0.10\t\t",rquantiles[4],"\n");
    cat("0.50\t\t",rquantiles[5],"\n");
    cat("0.90\t\t",rquantiles[6],"\n");
    cat("0.95\t\t",rquantiles[7],"\n");
    cat("0.975\t\t",rquantiles[8],"\n");
    cat("0.99\t\t",rquantiles[9],"\n");
    cat("---------------------------------------\n");
}