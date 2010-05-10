multiPC <- function(P,cop=c("gauss","frank","clayton","nelsen2","genest"),loss=c("MH","SS"),fitLoss="modelLoss",subdomains=1){

    output <- list();
    output$P <- P;
    output$models <- list();
    
    for(thisCop in cop){
    
        output$models[[thisCop]] <- list();

        for(thisLoss in loss){
        
            thisModel <- PC(P,thisCop,thisLoss,subdomains=subdomains);
            thisGF <- GF(thisModel,fitLoss=fitLoss);
            output$models[[thisCop]][[thisLoss]] <- thisGF;
        }
    }
    
    output$P <- P;
    output$cop <- cop;
    output$loss <- loss;
    
    class(output) <- "multiPC";
    return(output);
}

print.multiPC <- function(x, ...){

    cop <- x$cop;
    loss <- x$loss;

    cat(rep("---------- multiPC Report ------------\t\t",length(loss)),"\n",sep="");
    cat(paste("Loss Function: ",loss,"\t\t\t\t",sep=""),"\n",sep="");
    cat(rep("\t\t r \t\t p(fit)\t\t",length(loss)),"\n",sep="");
    
    for(thisCop in cop){
        outputString <- "";
        for(thisLoss in loss){
            outputString <- paste(outputString, substr(paste(thisCop,"        ",sep=""),0,10),"\t",sep="");
            outputString <- paste(outputString, round(attr(x$models[[thisCop]][[thisLoss]]$model,"r"),2),"\t\t");
            outputString <- paste(outputString, round(x$models[[thisCop]][[thisLoss]]$p,2),"\t\t");
        }
        cat(outputString,"\n");
    }

    cat(rep("---------------------------------------\t\t",length(loss)),"\n",sep="");
}



