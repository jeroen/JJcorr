########################## H_Vol #############################
#                                                            #
# H_Vol calculates the volume between the given coordinates  #
# for x and y, given a certain copula family and theta       #
# parameter                                                  #
#                                                            #
##############################################################


H_Vol <- function(cop,xlim,ylim,thet){

    if ((length(xlim)!=2) || (length(ylim)!=2)){
        stop('Both x and y should be vectors of length 2.');
    }   
    
    if ((xlim[1] >= xlim[2]) || (ylim[1] >=ylim[2])){
        stop('Upper limits should be higher then lower limits for x and y.');
    }

    copulaFunction <- get(paste(cop,"c",sep=""));
    out <- copulaFunction(c(xlim[2],ylim[2]),thet) - copulaFunction(c(xlim[1],ylim[2]),thet) - copulaFunction(c(xlim[2],ylim[1]),thet) + copulaFunction(c(xlim[1],ylim[1]),thet)

    return(out);
}

