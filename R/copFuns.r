################## Copula Family functions ################
#                                                         #
# Functions compute the CDF at some point x and y         #
# for a certain family given parameter value theta.       #
#                                                         #
###########################################################

claytonc <- function(inputVector=NULL,thet=NULL){

    u <- inputVector[1];
    v <- inputVector[2];
    domain = c(-1,40);

    if(is.null(thet)) return(list(domain=domain));

    if((thet < min(domain)) || (thet > max(domain))){
        warning("Theta parameter seems to be out of range");
    }

    eps <- 0.0000001;
    if (abs(thet) < eps){
        out <- u*v;
    } else {
        out <- max(c(u^(-thet) + v^(-thet)-1, matrix(c(0,dim(u)))))^(-1/thet);
    }

    attr(out,"domain") <- domain;
    return(out);
}

frankc <- function(inputVector=NULL,thet=NULL){

    u <- inputVector[1];
    v <- inputVector[2];
    domain <- c(-100,37);
    
    if(is.null(thet)) return(list(domain=domain));

    if((thet < min(domain)) || (thet > max(domain))){
        warning("Theta parameter seems to be out of range");
    }

    eps = 0.0000001;

    if (abs(thet) < eps){
        out=u*v;
    } else {
        out = as.vector(-(thet^(-1))%*%log(1+(exp(-thet*u)-1)*(exp(-thet*v)-1)/(exp(-thet)-1)));
    }

    attr(out,"domain") <- domain;
    return(out);
}

gaussc <- function(inputVector=NULL,thet=NULL){

    u <- inputVector[1];
    v <- inputVector[2];
    domain = c(-1,1);
    
    if(is.null(thet)) return(list(domain=domain));

    if((thet < min(domain)) || (thet > max(domain))){
        warning("Theta parameter seems to be out of range");
    }

    eps = 0.001;

    if(thet <= -1+eps){
        out = max(u+v-1,0);
        #when correlation approaches -1, the formula simplifies
    } else if (thet >= 1-eps){
        out = min(u,v);
        #when correlation approaches  1, the formula simplifies
    } else if ((u < eps) || (v < eps)){
        out = 0;
        # if values for u or v are close to the boundary of I^2, there is no volume
    } else if (u > 1-eps){
        out = v;
    } else if (v > 1-eps){
        out = u;
    } else {
       qu <- qnorm(u);
       qv <- qnorm(v);
       out = pmvnorm(upper=c(qu,qv),corr=cbind(c(1,thet),c(thet,1)))
    }

    attr(out,"domain") <- domain;
    return(out);
}

genestc <- function(inputVector=NULL,thet=NULL){

    u <- inputVector[1];
    v <- inputVector[2];

    domain = c(1,10);
    
    if(is.null(thet)) return(list(domain=domain));

    if((thet < min(domain)) || (thet > max(domain))){
        warning("Theta parameter seems to be out of range");
    }

    out = (max(1-((1-u^(1/thet))^thet+(1-v^(1/thet))^thet)^(1/thet),0))^thet;

    attr(out,"domain") <- domain;
    return(out);
}

nelsen2c <- function(inputVector=NULL,thet=NULL){

    u <- inputVector[1];
    v <- inputVector[2];

    domain = c(1,10);
    
    if(is.null(thet)) return(list(domain=domain));

    if((thet < min(domain)) || (thet > max(domain))){
        warning("Theta parameter seems to be out of range");
    }

    out = max(1-((1-u)^thet+(1-v)^thet)^(1/thet),matrix(c(0,dim(u))));

    attr(out,"domain") <- domain;
    return(out);
}


