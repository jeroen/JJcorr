spearmanR <- function(P){

    pNormalized <- P/sum(P);

    marg_x <- apply(pNormalized,1,sum);
    marg_y <- apply(pNormalized,2,sum);

    cum_x <- c(0,cumsum(marg_x));
    cum_y <- c(0,cumsum(marg_y));

    dim_x <- length(cum_x);
    dim_y <- length(cum_y);

    cum_x[cum_x > 1] <- 1;
    cum_y[cum_y > 1] <- 1;

    T <- matrix(NA,dim_x,dim_y);
    T[1,] <- 0;
    T[,1] <- 0;

    for(i in 1:(dim_x-1)){
        for(j in 1:(dim_y-1)){
            T[i+1,j+1] <- sum(pNormalized[1:i,1:j])
        }
    }

    H <- matrix(NA,dim_x-1,dim_y-1);

    for(i in 1:(dim_x-1)){
        for(j in 1:(dim_y-1)){
            a <- c(cum_x[i], cum_y[j], T[i,j]);
            b <- c(cum_x[i], cum_y[j+1], T[i,j+1]);
            c <- c(cum_x[i+1], cum_y[j], T[i+1,j]);
            d <- c(cum_x[i+1], cum_y[j+1], T[i+1,j+1]);

            H[i,j] <- (c[1]-a[1])*(b[2]-a[2])*(a[3]+b[3]+c[3]+d[3])/4
        }
    }

    rho_S = 12*sum(H)-3;
    out = rho_S * ((1-sum(marg_x^3))*(1-sum(marg_y^3)))^(-.5);

    attr(out,"method") <- "spearman";
    attr(out,"P") <- P;
    class(out) <- "spearmanR";
    return(out);
}
