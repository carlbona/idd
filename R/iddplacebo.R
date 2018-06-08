#' @title Run placebo studies on untreated units and estimate placebo p-values.
#' @description \code{iddplacebo} estimates placebo studies on untreated units and produces pseudo p-values based on the empirical distribution of post/pre-RMSE ratios.
#'
#' @param eventvar Name of the event count variable.
#' @param popvar The person-time (or another exposure) variable.
#' @param treatvar The treatment group indicator dummy (0 for controls, 1 for treated units).
#' @param postvar The post-period indicator dummy (0 for all time points in the pre-intervention period, 1 for all time points in the post-period)
#' @param timevar The time variable (can be coded on any scale).
#' @param idvar The panel id variable.
#' @param data A long-form panel dataset containing the supplied variables.
#' @param iter The number of subsampling iterations. Defaults to 10, but 500-1000 are usually recommended.
#'
#' @return Returns a list containing the following elements:
#' \code{$Resdata}: a data frame containing the results.
#' \code{$ECDF}: the empirical cumulative distribution function.
#' \code{$Effects}: a long-form data frame containing the estimated placebo effects.
#' \code{$Treat.ratio}: a data frame containing the post/pre-RMSE ratio for the treated unit.
#' \code{$supp_stats}: a data frame containing supplementary statistics (p-values etc).
#' @examples
#' \dontrun{
#' data(simpanel)
#' placebo.out <- iddplacebo(eventvar="y",
#'              popvar="pop",
#'              idvar="age",
#'              timevar="time",
#'              postvar="post",
#'              treatvar="treat",
#'              data=simpanel,
#'              iter=50)
#' }
#' @export



iddplacebo <- function(eventvar, popvar, treatvar, postvar, timevar, idvar, data, iter=10) {
        realtr=NULL
        #Reslist, for storage
        resmat <- as.list(rep(NA),iter)
        effectmat <- as.list(rep(NA), iter)
        #Gen data
        E <- data[,eventvar]
        P <- data[,popvar]
        TR <- data[,treatvar]
        PO <- data[,postvar]
        TI <- data[,timevar]
        id <- data[,idvar]
        TI_c <- TI-(min(TI)-1)
        df.act <- as.data.frame(cbind(E, P, TR, PO, TI, TI_c, id))
        df <- subset(df.act, TR==0)

        T.orig <- length(unique(df$TI))
        cat("Running ", iter, " placebo studies, please wait...", sep="")
        pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
        for (i in 1:iter) {

                T.min <- 1
                T.max <- T.orig
                T1 <- min(subset(df.act, PO==1)$TI_c)

                #Sample units
                N <- sample(3:length(unique(df$id)), 1) #set N for subsample
                N.list <- sample(unique(df$id), N) #randomly pick N units

                #Gen dataset

                df.sub <- subset(df, TI_c>=T.min & TI_c<=T.max)

                #Set post period

                POST <- as.numeric(df.sub$TI_c>=T1)
                df.sub$POST <- POST

                #Set treatment units
                N.treat <- sample(1:(N-2), 1)
                Treat.list <- sample(N.list, N.treat)
                df.sub$TREAT <- as.numeric(df.sub$id %in% Treat.list)

                #Run idd

                mtch <- idd(eventvar="E", popvar="P", treatvar="TREAT", postvar="POST", timevar="TI", idvar="id", print=FALSE, data=df.sub)

                #Obtain POST/PRE-RMSPE

                predat <- subset(mtch$Resdat, post==0)
                postdat <- subset(mtch$Resdat, post==1)
                RMSE_T1 <- sqrt(mean(postdat$effect^2))
                RMSE_T0 <-sqrt(mean(predat$effect^2))
                RMSE_RATIO <- RMSE_T1/RMSE_T0

                #Obtain effect matrix (for spaghetti plot)

                effect <- mtch$Resdat$effect
                time <- mtch$Resdat$time
                post <- mtch$Resdat$post
                subsample <- rep(i, length(time)) #Ids for plotting

                #Obtain parametric p-value (for evaluation)

                param_p <- mtch$supp_stats$dd_pval

                #Get results
                T <- length(T.min:T.max)
                effectmat[[i]] <- cbind(effect, time, subsample, post, row.names=NULL)
                subsample <- i
                resmat[[i]] <- cbind(RMSE_RATIO, RMSE_T0, RMSE_T1, param_p, T.min, T.max, T, T1, N, N.treat, subsample, row.names=NULL)
                utils::setTxtProgressBar(pb, i)

        }
        resdat <- plyr::ldply(resmat, data.frame) #Flatten list to data frame
        effectmat <- plyr::ldply(effectmat, data.frame)

        #Run original
        orig <- idd(eventvar="E", popvar="P", treatvar="TR", postvar="PO", timevar="TI", idvar="id", print=FALSE, data=df.act)
        T.min <- min(df.act$TI_c)
        T.max <- max(df.act$TI_c)
        T <- length(T.min:T.max)
        T1 <- min(subset(df.act, PO==1)$TI_c)
        N <- length(unique(df.act$id))
        dft <- subset(df.act, TR==1)
        N.treat <- length(unique(dft$id))
        predat <- subset(orig$Resdat, post==0)
        postdat <- subset(orig$Resdat, post==1)
        RMSE_T1 <- sqrt(mean(postdat$effect^2))
        RMSE_T0 <-sqrt(mean(predat$effect^2))
        RMSE_RATIO <- RMSE_T1/RMSE_T0
        subsample <- as.numeric(0)
        dd_orig <- orig$supp_stats$dd
        param_p <- orig$supp_stats$dd_pval
        #Set the results data.frame

        resdat2 <- as.data.frame(cbind(RMSE_RATIO, RMSE_T0, RMSE_T1, param_p, T.min, T.max, T, T1, N, N.treat, subsample, row.names=NULL))
        resdat$realtr <- as.numeric(0)
        resdat2$realtr <- as.numeric(1)
        resdata <- rbind(resdat, resdat2)

        #Set the effects data.frame (for spaghetti plot)

        effect <- orig$Resdat$effect
        time <- orig$Resdat$time
        post <- orig$Resdat$post
        subsample <- rep(0, length(time))

        #Store test results

        effectmat2 <- as.data.frame(cbind(effect, time, subsample, post, row.names=NULL))
        effectmat$realtr <- as.numeric(0)
        effectmat2$realtr <- as.numeric(1)
        effectdata <- rbind(effectmat, effectmat2)


        #Get inverse ranks and probabilities based on the number of unique ranks (adjusts for the possibility that the subsample is picked more than once)
        rankdata <- as.data.frame(rep(NA, nrow(resdata)))
        rankdata$RMSE_RATIO <- resdata$RMSE_RATIO
        rankdata$realtr <- resdata$realtr
        rankdata <- unique(rankdata)

        rankdata$rank <- rank(-rankdata$RMSE_RATIO)
        rankdata$prob <- rankdata$rank/length(unique((rankdata$rank)))

        #Inverse rank for treatment unit and p-values based on this rank.

        rank <- subset(rankdata, realtr==1)$rank
        denom <- length(unique((rankdata$rank)))
        pval <- subset(rankdata, realtr==1)$prob

        #Get p-value based on empirical distribution function

        Fn <- stats::ecdf(unique(resdata$RMSE_RATIO))
        pval_ecdf <- 1-Fn(RMSE_RATIO)

        #Store p-value data

        supp_stats <- as.data.frame(cbind(rank, denom, pval, pval_ecdf, dd_orig))

        #Print

        cat("\n", "Results", "\n",
            "-------", "\n",
            "Rank: ", rank, "\n",
            "Denominator (unique subsamples): ", denom, "\n",
            "Placebo p-value (rank/unique subsamples): ", pval, "\n", "\n",
            "Notes", "\n",
            "-----", "\n",
            "Use the placeboci function to obtain interval estimates.", sep="")


        #Store and return the full results matrix and the empirical distribution function Fn(RMSE_RATIO).

        reslist <- as.list(NULL)
        reslist[[1]] <- resdata
        reslist[[2]] <- Fn
        reslist[[3]] <- effectdata
        reslist[[4]] <- RMSE_RATIO
        reslist[[5]] <- supp_stats
        names(reslist) <- c("Resdata", "ECDF", "Effects", "Treat.ratio", "supp_stats")
        base::close(pb)
        return(reslist)
}
