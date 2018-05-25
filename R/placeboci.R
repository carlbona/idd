#' @title Estimate placebo-based confidence intervals.
#' @description \code{placeboci} takes the results from a set of placebo studies obtained using \code{iddplacebo} and estimates non-parametric placebo-based confidence intervals for a selected alpha level.
#'
#' @param object An R object containing the results from a set of placebo studies obtained using \code{iddplacebo}.
#' @param alpha The significance level, defaults to 0.05 for 95 percent CI.
#' @param mult Multiplier for the rate, defaults to 100000.
#'
#' @return The code prints the CI for the average post-intervention effect, and returns a data frame containing the time-varying CI. To be used with \code{iddplacebo.ciplot} or \code{iddplacebo.pplot}.
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
#' ci.out <- placeboci(placebo.out, alpha=0.05)
#' }
#' @export


placeboci <- function(object, alpha=0.05, mult=100000) {
        subsample = NULL
        time_c = NULL
        realtr.x = NULL
        realtr = NULL
        att <- object$supp_stats$dd
        att_R <- object$Treat.ratio
        att_pval <- object$supp_stats$pval
        df <- object$Effects
        df$RMSE_t <- sqrt(df$effect^2)
        #merge RMSE_T0 with effectdata by subsample
        df2 <- object$Resdata
        dnew <- merge(df, df2, by=c("subsample"))
        dnew$T_RATIO <- dnew$RMSE_t/dnew$RMSE_T0
        #store for treated unit
        xtemp <- subset(dnew, subsample==0)
        trat_tr <- xtemp$T_RATIO
        ti <- xtemp$time
        treat_rat <- as.data.frame(cbind(trat_tr,time))

        #find time range

        t.adjust <- min(dnew$time)-1
        dnew$time_c <- dnew$time-t.adjust
        Tmax <- max(dnew$time_c)


        #subset data by time point
        tsets <- as.list(rep(NA), Tmax)
        tpval <- as.list(rep(NA), Tmax)
        for (t in 1:Tmax) {
                tsets[[t]] <- subset(dnew, time_c==t)
                #Get p-values using inverse ranks
                rankdata <- as.data.frame(rep(NA, nrow(tsets[[t]]))) #Vectorize
                rankdata$T_RATIO <- tsets[[t]]$T_RATIO #Create new dataset for ranking
                rankdata$realtr.x <- tsets[[t]]$realtr.x #Obtain treatment indicator variable
                rankdata <- unique(rankdata) #Remove duplicates
                rankdata$rank <- rank(-rankdata$T_RATIO) #Find the inverse rank for R_t
                rankdata$prob <- rankdata$rank/length(unique((rankdata$rank))) #Find the p-value for each R_t
                rank <- subset(rankdata, realtr.x==1)$rank #Find the rank for the treated
                denom <- length(unique((rankdata$rank))) #Find the denominator (unique ranks)
                pval <- subset(rankdata, realtr.x==1)$prob #Find p-value
                time <- t+t.adjust
                effect <- subset(tsets[[t]], realtr.x==1)$effect

                #Find placebo-based CI using empirical quantiles

                f1 <- effect+unique(tsets[[t]]$T_RATIO)*subset(tsets[[t]], realtr.x==1)$RMSE_T0
                f2 <- effect-unique(tsets[[t]]$T_RATIO)*subset(tsets[[t]], realtr.x==1)$RMSE_T0
                densf1 <- stats::density(f1)
                densf2 <- stats::density(f2)

                t_rec_upper <- spatstat::quantile.density(densf1, p=1-alpha, warn=FALSE)

                t_rec_lower <- spatstat::quantile.density(densf2, p=alpha, warn=FALSE)

                #Store results

                tpval[[t]] <- as.data.frame(cbind(time, rank, denom, pval, effect, t_rec_lower, t_rec_upper))
        }
        res <- plyr::ldply(tpval, data.frame)
        #Obtain confidence interval for DD

        #Construct the key (see paper)

        key <- abs(att)/att_R

        res$post <- subset(df, realtr==1)$post
        #f1 <- att+unique(df2$RMSE_RATIO)*subset(df2, realtr==1)$RMSE_T0
        #f2 <- att-unique(df2$RMSE_RATIO)*subset(df2, realtr==1)$RMSE_T0
        f1 <- att+unique(df2$RMSE_RATIO)*key
        f2 <- att-unique(df2$RMSE_RATIO)*key
        densf1 <- stats::density(f1)
        densf2 <- stats::density(f2)

        att_rec_upper <- spatstat::quantile.density(densf1, p=1-alpha, warn=FALSE)

        att_rec_lower <- spatstat::quantile.density(densf2, p=alpha, warn=FALSE)


        #Store results

        res$att <- rep(att, nrow(res))
        res$att_pval <- rep(att_pval, nrow(res))
        res$att_lower <- rep(att_rec_lower, nrow(res))
        res$att_upper <- rep(att_rec_upper, nrow(res))

        #Print results

        cat("Difference-in-differences estimate (average effect): ", round(att*mult, 5), "\n",
            "Lower placebo-based CI for ATT: ", round(att_rec_lower*mult, 5), "\n",
            "Upper placebo-based CI for ATT: ", round(att_rec_upper*mult, 5), "\n", "\n",
            "Notes", "\n",
            "-----", "\n",
            "Selected alpha-level: ", alpha,". Time-varying results are stored in the object.", sep="")

        return(res)
}
