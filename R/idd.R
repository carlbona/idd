#' @title Identify controls and estimate intervention effects.
#' @description \code{idd} takes a long-form panel dataset containing event count and exposure data, obtains the optimal controls and estimates an incidence difference-in-differences model.
#'
#' @param eventvar Name of the event count variable.
#' @param popvar The person-time (or another exposure) variable.
#' @param treatvar The treatment group indicator dummy (0 for controls, 1 for treated units).
#' @param postvar The post-period indicator dummy (0 for all time points in the pre-intervention period, 1 for all time points in the post-period)
#' @param timevar The time variable (can be coded on any scale).
#' @param idvar The panel id variable.
#' @param names A variable containing unit names (optional). Changes the output of id.selected to show unit names of selected controls.
#' @param print Print results? Default=TRUE.
#' @param data A long-form panel dataset containing the supplied variables.
#' @param mult Multiplier for the rates produced by print (default=100000).
#'
#' @return Returns a list containing the following elements:
#' \code{$Resdat}: a data frame containing the results.
#' \code{$cv_errors}: a data frame containing the cross-validation errors for the k-minimizing function.
#' \code{$supp_stats}: a data frame containing supplementary statistics (average effects, relative effects, p-values etc).
#' \code{$id_controls}: a data frame containing information on the selected controls.
#' @examples
#' \dontrun{
#' data(simpanel)
#' idd.out <- idd(eventvar="y",
#'              popvar="pop",
#'              idvar="age",
#'              timevar="time",
#'              postvar="post",
#'              treatvar="treat",
#'              data=simpanel)
#' }
#' @export



idd <- function(eventvar, popvar, treatvar, postvar, timevar, idvar, names=NULL, print=TRUE, data, mult=100000) {

        #Generate data
        E <- data[,eventvar]
        P <- data[,popvar]
        TR <- data[,treatvar]
        PO <- data[,postvar]
        TI <- data[,timevar]
        id <- data[,idvar]
        namevar <- as.character(data[,names])
        TI_c <- TI-(min(TI)-1)

        df <- as.data.frame(cbind(E, P, TR, PO, TI, TI_c, id))


        ##Prelimiary checks

        ##Missing data?

        #Missing event data?

        if (sum(is.na(df$E))>0) {
                stop("Event data is missing on at least one observation, please remove units with gaps.")
        }

        #Missing exposure data?

        if (sum(is.na(df$P))>0) {
                stop("Exposure data is missing on at least one observation, please remove units with gaps.")
        }

        #Balanced panel?

        tfirst <- doBy::summaryBy(TI~id, FUN=min, data=df, var.names=c("tmin"))
        first.test <- length(unique(tfirst$tmin))
        tlast <- doBy::summaryBy(TI~id, FUN=max, data=df, var.names=c("tmax"))
        last.test <- length(unique(tlast$tmax))

        if (first.test>1) {
                stop("Data is unbalanced, please supply a balanced panel dataset.")
        }
        if (last.test>1) {
                stop("Data is unbalanced, please supply a balanced panel dataset.")
        }


        #Subset pre-period.
        predf <- subset(df, PO==0)

        #More than 2 data points in the pre-period?

        if (max(predf$TI_c)<3) {
                stop("Pre-period is too short, the idd requires at least 3 data points to run.")
        }

        #At least 1 data point in the post-period?

        if (max(df$PO)<1) {
                stop("No post period found, need data before and after the intervention (code post period dummy as 0 in the before period and 1 after).")
        }


        #Treatment data, pre-period
        X1 <- subset(predf, TR==1)

        #Control data, pre-period.
        X0 <- subset(predf, TR==0)
        #Aggregate treatment group for use in matching later.
        sumtreat <- doBy::summaryBy(E + P ~ TI, FUN=sum, data=X1, var.names=c("y", "p"))
        sumtreat$inc <- sumtreat$y/sumtreat$p ##Compute incidence rate (not multiplied at this point)
        X1 <- cbind(sumtreat$TI, sumtreat$inc, row.names = NULL)
        X1 <- as.data.frame(X1)
        X1$id <- rep(-1337, nrow(X1))
        colnames(X1) <- c("time", "inc", "id")

        #Generate incidence rate for controls
        X0$inc <- X0$E/X0$P
        X0 <- cbind(X0$TI, X0$inc, X0$id, row.names = NULL)
        colnames(X0) <- c("time", "inc", "id")

        #Reshape wide X1 and X0 matrices wide.
        X0.w <- stats::reshape(as.data.frame(X0), idvar = "id", timevar = "time", direction = "wide")

        X1.w <- stats::reshape(as.data.frame(X1), timevar = "time", direction = "wide")

        #Bind data
        X <- rbind(X1.w, X0.w) ##First row is treatment group (we use this knowledge later.)


        #Take out fixed effects
        Xprim <- X
        for (i in 1:nrow(Xprim[,2:ncol(Xprim)])) {
                Xprim[i,2:ncol(Xprim)] <- Xprim[i,2:ncol(Xprim)]-rowMeans(Xprim[i,2:ncol(Xprim)])
        }
        cl <- scale(as.matrix(Xprim[,2:ncol(Xprim)]))

        #Calculate squared distances from trunit
        x_ctrl <- cl[2:nrow(cl),]
        x_tr <- cl[1,]
        sq_dist <- apply(x_ctrl, 1, function(x) (x_tr-x)^2)
        distvec <- rank(as.vector(sqrt(colSums(sq_dist))))

        #Gen treat output data

        treatout <- doBy::summaryBy(E + P ~ TI_c, FUN=sum, data=subset(df, TR==1), var.names=c("y1", "p1"), id="PO")

        ###################################################
        #Identify and cross-validate k nearest neighbours##
        ###################################################

        cv.error <- as.list(rep(NA),nrow(sq_dist)) ##Empty store list
        cv.error2 <- as.list(rep(NA), length(rankvec)) ##Empty store list
        for (t in 1:nrow(sq_dist)) { #Loop over each t in T0 (pre-period)
                rankvec <- rank(as.vector(sqrt(colSums(sq_dist[-t,])))) #Remove 1 obs (in time) and rank and calculate Euc distance
                idrank <- as.data.frame(cbind(X[-1,1], rankvec, row.names = NULL)) #Match ranking to id numbers
                colnames(idrank) <- c("id", "rank") #Change column names
                kneardat <- merge(df, idrank, by=c("id")) #set rank for t-training data
                for (k in 1:length(rankvec)) { #Loop over all possible selections of k (1 nearest through N nearest - i.e. entire sample)
                        ##Summarize control data for rank<=k
                        ctrlout <- doBy::summaryBy(E + P ~ TI_c, FUN=sum, data=subset(kneardat, rank<=k), var.names=c("y0", "p0"), id="PO")
                        #Compute cv error using estimation formula
                        outdat <- cbind(ctrlout, treatout, row.names = NULL)
                        predat <- subset(outdat, PO==0)
                        t0 <- sum(predat$y1.sum)/sum(as.numeric(predat$p1.sum))
                        c0 <- sum(predat$y0.sum)/sum(as.numeric(predat$p0.sum))
                        y0 <- outdat$y0.sum/outdat$p0.sum
                        y1 <- outdat$y1.sum/outdat$p1.sum
                        time <- outdat$TI_c
                        cf <- (t0-c0)+y0
                        effect <- y1-cf
                        temp <- as.data.frame(cbind(effect, time, row.names = NULL))
                        cv.error[[k]] <- subset(temp$effect, time==t)^2 #Store squared cv error for k
                }
                cv.error2[[t]] <- plyr::ldply(cv.error, data.frame) #Store for t, k
                cv.error2[[t]][,2] <- 1:length(rankvec) #add reference number for k (for summary)
        }
        errordat <- plyr::ldply(cv.error2, data.frame) #Flatten list to data frame
        colnames(errordat) <- c("cv", "k")

        k.mean <- doBy::summaryBy(cv~k, FUN=mean, data=errordat) #Compute mean cross-validation error for each k
        k.mean$cv.mean <- sqrt(k.mean$cv.mean)
        k.min <- which.min(k.mean$cv.mean) #Extract best selection for k according to leave-one-out cross validation on RMSE
        cv.rmse <- k.mean$cv.mean[k.min]

        #######################################
        ##Obtain effect estimates for best k ##
        #######################################

        rankvec <- rank(as.vector(sqrt(colSums(sq_dist)))) #Rank using all T0 obs
        distance <- as.vector(sqrt(colSums(sq_dist)))
        idrank <- as.data.frame(cbind(X[-1,1], rankvec, distance, row.names = NULL)) #Match ranking to id numbers
        colnames(idrank) <- c("id", "rank", "euc_dist") #Change column names
        ctrldat <- merge(df, idrank, by=c("id")) #set ranks for matrix of control data

        #Obtain ids of k*-controls:

        if (is.null(names)) {
                id.selected <- subset(as.data.frame(idrank), idrank$rank<=k.min)

        }
        else if (!is.null(names)) { #If the names variable has been specified, supply names as well.
                idx <- subset(as.data.frame(idrank), idrank$rank<=k.min)
                namedat <- cbind(id, namevar, row.names = NULL)
                colnames(namedat) <- c("id", "name")
                id.selected <- merge(idx, namedat, by=c("id"))
                colnames(id.selected) <- c("id", "rank", "euc_dist","name")
                id.selected <- stats::aggregate(cbind(id, rank, euc_dist, row.names=NULL) ~ name, data=id.selected, FUN=mean)
        }


        #Summarize control data
        ctrlout <- doBy::summaryBy(E + P ~ TI, FUN=sum, data=subset(ctrldat, rank<=k.min), var.names=c("y0", "p0"), id="PO")
        donorpool <- doBy::summaryBy(E + P ~ TI, FUN=sum, data=ctrldat, var.names=c("yx", "px"), id="PO")

        #Summarize treat data
        treatout <- doBy::summaryBy(E + P ~ TI, FUN=sum, data=subset(df, TR==1), var.names=c("y1", "p1"), id="PO")
        #Combine
        outdat <- cbind(ctrlout, treatout, donorpool, row.names = NULL)
        #Get pre-period for means
        predat <- subset(outdat, PO==0)

        ##compute stats: prep

        t0 <- sum(as.numeric(predat$y1.sum))/sum(as.numeric(predat$p1.sum))
        c0 <- sum(as.numeric(predat$y0.sum))/sum(as.numeric(predat$p0.sum))
        y0 <- outdat$y0.sum/outdat$p0.sum
        y1 <- outdat$y1.sum/outdat$p1.sum
        time <- outdat$TI

        #Get original donor pool data for comparison to matched controls
        cx <- sum(as.numeric(predat$yx.sum))/sum(as.numeric(predat$px.sum))
        yx <- outdat$yx.sum/outdat$px.sum

        #Compute time-varying effects
        cf <- (t0-c0)+y0
        effect <- y1-cf
        cf.donor <- (t0-cx)+yx

        se <- sqrt(sum(predat$y1.sum)/(sum(as.numeric(predat$p1.sum))^2)+sum(predat$y0.sum)/(sum(as.numeric(predat$p0.sum))^2)+outdat$y0.sum/(as.numeric(outdat$p0.sum^2))+outdat$y1.sum/(as.numeric(outdat$p1.sum^2)))
        z <- abs(effect/se)
        pval <- 2*stats::pnorm(-z)

        #Store time-varying effects matrix
        post <- outdat$PO
        y0.ctrl <- y0
        y0.donors <- yx
        res <- as.data.frame(cbind(y1,y0.ctrl,cf,effect,se,pval,time,post,cf.donor,y0.donors, row.names = NULL))

        #Compute DD estimates for average effects
        postdat <- subset(outdat, PO==1)
        t1 <- sum(as.numeric(postdat$y1.sum))/sum(as.numeric(postdat$p1.sum))
        c1 <- sum(as.numeric(postdat$y0.sum))/sum(as.numeric(postdat$p0.sum))
        c1x <- sum(as.numeric(postdat$yx.sum))/sum(as.numeric(postdat$px.sum))
        dd <- (t1-t0)-(c1-c0)
        dd_se <- sqrt(sum(predat$y1.sum)/(sum(as.numeric(predat$p1.sum))^2)+sum(predat$y0.sum)/(sum(as.numeric(predat$p0.sum))^2)+sum(postdat$y1.sum)/(sum(as.numeric(postdat$p1.sum))^2)+sum(postdat$y0.sum)/(sum(as.numeric(postdat$p0.sum))^2))
        dd_z <- abs(dd/dd_se)
        dd_pval <- 2*stats::pnorm(-dd_z)
        dd_donor <- (t1-t0)-(c1x-cx)

        #Compute cumulative effects (in events, entire T1-period - based on time-varying effect)
        posteffect <- subset(res, post==1)$effect
        postse <- subset(res, post==1)$se
        cm.post <- sum(posteffect*as.numeric(postdat$p1.sum))
        cm.se <- abs(cm.post/dd_z)
        cm.events.lower <- (cm.post-cm.se*stats::qnorm(0.975))
        cm.events.upper <- (cm.post+cm.se*stats::qnorm(0.975))
        cm.events <- cm.post

        #Compute relative effect

        post.events <- sum(postdat$y1.sum)
        ratio <- (post.events)/(post.events-cm.events)
        lnratio <- log(ratio)
        rrse <- sqrt((1/sum(predat$y1))+(1/sum(predat$y0))+(1/sum(postdat$y1))+(1/sum(postdat$y0)))
        ratio.up <- exp(lnratio+stats::qnorm(0.975)*rrse)
        ratio.lo <- exp(lnratio-stats::qnorm(0.975)*rrse)


        #Store results

        #Spit out some basic info first. See stored object for time-varying results.

        if (isTRUE(print)) {

                print(paste0("k*: ", k.min))
                print(paste0("DD est (avg effect, 100.000 pop): ", round(dd*mult, 4)))
                print(paste0("DD est (se): ", round(dd_se*mult, 4)))
                print(paste0("DD est (p-value): ", round(dd_pval, 5)))
                print(paste0("Cumulative effect (in events): ", round(cm.events, 1)))
                print(paste0("Cumulative effect (lower, 95% CI): ", round(cm.events.lower, 1)))
                print(paste0("Cumulative effect (upper, 95% CI): ", round(cm.events.upper, 1)))
                print(paste0("Ratio effect (compared to counterfactual): ", round(ratio, 4)))
                print(paste0("Ratio effect (lower, 95% CI): ", round(ratio.lo, 4)))
                print(paste0("Ratio effect (upper, 95% CI): ", round(ratio.up, 4)))

        }

        results <- as.list(NULL)
        results[[1]] <- res
        results[[2]] <- k.mean
        results[[3]] <- as.data.frame(cbind(k.min,dd,dd_se,dd_pval,cm.events, cm.events.lower, cm.events.upper, ratio, ratio.lo, ratio.up, cv.rmse, dd_donor, row.names = NULL))
        results[[4]] <- id.selected
        names(results) <- c("Resdat", "cv_errors", "supp_stats", "id_controls")

        return(results)
}
