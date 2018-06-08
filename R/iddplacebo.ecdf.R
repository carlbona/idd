#' @title Plot the inverse empirical cumulative distribution function from the placebo studies.
#' @description \code{iddplacebo.ecdf} takes an output object after \code{iddplacebo} and plots the E.C.D.F. using \code{ggplot2}.
#'
#' @param x Name of the \code{iddplacebo} output object.
#' @param alpha controls the shaded significance region on the plot. Defaults to 0.05.
#' @param convert Convert x-axis to effect size? Defaults to FALSE.
#' @param mult Multiplier for the effect size (only used with convert=TRUE).
#' @param quantile Add reference line for the acceptance region? Defaults to FALSE.
#'
#' @return Returns a ggplot containing the E.C.D.F. of the post/pre-RMSE ratios from the placebo studies.
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
#' plot.out <- iddplacebo.ecdf(placebo.out)
#' }
#' @export

iddplacebo.ecdf <- function(x, alpha=0.05, convert=FALSE, mult=100000, quantile=FALSE) {
        invp=NULL
        RMSE_RATIO=NULL
        x <- x
        Fn <- x$ECDF
        data <- x$Resdata
        ratio <- x$Treat.ratio
        att <- x$supp_stats$dd
        key <- abs(att)/ratio
        data$invp <- 1-Fn(data$RMSE_RATIO)
        pval <- x$supp_stats$pval

        if (isTRUE(convert)==TRUE) {
        data$RMSE_RATIO <- data$RMSE_RATIO*key*mult
        densf1 <- stats::density(unique(data$RMSE_RATIO))
        AR <- spatstat::quantile.density(densf1, p=1-alpha, warn=FALSE) #acceptance region
        Fn <- stats::ecdf(unique(data$RMSE_RATIO))
        data$invp <- 1-Fn(data$RMSE_RATIO)
        p <- ggplot2::ggplot(data, ggplot2::aes(y=invp, x=RMSE_RATIO)) + ggplot2::geom_step() + ggplot2::geom_vline(xintercept=abs(att)*mult, linetype=2) +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                ggplot2::geom_hline(yintercept=pval, linetype=2) + ggplot2::xlab("Effect size (absolute)") + ggplot2::ylab("Inverse Cumulative Probability") + ggplot2::geom_ribbon(ggplot2::aes(ymin=0, ymax=alpha), alpha=0.2) + ggplot2::ggtitle("Inverse E.C.D.F (Absolute effect size)") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))

        if (isTRUE(quantile)==TRUE) {
                p <- p + ggplot2::geom_text(ggplot2::aes(x=AR, label=base::paste((1-alpha)*100,"th quantile = ", round(AR, 2), sep=""), y=0.6), angle=90, vjust=-1) + ggplot2::geom_vline(xintercept=AR, linetype=3)
        }
        }

        if (isTRUE(convert)==FALSE) {
        densf1 <- stats::density(unique(data$RMSE_RATIO))
        AR <- spatstat::quantile.density(densf1, p=1-alpha, warn=FALSE) #acceptance region
        p <- ggplot2::ggplot(data, ggplot2::aes(y=invp, x=RMSE_RATIO)) + ggplot2::geom_step() + ggplot2::geom_vline(xintercept=ratio, linetype=2) +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                ggplot2::geom_hline(yintercept=pval, linetype=2) + ggplot2::xlab("Post/Pre-RMSE Ratio") + ggplot2::ylab("Inverse Cumulative Probability") + ggplot2::geom_ribbon(ggplot2::aes(ymin=0, ymax=alpha), alpha=0.2) + ggplot2::ggtitle("Inverse E.C.D.F (RMSE Ratio)") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))

        if (isTRUE(quantile)==TRUE) {
                p <- p + ggplot2::geom_text(ggplot2::aes(x=AR, label=base::paste((1-alpha)*100,"th quantile = ", round(AR, 2), sep=""), y=0.6), angle=90, vjust=-1) + ggplot2::geom_vline(xintercept=AR, linetype=3)
        }
        }

        return(p)
}
