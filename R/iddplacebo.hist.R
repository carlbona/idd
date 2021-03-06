#' @title Plot a histogram of placebo effects.
#' @description \code{iddplacebo.spaghetti} takes an output object after \code{iddplacebo} and plots a histogram of post/pre-RMSE ratios using \code{ggplot2}.
#'
#' @param x Name of the \code{iddplacebo} output object.
#' @param binwidth Control the binwidth (optional). Defaults to 2*IQR(RMSE)/(length(RMSE)^(1/3)).
#' @param convert Convert to effect sizes? Defaults to FALSE, presenting RMSE Ratios.
#' @param mult Multiplier for the effect size (only used with convert=TRUE).
#' @param alpha Alpha for the acceptance region reference line.
#' @param quantile Add reference line for the acceptance region? Defaults to FALSE.
#'
#' @return Returns a ggplot containing a histogram of the post/pre-RMSE ratios from the placebo studies.
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
#' plot.out <- iddplacebo.hist(placebo.out)
#' }
#' @export
#'

iddplacebo.hist <- function(x, binwidth=NULL, convert=FALSE, mult=100000, alpha=0.05, quantile=FALSE) {
        ..count..=NULL
        time=NULL
        ratio <- x$Treat.ratio
        att <- x$supp_stats$dd
        key <- abs(att)/ratio
        df <- x$Resdata
        df <- df[!duplicated(df),]
        RMSE <- unique(df$RMSE_RATIO)

        if (isTRUE(convert)==TRUE) {
                RMSE <- (RMSE*key)*mult
                densf1 <- stats::density(RMSE)
                AR <- spatstat::quantile.density(densf1, p=1-alpha, warn=FALSE) #acceptance region
                if (is.null(binwidth)) {
                        binwidth <- 2*stats::IQR(RMSE)/(length(RMSE)^(1/3))
                }
                ratio <- abs(att)*mult
                p <- ggplot2::ggplot() + ggplot2::geom_histogram(ggplot2::aes(x=RMSE, ..count../sum(..count..)), binwidth=binwidth, boundary = 0, color="black", fill="gray50") +
                        ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                        ggplot2::geom_hline(yintercept=0, linetype=2) + ggplot2::geom_vline(xintercept=ratio, linetype=2) + ggplot2::xlab("Effect size (absolute)") + ggplot2::ylab("Density") + ggplot2::ggtitle("Histogram (converted to effect size)") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
                if (isTRUE(quantile)==TRUE) {
                plotinfo <- ggplot2::ggplot_build(p)
                ymax <- max(plotinfo$data[[1]]$ymax)
                p <- p + ggplot2::geom_text(ggplot2::aes(x=AR, label=base::paste((1-alpha)*100,"th quantile = ", round(AR, 2), sep=""), y=ymax/1.5), angle=90, vjust=-1) + ggplot2::geom_vline(xintercept=AR, linetype=3)
                }
        }
        if (isTRUE(convert)==FALSE){
                densf1 <- stats::density(RMSE)
                AR <- spatstat::quantile.density(densf1, p=1-alpha, warn=FALSE) #acceptance region
        if (is.null(binwidth)) {
                binwidth <- 2*stats::IQR(RMSE)/(length(RMSE)^(1/3))
        }
        p <- ggplot2::ggplot() + ggplot2::geom_histogram(ggplot2::aes(x=RMSE, ..count../sum(..count..)), binwidth=binwidth, boundary = 0, color="black", fill="gray50") +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                ggplot2::geom_hline(yintercept=0, linetype=2) + ggplot2::geom_vline(xintercept=ratio, linetype=2) + ggplot2::xlab("Post/Pre-RMSE Ratio") + ggplot2::ylab("Density") + ggplot2::ggtitle("Histogram (RMSE Ratio)") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        if (isTRUE(quantile)==TRUE) {
                plotinfo <- ggplot2::ggplot_build(p)
                ymax <- max(plotinfo$data[[1]]$ymax)
                p <- p + ggplot2::geom_text(ggplot2::aes(x=AR, label=base::paste((1-alpha)*100,"th quantile = ", round(AR, 2), sep=""), y=ymax/1.5), angle=90, vjust=-1) + ggplot2::geom_vline(xintercept=AR, linetype=3)
        }
        }
        return(p)
}
