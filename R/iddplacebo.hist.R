#' @title Plot a histogram of placebo effects.
#' @description \code{iddplacebo.spaghetti} takes an output object after \code{iddplacebo} and plots a histogram of post/pre-RMSE ratios using \code{ggplot2}.
#'
#' @param x Name of the \code{iddplacebo} output object.
#' @param binwidth Control the binwidth (optional). Defaults to 2*IQR(RMSE)/(length(RMSE)^(1/3)).
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

iddplacebo.hist <- function(x, binwidth=NULL) {
        ..density..=NULL
        time=NULL
        ratio <- x$Treat.ratio
        df <- x$Resdata
        df <- df[!duplicated(df),]
        RMSE <- unique(df$RMSE_RATIO)
        if (is.null(binwidth)) {
                binwidth <- 2*stats::IQR(RMSE)/(length(RMSE)^(1/3))
        }
        p <- ggplot2::ggplot() + ggplot2::geom_histogram(ggplot2::aes(x=RMSE, ..density..), binwidth=binwidth, color="black", fill="gray50") +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                ggplot2::geom_hline(yintercept=0, linetype=2) + ggplot2::geom_vline(xintercept=ratio, linetype=2) + ggplot2::xlab("Post/Pre-RMSE Ratio") + ggplot2::ylab("Density") + ggplot2::ggtitle("Histogram (RMSE Ratio)") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        return(p)
}