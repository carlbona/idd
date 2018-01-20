#' @title Plot time-varying placebo-based p-values.
#' @description \code{placebo.pplot} takes the results from \code{placeboci} and plots the time-varying p-values.
#'
#' @param object An R object obtained using \code{placeboci}.
#' @param alpha Alpha value for the significance region on the plot, defaults to 0.05.
#'
#' @return Returns a plot containing time-specific placebo-based confidence intervals.
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
#' placebo.pplot(ci.out)
#' }
#' @export

placebo.pplot <- function(object, alpha=0.05) {
        post=NULL
        pval=NULL
        time=NULL
        x <- object
        hl <- alpha
        vl <- min(subset(x, post==1)$time)
        p <- ggplot2::ggplot(x, ggplot2::aes(y=pval, x=time)) + ggplot2::geom_point() +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                ggplot2::geom_hline(yintercept=hl, linetype=2) + ggplot2::xlab("Time") + ggplot2::ylab("P-value") + ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::coord_cartesian(ylim=c(0, 1)) + ggplot2::ggtitle("Time-specific placebo-based p-values") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        return(p)
}
