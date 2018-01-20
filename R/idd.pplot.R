#' @title Plot the p-values.
#' @description \code{idd.gplot} takes an output object after \code{idd} and plots time-varying p-values for the effect using \code{ggplot2}.
#'
#' @param x Name of the \code{idd} output object.
#' @param alpha Set alpha value for the significance region (defaults to 0.05).
#'
#' @return Returns a ggplot of the estimated time-varying p-values for the effect based on the selected control group.
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
#' plot.out <- idd.pplot(idd.out)
#' }
#' @export

idd.pplot <- function(x, alpha=0.05) {
        post=NULL
        pval=NULL
        time=NULL
        x <- x$Resdat
        hl <- alpha
        vl <- min(subset(x, post==1)$time)
        p <- ggplot2::ggplot(x, ggplot2::aes(y=pval, x=time)) + ggplot2::geom_point() +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                ggplot2::geom_hline(yintercept=hl, linetype=2) + ggplot2::xlab("Time") + ggplot2::ylab("P-value") + ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::coord_cartesian(ylim=c(0, 1)) + ggplot2::ggtitle("Time-specific p-values") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        p
        return(p)
}
