#' @title Plot the estimated effect.
#' @description \code{idd.gplot} takes an output object after \code{idd} and plots the estimated effect using \code{ggplot2}.
#'
#' @param x Name of the \code{idd} output object.
#' @param mult Multiplier for the rates (default=100000).
#'
#' @return Returns a ggplot of the estimated time-varying effects based on the selected control group.
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
#' plot.out <- idd.gplot(idd.out)
#' }
#' @export

idd.gplot <- function(x, mult=100000) {
        post=NULL
        effect=NULL
        se=NULL
        time=NULL
        x <- x$Resdat
        x$effect <- x$effect*mult
        x$se <- x$se*mult
        vl <- min(subset(x, post==1)$time)
        p <- ggplot2::ggplot(x, ggplot2::aes(y=effect, x=time)) + ggplot2::geom_line(linetype=1) + ggplot2::geom_ribbon(ggplot2::aes(ymin=effect-se*1.96, ymax=effect+se*1.96), fill="black", alpha=0.2) +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                ggplot2::geom_hline(yintercept=0, linetype=2) + ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::xlab("Time") + ggplot2::ylab("Effect estimate") + ggplot2::ggtitle("Time-specific effect estimates (95% CI)") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        p
        return(p)
}
