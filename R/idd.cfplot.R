#' @title Plot the counterfactual.
#' @description \code{idd.cfplot} takes an output object after \code{idd} and plots the counterfactual using \code{ggplot2}.
#'
#' @param x Name of the \code{idd} output object.
#' @param mult Multiplier for the rates (default=100000).
#' @param donor Plot a counterfactual based on the entire donor pool as well? Defaults to \code{FALSE}.
#'
#' @return Returns a ggplot of the estimated counterfactual based on the selected control group.
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
#' plot.out <- idd.cfplot(idd.out)
#' }
#' @export

idd.cfplot <- function(x, mult=100000, donor=FALSE) {
        post=NULL
        y1=NULL
        cf=NULL
        se=NULL
        cf.donor=NULL
        time=NULL
        x <- x$Resdat
        x$y1 <- x$y1*mult
        x$cf <- x$cf*mult
        x$se <- x$se*mult
        x$cf.donor <- x$cf.donor*mult
        vl <- min(subset(x, post==1)$time)
        if (donor == FALSE) {
                p <- ggplot2::ggplot(x, ggplot2::aes(y=y1, x=time)) + ggplot2::geom_line() + ggplot2::geom_line(ggplot2::aes(y=cf), linetype=2, color="grey10") + ggplot2::geom_ribbon(ggplot2::aes(ymin=cf-se*1.96, ymax=cf+se*1.96), fill="black", alpha=0.2) +
                        ggplot2::theme_bw() + ggplot2::theme(panel.grid.major =  ggplot2::element_blank(), panel.grid.minor =  ggplot2::element_blank(), axis.line =  ggplot2::element_line(colour = "black")) + ggplot2::xlab("Time") + ggplot2::ylab("Outcome") +
                        ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::ggtitle("Estimated counterfactual (95% CI)") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        }
        if (donor == TRUE) {
                p <- ggplot2::ggplot(x, ggplot2::aes(y=y1, x=time)) + ggplot2::geom_line() + ggplot2::geom_line(ggplot2::aes(y=cf), linetype=2, color="grey10") + ggplot2::geom_ribbon(ggplot2::aes(ymin=cf-se*1.96, ymax=cf+se*1.96), fill="black", alpha=0.2) +
                        ggplot2::theme_bw() + ggplot2::theme(panel.grid.major =  ggplot2::element_blank(), panel.grid.minor =  ggplot2::element_blank(), axis.line =  ggplot2::element_line(colour = "black")) +  ggplot2::xlab("Time") + ggplot2::ylab("Outcome") +
                        ggplot2::geom_vline(xintercept=vl, linetype=2) +  ggplot2::ggtitle("Estimated counterfactual (95% CI)") + ggplot2::theme(plot.title = ggplot2::element_text(size=12)) + ggplot2::geom_line(ggplot2::aes(y=cf.donor), linetype=3, color="red")
        }
        p
        return(p)
}
