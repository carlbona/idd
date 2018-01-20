#' @title Plot the raw data.
#' @description \code{idd.dplot} takes an output object after \code{idd} and plots the raw data using \code{ggplot2}.
#'
#' @param x Name of the \code{idd} output object.
#' @param mult Multiplier for the rates (default=100000).
#' @param donor Plot data from the entire donor pool as well? Defaults to \code{FALSE}.
#'
#' @return Returns a ggplot of raw data (without standardization and removal of fixed effects).
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
#' plot.out <- idd.dplot(idd.out)
#' }
#' @export


idd.dplot <- function(x, mult=100000, donor=FALSE) {
        post=NULL
        y1=NULL
        y0.ctrl=NULL
        y0.donors=NULL
        time=NULL
        x <- x$Resdat
        x$y1 <- x$y1*mult
        x$y0.donors <- x$y0.donors*mult
        x$y0.ctrl <- x$y0.ctrl*mult
        vl <- min(subset(x, post==1)$time)
        if (donor == FALSE) {
                p <- ggplot2::ggplot(x, ggplot2::aes(y=y1, x=time)) + ggplot2::geom_line() + ggplot2::geom_line(ggplot2::aes(y=y0.ctrl), linetype=2, color="grey10") +
                        ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) + ggplot2::xlab("Time") + ggplot2::ylab("Outcome") +
                        ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::ggtitle("Data, treated vs. control") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        }
        if (donor == TRUE) {
                p <- ggplot2::ggplot(x, ggplot2::aes(y=y1, x=time)) + ggplot2::geom_line() + ggplot2::geom_line(ggplot2::aes(y=y0.ctrl), linetype=2, color="grey10") + ggplot2::geom_line(ggplot2::aes(y=y0.donors), color="red", linetype=3) +
                        ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) + ggplot2::xlab("Time") + ggplot2::ylab("Outcome") +
                        ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::ggtitle("Data, treated vs. control/donors") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        }
        p
        return(p)
}
