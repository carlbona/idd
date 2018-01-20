#' @title Plot the estimated effects from the placebo studies.
#' @description \code{iddplacebo.spaghetti} takes an output object after \code{iddplacebo} and plots the estimated placebo effects using \code{ggplot2}.
#'
#' @param x Name of the \code{iddplacebo} output object.
#' @param mult Multiplier for the rates. Defaults to 100000.
#' @param rm Remove placebo studies with poor pre-intervention fit? Default=TRUE.
#' @param alpha controls the shaded significance region on the plot. Defaults to 0.05.
#'
#' @return Returns a ggplot containing the estimated effects from the placebo studies along with the results from the main analysis.
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
#' plot.out <- iddplacebo.spaghetti(placebo.out)
#' }
#' @export
#'
iddplacebo.spaghetti <- function(x, mult=100000, rm=TRUE, alpha=0.05) {

        post=NULL
        realtr=NULL
        effect=NULL
        subsample=NULL
        time=NULL

        z <- x$Resdata
        x <- x$Effects
        x$effect <- x$effect*mult
        vl <- min(subset(x, post==1)$time)
        x0 <- subset(x, realtr==0)
        x1 <- subset(x, realtr==1)
        z0 <- subset(z, realtr==0)
        fun <- stats::ecdf(z$RMSE_T0)
        z$pvec <- fun(z$RMSE_T0)
        x0 <- merge(z, x0, by=c("subsample"))
        if (rm == TRUE) {
                x0 <- x0[which(x0$pvec<=(1-alpha)),]
                p <- ggplot2::ggplot(x0, ggplot2::aes(y=effect, x=time)) + ggplot2::geom_line(ggplot2::aes(group=subsample), colour="gray80") + ggplot2::geom_line(data=x1, ggplot2::aes(y=effect, x=time)) +
                        ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                        ggplot2::geom_hline(yintercept=0, linetype=2) + ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::xlab("Time") + ggplot2::ylab("Effect estimate") + ggplot2::ggtitle("Placebo studies") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))

        }


        if (rm == FALSE) {
                p <- ggplot2::ggplot(x0, ggplot2::aes(y=effect, x=time)) + ggplot2::geom_line(ggplot2::aes(group=subsample), colour="gray80") + ggplot2::geom_line(data=x1, ggplot2::aes(y=effect, x=time)) +
                        ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                        ggplot2::geom_hline(yintercept=0, linetype=2) + ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::xlab("Time") + ggplot2::ylab("Effect estimate") + ggplot2::ggtitle("Placebo studies") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        }
        return(p)
}
