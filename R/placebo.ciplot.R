#' @title Plot time-varying placebo-based confidence intervals.
#' @description \code{placebo.ciplot} takes the results from \code{placeboci} and plots the time-varying confidence intervals for the effect.
#'
#' @param object An R object obtained using \code{placeboci}.
#' @param mult A multiplier for the rates, defaults to 100000.
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
#' ci.out <- placebo.ci(placebo.out, alpha=0.05)
#' placebo.ciplot(ci.out)
#' }
#' @export

placebo.ciplot <- function(object, mult=100000) {
        post=NULL
        effect=NULL
        placebo_lower=NULL
        placebo_upper=NULL
        pinv_lower=NULL
        pinv_upper=NULL
        time=NULL
        x <- object
        x$effect <- x$effect*mult
        vl <- min(subset(x, post==1)$time)
        x$placebo_upper <- x$t_rec_upper*mult
        x$placebo_lower <- x$t_rec_lower*mult

        p <- ggplot2::ggplot(x, ggplot2::aes(y=effect, x=time)) + ggplot2::geom_ribbon(ggplot2::aes(ymin=placebo_lower, ymax=placebo_upper), fill="black", alpha=0.2) + ggplot2::geom_line(ggplot2::aes(y=effect, x=time)) +
                        ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                        ggplot2::geom_hline(yintercept=0, linetype=2) + ggplot2::geom_vline(xintercept=vl, linetype=2) + ggplot2::xlab("Time") + ggplot2::ylab("Effect estimate") + ggplot2::ggtitle("Time-specific placebo-based 95% CI") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        return(p)

}
