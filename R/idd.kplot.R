#' @title Plot the k-minimization function.
#' @description \code{idd.cfplot} takes an output object after \code{idd} and plots the k-minimization function using \code{ggplot2}.
#'
#' @param x Name of the \code{idd} output object.
#'
#' @return Returns a ggplot of the k-minimization function used to select the optimal value of k.
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
#' plot.out <- idd.kplot(idd.out)
#' }
#' @export

idd.kplot <- function(x) {
        cv.mean=NULL
        k=NULL
        x <- x$cv_errors
        min.k <- which.min(x$cv.mean)
        p <- ggplot2::ggplot(x, ggplot2::aes(y=cv.mean, x=k)) + ggplot2::geom_point() +
                ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                ggplot2::geom_vline(xintercept=min.k, linetype=2) + ggplot2::xlab("k") + ggplot2::ylab("Cross-validation error (RMSE)") + ggplot2::ggtitle("k-optimization function") + ggplot2::theme(plot.title = ggplot2::element_text(size=12))
        p
        return(p)
}
