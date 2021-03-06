% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iddplacebo.hist.R
\name{iddplacebo.hist}
\alias{iddplacebo.hist}
\title{Plot a histogram of placebo effects.}
\usage{
iddplacebo.hist(x, binwidth = NULL, convert = FALSE, mult = 100000,
alpha = 0.05, quantile = FALSE)
}
\arguments{
\item{x}{Name of the \code{iddplacebo} output object.}

\item{binwidth}{Control the binwidth (optional). Defaults to 2*IQR(RMSE)/(length(RMSE)^(1/3)).}

\item{convert}{Convert to effect sizes? Defaults to FALSE, presenting RMSE Ratios.}

\item{mult}{Multiplier for the effect size (only used with convert=TRUE).}

\item{alpha}{Alpha for the acceptance region reference line.}

\item{quantile}{Add reference line for the acceptance region? Defaults to FALSE.}
}
\value{
Returns a ggplot containing a histogram of the post/pre-RMSE ratios from the placebo studies.
}
\description{
\code{iddplacebo.spaghetti} takes an output object after \code{iddplacebo} and plots a histogram of post/pre-RMSE ratios using \code{ggplot2}.
}
\examples{
\dontrun{
data(simpanel)
placebo.out <- iddplacebo(eventvar="y",
             popvar="pop",
             idvar="age",
             timevar="time",
             postvar="post",
             treatvar="treat",
             data=simpanel,
             iter=50)
plot.out <- iddplacebo.hist(placebo.out)
}
}
