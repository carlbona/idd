#' A panel dataset of residential fires in 153 rescue service confederations in Sweden (2000-2015).
#'
#'
#' @format A data frame with N=153 units observed over T=17 years, containing the following variables:
#' \describe{
#'   \item{year}{actual year}
#'   \item{id}{unit id}
#'   \item{rtjf_15}{unit name}
#'   \item{pop}{population size}
#'   \item{fires}{number of residential fire events}
#'   \item{treat}{treatment unit dummy variable}
#'   \item{post}{post-intervention dummy, coded as 1 if year>=2010}
#' }
#' @examples
#' \dontrun{
#' #This code replicates the results from the paper.
#'
#' data(homevisits) #Open the data
#'
#' #Find the optimal value for k and estimate effects
#'
#' firmod <- idd(eventvar="fires",
#'              popvar="pop",
#'              idvar="id",
#'              timevar="year",
#'              postvar="post",
#'              treatvar="treat",
#'              names="rtjf_15",
#'              data=homevisits)
#'
#' firmod$id_controls #Show the selected controls
#'
#' idd.kplot(firmod) #Plot the k-minimization function (figure 1)
#'
#' #Reproduce figure 2 using grid and gridExtra
#'
#' A <- idd.dplot(firmod, donor=TRUE)
#' B <- idd.cfplot(firmod, mult=100000, donor=TRUE)
#' C <- idd.gplot(firmod, mult=100000)
#' D <- idd.pplot(firmod)
#'
#' require(grid)
#' require(gridExtra)
#'
#' grid.arrange(A, B, C, D, ncol=2)
#'
#' #Perform the placebo studies (fair warning: 1000 iterations will take a while; expect 15-30 min)
#'
#' set.seed(12049135) ##Set seed to replicate the results in the paper
#'
#' placebt <- iddplacebo(eventvar="fires",
#'                      popvar="pop",
#'                      idvar="id",
#'                      timevar="year",
#'                      postvar="post",
#'                      treatvar="treat",
#'                      data=homevisits,
#'                      iter=1000)
#'
#' #Obtain placebo-based confidence intervals and p-values
#'
#' pci <- placeboci(placebt, alpha=0.05)
#'
#' #Reproduce figure 3
#'
#' A1 <- iddplacebo.ecdf(placebt)
#' B1 <- iddplacebo.hist(placebt)
#' C1 <- placebo.ciplot(pci)
#' D1 <- placebo.pplot(pci)
#' grid.arrange(A1, B1, C1, D1, ncol=2)
#' }
"homevisits"

#' A toy panel dataset of five age groups, with two poorly fitting and two good fitting controls.
#'
#'
#' @format A small toy dataset containing the following variables:
#' \describe{
#'   \item{age}{age group, used as id variable}
#'   \item{treat}{treatment unit dummy variable}
#'   \item{y}{outcome variable, n events}
#'   \item{time}{time variable}
#'   \item{post}{post-intervention dummy}
#'   \item{pop}{population denominator}
#'   \item{inc}{incidence rate (y/pop), not used}
#' }
"simpanel"
