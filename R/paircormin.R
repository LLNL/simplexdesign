#' Minimize Pairwise Correlation
#'
#' @description Optimizes run order within classes of homogeneous agents to minimize pairwise correlation using Simulated Annealing.
#'
#' @usage paircormin(Design, kvec, ti = 10, tf = 0.0001, alph = 0.99, iter = 2000)
#'
#'
#' @param Design The unoptimized design
#' @param kvec A vector containing the number of columns in each block of the design
#' @param ti Initial temperature value for SA
#' @param tf Final temperature value for SA
#' @param alph Decay parameter for SA (should be less than 1)
#' @param iter Number of iterations at each temperature value, should increase with design size
#'
#'
#' @return \describe{
#' \item{Design}{The design optimized to minimize pairwise correlations between input classes.}
#' \item{CritVals}{The pairwise correlations at each temperature change, useful for checking if the SA algorithm has converged.}
#' }
#' @details This function optimizes run order within each block of a design made of multiple simplexes to minimize pairwise correlations using a Simulated Annealing algorithm.  The criterion of interest is the sum of squared cor(\eqn{xi},\eqn{xj}) for all pairs of columns between classes.  It is recommended that the parameters (such as ti, tf, and iter) are scaled with the size of the design size, and that the criterion values should settle to some local optimal value.
#'
#' @export
#'
#' @examples
#'
#' # Generate Unoptimized Design
#' D1 <- MmSimplex(3,30,10,cords = 1,randst = 1,phival = 50)
#' D2 <- MmSimplex(3,30,10,cords = 1,randst = 1,phival = 50)
#' D3 <- MmSimplex(3,30,10,cords = 1,randst = 1,phival = 50)
#' D <- cbind(D1[[1]],D2[[1]],D3[[1]])
#' c1 <- sum(cor(D1[[1]])[upper.tri(cor(D1[[1]]))]^2)
#' c2 <- sum(cor(D2[[1]])[upper.tri(cor(D2[[1]]))]^2)
#' c3 <- sum(cor(D3[[1]])[upper.tri(cor(D3[[1]]))]^2)
#' ## Sum all pariwise correlations and
#' sum(cor(D)[upper.tri(cor(D))]^2) - (c1 + c2 + c3)
#' ## Optimize run order
#' Dopt <- paircormin(D, c(3,3,3), ti = 0.01, tf = 0.0005, iter = 25)
#' sum(cor(Dopt[[1]])[upper.tri(cor(Dopt[[1]]))]^2) - (c1 + c2 + c3)
#'
#'
paircormin <- function(Design, kvec, ti = 10, tf = 1e-04, alph = 0.99, iter = 2000) {

  # internal function for swapping to items in an array
  swap2 <- function(inpvect, i, j, cols) {
    temp <- inpvect[i, cols]
    inpvect[i, cols] <- inpvect[j, cols]
    inpvect[j, cols] <- temp
    return(inpvect)
  }



  paircor2 <- function(Design, kvec) {
    cormat <- stats::cor(Design)
    for (i in 1:length(kvec)) {
      values <- (sum(kvec[0:(i - 1)]) + 1):sum(kvec[1:i])
      cormat[values, values] <- matrix(rep(0, kvec[i]^2), ncol = kvec[i])
    }
    return(sum(cormat^2))
  }
  paircormat2 <- function(Design, kvec) {
    cormat <- stats::cor(Design)
    for (i in 1:length(kvec)) {
      values <- (sum(kvec[0:(i - 1)]) + 1):sum(kvec[1:i])
      cormat[values, values] <- matrix(rep(0, kvec[i]^2), ncol = kvec[i])
    }
    return(cormat)
  }

  pb <- progress::progress_bar$new(format = "Optimization in Progress [:bar] :percent eta: :eta", total = ceiling(log(tf/ti)/log(alph)),
    clear = FALSE)

  phibest <- NULL
  phistore <- NULL
  phibest <- 1e+05
  testd2 <- Design
  cols <- list()
  for (i in 1:length(kvec)) {
    cols[[i]] <- (sum(kvec[0:(i - 1)]) + 1):sum(kvec[0:i])
  }
  kl <- length(kvec)
  N <- length(Design[, 1])

  while (ti > tf) {
    for (i in 1:iter) {
      runsample <- sample(N, 2)
      colv <- sample(1:kl, 1)
      temprun <- swap2(testd2, runsample[1], runsample[2], cols[[colv]])
      phivals <- paircor2(temprun, kvec)
      if (phibest > phivals) {
        phibest <- phivals
        testd2 <- temprun
      }

      if (phibest < phivals) {
        if (exp(-(phivals - phibest)/ti) > stats::runif(1)) {
          phibest <- phivals
          testd2 <- temprun
        }
      }
    }
    ti <- ti * alph
    pb$tick()
    phistore <- c(phistore, phibest)
  }
  return(list(Design <- testd2, CritVals <- phistore/2))
}



