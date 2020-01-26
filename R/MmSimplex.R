#' Maxmin Simplex Design
#'
#' @description Constructs a Maximin space-filling design on a k-dimensional simplex.
#'
#' @usage MmSimplex(k,N,l,cords = 1, randst = 1, phival = 50, tol = 0.0001)
#'
#' @references  \itemize{
#'  \item Johnson, M. E., L. M. Moore, and D. Ylvisaker (1990).  Minimax and maximin distance designs. Journal of statistical planning and inference 26 (2). 131-148.
#'  \item Meyer, R. K., and C.J. Nachtsheim (1995). The coordinate-exchange algorithm for constructing exact optimal experimental designs. Technometrics 37 (1). 60-69.
#' }
#'
#' @param k Number of factors in the design
#' @param N Number of experimental runs
#' @param l Number of levels for the search grid
#' @param cords Number of coordinates to exchange simultaneously using an exchange algorithm
#' @param randst Number of random starts to the design
#' @param phival Value of p in the PhiP criterion
#' @param tol Tolerance of the optimization, the value for which an improvement smaller than this ends the optimization
#'
#' @return \describe{
#' \item{D1}{The optimal design among the random starts standardized to [0,1]}
#' \item{phival}{The phiP value of the design}
#' \item{phistore}{A vector of phiP values from each random start}
#' }
#' @details This function applies a coordinate-exchange algorithm to optimize the Maximin distance criterion for a simplex.
#' A maximin design maximizes the minimum interpoint distance and is one commonly used space-filling criterion.  We do not optimize this criterion directly, but rather optimizes \deqn{ \phi_p = ( \sum\sum d_{ij}^{-p})^{1/p} }
#' This is done for optimization purposes since it is a smoother criteria.
#' @export
#'
#' @examples
#' ## Generate Design
#' time1 <- Sys.time()
#' D1 <- MmSimplex(5,40,10,cords = 1,randst = 2,phival = 50)
#' Sys.time() - time1
#'
#' ##View best design scaled to [0,1]^5
#' D1[[1]]
#'
#' ##View PhiP criterion for best design among 2 random starts
#' D1[[2]]
#'
MmSimplex <- function(k, N, l, cords = 1, randst = 1, phival = 50, tol = 1e-04) {

  if (randst > 1) {
    pb <- progress::progress_bar$new(format = "Checking Random Starts [:bar] :percent eta: :eta", total = randst, clear = FALSE)
  }

  if (cords > k) {
    stop("cords must be less than or equal to k")
  }


  # generates all ordered settings of length q
  gengenords <- function(minv, maxv, q) {
    if (q == 1) {
      testmatnew <- minv:maxv
      testmatnew <- matrix(testmatnew, ncol = 1)
      return(testmatnew)
    }
    testm <- NULL
    testm2 <- NULL
    for (m in minv:maxv) {
      testm <- c(testm, rep(m, maxv + 1 - m))
    }
    tvals <- unique(testm)
    for (m in 1:length(tvals)) {
      testm2 <- c(testm2, (minv + (m - 1)):maxv)
    }
    testmat <- cbind(testm, testm2)
    if (q == 2) {
      return(testmat)
    }
    if (q > 2) {
      testmatnew <- NULL
      for (m in 1:length(testmat[, 1])) {
        repval <- testmat[m, 2]
        for (j in repval:maxv) {
          testmatnew <- rbind(testmatnew, c(testmat[m, ], j))
        }
      }
    }
    if (q > 3) {
      for (i in 4:q) {
        testmat <- testmatnew
        testmatnew <- NULL
        for (m in 1:length(testmat[, 1])) {
          repval <- testmat[m, (i - 1)]
          for (j in repval:maxv) {
          testmatnew <- rbind(testmatnew, c(testmat[m, ], j))
          }
        }

      }
    }

    return(testmatnew)
  }

  fast_phiP <- function(distmat, desmat, new, col, p = 50) {
    desmat[col, ] <- new
    newvs <- Rfast::dista(new, desmat)
    distmat[, col] <- newvs^-p
    distmat[col, ] <- newvs^-p
    fi_p <- sum(stats::as.dist(distmat))^(1/p)
    return(list(fi_p, distmat))
  }

  # phistore will contain all of the phiP values of the designs
  x <- matrix(rep(1, N), ncol = 1)
  x <- cbind(x, matrix(rep(0, N * (k - 1)), ncol = k - 1))

  for (j in 1:(N)) {
    x[j, ] <- sort(sample(1:l, k, replace = TRUE))
  }

  # check there are no repeats, if there are add more random points until there are none
  x <- unique(x)
  counter <- 0
  while (length(x[, 1]) < N) {
    x2 <- rep(1, N)
    x2 <- cbind(x2, matrix(rep(0, N * (k - 1)), ncol = k - 1))
    for (j in 1:(N)) {
      x2[j, ] <- sort(sample(1:l, k, replace = TRUE))
    }
    x <- rbind(x, x2)
    x <- unique(x)
    counter <- counter + 1
    if (counter > 10) {
      stop("Difficulty finding unique random starts.  Try increasing grid size.")
    }
  }

  # make sure right number of runs
  x <- x[1:N, ]

  xstore <- x
  phistore <- NULL

  distmat <- as.matrix(stats::dist(x))^-phival
  diststore <- distmat
  distx <- diststore
  phimin <- phicomp <- 1e+08
  for (u in 1:randst) {
    # Generate Random Start
    x <- matrix(rep(1, N), ncol = 1)
    x <- cbind(x, matrix(rep(0, N * (k - 1)), ncol = k - 1))

    for (j in 1:(N)) {
      x[j, ] <- sort(sample(1:l, k, replace = TRUE))
    }

    x <- unique(x)
    while (length(x[, 1]) < N) {
      x2 <- rep(1, N)
      x2 <- cbind(x2, matrix(rep(0, N * (k - 1)), ncol = k - 1))
      for (j in 1:(N)) {
        x2[j, ] <- sort(sample(1:l, k, replace = TRUE))
      }
      x <- rbind(x, x2)
      x <- unique(x)
    }


    x <- x[1:N, ]
    imp <- 1
    xbest <- x
    distbest <- as.matrix(stats::dist(xbest))^-phival
    distx <- as.matrix(stats::dist(xbest))^-phival
    phibest <- DiceDesign::phiP(xbest, p = phival)

    # Coordinate Exchange Algorithm
    fwd <- 1
    while (imp == 1) {

      upd <- 0
      for (j in 1:N) {
        if (fwd == 1) {
          for (i in 1:(k - (cords - 1))) {
          x <- xbest
          maxv <- l
          minv <- 1
          if (i > 1) {
            minv <- xbest[j, i - 1]

          }
          if (i < (k - (cords - 1))) {
            maxv <- xbest[j, i + cords]
          }

          # generate orders
          testvs <- gengenords(minv, maxv, cords)
          currentvalue <- x[j, c(i:(i + (cords - 1)))]
          if (maxv > minv) {
            for (q in 1:length(testvs[, 1])) {
            if (any(testvs[q, ] != currentvalue)) {
              x <- xbest
              distx <- distbest
              x[j, c(i:(i + (cords - 1)))] <- testvs[q, ]
              eval <- fast_phiP(distx, x, matrix(x[j, ], ncol = k), j, p = phival)
              phix <- eval[[1]]
              distx <- eval[[2]]
              # phix <- DiceDesign::phiP(x)
              if (phix < phibest) {
              xbest <- x
              distbest <- distx
              phibest <- phix
              upd <- 1
              }
            }
            }
          }




          }
        }
        if (fwd == 0) {
          for (i in (k - (cords - 1)):1) {
          x <- xbest
          maxv <- l
          minv <- 1
          if (i > 1) {
            minv <- xbest[j, i - 1]

          }
          if (i < (k - (cords - 1))) {
            maxv <- xbest[j, i + cords]
          }

          # generate triples
          testvs <- gengenords(minv, maxv, cords)
          currentvalue <- x[j, c(i:(i + (cords - 1)))]
          if (maxv > minv) {
            for (q in 1:length(testvs[, 1])) {
            if (any(testvs[q, ] != currentvalue)) {
              x <- xbest
              distx <- distbest
              x[j, c(i:(i + (cords - 1)))] <- testvs[q, ]
              eval <- fast_phiP(distx, x, matrix(x[j, ], ncol = k), j, p = phival)
              phix <- eval[[1]]
              distx <- eval[[2]]
              # phix <- DiceDesign::phiP(x)
              if (phix < phibest) {
              xbest <- x
              distbest <- distx
              phibest <- phix
              upd <- 1
              }
            }
            }
          }




          }
        }
      }


      fwd <- -(fwd - 0.5) + 0.5


      if (upd == 0) {
        imp <- 0
      }

      if (abs(phicomp - phibest) < tol) {
        imp <- 0
      }
      phicomp <- phibest
    }

    if (phimin > phibest) {
      xstore <- xbest
      diststore <- distbest
      phimin <- phibest
    }
    phistore <- c(phistore, DiceDesign::phiP((xbest - 1)/(l - 1)))
    if (randst > 1) {
      pb$tick()
    }
  }

  features <- c(sprintf("X%02d", seq(1, k)))
  colnames(xstore) <- features

  return(list(((xstore - 1)/(l - 1)), min(phistore), phistore))

}

