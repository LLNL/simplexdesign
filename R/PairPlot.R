#' Design Projection Plot
#'
#' Plots the projection of all pairs of inputs from a design from MmSimplex.  Can actually be used for any design on a Simplex.
#'
#' @usage PairPlot(Design, greys = TRUE)
#'
#' @param Design List object from MmSimplex, also works with a matrix input
#' @param greys Draws grey boxes to indicate the region of valid design points
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' ## Generate Design
#' D1 <- MmSimplex(3,75,10,cords = 1,randst = 1,phival = 50)
#'
#' ##Plot Design
#' PairPlot(D1)
PairPlot <- function(Design, greys = TRUE) {
  if (is.list(Design) == FALSE) {
    Design <- list(Design)
  }
  xstore <- Design[[1]]
  l <- max(Design[[1]])
  pm <- GGally::ggpairs(data.frame(xstore)) + ggplot2::theme_bw()
  pm2 <- pm
  if (greys == TRUE) {
    tmp <- data.frame(x = c(-0.05, -0.05, l + 0.05), y = c(-0.05, l + 0.05, l + 0.05))
    for (i in 2:pm$nrow) {
      for (j in 1:(i - 1)) {
        pm2[i, j] <- pm[i, j] + ggplot2::geom_polygon(data = tmp, ggplot2::aes(tmp$x, tmp$y), fill = "gray45", alpha = 0.3)
      }
    }
  }
  return(pm2)
}

