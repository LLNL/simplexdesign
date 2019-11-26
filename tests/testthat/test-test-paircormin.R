test_that("paircormin Example", {
  set.seed(1)
  D1 <- MmSimplex(5,40,10,cords = 1,randst = 1,phival = 50)
  D2 <- MmSimplex(5,40,10,cords = 1,randst = 1,phival = 50)
  D3 <- MmSimplex(5,40,10,cords = 1,randst = 1,phival = 50)
  D <- cbind(D1[[1]],D2[[1]],D3[[1]])
  c1 <- sum(cor(D1[[1]])[upper.tri(cor(D1[[1]]))]^2)
  c2 <- sum(cor(D2[[1]])[upper.tri(cor(D2[[1]]))]^2)
  c3 <- sum(cor(D3[[1]])[upper.tri(cor(D3[[1]]))]^2)
  ## Sum all pariwise correlations and
  sum(cor(D)[upper.tri(cor(D))]^2) - (c1 + c2 + c3)
  ## Optimize run order
  Dopt <- paircormin(D, c(5,5,5), ti = 0.1, tf = 0.0005, iter = 50)
  sum(cor(Dopt[[1]])[upper.tri(cor(Dopt[[1]]))]^2) - (c1 + c2 + c3)
  expect_known_value(Dopt, "paircormin.test.output")
})
