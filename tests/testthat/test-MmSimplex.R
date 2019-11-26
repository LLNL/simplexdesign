test_that("MmSimple Example", {
  set.seed(1)
  D1 <- MmSimplex(10,40,10,cords = 1,randst = 1,phival = 50)
  expect_known_value(D1, "MmSimplex.test.output")
})
