test_that("simdat() generates a data frame", {
  expect_equal(is.data.frame(simdat(c(100, 100, 100), 1)), TRUE)
})
