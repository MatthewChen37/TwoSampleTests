test_that("HTT Gives Right Statistic", {
  set.seed(1)

  pros.data = read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
  X=pros.data[pros.data$age>65, c('lpsa','lweight','lcp','lbph')]
  Y=pros.data[pros.data$age<=65, c('lpsa','lweight','lcp','lbph')]
  Z = rbind(X, Y)
  dst <- as.matrix(dist(Z))
  HTT_statistic <- HTT(Z, 1:nrow(Z), c(nrow(X), nrow(Y)))
  expect_equal(11.6763633, HTT_statistic)
})

test_that("HTT Perm Gives Decision", {
  set.seed(1)

  pros.data = read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
  X=pros.data[pros.data$age>65, c('lpsa','lweight','lcp','lbph')]
  Y=pros.data[pros.data$age<=65, c('lpsa','lweight','lcp','lbph')]
  Z = rbind(X, Y)
  decision <- HTT.perm(Z, c(nrow(X), nrow(Y)), 0.05, 1000)
  expect_equal(TRUE, decision)
})
