test_that("NNT Gives Right Statistic", {
  set.seed(1)

  pros.data = read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
  X=pros.data[pros.data$age>65, c('lpsa','lweight','lcp','lbph')]
  Y=pros.data[pros.data$age<=65, c('lpsa','lweight','lcp','lbph')]
  Z = rbind(X, Y)
  NNT_statistic <- NNT(Z, 1:nrow(Z), c(nrow(X), nrow(Y)), 2)
  expect_equal(0.53608247, NNT_statistic)
})

test_that("NNT Perm Gives Decision", {
  set.seed(1)

  pros.data = read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
  X=pros.data[pros.data$age>65, c('lpsa','lweight','lcp','lbph')]
  Y=pros.data[pros.data$age<=65, c('lpsa','lweight','lcp','lbph')]
  Z = rbind(X, Y)
  decision <- NNT.perm(Z, c(nrow(X), nrow(Y)), R = 2, 0.05, 1000)
  expect_equal(FALSE, decision)
})
