test_that("EBT Gives Right Statistic", {
  set.seed(1)

  pros.data = read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
  X=pros.data[pros.data$age>65, c('lpsa','lweight','lcp','lbph')]
  Y=pros.data[pros.data$age<=65, c('lpsa','lweight','lcp','lbph')]
  Z = rbind(X, Y)
  dst <- as.matrix(dist(Z))
  EBT_statistic <- EBT(dst, 1:nrow(dst), c(nrow(X), nrow(Y)))
  expect_equal(7.9336746, EBT_statistic)
})

test_that("EBT Perm Gives Decision", {
  set.seed(1)

  pros.data = read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
  X=pros.data[pros.data$age>65, c('lpsa','lweight','lcp','lbph')]
  Y=pros.data[pros.data$age<=65, c('lpsa','lweight','lcp','lbph')]
  Z = rbind(X, Y)
  dst <- as.matrix(dist(Z))
  decision <- EBT.perm(dst, c(nrow(X), nrow(Y)), 0.05, 1000)
  expect_equal(TRUE, decision)
})
