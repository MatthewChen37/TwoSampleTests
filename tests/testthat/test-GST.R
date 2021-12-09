test_that("GST Gives Right Statistic", {
  set.seed(1)
  pros.data = read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
  X=pros.data[pros.data$age>65, c('lpsa','lweight','lcp','lbph')]
  Y=pros.data[pros.data$age<=65, c('lpsa','lweight','lcp','lbph')]
  Z = rbind(X, Y)
  dst <- as.matrix(dist(Z))
  GST_statistic <- GST(dst, 1:nrow(dst), c(nrow(X), nrow(Y)), 2)
  expect_equal(0.51609848, GST_statistic)
})

test_that("GST Perm Gives Decision", {
  set.seed(1)
  pros.data = read.table("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
  X=pros.data[pros.data$age>65, c('lpsa','lweight','lcp','lbph')]
  Y=pros.data[pros.data$age<=65, c('lpsa','lweight','lcp','lbph')]
  Z = rbind(X, Y)
  dst <- as.matrix(dist(Z))
  decision <- GST.perm(dst, c(nrow(X), nrow(Y)), Q=2, 0.05, 1000)
  expect_equal(FALSE, decision)
})
