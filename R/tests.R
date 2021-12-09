#' NNT Function
#'
#' Calculates the test statistic for the two-sampled Nearest Neighbor Test
#' @param z is the matrix of the data set (each row is an observation)
#' @param ix is a permutation of row indices of z
#' @param sizes is a vector of the sample sizes
#' @param R is the number of neighbors to look at
#' @export
#' @return The NNT test statistic
#' @examples
#' x = matrix(rexp(100 * 2), 100, 2)
#' y = matrix(rexp(200 * 2), 200, 2)
#' z <- rbind(x,y)
#' NNT_statistic <- NNT(z, 1:nrow(z), c(100, 200), 4)

NNT <- function(z, ix, sizes, R) {
  n1 <- sizes[1]
  n2 <- sizes[2]
  n <- n1 + n2
  z <- z[ix, ]
  NN <- RANN::nn2(z, z, k=R+1)
  block1 <- NN$nn.idx[1:n1, -1]
  block2 <- NN$nn.idx[(n1+1):n, -1]
  i1 <- sum(block1 < n1 + .5)
  i2 <- sum(block2 > n1 + .5)
  return((i1 + i2) / (R * n))
}

#' NNT.perm Function
#'
#' Runs a permutation test using the Nearest Neighbors test statistic
#' @param z is the matrix of the data set (each row is an observation)
#' @param sizes is a vector of the sample sizes
#' @param R is the number of neighbors to look at
#' @param alpha is the significance level
#' @param B is the number of replicates in the permutation test
#' @export
#' @return \code{TRUE} if it decides to reject the NNT null hypothesis or \code{FALSE} if not
#' @examples
#' x = matrix(rexp(100 * 2), 100, 2)
#' y = matrix(rexp(200 * 2), 200, 2)
#' z <- rbind(x,y)
#' decision <- NNT.perm(z, c(100, 200), 4, 0.05, 1000)

NNT.perm <- function(z, sizes, R,alpha,B) {
  T0=NNT(z,1:nrow(z),sizes, R)
  TPerm=rep(0,B)
  for (b in 1:B) {
    permindex=sample(1:nrow(z))
    TPerm[b]=NNT(z,permindex,sizes, R)
  }
  P <- mean(c(T0, TPerm) >= T0)
  return(P<alpha)
}

#' EBT Function
#'
#' Calculates the test statistic for the two-sampled Energy Distance Test
#' @param dst is the distance matrix of the data set
#' @param ix is a permutation of row indices of z
#' @param sizes is a vector of the sample sizes
#' @export
#' @return The EBT test statistic
#' @examples
#' x = matrix(rexp(100 * 2), 100, 2)
#' y = matrix(rexp(200 * 2), 200, 2)
#' z <- rbind(x,y)
#' dst <- as.matrix(dist(z))
#' EBT_statistic <- EBT(dst, 1:nrow(dst), c(100, 200))

EBT <- function(dst, ix, sizes) {
  n1 <- sizes[1]
  n2 <- sizes[2]
  ii <- ix[1:n1]
  jj <- ix[(n1+1):(n1+n2)]
  w <- n1 * n2 / (n1 + n2)
  m11 <- sum(dst[ii, ii]) / (n1 * n1)
  m22 <- sum(dst[jj, jj]) / (n2 * n2)
  m12 <- sum(dst[ii, jj]) / (n1 * n2)
  e <- w * ((m12 + m12) - (m11 + m22))
  return (e)
}

#' EBT.perm Function
#'
#' Runs a permutation test using the Energy Distance test statistic
#' @param dst is the distance matrix of the data set
#' @param sizes is a vector of the sample sizes
#' @param alpha is the significance level
#' @param B is the number of replicates in the permutation test
#' @export
#' @return \code{TRUE} if it decides to reject the EBT null hypothesis or \code{FALSE} if not
#' @examples
#' x = matrix(rexp(100 * 2), 100, 2)
#' y = matrix(rexp(200 * 2), 200, 2)
#' z <- rbind(x,y)
#' dst <- as.matrix(dist(z))
#' decision <- EBT.perm(dst, c(100, 200), 0.05, 1000)

EBT.perm <- function(dst, sizes,alpha,B) {
  T0=EBT(dst,1:nrow(dst),sizes)
  TPerm=rep(0,B)
  for (b in 1:B) {
    permindex=sample(1:nrow(dst))
    TPerm[b]=EBT(dst,permindex,sizes)
  }
  P <- mean(c(T0, TPerm) >= T0)
  return(P<alpha)
}

#' HTT Function
#'
#' Calculates the test statistic for the two-sampled Hotelling's T-squared Test
#' @param z is the matrix of the data set (each row is an observation)
#' @param ix is a permutation of row indices of z
#' @param sizes is a vector of the sample sizes
#' @export
#' @return The HTT test statistic
#' @examples
#' x = matrix(rexp(100 * 2), 100, 2)
#' y = matrix(rexp(200 * 2), 200, 2)
#' z <- rbind(x,y)
#' HTT_statistic <- HTT(z, 1:nrow(z), c(100, 200))
HTT <- function(z, ix, sizes){
  n <- sizes[1]
  m <- sizes[2]
  C <- (n * m) / (n + m)
  ii <- ix[1:n]
  jj <- ix[(n+1):(n+m)]

  X <- z[ii,]
  Y <- z[jj,]

  SampleVarX <- stats::var(X)
  SampleVarY <- stats::var(Y)

  Sigma_hat <- (1 / (n + m - 2)) * ((n - 1) * SampleVarX + (m - 1) * SampleVarY)
  t1 <- C * t(colMeans(X) - colMeans(Y))
  t2 <- solve(Sigma_hat) %*% (colMeans(X) - colMeans(Y))

  t_squared <- t1 %*% t2
  return(t_squared[1, 1]) # number is a 1x1 matrix so we convert it into a scalar
}

#' HTT.perm Function
#'
#' Runs a permutation test using the Hotelling's T-squared statistic
#' @param z is the matrix of the data set (each row is an observation)
#' @param sizes is a vector of the sample sizes
#' @param alpha is the significance level
#' @param B is the number of replicates in the permutation test
#' @export
#' @return \code{TRUE} if it decides to reject the HTT null hypothesis or \code{FALSE} if not
#' @examples
#' x = matrix(rexp(100 * 2), 100, 2)
#' y = matrix(rexp(200 * 2), 200, 2)
#' z <- rbind(x,y)
#' decision <- HTT.perm(z, c(100, 200), 0.05, 1000)
HTT.perm <- function(z, sizes, alpha, B){
  T0=HTT(z,1:nrow(z),sizes)
  TPerm=rep(0,B)
  for (b in 1:B) {
    permindex=sample(1:nrow(z))
    TPerm[b]=HTT(z,permindex,sizes)
  }
  P <- mean(c(T0, TPerm) >= T0)
  return(P<alpha)
}

#' GST Function
#'
#' Calculates the test statistic for the two-sampled Graph Based Test
#' @param dst is the distance matrix of the data set
#' @param ix is a permutation of row indices of z
#' @param sizes is a vector of the sample sizes
#' @param Q is the distance threshold parameter
#' @export
#' @return The GST test statistic
#' @examples
#' x = matrix(rexp(100 * 2), 100, 2)
#' y = matrix(rexp(200 * 2), 200, 2)
#' z <- rbind(x,y)
#' dst <- as.matrix(dist(z))
#' GST_statistic <- GST(dst, 1:nrow(dst), c(100, 200), 2)
GST <- function(dst, ix, sizes, Q){
  n <- sizes[1]
  m <- sizes[2]

  ii <- ix[1:n]
  jj <- ix[(n+1):(n + m)]

  num_of_edges <- (sum(dst < Q) - (n + m)) %/% 2

  # intepret counting the sum as counting the number of
  # edges where the two vertices are in the same sample

  sum_block1 <- (sum(dst[ii, ii] < Q) - n) / 2
  sum_block2 <- (sum(dst[jj, jj] < Q) - m) / 2

  if (is.na(sum_block1)) {
    sum_block1 <- 0
  }

  if (is.na(sum_block2)) {
    sum_block2 <- 0
  }


  if (num_of_edges == 0){
    t_stat <- 0
  }
  else{
    t_stat <- (sum_block1 + sum_block2) / num_of_edges
  }
  return(t_stat)
}

#' GST.perm Function
#'
#' Runs a permutation test using the Energy Distance test statistic
#' @param dst is the distance matrix of the data set
#' @param sizes is a vector of the sample sizes
#' @param Q is the distance threshold parameter
#' @param alpha is the significance level
#' @param B is the number of replicates in the permutation test
#' @export
#' @return \code{TRUE} if it decides to reject the EBT null hypothesis or \code{FALSE} if not
#' @examples
#' x = matrix(rexp(100 * 2), 100, 2)
#' y = matrix(rexp(200 * 2), 200, 2)
#' z <- rbind(x,y)
#' dst <- as.matrix(dist(z))
#' decision <- GST.perm(dst, c(100, 200),2, 0.05, 1000)
#'
GST.perm <- function(dst, sizes, Q, alpha,B){
  T0=GST(dst,1:nrow(dst), sizes, Q)
  TPerm=rep(0,B)
  for (b in 1:B) {
    permindex=sample(1:nrow(dst))
    TPerm[b]=GST(dst,permindex,sizes, Q)
  }
  P <- mean(c(T0, TPerm) >= T0, na.rm=TRUE)
  return(P<alpha)
}



