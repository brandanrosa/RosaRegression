#' myLAreg
#'
#' A function which takes an `lm()` object and produces a list of fun and interesting values which Regression covets
#'
#' @param obj an `lm()` object
#'
#' @return a named list with the following: model summary, n (num observations), k (num beta parameters), beta estimates, RSS, s-squared, s, and the Design & Response Matrices.
#' @export
#'
#' @examples
#' ylm <- lm(LENGTH ~ WEIGHT, data = ddt)
#' myLAreg(ylm)
myLAreg <- function(obj) { # obj = lm() object

  # summary
  sm <- summary(obj)

  # data frame
  dat <- stats::model.frame(obj)

  # n and k
  n <- dim(dat)[1]
  k <- dim(dat)[2]

  # matrices
  X <- stats::model.matrix(obj)
  Y <- dat[, 1]
  Y <- as.matrix(Y)

  # beta estimates
  betahat <- solve(t(X) %*% X) %*% t(X) %*% Y

  # RSS
  RSS <- t(Y) %*% Y - t(betahat) %*% t(X) %*% Y

  # ssq & s
  ssq <- RSS / (n - (k + 1))
  s <- sqrt(ssq)

  # List
  list(Summary=sm,
       n=n,
       k=k,
       beta.estimates=betahat,
       RSS=RSS,
       s.squared=ssq,
       s=s,
       Design.Matrix=X,
       Response.Matrix=Y,
       dat=dat)
}
