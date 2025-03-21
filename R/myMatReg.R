#' myMatReg
#'
#' A function which returns many great things from matrix algebra for regression
#'
#' @param obj an object from an `lm()`
#'
#' @return an invisible list to call including: summary, data frame, characteristic & response matrices, beta estimates, RSS, ssq, s, fitted values, residuals, leverage values, and a resiudal plot.
#' @export
#'
#' @importFrom ggplot2 ggplot geom_point geom_hline scale_color_manual labs aes
#'
#' @examples \dontrun{myMatReg(obj = ylm)}
myMatReg <- function(obj) { # obj = lm() object

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

  # fitted and residuals
  yhat <- X %*% betahat
  res <- Y - yhat

  # leverage
  P <- X %*% solve(t(X) %*% X) %*% t(X)
  leverage <- diag(P)

  # diagnostic residual plot
  dat2 <- data.frame(yhat, res)

  out <- c()

  dat2$out <- factor(
    ifelse(abs(res) > 3*as.numeric(s), "out",
           ifelse(abs(res) > 2*as.numeric(s) & abs(res) <= 3*as.numeric(s), "maybe", "in")))

  res.plot <- ggplot(data = dat2, aes(x = yhat, y = res, color = out)) +
    geom_point(cex = 3, aes(x = yhat, y = res)) +
    geom_hline(yintercept = c(0, 2*s, -2*s, 3*s, -3*s), color = "black", linetype = "dashed") +
    scale_color_manual(name = "Outliers & Friends",
                       values = stats::setNames(c("darkred", "navy", "darkgreen"), c("out", "maybe", "in")),
                       labels = c("Not Outlier", "Potential", "Outlier")) +
    labs(title = "Residuals ~ Fitted",
         subtitle = "with Outlier ID",
         x = "Fitted Values",
         y = "Residuals")

  # list
  invisible(list(summary=sm,
                 dat=dat,
                 X=X,
                 Y=Y,
                 betahat=betahat,
                 RSS=RSS,
                 ssq=ssq,
                 s=s,
                 fit=yhat,
                 resid=res,
                 leverage=leverage,
                 res.plot=res.plot,
                 n=n,
                 k=k))
}
