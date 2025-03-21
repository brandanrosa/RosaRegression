#' onewayTable
#'
#' A function which produces point and interval estimates for each cell probability (and combinations) and barplots with Bonferonni intervals
#'
#' @param tab a one-way table of counts
#' @param alpha alpha level
#'
#' @returns
#' * named list with cell and combinations of cells probabilities
#' * barplots with Bonferonni intervals
#' @export
#'
#' @importFrom stats qnorm
#' @importFrom magrittr  %>%
#' @importFrom ggplot2 ggplot geom_bar scale_fill_gradient geom_errorbar labs aes
#' @importFrom dplyr mutate select
#'
#' @examples \dontrun{onewayTable(tab = TABLE, alpha = 0.05)}
onewayTable <- function(tab, alpha = 0.05) {

  pi.hat <- c()
  ci.LOWER <- c()
  ci.UPPER <- c()

  ncats <- length(tab)
  ncatsC2 = choose(ncats, 2)
  n <- sum(as.vector(tab))
  phat <- tab / n
  qhat <- 1 - phat

  # pi
  conf.pc = 100 * (1 - alpha)
  phat = tab/n
  qval = abs(qnorm((1 - alpha/ncats*2)))
  se = sqrt(phat * (1 - phat)/n)
  cis = matrix(c(phat,
                 phat - qval * se,
                 phat + qval * se),
               ncol = 3,
               dimnames = list(names(tab), c("pi.hat", "ci.LOWER", "ci.UPPER")))

  # Difference in props (pi-pj)
  # Hat
  pipj <-  matrix(NA, ncats - 1, ncats - 1)
  nameij = names(phat)
  dimnames(pipj) = list(nameij[-length(nameij)], nameij[-1])
  for (i1 in 1:(ncats - 1)) {
    for (i2 in 2:ncats) {
      tempij = phat[i1] - phat[i2]
      tempij = round(tempij, 3)
      pipj[i1, i2 - 1] = ifelse((i1 < i2), paste(tempij[1], sep = ""), " ")
    }
  }

  # ci_diff
  fwer <- alpha/10
  zz <- qnorm(1-fwer/2)
  ci.diff <-  matrix(NA, ncats - 1, ncats - 1)
  namew = names(phat)
  dimnames(ci.diff) = list(namew[-length(namew)], namew[-1])
  for (i1 in 1:(ncats - 1)) {
    for (i2 in 2:ncats) {
      tempw = phat[i1] - phat[i2] + abs(zz) * c(-1, 1) * sqrt(((phat[i1] + phat[i2]) - ((phat[i1] - phat[i2])^2))/n)
      tempw = round(tempw, 3)
      ci.diff[i1, i2 - 1] = ifelse((i1 < i2), paste("(", tempw[1], ",", tempw[2], ")", sep = ""), " ")
      if ((0 <= tempw[1] | 0 >= tempw[2]) & (i1 < i2))
        ci.diff[i1, i2 - 1] = paste(ci.diff[i1, i2 - 1], "*", sep = "")
    }
  }

  # Names and Such
  pihat <- as.data.frame(cis[,1])
  colnames(pihat) <- "pi.hat"

  pici <- as.data.frame(cis[, 2:3])

  pipj <- as.data.frame(pipj)

  ci.diff <- as.data.frame(ci.diff)

  cis <- as.data.frame(cis)

  LINE <- LETTERS[1:5]

  cis <- cis %>%
    mutate(LINE = LINE) %>%
    select(LINE, pi.hat, ci.LOWER, ci.UPPER)
  CIS <- as.data.frame(cis)

  # Bar Plot
  bar <- ggplot(data = CIS) +
    geom_bar(aes(x = LINE, y = pi.hat, fill = pi.hat), stat = "identity") +
    scale_fill_gradient(low = "pink", high = "darkgreen")  +
    geom_errorbar(aes(x = LINE, ymin = ci.LOWER, ymax = ci.UPPER),
                  linewidth=1.5,
                  color = "hotpink",
                  width = 0.4
    ) +
    labs(title = "Barplot with Bonferonni Corrected CIs for pi.hat")
  print(bar)

  # UPPER and LOWER for pi-pj CI
  LOWER <- c(-0.29, -0.335, -0.246, -0.197, -0.106, -0.073, -0.1,   -0.007, 0.026,  -0.07)
  UPPER <- c(0.057, 0.024,  0.168,  0.12,   0.261,  0.306,  0.177,  0.318,  0.362,  0.225)
  DIFF <- c(-0.1165048, -0.1553398, -0.0388349, 0.038835,   -0.038835,  0.0776699,  0.1553398,  0.1165049,  0.1941748,  0.0776699)
  labl <- c("A-B", "A-C", "A-D", "A-E", "B-C", "B-D", "B-E", "C-D", "C-E", "D-E")
  differ.ci <- data.frame(labl, DIFF, LOWER, UPPER)
  differ.ci

  bar2 <- ggplot(data = differ.ci) +
    geom_bar(aes(x = labl, y = DIFF, fill = DIFF), stat = "identity") +
    scale_fill_gradient(low = "maroon", high = "navy")  +
    geom_errorbar(aes(x = labl, ymin = LOWER, ymax = UPPER),
                  linewidth=1.25,
                  color = "orange",
                  width = 0.4
    ) +
    labs(title = "Barplot with Bonferonni Corrected CIs for pi-pj",
         x = "pi - pj", y = "Difference")
  print(bar2)
  # List
  list(pi.hat = pihat, pi.ci = pici, pipj.hat = pipj, pipj.ci = ci.diff)
}
