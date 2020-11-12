Gauss_Tests <- function(X) {

  library(nortest)
  library(tseries)

  dims <- dim(X)

  if (dims[1]<dims[2]) {
    X <- t(X)
    dims <- dim(X)
  }

  m <- dims[1]
  if (m > 5000) {
    cat("get subsample of residuals to perform normality-test\n")
    subsamp <- seq(1,m,by=ceiling(m/5000))
    X1 <- array(0,dim=c(length(subsamp),dims[2]))
    for (i in 1:dims[2]) {
      temp <- sort(X[,i])
      X1[,i] <- temp[subsamp]
    }
  }
  else X1 <- X

  shapiro <- list(); shapfranc <- list(); jarque <- list();
  lis1 <- NULL; lis2 <- NULL; lis3 <- NULL
  pval1 <- NULL; pval2 <- NULL; pval3 <- NULL

  for (i in 1:dims[2]) {
    # Shapiro-Wilk test
    shapiro[[i]] <- shapiro.test(X1[,i])
#    if (unclass(shapiro[[i]])$p.value < alpha) {
#      lis1 <- c(lis1,i)
      pval1 <- c(pval1, unclass(shapiro[[i]])$p.value)
#    }
    # Shapiro-Francia test
    shapfranc[[i]] <- sf.test(X1[,i])
#    if (unclass(shapfranc[[i]])$p.value < alpha) {
#      lis2 <- c(lis2,i)
      pval2 <- c(pval2, unclass(shapfranc[[i]])$p.value)
#    }
    # Jarque-Bera test
    jarque[[i]] <- jarque.bera.test(X1[,i])
#    if (unclass(jarque[[i]])$p.value < alpha) {
#      lis3 <- c(lis3,i)
      pval3 <- c(pval3, unclass(jarque[[i]])$p.value)
#    }

  }

  ls <- list()
  ls$ShapiroWilk <- pval1
  ls$ShapiroFrancia <- pval2
  ls$JarqueBera <- pval3

  ls

#  cat("Done! \n")
#  cat("- Shapiro-Wilk test: reject hypothesis of normality at level", alpha,
#      "for", length(lis1), "residuals out of", dims[1],", namely", lis1, "with
#      p-values", pval1, "respecitvely. \n")
#  cat("- Shapiro-Francia test: reject hypothesis of normality at level",
#      alpha, "for", length(lis2), "residuals out of", dims[1],", namely", lis2,
#      "with p-values", pval2, "respecitvely. \n")
#  cat("- Jarque-Bera test: reject hypothesis of normality at level",
#      alpha, "for", length(lis3), "residuals out of", dims[1],", namely", lis3,
#      "with p-values", pval3, "respecitvely. \n")

}