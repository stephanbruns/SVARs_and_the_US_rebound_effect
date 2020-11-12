nzdiaghungarian <- function( W ) {

  # ---------------------------------------------------------------------------
  # Find best row permutation by hungarian algorithm
  # [Uses function solve_LSAP (solves linear sum assignment problem using the
  # hungarian algoritm) from R-package clue]
  # ---------------------------------------------------------------------------

  library('clue')

  n <- dim(W)[1]

  # Create score matrix
  S <- 1./abs(W)

  # Call the hungarian algorithm
  temp <- solve_LSAP(t(S),maximum=FALSE)
  rowp <- unclass(temp) # this is the optimal permuatation
  dim(rowp) <- c(1,n)

  # Permute W to get Wopt
  Wopt <- W[rowp,]

  res <- list()
  res$Wopt <- Wopt
  res$rowp <- rowp

  res

}