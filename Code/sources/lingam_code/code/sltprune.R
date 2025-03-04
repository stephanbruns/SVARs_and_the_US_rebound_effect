sltprune <- function( B ) {

  # Finds an identical row and column permutation so that the permuted
  # B is approximately strictly lower triangular. This algo also works
  # in high dimensions!

  # The basic idea is the following: First, set to zero the n(n+1) smallest
  # (in absolute value) entries of B, test if we can permute to strict
  # lower triangularity. If not, then iteratively set to zero the next
  # smallest entry and test again. Continue this process until we succeed.

  # Dimensionality of problem
  n <- dim(B)[1]

  # First, order the indexes in terms of absolute value
  # (Here, a better strategy would be to use some sort of p-values,
  # but this first implementation will use this simple approach.)
  temp <- sort(abs(B),index.return=TRUE);
  y <- temp$x
  ind <- temp$ix

  for (i in (n*(n+1)/2):(n*n)) {

    # Copy original B into Bi
    Bi <- B

    # Set 'i' smallest (in absolute value) coefficients to zero
    Bi[ind[1:i]] <- 0

    # Try to do permutation
    p <- slttestperm( Bi )

    # If we succeeded, then we're done!
    if (any(p!=0)) {
      res <- list()
      res$Bopt <- B[p,p]
      res$optperm <- p
      return(res)
    }

    # ...else we continue, setting one more to zero!
  }

  res

}