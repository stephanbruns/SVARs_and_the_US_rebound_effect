slttestperm <- function( B ) {

  # slttestperm - tests if we can permute B to strict lower triangularity
  # If we can, then we return the permutation in p, otherwise p=0.

  # Dimensionality of the problem
  n <- dim(B)[1]

  # This will hold the permutation
  p <- NULL;

  # Remaining nodes
  remnodes <- 1:n

  # Remaining B, take absolute value now for convenience
  Brem = abs(B)

  # Select nodes one-by-one
  for (i in 1:n) {
    # Find the row with all zeros
    if (!is.null(dim(Brem))) {
      therow <- which(rowSums(Brem)<1e-12)
    }
    else {
      if (Brem<1e-12) therow <- 1
      else therow <- integer(0)
    }

    # If empty, return 0
    if (length(therow)==0) {
      p <- 0
      return(p)
    }

    # If more than one, arbitrarily select the first 
    therow = therow[1]

    # If we made it to the end, then great!
    if (i==n) {
      p <- c(p,remnodes)
      return(p)
    }

    # Take out that row and that column
    inds <- which((1:(n-i+1)) != therow)
    Brem <- Brem[inds,inds]
    if (is.null(dim(Brem))) dim(Brem) <- c(1,1)

    # Update remaining nodes
    p <- c(p,remnodes[therow])
    remnodes <- remnodes[inds]

  }

  p

}
