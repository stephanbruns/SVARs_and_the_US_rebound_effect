# PRUNE Prune weight matrix B.
#
# [Bpruned stde ci] = prune(X, k, sources, sinks)
#     Prunes and re-estimates the weight matrix using a simple 
#     resampling method. Also re-estimates the standard deviations of error
#     terms and constant terms.
#
#
# Input arguments:
#     X        - the original data matrix used in model estimation.
#     k        - the estimated causal ordering of the variables.
#     sources  - vector of indices guaranteed to have no parents (optional) 
#     sinks    - vector of indices guaranteed to have no children (optional) 
#
# Output:
#     res$Bpruned  - the pruned matrix of connection strenghts.
#     res$stde     - the re-estimated standard deviations of error terms.
#     res$ci       - the re-estimated constant terms.
#
# See also ESTIMATE, LINGAM.

prune_mute <- function( X, k, sources=c(), sinks=c() ) {

  # ---------------------------------------------------------------------------
  # Default values for parameters
  # ---------------------------------------------------------------------------

  method <- 'resampling'  # the pruning method

  # 'prunefactor' determines how easily weak connections are pruned
  # in the simple resampling based pruning. Zero gives no pruning 
  # whatsoever. We use a default of 1, but remember that it is quite arbitrary
  # (similar to setting a significance threshold). It should probably
  # become an input parameter along with the data, but this is left
  # to some future version.
  #prunefactor <- 1
  prunefactor <- 0

  # ---------------------------------------------------------------------------
  # Pruning
  # ---------------------------------------------------------------------------

 # cat('Pruning the network connections...\n')
  dims <- nrow(X)
  ndata <- ncol(X)

  if (method == 'resampling') {

    # -------------------------------------------------------------------------
    # Pruning based on resampling: divide the data into several equally
    # sized pieces, calculate B using covariance and QR for each piece
    # (using the estimated causal order), and then use these multiple
    # estimates of B to determine the mean and variance of each element.
    # Prune the network using these.
  
    npieces <- 10
    piecesize <- floor(ndata/npieces)

    Bpieces <- array(0,c(dims,dims,npieces))
    diststdpieces <- array(0,c(dims,npieces))
    cpieces <- array(0,c(dims,npieces))
    Bfinal <- matrix(0,dims,dims)
    
    for (i in 1:npieces) {
    
      # Select the appropriate subset of data
      Xs <- X[,((i-1)*piecesize+1):(i*piecesize)]
    
      # Step through all of the variables in the given ordering
      for (j in 1:length(k)) {
      
        # Current variable number
        v <- k[j]
      
        #-------- Find the set of possible parents -------
      
        # All previous ones are possible parents...
        if (j>1) parents <- k[1:(j-1)] else parents <- c()
        
        # ...but not any previous sinks...
        parents <- setdiff(parents,sinks)
        
        # ...and if the current node is a source it has no parents!
        if (is.element(v,sources)) parents <- c()
        
        # Number of parents
        np <- length(parents)
        
        #----- Regress current variable on set of possible parents ----

		# Xp holds the parents, Xv the current variable
        if (np>0) Xp <- Xs[parents,]
        Xv <- Xs[v,]

        # Make sure we have appropriate dimensions 
        if (np>0) {
          dim(Xp) <- c(np,piecesize)
          dim(Xv) <- c(1,piecesize)
        }
        
        # Take out means
        if (np>1) Xpzm <- Xp - rowMeans(Xp)
        else if (np>0) Xpzm <- Xp - mean(Xp)
        Xvzm <- Xv - mean(Xv)
                
        # Compute regression coefficient(s)
        if (np>0) {
          r <- (Xvzm %*% t(Xpzm)) %*% solve(Xpzm %*% t(Xpzm))
  		}      
  		
  		# Compute residual
  		if (np>0) E <- Xv - (r %*% Xp)
  		else E <- Xv
  		
  		# Compute mean and std of residual
  		mE <- mean(E)
  		sdE <- sqrt(mean((E-mE)^2))
  		
        #----- Fill into appropriate results matrices ----

		if (np>0) Bpieces[v,parents,i] <- r
    	diststdpieces[v,i] <- sdE
    	cpieces[v,i] <- mE
  		
      }

    }

	# Pruning matrix B
    for (i in 1:dims) {
      for (j in 1:dims) {
      
        themean <- mean(Bpieces[i,j,])
        thestd <- sd(Bpieces[i,j,])
        if (abs(themean)<prunefactor*thestd) {	    
          Bfinal[i,j] <- 0
        }
        else {
          Bfinal[i,j] <- themean
        }
        
      }
    }

    # Estimating std and c
    diststdfinal <- rowMeans(diststdpieces)
    cfinal <- rowMeans(cpieces)

    # Finally, rename all the variables to the way we defined them
    # in the function definition  
    Bpruned <- Bfinal
    stde <- diststdfinal
    ci <- cfinal

  }
  
  if (method == 'olsboot') {
    stop('Not implemented yet!')
  }
  
  if (method == 'wald') {
    stop('Not implemented yet!')
  }

  if (method == 'bonferroni') {
    stop('Not implemented yet!')
  }

  if (method == 'hochberg') {
    stop('Not implemented yet!')
  }

  if (method == 'modelfit') {
    stop('Not implemented yet!')
  }

  #cat('Done!\n')

  # Return the result
  res <- list()
  res$Bpruned <- Bpruned
  res$stde <- stde
  res$ci <- ci
  res

}

