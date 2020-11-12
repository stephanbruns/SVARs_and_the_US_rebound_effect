# ESTIMATE Estimate LiNGAM model from data.
#
# Estimates a LiNGAM model from data. The returned weight matrix B
# is strictly lower triangular, but the weights are not pruned. To
# prune the matrix, call PRUNE. (If sources and sinks are specified,
# then these appropriately have no parents or children, respectively,
# but otherwise the network is fully connected.)
#
# SYNTAX:
# res <- estimate(X, sources=c(), sinks=c())
#
# INPUT:
# X        - Data matrix: each row is an observed variable, each
#            column one observed sample. The number of columns should
#            be far larger than the number of rows.
# sources  - vector of indices guaranteed to have no parents (optional) 
# sinks    - vector of indices guaranteed to have no children (optional) 
# 
# OUTPUT:
# B     - Matrix of estimated connection strenghts
# stde  - Standard deviations of disturbance variables
# ci    - Constants
# k     - An estimated causal ordering
#
# (Note that B, stde, and ci are all ordered according to the same 
# variable ordering as X.)
#
# See also PRUNE, LINGAM.

estimate_mute <- function( X, sources=c(), sinks=c() ) {

  # Using the fastICA R package so make sure it is loaded
  library('fastICA')
  
  # Store number of samples
  nsamples <- dim(X)[2]
    
  # If we have sources and/or sinks, take them out of the 
  # analysis completely at this stage. Note: if less than two
  # variables left after taking out sources or sinks then there
  # really is no reason to run ICA or LiNGAM at all, so at 
  # least for now we do not gracefully handle this case, sorry!
  if ((length(sources)>0) || (length(sinks)>0)) {
  	Xorig <- X
  	Xorigzm <- Xorig - rowMeans(Xorig)
  	X <- X[-c(sources,sinks),]
  	if (nrow(X)<=1) { stop('No ICA left to do!') }
  }

  # Remove means (for easily computing regression coefficients)
  Xzm <- X-rowMeans(X)
  
  # If we have sources, regress the remaining variables on them,
  # then continue with residuals only
  if (length(sources)>0) {
    Xs <- Xorig[sources,]
    if (length(sources)==1) { Xszm <- Xs-mean(Xs) }
    else { Xszm <- Xs-rowMeans(Xs) }
    dim(Xszm) <- c(length(sources),nsamples)
  	Rs <- Xzm %*% t(Xszm) %*% solve(Xszm %*% t(Xszm))
  	X <- X - Rs %*% Xs
  }
  
  # Call the fastICA algorithm
  icares <- fastICA( t(X), nrow(X), tol=1e-14 ) # Doris: tol=1e-14
  W <- t((icares$K) %*% (icares$W))  
  
  # [Here, we really should perform some tests to see if the 'icasig' 
  # really are independent. If they are not very independent, we should
  # issue a warning. This is not yet implemented.]

  # Try to permute the rows of W so that sum(1./abs(diag(W))) is minimized
  #cat('Performing row permutation...\n');
  dims <- nrow(X)
  if (dims <= 8) {  
   # cat('(Small dimensionality, using brute-force method.)\n')
    temp <- nzdiagbruteforce( W )
    Wp <- temp$Wopt
    rowp <- temp$rowp
  }
  else {
    cat('(Using the Hungarian algorithm.)\n')
#    stop('Not implemented yet, sorry!')
    temp <- nzdiaghungarian( W )
    Wp <- temp$Wopt
    rowp <- temp$rowp
  }
  #cat('Done!\n')

  # Divide each row of Wp by the diagonal element
  estdisturbancestd <- 1/diag(abs(Wp))
  Wp <- Wp/diag(Wp)

  # Compute corresponding B
  Best <- diag(dims)-Wp

  # Estimate the constants c
  m <- rowMeans(X)
  dim(m) <- c(dims,1)
  cest <- Wp %*% m

  # Next, identically permute the rows and columns of B so as to get an
  # approximately strictly lower triangular matrix
  #cat('Performing permutation for causal order...\n');

  if (dims <= 8) {  
   # cat('(Small dimensionality, using brute-force method.)\n');
    temp <- sltbruteforce( Best )
    Bestcausal <- temp$Bopt
    causalperm <- temp$optperm
  }
  else {
    cat('(Using pruning algorithm.)\n')
    temp <- sltprune( Best )
    Bestcausal <- temp$Bopt
    causalperm <- temp$optperm
  }
  #cat('Done!\n');

#  print(Bestcausal)
  
  # Here, we report how lower triangular the result was, and in 
  # particular we issue a warning if it was not so good!
  percentinupper <- sltscore(Bestcausal)/sum(Bestcausal^2)
  if (percentinupper>0.2) 
    cat('WARNING: Causal B not really triangular at all!!\n')
  else if (percentinupper>0.05)
    cat('WARNING: Causal B only somewhat triangular!\n')
  else
   # cat('Causal B nicely triangular. No problems to report here.\n')

  # Set the upper triangular to zero
  Bestcausal[upper.tri(Bestcausal,diag=FALSE)] <- 0

  # Finally, permute 'Bestcausal' back to the original variable
  # ordering and rename all the variables to the way we defined them
  # in the function definition
  icausal <- iperm(causalperm);
  res <- list()
  res$B <- Bestcausal[icausal, icausal];
  res$stde <- estdisturbancestd
  res$ci <- cest
  res$k <- causalperm
  
  # If we have sources and/or sinks, we need to add them in
  if ((length(sources)>0) || (length(sinks)>0)) {
	fullB <- matrix(0,nrow(Xorig),nrow(Xorig))
	fullB[-c(sources,sinks),-c(sources,sinks)] <- res$B
	if (length(sources)>0) {
		fullB[-c(sources,sinks),sources] <- 
			(diag(nrow(res$B))-res$B) %*% Rs
	}
	if (length(sinks)>0) {		
		Rsi <- Xorigzm[sinks,] %*% t(Xorigzm[-sinks,]) %*% 
			solve(Xorigzm[-sinks,] %*% t(Xorigzm[-sinks,]))
		fullB[sinks,-sinks] <- Rsi
	}
	res$B <- fullB
	fullstde <- rep(0,nrow(Xorig))
	if (length(sources)>0) {
		if (length(sources)==1) {
			fullstde[sources] <- sqrt(mean(Xorigzm[sources,]^2))
		} else {
    		fullstde[sources] <- sqrt(rowMeans(Xorigzm[sources,]^2))
		}
	}
	if (length(sinks)>0) {		
		if (length(sinks)==1) {
			fullstde[sinks] <- sqrt(mean((
				Xorigzm[sinks,] - (Rsi %*% Xorigzm[-sinks,]))^2))		
		} else {
			fullstde[sinks] <- sqrt(rowMeans((
				Xorigzm[sinks,] - (Rsi %*% Xorigzm[-sinks,]))^2))
		}
	}
	fullstde[-c(sources,sinks)] <- res$stde
	res$stde <- fullstde
	fullci <- rep(0,nrow(Xorig))
	fullci[-c(sources,sinks)] <- res$ci
	if (length(sources)>0) {
		if (length(sources)==1) {
			fullci[sources] <- mean(Xorig[sources,])
		} else {
    		fullci[sources] <- rowMeans(Xorig[sources,])
		}
	}
	if (length(sinks)>0) {		
		if (length(sinks)==1) {
			fullci[sinks] <- mean(
				Xorig[sinks,] - (Rsi %*% Xorig[-sinks,]))		
		} else {
			fullci[sinks] <- rowMeans(
				Xorig[sinks,] - (Rsi %*% Xorig[-sinks,]))
		}
	}	
	res$ci <- fullci
	origind <- 1:nrow(Xorig)
	origind <- origind[-c(sources,sinks)]
	res$k <- c(sources,origind[res$k],sinks)
  }

  # Return the result
  res

}
