# lingam - perform the LiNGAM analysis
#
# SYNTAX:
# res <- lingam( X, sources=c(), sinks=c() )
#
# INPUT:
# X        - Data matrix: each row is an observed variable, each
#            column one observed sample. The number of columns should
#            be far larger than the number of rows.
# sources  - Vector of indices of guaranteed sources (no parents)
# sinks    - Vector of indices of guaranteed sinks (no children)
#
# OUTPUT:
# res$B     - Matrix of estimated connection strenghts
# res$stde  - Standard deviations of disturbance variables
# res$ci    - Constants
# res$k     - An estimated causal ordering
#
# (Note that B, stde, and ci are all ordered according to the same
#  variable ordering as X.)
#
# VERSION: 
# 0.1 (8 Jan 2008) - first version for R
# 0.2 (23 Jun 2009) - possible to define sources and sinks information
#
# EXAMPLES:
# res <- lingam( X )
# res <- lingam( X, sources=c(3,7) )
# res <- lingam( X, sinks=c(1,2) )
# res <- lingam( X, sources=c(3,7), sinks=c(1,2) )
# 
# See also ESTIMATE, PRUNE.

lingam <- function( X, sources=c(), sinks=c() ) {

  temp <- estimate( X, sources, sinks )
  res <- prune( X, temp$k, sources, sinks )
  res$k<- temp$k
  res$Bnopruned<-temp$B #added oct 2020 by a.m.
  res  
}
