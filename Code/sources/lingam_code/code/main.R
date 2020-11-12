
# STARTING UP:
#
# > source('loadmain.R')
# > options(error=recover)
# > install.packages('fastICA') # only needed if not already installed
# > install.packages('clue') # only needed if data has more than 8
#        variables (for hungarian algorithm) and if not already installed
# > main(1)

# this gets called by main() from prompt
runmain <- function( task, theseed=-1 ) {

  cat('--------------------------------------------------------\n')

  # Randomly/manually set seed for reproducibility
  # (omit 'theseed' for random seed; specify 'theseed' to start
  # from a specific seed)
  
  if (theseed != -1) {
    set.seed( theseed )
    cat('Using seed: ',theseed,'\n')
  }
  else {
    theseed <- floor(runif(1,0,1)*100000)
    cat('Using seed: ',theseed,'\n')
    set.seed( theseed )
  }

  # Use the fastICA R package
  library('fastICA')

  #------------------------------------------------------------------------
  # Task 1: Simple three-variable network (a chain) example
  #------------------------------------------------------------------------

  if (task==1) {

    x <- matrix(rnorm(10000)^3,1,10000)
    y <- 0.3*x + matrix(rnorm(10000)^3,1,10000)
    z <- -2*y + matrix(rnorm(10000)^3,1,10000)
    D <- rbind(z,y,x)
    res <- lingam(D)    
    print(res)
    plotgraph(res$Bpruned)
    
  } 
  
  #------------------------------------------------------------------------
  # Task 2: randomly choose DAG and params, simulate data, and estimate
  #------------------------------------------------------------------------
    
  if (task==2) {

    testlingam()
    
  }

  #------------------------------------------------------------------------
  # Task 3: Test the sources/sinks restrictions in a trivial example
  #------------------------------------------------------------------------

  if (task==3) {

    rsuperg <- function(N) { 
    	a <- rnorm(N)^3
    	a <- a - mean(a)
    	a <- a/sd(a)
    	a <- matrix(a,1,N)
		a
    }

    N <- 10000
    x <- rsuperg(N) + 1
    y <- x + rsuperg(N) + 3
    z <- y + rsuperg(N) + 3
    u <- z + rsuperg(N) + 4
    v <- z + 0.1*u+ rsuperg(N) + 5
    D <- rbind(x,y,z,u,v)
    sources <- c()
    sinks <- c(3,4,5)
    res <- lingam(D, sources, sinks)    
    print(res)
    plotgraph(res$Bpruned)
    
  } 

}

