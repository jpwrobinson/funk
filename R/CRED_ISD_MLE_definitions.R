#!bin/env Rscript

### R function definitions ###
# Negative log-likelihood function for power law distribution over
#  infinite range (PL model). From equation (5) in Edwards et al.
#  (2007), Nature 449:1044-1048. Note here using the negative
#  log-likelihood to be minimised. [The 'raw' in the function name
#  corresponds to using the raw data - earlier studies have required
#  explicit likelihood functions that took into account the data
#  collection methods (e.g. albatross data in 2007 paper), or
#  analysing data that were already binned and not available 'raw'
#  (e.g. bees data in 2007 paper)].

powlaw = function(mu)
  {     
 -n_obs * log(mu - 1) - n_obs *(mu - 1) * log(a) + mu * sumlogx
  }

# Negative log-likelihood function for exponential distribution
#  over infinite range (Exp model). From equation (6) in Edwards
#  et al. (2007). 
explaw = function(lambda)   
  {    # Negative log-likelihood is just:                         
     - n_obs * log(lambda) - n_obs * lambda * a + lambda * sumx
  }


# Negative log-likelihood function for bounded power law
#  distribution (PLB model). From equation (A.23) in Edwards (2011),
#  Ecology 92:1247-1257.


powlawbound = function(mu)
 

  {              #  Negative log-likelihood is:
  if(mu == 1)    # See equation (A.25) in Ecology paper.
    { n_obs * log( log( b/a)) + sumlogx}
  else
    {  -n_obs * log( (mu - 1) / (a^(1 - mu) - b^(1 - mu))) + mu * sumlogx}
                 # This works for mu<1 as well as mu>1
  }

# Negative log-likelihood function for bounded exponential
#  distribution (ExpB model). From equation (A.27) in
#  Edwards (2011).
explawbound = function(lambda)
  {               # Negative log-likelihood is just:  
     - n_obs * log(abs(lambda)) + n_obs * log( abs( exp(- lambda * a) - exp( -lambda * b) ) ) + lambda * sumx    # abs allows lambda < 0
  }

### End R function definitions ###



### R parameter definitions ###




numparinf = 2               # if a is min of data then
#  estimating a and mu/lambda for
#  unbounded models; numparinf is
#  number of parameters estimated
#  for unbounded (inf) models.
#  Else if prescribing a then just
#  mu is estimated.

numpar = numparinf + 1           # num of pars for bounded models
# bounded model has 1 extra parameter


#bround = TRUE                    # Whether to round b up. If TRUE 
#  then b=ceiling(max(rawvalsall))power
#  FALSE b= max(rawvalsall) [below].


# These are for numerical optimizations.


# Similarly, these are vectors for generality.
#  They are the starting value for optimisation
#  and the start, end, and increments to use
#  to generate values to calculate 95%
#  confidence intervals.
mu.startinfvec = vector()        # inf corresponds to unbounded 
lambda.startinfvec = vector()    #  (infinite) models.
muvecinfstartvec = vector()
muvecinfendvec   = vector()
muvecinfincvec   = vector()
muvecstartvec    = vector()
muvecendvec      = vector()
muvecincvec      = vector()
lambdaincvec     = vector()
mu.start2vec     = vector()
lambda.start2vec = vector()
lambdainfstartvec=vector()
lambdastartvec=vector()
lambdainfendvec=vector()
lambdaendvec=vector()
lambdainfincvec=vector()
lambdaincvec=vector()



B0.startvec<-vector()
B1.startvec<-vector()
B2.startvec<-vector()
B3.startvec<-vector()
B4.startvec<-vector()
B5.startvec<-vector()



# Initial starting values for finding MLE using nlm. Setting to
#  close to what it eventually finds to avoid NA/Inf warnings, which
#  can occur when initial gradient is too steep.

mu.startinfvec[1] = c(1.001)
lambda.startinfvec[1] = c(0.01)
mu.start2vec[1] = c(0.5) 
lambda.start2vec[1] = c(0.01)

#B0.startvec<-c(-0.1)
#B1.startvec<-c(-0.1)
#B2.startvec<-c(-0.1)
#B3.startvec<-c(-0.1)
#B4.startvec<-c(-0.1)
#B5.startvec<-c(-0.1)


#par<-c(B0.startvec, B1.startvec, B2.startvec, B3.startvec, B4.startvec, B5.startvec)

# To define range of mu to try for inf, must be >1
muvecinfstartvec[1] = c(1.001)
muvecinfendvec[1] = c(2)
muvecinfincvec[1] = c(0.001)

# for mu bounded.
muvecstartvec[1]  = c(0.5)
muvecendvec[1]  = c(2)
muvecincvec[1]  = c(0.001)

# increment for lambda (0.5 and 2* mle)
lambdainfstartvec[1]  = c(0.0001)
lambdainfendvec[1]  = c(0.5)
lambdainfincvec[1] = c(0.001)
# increment for lambda (0.5 and 2* mle)
lambdastartvec[1]  = c(0.0001)
lambdaendvec[1]  = c(0.5)
lambdaincvec[1] = c(0.001)



# Likelihood calculations.
#  Range of mu to try for PL model, can't be <=1
muvecinf = seq(muvecinfstartvec[1], muvecinfendvec[1], muvecinfincvec[1])                    

#  Range of mu to try for PLB model, can't be <=1
muvec = seq(muvecstartvec[1], muvecendvec[1], muvecincvec[1])

lambdavecinf = seq(lambdainfstartvec[1], lambdainfendvec[1], lambdainfincvec[1])
lambdavec = seq(lambdastartvec[1], lambdaendvec[1], lambdaincvec[1])


### End of R parameter definitions ###

