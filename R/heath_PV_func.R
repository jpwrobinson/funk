

## Adapting matlab code from Appendix 1 of Heath (2006) 
## Quantifying temporal variability in population abundances. Oikos

library(caTools)

pv<-function(P){

	# P = vector of population abundances

	# estimate all pairwise combinations of abundances
	C  <- combs(P, k = 2)

	# how many combinations?
	y  <- length(P)
	## empty diff vector
	Diff<-numeric()

	# now difference function - for each combination,
	# calculate proportional difference between each pair
		for(m in 1:y) {
			Num  <- abs((C[m,1])-(C[m,2]))
			Denom <- max((C[m,1]),(C[m,2]))
			# if abundances are equal, set Diff to 0 
			# (this ensures no NaN when two zero counts are compared)
			if(Num == Denom){ Diff[m] <- 0 }
			if(Num != Denom){ Diff[m] <- Num/Denom }
		 }

	# estimate PV
	PV  <- mean(Diff)

	return(PV)
}

