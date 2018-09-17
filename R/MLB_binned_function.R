
## PLB_binned_function.R

## Function for calculating ISD exponent from binned data with likelihood methods

## Dependencies: countsFunctions.R (https://github.com/andrew-edwards/size-spectra-man2/blob/master/code/countsFunctions.r)
# 			       	 PLBfunctions.R (https://github.com/andrew-edwards/fitting-size-spectra/blob/master/code/PLBfunctions.r)


PLB_binned<-function(x=NULL, counts=NULL, binWidth=1, size.var='mass', count.var='count',integer=TRUE,...)
	{

	if(is.null(x)){
		# subset data to counts and masses
		counts<-subset(counts, select=c(size.var, count.var))
		samp.bins<-binData(counts=counts, binWidth=binWidth, integer=integer)
	}

	if(is.null(counts)){
		# subset data to vector of individual values
		# x<-subset(x, select=c(size.var))
    counts<-x
		samp.bins<-binData(x=counts, binWidth=binWidth, integer=integer)
	}

## Call binData function from countsFunctions.r
		 # Args:
    #  x: vector of individual values (e.g. body masses).
    #   ORto
    #  counts: dataframe (or array) with first column being an x value
    #  (e.g. body mass), and second column being the counts of the
    #   number of individuals for that value.
    #   Only x or counts can be specified.
	#  binWidth = type of bin to use (SPC/TOW processing uses 1g rounded)  

	## now extract bin information for PLB likelihood
	    num.bins = dim(samp.bins$binVals)[1]
      binBreaks = samp.bins$binVals[,"binMin"]$binMin   # Loses the column names
      maxOfMaxBin = samp.bins$binVals[num.bins, "binMax"]$binMax
      binBreaks = c(binBreaks, maxOfMaxBin) # Append endpoint of final bin
    
      binCounts = samp.bins$binVals[,"binCount"]$binCount
      binMids = samp.bins$binVals[,"binMid"]$binMid  # Midpoints of bins
    

      sumCntLogMids = sum(binCounts * log(binMids))
    
      # MLEmid (maximum likelihood using midpoints) calculations.
    
      # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
      #  as a starting point for nlm for MLE of b for PLB model.
      PL.bMLE = 1/( log(min(binBreaks)) - sumCntLogMids/sum(binCounts)) - 1

      # MLEbin (maximum likelihood on binned data) calculations.
      # using calcLike from countsFunctions.r

      MLEbin.res = calcLike(negLL.PLB.binned, p=PL.bMLE,
          w = binBreaks, d = binCounts, J = length(binCounts))

      MLE.result<-data.frame(b = MLEbin.res$MLE, conf.upper = MLEbin.res$conf[2], conf.lower = MLEbin.res$conf[1])

      return(MLE.result)
  }