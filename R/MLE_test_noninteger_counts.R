


# Testing PLB likelihood (adapted for counts by Andy in github.com/andrew-edwards/size-spectra-man2/fitMLEmidMLEbin.R)

setwd("/Users/jpwrobinson/Documents/git-repos/reef-indicators2")
setwd("/Users/IMAC3/Documents/git-jpwrobinson/reef-indicators2")
rm(list=ls())


source('scripts/functions/countsFunctions.r')
source('scripts/functions/MLB_binned_function.r')
source('scripts/functions/PLBfunctions.r')

## simulate power law data, 
sim<-rPLB(n=10000, b=-2, xmin=16, xmax=10000)
## bin data to 1g (i.e. rounding masses to nearest gram, as for reef data)
sim<-round(sim, digits=0)



#------------------------------------------------------------
#------------ Now estimate exponent ------------------------#
#------------------------------------------------------------
# ## params
#   log.x = log(sim)                     
#   sum.log.x = sum( log.x ) 
#   xmin = min(sim)
#   xmax = max(sim)

# #-------------------------------------------------------------------------
# #------------ PLB for binned (1g), integer frequency data ------------------------#
# #-------------------------------------------------------------------------

#   # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
#   #  as a starting point for nlm for MLE of b for PLB model.
#   PL.bMLE = 1/( log(min(sim)) - sum.log.x/length(sim)) - 1
      
#   PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=sim, n=length(sim),
#       xmin=xmin, xmax=xmax, sumlogx=sum.log.x) #, print.level=2 )
  
#  PLB.bMLE = PLB.minLL$estimate
#    # 95% confidence intervals for MLE method.
  
#   PLB.minNegLL = PLB.minLL$minimum
  
#   # Values of b to test to obtain confidence interval. For the real movement data
#   #  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
#   #  symmetric interval here.
  
#   bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 0.001)  # If make 0.0001 then do 
#                               # get an interval for raw 1980 data 
      
#   PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
#   for(i in 1:length(bvec))
#       {
#           PLB.LLvals[i] = negLL.PLB(bvec[i], x=sim, n=length(sim), xmin=xmin,
#               xmax=xmax, sumlogx=sum.log.x)   
#       }
#   critVal = PLB.minNegLL  + qchisq(0.95,1)/2
#                       # 1 degree of freedom, Hilborn and Mangel (1997) p162.
#   bIn95 = bvec[ PLB.LLvals < critVal ]
#                       # b values in 95% confidence interval
#   PLB.MLE.bConf = c(min(bIn95), max(bIn95))
 
# PLB_freq<-data.frame(b=PLB.bMLE, conf.upper=PLB.MLE.bConf[2],conf.lower=PLB.MLE.bConf[1])

#-------------------------------------------------------------------------
#------------ PLB for binned (1g), integer count data ------------------------#
#-------------------------------------------------------------------------

## MLE for integer values of frequency mass data
PLB_x_integer<-PLB_binned(x = sim, binWidth=1,integer=TRUE)


## MLE for integer values of count ~ mass data
# ## aggregate masses to get counts
sim.counts<-data.frame(table(sim))
colnames(sim.counts)<-c('mass',  'count')
sim.counts$mass<-as.numeric(as.character(sim.counts$mass))
PLB_counts_integer<-PLB_binned(counts = sim.counts, binWidth=1,integer=TRUE)

#---------------------------------------------------------------------------------------------------
#------------ PLB for binned (1g), non-integer count data (divided by 2) ------------------------#
#---------------------------------------------------------------------------------------------------
sim.counts.2<-sim.counts 
sim.counts.2$count<-sim.counts.2$count/2

PLB_counts_noninteger_2<-PLB_binned(counts = sim.counts.2, binWidth=1, integer=FALSE)
PLB_counts_noninteger_2

# tt<-binData(counts=sim.counts.2, binWidth=1, integer=FALSE)
# tail(tt$binVals)

#---------------------------------------------------------------------------------------------------
#------------ PLB for binned (1g), non-integer count data (divided by 3) ------------------------#
#---------------------------------------------------------------------------------------------------
sim.counts.3<-sim.counts 
sim.counts.3$count<-sim.counts.3$count/3

PLB_counts_noninteger_3<-PLB_binned(counts = sim.counts.3, binWidth=1, integer=FALSE)
PLB_counts_noninteger_3

#---------------------------------------------------------------------------------------------------
#------------ PLB for binned (1g), non-integer count data (divided by 4) ------------------------#
#---------------------------------------------------------------------------------------------------
sim.counts.4<-sim.counts 
sim.counts.4$count<-sim.counts.4$count/4

PLB_counts_noninteger_4<-PLB_binned(counts = sim.counts.4, binWidth=1, integer=FALSE)
PLB_counts_noninteger_4

#---------------------------------------------------------------------------------------------------
#------------ PLB for binned (1g), non-integer count data (divided by 5) ------------------------#
#---------------------------------------------------------------------------------------------------
sim.counts.5<-sim.counts 
sim.counts.5$count<-sim.counts.5$count/5

PLB_counts_noninteger_5<-PLB_binned(counts = sim.counts.5, binWidth=1, integer=FALSE)
PLB_counts_noninteger_5

### visualise
plb_sens<-rbind(PLB_x_integer, PLB_counts_integer, PLB_counts_noninteger_2, PLB_counts_noninteger_3, PLB_counts_noninteger_4, PLB_counts_noninteger_5)
plb_sens$scale<-c('MLEbin_freq', 'MLEbin_counts', 'counts/2', 'counts/3', 'counts/4', 'counts/5')

theme_set(theme_bw())
pdf(file='figures/PLB_MLE_sizespec_noninteger_test.pdf', height=7, width=11)
ggplot(plb_sens, aes(scale, y= b, ymax=conf.upper, ymin=conf.lower)) + geom_errorbar() + geom_point()
dev.off()
