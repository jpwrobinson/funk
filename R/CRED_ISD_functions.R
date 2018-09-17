
theme_set(theme_bw())

### function to run ANY ISLAND
dataorder<-function(x, catch, size) data.frame(rep(x[,size], times=x[,catch] ))

mle_params<-function(df, size, catch){
		#df$mass[df$mass==0]<-0.001
	
df<-dataorder(df, size=size, catch=catch)
rawvalsall<-df[,1]	
x <<-df[,1]
a <<- min(rawvalsall)
b <<- max(rawvalsall)
n_obs <<- length(rawvalsall)

sumlogx <<- sum(log(rawvalsall)) # for powlaw function
sumx<<-sum(rawvalsall)
}


pl_mlefunc<-function(df, size, catch, dataorder=TRUE) {



if(dataorder==TRUE){
	mle_params(df, size, catch)
}

if(dataorder==FALSE){
	rawvalsall<-df[,size]	
x <<-df[,size]
a <<- min(rawvalsall)
b <<- max(rawvalsall)
n_obs <<- length(rawvalsall)

sumlogx <<- sum(log(rawvalsall)) # for powlaw function
sumx<<-sum(rawvalsall)
}

# Analytically (Box 1 of Edwards et al. 2007):
mumleinfanal<- 1/( - log(a) + sumlogx/n_obs) + 1
                        # maximum likelihood estimate for mu

# Numerically:
outpowlawinf<- nlm(powlaw, mu.startinfvec[1])  
                         # , print.level=2)
mumleinf<- outpowlawinf$estimate
negloglikmlepowlawinf<- outpowlawinf$minimum
                         # negative log likelihood

# Check that analytical and numerical agree:
#if( abs(mumleinf - mumleinfanal) >10^(-5))
#  {stop("Analytical and numerical disagree for PL model")}

# Normalisation constant:
Cpowlawinf<- (mumleinf - 1) * a^(mumleinf - 1)

AICpowlawinf<- 2 * negloglikmlepowlawinf + 2*numparinf 
                           # will be +4 if est mu and a



# Now for 95% CI of MLE:
muvarynegloglikinf<- vector()  # negative log lik for each mu value
  for(i in 1:length(muvecinf))
    {
      muvarynegloglikinf[i]<- powlaw(muvecinf[i])
    }
critval1inf<- negloglikmlepowlawinf + qchisq(0.95,1)/2  
                  # 1 dof. Hilborn and Mangel (1997) p163.
muin95inf<- muvecinf[ muvarynegloglikinf < critval1inf] 
                  # mu values definitely in 95% confidence interval
minmuin95inf<- min(muin95inf)
maxmuin95inf<- max(muin95inf)


coefs<<-data.frame(mumleinfanal, mumleinf, Cpowlawinf, 
	minmuin95inf, maxmuin95inf, AICpowlawinf)
colnames(coefs)<<-c("exponent_analytical", "exponent_num",
 "constant", "CImin", "CImax", "AIC")

coefs
}




plb_mlefunc<- function(df, size, catch, dataorder=TRUE) {
	
if(dataorder==TRUE){
	mle_params(df, size, catch)
}

if(dataorder==FALSE){
	rawvalsall<-df[,size]	
x <<-df[,size]
a <<- min(rawvalsall)
b <<- max(rawvalsall)
n_obs <<- length(rawvalsall)

sumlogx <<- sum(log(rawvalsall)) # for powlaw function
sumx<<-sum(rawvalsall)
}

outpowlaw<- nlm(powlawbound, mu.start2vec[1]) 
                  #, print.level=2)
mumle<- outpowlaw$estimate
Cpowlaw<- (mumle - 1) / (a^(1-mumle) - b^(1-mumle)) 
negloglikmlepowlaw<- outpowlaw$minimum

AICpowlaw<- 2 * negloglikmlepowlaw + 2*numpar

# Now for 95% CI of MLE:    
muvarynegloglik<- vector()       # negative log lik for each mu value
  for(i in 1:length(muvec))
    {
      muvarynegloglik[i] = powlawbound(muvec[i])
    }
critval1<- negloglikmlepowlaw + qchisq(0.95,1)/2
                    # 1 dof. Hilborn and Mangel (1997) p163.
muin95<- muvec[ muvarynegloglik < critval1] 
                    # mu values definitely in 95% confidence interval
minmuin95<- min(muin95)
maxmuin95<- max(muin95)


coefs<-data.frame(mumle, Cpowlaw, minmuin95, maxmuin95, AICpowlaw)
colnames(coefs)<-c("exponent_num", "constant", "CImin", "CImax", "AIC")
coefs
}


exp_mlefunc<- function(df, size, catch, dataorder=TRUE) {
	

if(dataorder==TRUE){
	mle_params(df, size, catch)
}

if(dataorder==FALSE){
	rawvalsall<-df[,size]	
x <<-df[,size]
a <<- min(rawvalsall)
b <<- max(rawvalsall)
n_obs <<- length(rawvalsall)


sumx<<-sum(rawvalsall)
}


# Analytically (Box 1, Edwards et al. 2007):
lambdamleinfanal = 1 / ( sumx/n_obs - a)    # MLE for lambda
outexpinf = nlm(explaw, lambda.startinfvec[1])
                                        # , print.level=2)
lambdamleinf = outexpinf$estimate
# if( abs(lambdamleinf - lambdamleinfanal) >10^(-1))
#    {stop("Analytical and numerical disagree for Exp model")}

negloglikmleexpinf = outexpinf$minimum

Cexpinf = lambdamleinf * exp(lambdamleinf * a)

AICexpinf = 2 * negloglikmleexpinf + 2*numparinf   
# Now for 95% CI of MLE:   

lambdavarynegloglikinf = vector()
    #                        negative log lik for each lambda value
for(i in 1:length(lambdavecinf))
 {
   lambdavarynegloglikinf[i] = explaw(lambdavecinf[i])
 }

critval2inf = negloglikmleexpinf + qchisq(0.95,1)/2 
      #            1 dof. Hilborn and Mangel (1997) p163.
lambdain95inf = lambdavecinf[ lambdavarynegloglikinf < critval2inf]
     #            lambda values definitely in 95% confidence interval
minexpin95<- min(lambdain95inf)
maxexpin95<- max(lambdain95inf)

coefs<<-data.frame(lambdamleinfanal,lambdamleinf, Cexpinf,minexpin95, maxexpin95,  AICexpinf)
colnames(coefs)<<-c("exponent_analytical", "exponent_num", "constant", "CImin", "CImax", "AIC")
coefs
}

expb_mlefunc<- function(df, size, catch, dataorder=TRUE) {
	


if(dataorder==TRUE){
	mle_params(df, size, catch)
}

if(dataorder==FALSE){
	rawvalsall<-df[,size]	
x <<-df[,size]
a <<- min(rawvalsall)
b <<- max(rawvalsall)
n_obs <<- length(rawvalsall)


sumx<<-sum(rawvalsall)
}

outexp = nlm(explawbound, lambda.start2vec[1]) 
                     #, print.level=2) 
lambdamle = outexp$estimate
negloglikmleexp = outexp$minimum

Cexp = lambdamle / (exp(- lambdamle * a) - exp( -lambdamle * b))
AICexp = 2 * negloglikmleexp + 2*numpar    

 #sNow for 95% CI of MLE:

lambdavarynegloglik = vector()
   #                         negative log lik for each lambda value
 for(i in 1:length(lambdavec))
  {
    lambdavarynegloglik[i] = explawbound(lambdavec[i])
  }
critval2 = negloglikmleexp + qchisq(0.95,1)/2
       #          1 dof. Hilborn and Mangel (1997) p163.
lambdain95 = lambdavec[ lambdavarynegloglik < critval2] 
plot(lambdavec, lambdavarynegloglik)


        #         lambda values definitely in 95% confidence interval
minexpbin95<- min(lambdain95)
maxexpbin95<- max(lambdain95)

coefs<<-data.frame(lambdamle, Cexp,minexpbin95, maxexpbin95,  AICexp)
colnames(coefs)<<-c("exponent_num", "constant","CImin", "CImax" , "AIC")
coefs
}



## lm function for comparisons
#lmfunction<-function(df) lm( log10(COUNT_SUM) ~ log10(MASS), data = df)

#### master function for AIC, estimates, all dists, all methods

isd_function<-function(df, vars, func, size, catch, dataorder=TRUE) {
	#data<-func(df)
	data<-df

	pl_coefs<-ddply(data, vars, pl_mlefunc, size, catch, dataorder)
	pl_coefs$method<-"PL"
	
	plb_coefs<-ddply(data, vars, plb_mlefunc, size, catch, dataorder)
	plb_coefs$method<-"PLB"
	plb_coefs$exponent_analytical<-"NA"

	exp_coefs<-ddply(data, vars, exp_mlefunc, size, catch, dataorder)
	exp_coefs$method<-"EXP"
	#exp_coefs$CImin<-"NA"
	#exp_coefs$CImax<-"NA"

	expb_coefs<-ddply(data, vars, expb_mlefunc, size, catch, dataorder)
	expb_coefs$method<-"EXPB"
	expb_coefs$exponent_analytical<-"NA"
	#	expb_coefs$CImin<-"NA"
	#expb_coefs$CImax<-"NA"


mle_est<-rbind(pl_coefs, plb_coefs, exp_coefs, expb_coefs)

AICall = data.frame(pl_coefs[,vars], pl_coefs$AIC, plb_coefs$AIC, exp_coefs$AIC, expb_coefs$AIC)
colnames(AICall)<-c(vars, "AICpl", "AICplb", "AICexp", "AICexpb")
var_l<-length(vars)
AICmin = data.frame(apply(AICall[,(var_l+1):(var_l+4)],1, min))	
Delta = data.frame(AICall[,(var_l+1)]-AICmin, AICall[,(var_l+2)]-AICmin,
 AICall[,(var_l+3)]-AICmin, AICall[,(var_l+4)]-AICmin)
colnames(Delta)<-c( "Deltapl", "Deltaplb", "Deltaexp", "Deltaexpb")
like = data.frame(exp( - 0.5 * Delta$Deltapl), exp( - 0.5 * Delta$Deltaplb),
	exp( - 0.5 * Delta$Deltaexp), exp( - 0.5 * Delta$Deltaexpb))
Aweights = data.frame(like/rowSums(like))

colnames(Aweights)<-c("PL", "PLB", "EXP", "EXPB")
Aweights[,vars]<-pl_coefs[,vars]

return(list(mle_est,AICall, Aweights))

}


# ## sample functions for standardising sample size

# sample.func<-function(df, sample_n=40, replicate="REPLICATEID") {
	
# 	sites<-unique(df[[replicate]])
# 	sites_sample<-data.frame(sample(sites, size=sample_n, replace=T))
	
# 	data<-df[df[[replicate]]%in%sites_sample[,1],]
# 	dups<-sites_sample[duplicated(sites_sample),]
# 	data<-rbind(data,df[df[[replicate]]%in%dups,] )

# 	#sample_agg<<-aggregate(COUNT ~ Region1 + ISLAND + MASS + STATE,data=data, sum)
	

# 	#sample_agg
#     data 
# }

#### rank frequency plot functions

## for ISD 
rf_plot_func<-function(df, size,catch,location="sample", cex=1, labs=FALSE){



size_data<-dataorder(df, size=size, catch=catch)



		plot(x=sort(size_data[,1], decreasing=T), y=1:dim(size_data)[1],
			  log="xy",axes=F,xlab="", ylab="",
			  main=unique(df[,location]), cex.main=cex) 
		axis(1, at=c(1,10, 20, 100, 1000, 10000, 100000, 1000000), labels=c(1,10, 20, 100, 1000, 10000, 100000, 1000000))
	axis(2, at=c(1, 10, 100, 1000, 10000))
		
		if (labs==TRUE){
	mtext( side=2, text="Number of body sizes > x", cex=cex, line=2)
	mtext( side=1, text="Body size, x (grams)", cex=cex, line=2)
#mtext(unique(df$location), cex=0.5)
}

	a<-min(size_data)
	b<-max(size_data)
	#mumle<-round(plb_all[which(plb_all$names==names(samples[i])),]$SLOPE, digits=2)
	n_obs<-dim(size_data)[1]



	isd<-plb_mlefunc(df, size=size, catch=catch)
	mumle<-isd$exponent_num


	xrankfreq<- seq(a, b, length=n_obs)
	yrankfreqpow<- (xrankfreq^( 1 - mumle) - b^(1 - mumle)) / (a^(1 - mumle) - b^(1 - mumle) ) * n_obs

	lines(x=xrankfreq, y=yrankfreqpow, lty=1 ) 
	text(100,10,  paste("ISDE =", round(-mumle, digits=2)))

}

### for div spectra
rf_div_plot_func<-function(df, size="MASS", location="sample", cex=1, labs=FALSE){


	size_data<-data.frame(df[,size] )

		plot(x=sort(size_data[,1], decreasing=T), y=1:dim(size_data)[1],
			  log="xy",axes=F,xlab="", ylab="",
			  main=unique(df[,location]),
			   cex.main=cex) 
		axis(1, at=c(1,10, 20, 100, 1000, 10000, 100000, 1000000), labels=c(1,10, 20, 100, 1000, 10000, 100000, 1000000))
	axis(2, at=c(1, 10, 100, 1000, 10000))
		
		if (labs==TRUE){
	mtext( side=2, text="Number of body masses > x", cex=cex, line=2)
	mtext( side=1, text="Asymptotic body mass, x (grams)", cex=cex, line=2)
#mtext(unique(df$location), cex=0.5)
}

	a<-min(size_data[,1])
	b<-max(size_data[,1])
	#mumle<-round(plb_all[which(plb_all$names==names(samples[i])),]$SLOPE, digits=2)
	n_obs<-length(size_data[,1])



	isd<-plb_mlefunc(size_data,size=colnames(size_data)[1],  dataorder=FALSE)
	mumle<-isd$exponent_num


	xrankfreq<- seq(a, b, length=n_obs)
	yrankfreqpow<- (xrankfreq^( 1 - mumle) - b^(1 - mumle)) / (a^(1 - mumle) - b^(1 - mumle) ) * n_obs

	lines(x=xrankfreq, y=yrankfreqpow, lty=1 ) 
	text(100,10,  paste("Exponent =", round(-mumle, digits=2)))

}

lbn_func<-function(data, size="MASS", catch){

print(paste("ISD estimated using", size))

		a<<-min(data[,size])
	b<<-max(data[,size])


breaks=vector()
x=1
repeat {
	x<-x+1
	breaks[1]=a
	breaks[x]=2*breaks[x-1]
	if (breaks[x] > b) {break}
}


widths=breaks[1:length(breaks)-1]




### to get midpoints
freqs<-dataorder(data, size=size, catch=catch)

### identify size classes that double in size
histogram<-hist(freqs[,1],breaks=breaks, plot=F)




### break each size into a size class
cuts<-cut(freqs[,1], breaks, right=F)

bios<-cbind(freqs, cuts)
colnames(bios)<-c("freqs", "cuts")


### find sum of sizes within each size class
biomass<-tapply(bios$freqs, bios$cuts, sum)

normbiomass<-biomass/widths

### take midpoints of size classes
mids<-histogram$mids

## add counts in each size class for ISD
counts<-histogram$counts
## normalise counts
normcounts<-counts/widths

# 	slopes<-data.frame(coef(lm(log10(counts[counts>0]) ~ log10(mids[which(counts>0)]))))
# 	slopes[1,2]<-"intercept"
# 	slopes[2,2]<-"slope"


# colnames(slopes)<-c("value", "stat")


	norm_slopes<-data.frame(coef(lm(log10(normcounts[normcounts>0]) ~ log10(mids[which(normcounts>0)]))))
	norm_slopes[1,2]<-"intercept"
	norm_slopes[2,2]<-"slope"
colnames(norm_slopes)<-c("value", "stat")

	
# slopes<-cbind(slopes, norm_slopes)
norm_slopes$type<-"ISD"

	return(norm_slopes)
}

