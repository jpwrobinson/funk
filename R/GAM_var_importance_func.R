## Function builds GAM and calculates deviance explained by each covariate. 

## REFERENCE:
## Simon Wood; http://r.789695.n4.nabble.com/variance-explained-by-each-term-in-a-GAM-td836513.html
## Litzow et al. (2013) 
## Reassessing regime shifts in the North Pacific: incremental climate change and commercial fishing are necessary for explaining decadal-scale biological variability. 
## Global Change Biology.

gam.dev.exp<-function(dataset, exp.names, indicator, family, base_k=-1, smoother='cr'){

  # library(gamm4)
  library(mgcv)


  # Arguments:
  # dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
  # exp.names = explanatory covariate names, passed as vector of characters
  # indicator = y variable of interest (character)
  # family = GAM family distribution, takes 'gaussian' or 'Gamma'
  # base_k = number of knots for each smoother term. Default = -1, mgcv  uses generalized cross-validation
  # smoother = cubic regression spline method (cr is default)

# Returns: 
  # proportion deviance explained for each exp. covariate.

  #assigning dataset to global env. to make function general. had issues using local objects with MuMIn 
  dataset<<-dataset

#--------------------------------------------------------------------------------#
                          # Build global model 
#--------------------------------------------------------------------------------#

## rearrange exp.names for gam notation, specifying knots ('base_k') and smoother term ('smoother')
t<-unlist(strsplit(exp.names, ' + ', fixed=TRUE))
smooth.exp.names<-paste('s(',t, ', k =', base_k, ', bs = ',shQuote(smoother), ')', sep='')
smooth.exp.names<-paste(smooth.exp.names,  collapse=' + ')

# create formula for GAM
f.smooth <- as.formula(paste(indicator, smooth.exp.names, sep="~"))
f.null <- as.formula(paste(indicator, 1, sep='~'))

## build global GAM and null GAM
if(family == 'Gamma'){
M_FULL<-gam(f.smooth,   data=dataset, family=Gamma(link=log))
M_NULL<-gam(f.null, data=dataset, family=Gamma(link=log))
}

if(family == 'gaussian'){
M_FULL<-gam(f.smooth,  data=dataset, family='gaussian')
M_NULL<-gam(f.null, data=dataset, family='gaussian')
}

## save global deviance
dev.global<-deviance(M_FULL)
## save null deviance
dev.null<-deviance(M_NULL)


#--------------------------------------------------------------------------------#
                          # Build saturated models, save deviance for each 
#--------------------------------------------------------------------------------#

dev<-data.frame(EXP = t, deviance.sat=NA)


for(i in 1:length(t)){

  # recreate GAM formula dropping one exp. covariate
  t.temp<-unlist(strsplit(exp.names, ' + ', fixed=TRUE))
  t.temp<-t.temp[-i]
  smooth.exp.temp<-paste('s(',t.temp, ', k =', base_k, ', bs = ',shQuote(smoother), ')', sep='')
  smooth.exp.temp<-paste(smooth.exp.temp,  collapse=' + ')

  f.smooth.temp<-as.formula(paste(indicator, smooth.exp.temp, sep="~"))


  ## build saturated GAM
  if(family == 'Gamma'){
  M_temp<-gam(f.smooth.temp, sp=M_FULL$sp[-i],  data=ind.scaled, family=Gamma(link=log))}

  if(family == 'gaussian'){
  M_temp<-gam(f.smooth.temp,  sp=M_FULL$sp[-i], data=ind.scaled, family='gaussian')}

  dev.temp<-deviance(M_temp)
  dev$deviance.sat[i]<-dev.temp

}

#--------------------------------------------------------------------------------#
                  # Calculate proportion deviance explained  
#--------------------------------------------------------------------------------#
dev$deviance.exp<-(dev$deviance.sat - dev.global)/dev.null
dev$deviance.sat<-NULL
dev$indicator<-indicator


return(dev)

}