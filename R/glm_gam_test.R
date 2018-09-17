#' Comparison of GAM and GLM models
#'
#' Function builds GLM or GAM, compares AIC and outputs the winner.
#' @param dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
#' @param exp.names = explanatory covariate names, passed as vector of characters
#' @param indicator = y variable of interest (character)
#' @param family = GLM family distribution, takes 'gaussian' or 'Gamma'
#' @keywords model comparison
#' @export
#' @examples
#' glm.gam.test

glm_gam_test<-function(dataset, exp.names, indicator, family){

  # library(gamm4)
  library(mgcv)
  library(MuMIn)

  # Arguments:
  # dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
  # exp.names = explanatory covariate names, passed as vector of characters
  # indicator = y variable of interest (character)
  # family = GLM family distribution, takes 'gaussian' or 'Gamma'

#--------------------------------------------------------------------------------#
  # Build global model and check for collinearity, and print model diagnostic plots#
#--------------------------------------------------------------------------------#
## rearrange exp.names for gam notation
t<-strsplit(exp.names, ',', fixed=TRUE)

linear.exp.names<-paste(t, collapse=' + ')
smooth.exp.names<-paste('s(', t, ')', sep='')
smooth.exp.names<-paste(smooth.exp.names,  collapse=' + ')

# create formula for GLM and GAM
f <- as.formula(paste(indicator, linear.exp.names, sep="~"))
f.smooth <- as.formula(paste(indicator, smooth.exp.names, sep="~"))

# build GLM
if(family == 'Gamma'){
M_FULL.glm<-glm(f,   data=dataset, family=Gamma(link=log))}

if(family == 'gaussian'){
M_FULL.glm<-glm(f,   data=dataset, family='gaussian')}



## build GAM
if(family == 'Gamma'){
M_FULL.gam<-gam(f.smooth,   data=dataset, family=Gamma(link=log))}

if(family == 'gaussian'){
M_FULL.gam<-gam(f.smooth,   data=dataset, family='gaussian')}

# extract AIC for GLM and GAM
aic.glm<-AICc(M_FULL.glm)
aic.gam<-AICc(M_FULL.gam)

## compare full gam with full linear model. Print results.
if(aic.glm>aic.gam){print(paste("GAM supported. AICc", round(aic.glm-aic.gam, digits=2), 'units lower'))}
if(aic.gam>aic.glm){print(paste("GLM supported. AICc", round(aic.gam-aic.glm, digits=2), 'units lower'))}

## check for collinearity
if(aic.glm>aic.gam){  print(concurvity(M_FULL.gam)) }
if(aic.gam>aic.glm){  print(vif(M_FULL.glm)) }



## Save best supported model type for remaining analyses
if(aic.glm>aic.gam){  M_FULL <- M_FULL.gam }
if(aic.gam>aic.glm){  M_FULL <- M_FULL.glm }

# pdf(file=paste('MMI_GLM_GAM_diagnostic_', indicator,'.pdf', sep=""))
par(mfrow=c(2,2))
hist(dataset[,indicator], main=paste('Hist of', indicator, sep=' '))
plot(fitted(M_FULL), dataset[,indicator], main=paste('Fitted vs. obs', indicator, sep=' '))
plot(resid(M_FULL), main=paste('Residuals of', indicator, sep=' '))
hist(resid(M_FULL), main=paste('Hist of residuals of', indicator, sep=' '))
# dev.off()

return(M_FULL)

}
