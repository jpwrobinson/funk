#' Fit jacknife GAM analysis.
#'
#' This function runs a drop-one jacknife to estimate relative variable importance.
#' @param dataset is a data frame of response and explanatory covariates
#' @param exp.names is a vector of explanatory covariate names, matching colnames of dataset
#' @param indicator is the name of the focal response variable
#' @param family specifies the fitted distribution (e.g. 'gaussian', 'poisson', 'binomial')
#' @param base_k is the number of knots for each smoother term. Defaults to -1, where mgcv uses generalized cross-validation
#' @param smoother is the spline basis, defaults to 'cr'
#' @param n.param.max is the maximum number of fitted explanatory covariates, defaults to 3
#' @keywords gam
#' @export
#' @examples
#' gam.jacknife



gam_jacknife<-function(dataset, exp.names, indicator, family, base_k=-1, smoother='cr',n.param.max=3){

  library(mgcv)

### WORKS FOR UVC BIOMASS GAMS ONLY
  # Arguments:
  # dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
  # exp.names = explanatory covariate names, passed as vector of characters
  # indicator = y variable of interest (character)
  # family = GAM family distribution, takes 'gaussian' or 'Gamma'
  # base_k = number of knots for each smoother term. Default = -1, mgcv  uses generalized cross-validation

# Returns data frame:
#                 1) jackknife stability percentage

#--------------------------------------------------------------------------------#
  # Build global model and check for collinearity, and print model diagnostic plots#
#--------------------------------------------------------------------------------#

## rearrange exp.names for gam notation, specifying knots ('base_k') and smoother term ('smoother')
t<-unlist(strsplit(exp.names, ' + ', fixed=TRUE))
smooth.exp.names<-paste('s(',t, ', k =', base_k, ', bs = ',shQuote(smoother), ')', sep='')
smooth.exp.names<-paste(smooth.exp.names,  collapse=' + ')

# create formula for GAM
f.smooth <- as.formula(paste(indicator, smooth.exp.names, sep="~"))

## build GAM
if(family == 'Gamma'){
M_FULL<-gam(f.smooth,   data=dataset, family=Gamma(link=log))}

if(family == 'gaussian'){
M_FULL<-gam(f.smooth,  data=dataset, family='gaussian')}


### now dredge global model 
## commented out for later - wave and SST can't be contained in the same model
## m.max determines maximum number of covariates included in each model

M.set<-dredge(M_FULL,beta=FALSE, rank="AICc", 
                # subset=!(`s(WV, k = base_k, bs = "cr")`  &&  `s(SSTL, k = base_k, bs = "cr")`), 
                # extra = alist(AIC, "R^2", "adjR^2"), 
                m.lim=c(1,n.param.max))

## extract models with AIC < 2
top.models<-get.models(M.set, delta<2)
top.models.call<-get.models(M.set, delta<2)[[1]]


#--------------------------------------------------------------------------------#
### IF single top model with AICc > 2 lower than other candidate models
#--------------------------------------------------------------------------------#
    # compare jacknifed models to single model
#--------------------------------------------------------------------------------#

if(length(top.models)==1) {

  print(paste(indicator, ': single top model supported by AICc'))
  onemodel<-TRUE

}

#--------------------------------------------------------------------------------#
### IF no single top model:
#--------------------------------------------------------------------------------#
    # Compare jacknifed models to top model set (< 7 AIC) #
#--------------------------------------------------------------------------------#

if(length(top.models)>1) {

  ## extract models with AIC < 7
  top.models<-get.models(M.set, delta<7)
  len<-length(top.models)

  print(paste(indicator, ': no single top model. Comparing models to set of', length(top.models), 'models'))
  onemodel<-FALSE

}


## Now run same process for jacknifed dataset, dropping one data point each loop
jack<-numeric()
n.models<-numeric()

for(i in 1:nrow(dataset)){
  #i<-1
  # assigning jack.dat to global env. to allow get.models to work (doesn't recognise local env. data frame)
  jack.dat <<- dataset[-i,]

  ## build GAM
  if(family == 'Gamma'){
  M_JACK<-gam(f.smooth,   data=jack.dat, family=Gamma(link=log))}

  if(family == 'gaussian'){
  M_JACK<-gam(f.smooth,  data=jack.dat, family='gaussian')}


  Mjack.set<-dredge(M_JACK,beta=FALSE, rank="AICc", 
                  # subset=!(`s(WV, k = base_k, bs = "cr")`  &&  `s(SSTL, k = base_k, bs = "cr")`), 
                  # extra = alist(AIC, "R^2", "adjR^2"), 
                  m.lim=c(1,n.param.max))

  ## extract models with AIC < 2
  top.jack.models<-get.models(Mjack.set, delta<2)
  top.jack.models.call<-get.models(Mjack.set, delta<2)[[1]]


  if(length(top.jack.models)==1) {

  print(paste(indicator, ' jacknife', i, ': single top model supported by AICc'))
  onemodel<-TRUE

}

  # for multiple top models
  if(length(top.jack.models)>1) {

    print(paste(indicator, ' jacknife', i, ': no single top model. Comparing models to set of', length(top.jack.models), 'models'))
    ## extract models with AIC < 7
    top.jack.models<-get.models(Mjack.set, delta<7)
  }

## save number of models in top set
  n.models[i]<-length(top.jack.models)

  ## Now compare jacked model structure with full model structure

  ## For single top model and single jacked model
    if(length(top.jack.models)==1 & onemodel==TRUE) {

        jack[i]<-identical(top.models.call$pred.formula, top.jack.models.call$pred.formula)

    }

  ## For single top model but multiple jacked models
    if(length(top.jack.models)==1 & onemodel==FALSE) {

        jack[i]<-FALSE
    }

  ## For multiple top models but jacked single model
    if(length(top.jack.models)>1 & onemodel==TRUE) {

        jack[i]<-FALSE
    }

  ## For multiple top models and multiple jacked models
    if(length(top.jack.models)>1 & onemodel==FALSE) {

      if(length(top.models) != length(top.jack.models)) {
        jack[i]<-FALSE

      } else {
    
        multi.jack.ticker<-numeric()

        for (a in 1:length(top.models)) {

          # extract calls from each top model in set, jacked and full
           temp.top<-get.models(M.set, delta<7)[[a]]
           temp.jack<-get.models(Mjack.set, delta<7)[[a]]
         
            multi.jack.ticker[a]<-identical(temp.top$pred.formula, temp.jack$pred.formula)
        }

      if(sum(multi.jack.ticker) == length(top.models)) {jack[i]<-TRUE} 

      }

        
    }

}

jacks<-data.frame(jack)
jacks$model.set<-n.models
jacks$indicator<-indicator


print(paste('Jacknife stability:', round(sum(jacks$jack)/length(jacks$jack)*100,digits=2), '% of jacknifed models identical to full model structure'))

return(jacks)



}
