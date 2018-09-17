#' Testing predictive power for GAMs
#'
#' Function builds GAM and dredges saturated models, fixed to maximum of 3 covariates to prevent overfitting. Then predicts response against new predictor covariates
#' 1. the top-ranked model (AICc > 2 lower than candidate models) OR
#' 2. across a top-ranking model set (AICc < 7 of top-ranked model)
#' @param dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
#' @param newdataset = generating predictions with 'newdata'. should be scaled + centered (mean = 0, sd = 1) 
#' @param exp.names = explanatory covariate names, passed as vector of characters
#' @param indicator = y variable of interest (character)
#' @param family = GAM family distribution, takes 'gaussian' or 'Gamma'
#' @param base_k = number of knots for each smoother term. Default = -1, mgcv  uses generalized cross-validation
#' @param smoother = cubic regression spline method (cr is default)
#' @param n.param.max = maximum number of parameters included in each reduced model
#' @keywords multimodel, GAM, prediction
#' @export
#' @examples
#' gam.predict.newdata


gam_predict_newdata<-function(dataset, newdataset, exp.names, indicator, family, base_k=-1, smoother='cr',n.param.max=3){

  # library(gamm4)
  library(mgcv)
  library(MuMIn)


  # Arguments:
  # dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
  # newdataset = generating predictions with 'newdata'. should be scaled + centered (mean = 0, sd = 1) (using here for HUMANS set = 0)
  # exp.names = explanatory covariate names, passed as vector of characters
  # indicator = y variable of interest (character)
  # family = GAM family distribution, takes 'gaussian' or 'Gamma'
  # base_k = number of knots for each smoother term. Default = -1, mgcv  uses generalized cross-validation
  # smoother = cubic regression spline method (cr is default)
  # n.param.max = maximum number of parameters included in each reduced model

# Returns list:
#                 1) predicted indicator estimates (across top 7 AICc models)
#                 2) confidence intervals around predictions

#assigning dataset to global env. to make function general. had issues using local objects with MuMIn 
  dataset<<-dataset

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

# print(importance(M.set))
print(head(M.set))
## extract models with AIC < 2
top.models<-get.models(M.set, delta<2)

#--------------------------------------------------------------------------------#
### IF single top model with AICc > 2 lower than other candidate models
#--------------------------------------------------------------------------------#
    # Save predicted effects of each covariate  #
#--------------------------------------------------------------------------------#

if(length(top.models)==1) {

print(paste(indicator, ': single top model supported by AICc'))

    ## predict response at each datapoint with newdataset predictors
      temp.mod <- get.models(M.set, delta<2)
      pred <- predict(temp.mod[[1]], newdata = newdataset, type='response') # Get the raw prediction for the given model

}

#--------------------------------------------------------------------------------#
### IF no single top model:
#--------------------------------------------------------------------------------#
    # Save predicted effects of each covariate ACROSS TOP MODELS (< 7 AIC) #
#--------------------------------------------------------------------------------#

if(length(top.models)>1) {

## extract models with AIC < 7
top.models<-get.models(M.set, delta<7)
len<-length(top.models)
#recalc model weights for the top model set
top.weights <- M.set$weight[1:len]/sum(M.set$weight[1:len])


print(paste(indicator, ': no single top model. Using weighted predictions across', length(top.models), 'models'))

    
    #Create matrices to hold the raw and weighted residuals for each top model for the given variable
    wt.pred <- matrix(NA, nrow = len , ncol = nrow(newdataset))
    raw.pred <-  matrix(NA, nrow = len , ncol = nrow(newdataset))
    var.wt.pred<-matrix(NA, nrow=1, ncol=nrow(newdataset))

    for (j in 1:len){
      temp.mod <- get.models(M.set, subset = j)
      pred <- predict(temp.mod[[1]], newdata = newdataset, type='response') # Get the raw prediction for the given model
      raw.pred[j,] <- pred
      pred.wt.temp <- pred*top.weights[j] # Weight the prediction by the corresponding model weight
      wt.pred[j,] <- pred.wt.temp
    }
    
    #Calculated the weighted mean prediction
    sum.pred.wt <- colSums(wt.pred)
    
    # calculate the weighted sample variance for each model
    # for (i in 1:length(pred)){
    # var.wt.pred[i]<-sum(top.weights*((raw.pred[,i]-sum.pred.wt[i])^2))
    #   }
    ## add to prediction dataframe
    pred<-sum.pred.wt
    # preds.var<-as.vector(var.wt.pred)

}

pred<-data.frame(pred)
pred$indicator<-indicator
pred$ISLAND<-newdataset$ISLAND
# preds.var$indicator<-indicator

# result<-list(preds, preds.var)

return(pred)


}
