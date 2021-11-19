#' GAM predictions
#'
#' Function builds GAM and dredges saturated models, fixed to maximum of 3 covariates to prevent overfitting, then predicts response against range of predictor covariates, for either 
#' 1. the top-ranked model (AICc > 2 lower than candidate models) OR 2. across a top-ranking model set (AICc < 7 of top-ranked model)
#' @param dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
#' @param exp.names = explanatory covariate names, passed as vector of characters
#' @param indicator = y variable of interest (character)
#' @param family = GLM family distribution, takes 'gaussian' or 'Gamma'
#' @param base_k = number of knots for each smoother term. Default = -1, mgcv  uses generalized cross-validation
#' @param smoother = cubic regression spline method (cr is default)
#' @param n.param.max = maximum number of parameters included in each reduced model
#' @keywords multimodel
#' @export
#' @examples
#' gam_predict


gam_predict<-function(dataset, exp.names, indicator, family, base_k=-1, smoother='cr',n.param.max=3){

  # library(gamm4)
  library(mgcv)
  library(MuMIn)

# if(!is.null(new.preddata)) {dataset <- new.preddata}


  # Arguments:
  # dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
  # exp.names = explanatory covariate names, passed as vector of characters
  # newdata = include for building models with 'dataset' but generating predictions with 'newdata' (e.g. HUMANS = 0)
  # indicator = y variable of interest (character)
  # family = GAM family distribution, takes 'gaussian' or 'Gamma'
  # base_k = number of knots for each smoother term. Default = -1, mgcv  uses generalized cross-validation
  # smoother = cubic regression spline method (cr is default)
  # n.param.max = maximum number of parameters included in each reduced model

# Returns list:
#                 1) predicted indicator estimates (across top 7 AICc models)
#                 2) confidence intervals around predictions
#                 3) deviance explained by top model(s)

  #assigning dataset to global env. to make function general. had issues using local objects with MuMIn 
  dataset<<-dataset

# pdf(file=paste('MMI_GAM_diagnostic_', indicator,'.pdf', sep=""))

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

# model diagnostic saved to pdf
gam.check(M_FULL)

## check for collinearity
# print(corvif(M_FULL))

# par(mfrow=c(2,2))
# hist(dataset[,indicator], main=paste('Hist of', indicator, sep=' '))
# plot(fitted(M_FULL), dataset[,indicator], main=paste('Fitted vs. obs', indicator, sep=' '))
# plot(resid(M_FULL), main=paste('Residuals of', indicator, sep=' '))
# hist(resid(M_FULL), main=paste('Hist of residuals of', indicator, sep=' '))
# dev.off()


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

t<-unlist(strsplit(exp.names, ' + ', fixed=TRUE))
vars <- as.data.frame(matrix(t, nrow = length(t), ncol = 1))
colnames(vars) <- c("Var")
vars$Var<-as.character(vars$Var)

# build empty predictor data frame, predicting indicator across 100 values of each exp.names
preds<-as.data.frame(matrix(NA, nrow=100, ncol=length(vars$Var)))
colnames(preds)<-vars$Var
preds.var<-as.data.frame(matrix(NA, nrow=100, ncol=length(vars$Var)))
colnames(preds.var)<-vars$Var




  #looping through all variables
for (y in 1:(nrow(vars))){
  #Pull out the variable of interest
  var.temp <-vars$Var[y]
    #Create a new data frame from prediction the focal predictor is sweep across the range of observed 
    #values and everything else is held at its mean
    newdata <- data.frame(matrix(0, nrow=100, ncol=length(vars$Var)))
    colnames(newdata)<-vars$Var
    newdata$test.var = seq(from = min(dataset[,var.temp]), to = max(dataset[,var.temp]), length.out = 100)
    

    #Find which is the given variable and update the column name
    newdata <- newdata[,!(names(newdata) %in% var.temp)]
    colnames(newdata)[which( colnames(newdata)== "test.var")] <- paste(var.temp)
    
    ## predict response across range of predictor covariate
      temp.mod <- get.models(M.set, delta<2)
      pred <- predict(temp.mod[[1]], newdata = newdata, type='response') # Get the raw prediction for the given model
      ## add to prediction dataframe
     preds[,y]<-pred
     

    }

#--------------------------------------------------------------------------------#
        # Get deviance explained for top model (< 2 ∆AIC) #
#--------------------------------------------------------------------------------#
  dev.exp<-data.frame(dev.exp=NA, indicator=NA)
  dev.exp$dev.exp<-1-(temp.mod[[1]]$dev/temp.mod[[1]]$null)

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

# check for error with top weights
# if(sum(top.weights)!= 1) stop("Top weights do not sum to 1") 

t<-unlist(strsplit(exp.names, ' + ', fixed=TRUE))
vars <- as.data.frame(matrix(t, nrow = length(t), ncol = 1))
colnames(vars) <- c("Var")
vars$Var<-as.character(vars$Var)

# build empty predictor data frame, predicting indicator across 100 values of each exp.names
preds<-as.data.frame(matrix(NA, nrow=100, ncol=length(vars$Var)))
colnames(preds)<-vars$Var
preds.var<-as.data.frame(matrix(NA, nrow=100, ncol=length(vars$Var)))
colnames(preds.var)<-vars$Var


  #looping through all variables
for (y in 1:(nrow(vars))){
  #Pull out the variable of interest
  var.temp <-vars$Var[y]
    #Create a new data frame from prediction the focal predictor is sweep across the range of observed 
    #values and everything else is held at its mean
    newdata <- data.frame(matrix(0, nrow=100, ncol=length(vars$Var)))
    colnames(newdata)<-vars$Var
    newdata$test.var = seq(from = min(dataset[,var.temp]), to = max(dataset[,var.temp]), length.out = 100)

    #Find which is the given variable and update the column name
    newdata <- newdata[,!(names(newdata) %in% var.temp)]
    colnames(newdata)[which( colnames(newdata)== "test.var")] <- paste(var.temp)
    
    
    #Create matrices to hold the raw and weighted residuals for each top model for the given variable
    wt.pred <- matrix(NA, nrow = len , ncol = nrow(newdata))
    raw.pred <-  matrix(NA, nrow = len , ncol = nrow(newdata))
    var.wt.pred<-matrix(NA, nrow=1, ncol=nrow(newdata))
    for (j in 1:len){
      temp.mod <- get.models(M.set, subset = j)
      pred <- predict(temp.mod[[1]], newdata = newdata, type='response') # Get the raw prediction for the given model
      raw.pred[j,] <- pred
      pred.wt.temp <- pred*top.weights[j] # Weight the prediction by the corresponding model weight
      wt.pred[j,] <- pred.wt.temp
    }
    
    #Calculated the weighted mean prediction
    sum.pred.wt <- colSums(wt.pred)
    # Caculate weighted sample variance
    
    # calculate the weighted sample variance for each model
    for (i in 1:length(pred)){
    var.wt.pred[i]<-sum(top.weights*((raw.pred[,i]-sum.pred.wt[i])^2))
      }
    ## add to prediction dataframe
    preds[,y]<-sum.pred.wt
    preds.var[,y]<-as.vector(var.wt.pred)

}
#--------------------------------------------------------------------------------#
        # Get deviance explained for each top model (< 7 AIC) #
#--------------------------------------------------------------------------------#
tt<-data.frame(subset(M.set, delta<7)) 
models<-rownames(tt)
dev.exp<-data.frame(models)

for (i in 1:len){
  # pseudo-R2: 1 - residualdeviance/nulldeviance; from Extending the linear model in R (Faraway)
  temp.mod <- get.models(M.set, subset = i)
  dev.exp$dev.exp[i]<-1-(temp.mod[[1]]$dev/temp.mod[[1]]$null)
  }

dev.exp$models<-NULL


}



preds$indicator<-indicator
preds.var$indicator<-indicator
dev.exp$indicator<-indicator

result<-list(preds, preds.var, dev.exp)
return(result)


}