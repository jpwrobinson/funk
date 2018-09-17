
## Function builds GLM, runs multi-model inference to produce:
# 1. Standardized t-values (and standard error) for each covariate across full model set (effect sizes)
# 2. R-squareds of each model < 7 AICc of top model
# 3. Predicted y across each covariate, weighted by AIC for all models < 7 AICc of top model (and variance in predictions)



glm.avg.tvalue<-function(dataset, exp.names, indicator, family){

  library(MuMIn)

  # Arguments:
  # dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
  # exp.names = explanatory covariate names, passed as vector of characters
  # indicator = y variable of interest (character)
  # family = GLM family distribution, takes 'gaussian' or 'Gamma'


pdf(file=paste('figures/models/MMI_GLM_diagnostic_', indicator,'.pdf', sep=""))

#--------------------------------------------------------------------------------#
	# Build global model and check for collinearity, and print model diagnostic plots#
#--------------------------------------------------------------------------------#
# create formula for 
f <- as.formula(paste(indicator, exp.names, sep="~"))

## compare full gam with full linear model
if(family == 'Gamma'){
M_FULL<-glm(f,   data=dataset, family=Gamma(link=log))}

if(family == 'gaussian'){
M_FULL<-glm(f,   data=dataset, family='gaussian')}

## check for collinearity
print(vif(M_FULL))

par(mfrow=c(2,2))
hist(dataset[,indicator], main=paste('Hist of', indicator, sep=' '))
plot(fitted(M_FULL), dataset[,indicator], main=paste('Fitted vs. obs', indicator, sep=' '))
plot(resid(M_FULL), main=paste('Residuals of', indicator, sep=' '))
hist(resid(M_FULL), main=paste('Hist of residuals of', indicator, sep=' '))
dev.off()
#--------------------------------------------------------------------------------#
			# Dredge for delta < 7 top models #
#--------------------------------------------------------------------------------#

M_FULL_SET<-dredge(M_FULL, rank="AICc")
print(head(M_FULL_SET))
## extract models with AIC < 7
top.models<-get.models(M_FULL_SET, delta<7)

#--------------------------------------------------------------------------------#
# Model average absolute value of t-statistics for measure of variable importance ACROSS ALL MODELS #
#--------------------------------------------------------------------------------#
var.imp <- as.data.frame(matrix(NA, nrow = (ncol(M_FULL_SET) - 6), ncol = 4))
colnames(var.imp) <- c("Var", "RI.t.abs", "RI.t.ratio",  "var.t")

#Loop through all the variables in the model
for (i in 1:(ncol(M_FULL_SET) - 6)){
  var.temp <- colnames(M_FULL_SET)[i+1]
  var.imp[i,"Var"] <- var.temp
  
  #Subset for only models for a given predictor, don't recalulate the model weights or deltaAICcs
  M_FULL_SET_x <- M_FULL_SET[i = !is.na(M_FULL_SET[,paste(var.temp)]), j = 1:ncol(M_FULL_SET), recalc.weights =FALSE, recalc.delta = FALSE]
  
  #Pull out the t-statistics for the given variables for each model in which it appears
  #Put all values in a data frame
  t.values <- as.data.frame(matrix(NA, nrow = nrow(M_FULL_SET_x), ncol = 4))
  colnames(t.values) <- c('varx',  "model.wt", 'imp.t.abs', "imp.t.ratio")
  for (l in 1:nrow(M_FULL_SET_x)){
    #put in variable name
    t.values[l,"varx"] <- var.temp
    
    #Pull out one model for the set that includes the variable of interest
    temp_mod <- get.models(M_FULL_SET, which(!is.na(M_FULL_SET[,paste(var.temp)]))[l])[[1]]
    
    #Pull out the model weight for this model
    wt.temp <- M_FULL_SET$weight[which(!is.na(M_FULL_SET[,paste(var.temp)]))][l]
    t.values[l,"model.wt"] <- wt.temp

    #abs value of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])
    t.values[l, "imp.t.abs"] <-  RI.x.temp 
    
    #ratios of t-statistics for variable of interest
    RI.x.temp <- abs(coef(summary(temp_mod))[paste(var.temp),3])/max(abs(coef(summary(temp_mod))[-1,3]))
    t.values[l, "imp.t.ratio"] <-  RI.x.temp 
  }
  #AICc weighted average of t-statistics for variable of interest
  avgRI.x.abs <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[2]*x[1]))
  var.imp[i,"RI.t.abs"] <- avgRI.x.abs
  
   #AICc weighted average ratio of t-statistics for variable of interest
  avgRI.x <- sum(apply(as.matrix(t.values[,2:4]), 1, function(x) x[3]*x[1]))
  var.imp[i,"RI.t.ratio"] <- avgRI.x

  ### variance
  varT.x<-sum(t.values[,2]*((t.values[,3]-avgRI.x.abs)^2))
  var.imp[i, "var.t"] <- varT.x
}

#--------------------------------------------------------------------------------#
				# Get R-squareds for each top model (< 7 AIC) #
#--------------------------------------------------------------------------------#
top<-get.models(M_FULL_SET, subset=delta<7)
tt<-data.frame(subset(M_FULL_SET, delta<7))	
models<-rownames(tt)
r2<-data.frame(ncol=1, models)

for (i in 1:length(top)){
  # pseudo-R2: 1 - residualdeviance/nulldeviance; from Extending the linear model in R (Faraway)
  r2$r2[i]<-1-top[[models[i]]]$deviance/top[[models[i]]]$null.deviance
  }

#--------------------------------------------------------------------------------#
		# Save predicted effects of each covariate ACROSS TOP MODELS (< 7 AIC) #
#--------------------------------------------------------------------------------#
len<-length(top.models)
#recalc model weights for the top model set
top.weights <- M_FULL_SET$weight[1:len]/sum(M_FULL_SET$weight[1:len])

## check for error with top weights
if(sum(top.weights)!= 1) stop("Top weights do not sum to 1") 

# build empty predictor data frame, predicting indicator across 100 values of each exp.names
preds<-as.data.frame(matrix(NA, nrow=100, ncol=length(var.imp$Var)))
colnames(preds)<-var.imp$Var
preds.var<-as.data.frame(matrix(NA, nrow=100, ncol=length(var.imp$Var)))
colnames(preds.var)<-var.imp$Var

	#looping through all variables
for (y in 1:(nrow(var.imp))){
  #Pull out the variable of interest
  var.temp <-var.imp$Var[y]
    #Create a new data frame from prediction the focal predictor is sweep across the range of observed 
    #values and everything else is held at its mean
    newdata <- data.frame(matrix(0, nrow=100, ncol=length(var.imp$Var)))
    colnames(newdata)<-var.imp$Var
  	newdata$test.var = seq(from = min(dataset[,var.temp]), to = max(dataset[,var.temp]), length.out = 100)

    #Find which is the given variable and update the column name
    newdata <- newdata[,!(names(newdata) %in% var.temp)]
    colnames(newdata)[which( colnames(newdata)== "test.var")] <- paste(var.temp)
    
    #Create matrices to hold the raw and weighted residuals for each top model for the given variable
    wt.pred <- matrix(NA, nrow = len , ncol = nrow(newdata))
    raw.pred <-  matrix(NA, nrow = len , ncol = nrow(newdata))
    var.wt.pred<-matrix(NA, nrow=1, ncol=nrow(newdata))
    for (j in 1:len){
      temp.mod <- get.models(M_FULL_SET, subset = j)
      pred <- predict(temp.mod[[1]], newdata = newdata, type='response') # Get the raw prediction for the given model
      raw.pred[j,] <- pred
      pred.wt.temp <- pred* top.weights[j] # Weight the prediction by the corresponding model weight
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

# save each data frame within a list
var.imp$indicator<-indicator
r2$indicator<-indicator
preds$indicator<-indicator
preds.var$indicator<-indicator

mmi<-list(var.imp, r2, preds, preds.var)

# return list of data frames renamed with indicator of interest
return(assign(paste('mmi', indicator, sep="."), mmi))
}



## Function builds GLM or GAM, compares AIC and outputs the winner.

glm.gam.test<-function(dataset, exp.names, indicator, family){

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
t<-strsplit(exp.names, ' + ', fixed=TRUE)
smooth.exp.names<-paste('s(', t[[1]], ')', sep='')
smooth.exp.names<-paste(smooth.exp.names,  collapse=' + ')

# create formula for GLM and GAM
f <- as.formula(paste(indicator, exp.names, sep="~"))
f.smooth <- as.formula(paste(indicator, smooth.exp.names, sep="~"))
print(f.smooth)
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

return(M_FULL)

pdf(file=paste('figures/models/MMI_GLM_GAM_diagnostic_', indicator,'.pdf', sep=""))
par(mfrow=c(2,2))
hist(dataset[,indicator], main=paste('Hist of', indicator, sep=' '))
plot(fitted(M_FULL), dataset[,indicator], main=paste('Fitted vs. obs', indicator, sep=' '))
plot(resid(M_FULL), main=paste('Residuals of', indicator, sep=' '))
hist(resid(M_FULL), main=paste('Hist of residuals of', indicator, sep=' '))
dev.off()

}



## Function builds GAM and dredges saturated models, fixed to maximum of 3 covariates to prevent overfitting 
## Then predicts response against range of predictor covariates, for either 
## 1. the top-ranked model (AICc > 2 lower than candidate models) OR
## 2. across a top-ranking model set (AICc < 7 of top-ranked model)

gam.predict<-function(dataset, exp.names, indicator, family, base_k=-1, smoother='cr',n.param.max=3){

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

pdf(file=paste('figures/models/MMI_GAM_diagnostic_', indicator,'.pdf', sep=""))

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
dev.off()


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


## Function builds GAM and dredges saturated models, fixed to maximum of 3 covariates to prevent overfitting 
## Then predicts response against new predictor covariates
## 1. the top-ranked model (AICc > 2 lower than candidate models) OR
## 2. across a top-ranking model set (AICc < 7 of top-ranked model)

gam.predict.newdata<-function(dataset, newdataset, exp.names, indicator, family, base_k=-1, smoother='cr',n.param.max=3){

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


## Function builds GAM and dredges saturated models, fixed to maximum of 3 covariates to prevent overfitting 
## Then compares top model structure of full dataset vs. dataset reduced by one data point, for each data point removed
## 1. output = % of times that model terms are the same as the main results (top model or top model set)

gam.jacknife<-function(dataset, exp.names, indicator, family, base_k=-1, smoother='cr',n.param.max=3){

  # library(gamm4)
  library(mgcv)
  library(MuMIn)


  #assigning dataset to global env. to make function general. had issues using local objects with MuMIn 
  dataset<<-dataset

# if(!is.null(new.preddata)) {dataset <- new.preddata}

### CURRENTLY ONLY WORKS FOR IND.SCALED DATASET. FUNCTION IS NOT GENERAL. 
  # Arguments:
  # dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
  # exp.names = explanatory covariate names, passed as vector of characters
  # newdata = include for building models with 'dataset' but generating predictions with 'newdata' (e.g. HUMANS = 0)
  # indicator = y variable of interest (character)
  # family = GAM family distribution, takes 'gaussian' or 'Gamma'
  # base_k = number of knots for each smoother term. Default = -1, mgcv  uses generalized cross-validation
  # smoother = cubic regression spline method (cr is default)
  # n.param.max = maximum number of parameters included in each reduced model

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
