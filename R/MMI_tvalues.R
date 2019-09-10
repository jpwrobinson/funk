#' Multimodel averaging with t-values
#'
#' Function builds GLM, runs multi-model inference to produce: 
#' 1. Standardized t-values (and standard error) for each covariate across full model set (effect sizes)
#' 2. R-squareds of each model < 7 AICc of top model
#' 3. Predicted y across each covariate, weighted by AIC for all models < 7 AICc of top model (and variance in predictions)
#' @param dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
#' @param exp.names = explanatory covariate names, passed as vector of characters
#' @param indicator = y variable of interest (character)
#' @param family = GLM family distribution, takes 'gaussian' or 'Gamma'
#' @keywords multimodel
#' @export
#' @examples
#' glm.avg.tvalue



mmi_tvalue<-function(M_FULL, dataset, exp.names, ranef, indicator, family, t.subset=FALSE){

  library(MuMIn)
  library(piecewiseSEM)
  library(car)

  # Arguments:
  # dataset = dataset containing y and all x covariates. should be scaled and centered (mean = 0, sd = 1)
  # exp.names = explanatory covariate names, passed as vector of characters
  # ranef = random effects, passed as vector of characters
  # indicator = y variable of interest (character)
  # family = GLM family distribution, takes 'gaussian' or 'Gamma'


# pdf(file=paste('MMI_GLM_diagnostic_', indicator,'.pdf', sep=""))

#--------------------------------------------------------------------------------#
	# Build global model and check for collinearity, and print model diagnostic plots#
#--------------------------------------------------------------------------------#

## check for collinearity
print('Checking collinearity. VIF values are:')
print(vif(M_FULL))

par(mfrow=c(2,2))
hist(dataset[,indicator], main=paste('Hist of', indicator, sep=' '))
plot(fitted(M_FULL), dataset[,indicator], main=paste('Fitted vs. obs', indicator, sep=' '))
plot(resid(M_FULL), main=paste('Residuals of', indicator, sep=' '))
hist(resid(M_FULL), main=paste('Hist of residuals of', indicator, sep=' '))
# dev.off()
#--------------------------------------------------------------------------------#
			# Dredge for delta < 7 top models #
#--------------------------------------------------------------------------------#

M_FULL_SET<-dredge(M_FULL, rank="AICc")
print(head(M_FULL_SET))

## if variance estimates are problematic (i.e. big effects scoop up variation in models with low numbers of covariates), then
## extract models with AIC < 7
if(t.subset == TRUE){
  M_FULL_SET<-subset(M_FULL_SET, delta<7)}

#--------------------------------------------------------------------------------#
# Model average absolute value of t-statistics for measure of variable importance ACROSS ALL MODELS #
#--------------------------------------------------------------------------------#
var.imp <- as.data.frame(matrix(NA, nrow = (ncol(M_FULL_SET) - 6), ncol = 4))
colnames(var.imp) <- c("Var", "RI.t.abs", "RI.t.ratio",  "var.t")

# Loop through all the variables in the model
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
top.models<-get.models(M_FULL_SET, subset=delta<7)

tt<-data.frame(subset(M_FULL_SET, delta<7))	
models<-rownames(tt)
r2<-data.frame(models)

for (i in 1:length(top.models)){
  # pseudo-R2: 1 - residualdeviance/nulldeviance; from Extending the linear model in R (Faraway)
  r2$r2.marg[i]<-piecewiseSEM::rsquared(top.models[[models[i]]])[1,5]
  r2$r2.cond[i]<-piecewiseSEM::rsquared(top.models[[models[i]]])[1,6]
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

    ## add first level of ran effects - these get cancelled
    ran.cancel<-data.frame(unique(dataset[,ranef[1]])[1], unique(dataset[,ranef[1]])[1])
    colnames(ran.cancel)<-ranef
    newdata<- cbind(newdata, ran.cancel)

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
      pred <- predict(temp.mod[[1]], newdata = newdata, type='response', re.form=NA) # Get the raw prediction for the given model
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
system("say YER DONE YA BISH")
}








