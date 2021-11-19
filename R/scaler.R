#' Scaling and centering explanatory covariates
#'
#' function to scale numeric and categorical variables for multi-model approaches. Categorical covariates must be provided as factors and will be converted to 0,1 dummy variables.
#' @param df = dataset containing all response, explanatory, and grouping factors
#' @param ID = names of grouping variables that should not be scaled; provide as c("Var1", "Var2")
#' @param centered = should the function center the mean at 0? (defaults to TRUE)
#' @param scaled = should the function scale each variable by its standard deviation? (defaults to TRUE)
#' @keywords multimodel, demean, scale
#' @export
#' @examples
#' scaler

scaler<-function(df, ID, centered=TRUE, scaled=TRUE, cats = TRUE, ...){

  ## Function parameters
  #	df = dataframe containing all response, explanatory, and grouping factors
  #	ID = names of grouping variables that should not be scaled; provide as c("Var1", "Var2")
  #	centered = should the function center the mean at 0?
  #	scaled = should the function scale each variable by its standard deviation?

  ###--------------------Begin function--------------------###


  # strip NAs and drop NA factor levels
  df<-df[complete.cases(df),]
  df<-droplevels(df)

  ## extract ID variables
  ID.vars<-df[,colnames(df)%in%ID]

  #--------------------scale the numeric variables--------------------#

  numerics<-sapply(df, is.numeric)
  dat_cont<-df[, numerics]
  dat_cont<-dat_cont[,!colnames(dat_cont)%in%ID, drop=F] ## drop any numeric ID variables
  scaled_cont<-scale(dat_cont, center=TRUE)
  colnames(scaled_cont)<-colnames(dat_cont)

  scaled.df<-cbind(ID.vars, data.frame(scaled_cont))
  return(scaled.df)
  }

  ## END

