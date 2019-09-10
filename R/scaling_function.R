#' Scaling and centering explanatory covariates
#'
#' function to scale numeric and categorical variables for multi-model approaches. Categorical covariates must be provided as factors and will be converted to 0,1 dummy variables.
#' @param df = dataset containing all response, explanatory, and grouping factors 
#' @param ID = names of grouping variables that should not be scaled; provide as c("Var1", "Var2")
#' @param centered = should the function center the mean at 0? (defaults to TRUE)
#' @param scaled = should the function scale each variable by its standard deviation? (defaults to TRUE)
#' @param family = GLM family distribution, takes 'gaussian' or 'Gamma'
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
#   cats = do categorical variables exist, and need converting to dummy?

###--------------------Begin function--------------------###


	# strip NAs and drop NA factor levels
	df<-df[complete.cases(df),]
	df<-droplevels(df)

	## extract ID variables
	ID.vars<-df[,colnames(df)%in%ID] 

#--------------------scale the numeric variables--------------------#
	
	numerics<-sapply(df, is.numeric)
	dat_cont<-df[, numerics]
	# dat_cont<-dat_cont[,!colnames(dat_cont)%in%ID, drop=F] ## drop any numeric ID variables
	scaled_cont<-scale(dat_cont, center=TRUE)
	colnames(scaled_cont)<-colnames(dat_cont)

#--------------------scale the categorical variables-----------------#
	if(cats == TRUE){

	cats<-df[,!numerics] ## drop numerics
	cats<-cats[,!colnames(cats)%in%ID] ## drop ID variables
	cats<-as.data.frame(cats) ## convert to dataframe
	cat.names<-colnames(df[,!numerics])
	cat.names<-cat.names[!cat.names %in% ID]
	
	# ## if you only have 1 categorical variable, do this...
	# if(dim(cats)[2]==1){

	# 	i.levels<-levels(cats[,1])
	# 	cats[, 2]<-0
	# 	cats[, 2][cats[,1]==i.levels[2]]<-1
	# 	colnames(cats)[2]<-paste(i.levels[1],i.levels[2],"dummy", sep=".")

	# } else if (dim(cats)[2]>1){
	
	nvars <- dim(cats)[2]


	scaled_cat <- as.data.frame(matrix(0, nrow=nrow(cats), ncol=1))

	for(a in 1:nvars){ 

		nd <- as.data.frame(matrix(0, nrow=nrow(cats), ncol=1))
		counter = 0
		level.max<-uniques(cats[,a])

		repeat {

		counter <- counter + 1

		i.levels<-levels(cats[,a])
		nd[, counter]<-cats[,a]
		colnames(nd)[counter]<-cat.names[a]

		# if the first categorical variable has 2 levels, do this...
			if(length(i.levels)==2){
			
			nd[, 1+counter]<-0
			nd[, 1+counter][nd[,counter]==i.levels[2]]<-1
			colnames(nd)[1+counter]<-paste(i.levels[1],i.levels[2],"dummy", sep=".")
			} else {
			
		## for variables with more than 2 levels we need to add more than 2 dummy variables
				if(length(i.levels)>2){

				for(j in 2:length(i.levels)){
					nd[, 1+ counter + (j-2)]<-0
					nd[, 1+ counter + (j-2)][cats[,a]==i.levels[j]]<-1
					colnames(nd)[1 + counter + (j-2)]<-paste(i.levels[1],i.levels[j],"dummy", sep=".")
					
				}
			}
		}
		counter<-dim(nd)[2]
		if(counter == level.max){break}

	}

	scaled_cat<-cbind(scaled_cat, nd)

	}




#--------------------bind numeric and categorical together with ID.vars-----------------#
	scaled.df<-cbind(ID.vars, scaled_cont, scaled_cat[,-1])
	return(scaled.df)
}

else {
	scaled.df<-cbind(ID.vars, scaled_cont)
	return(scaled.df)
}

	## END
	}

