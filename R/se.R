#' Estimating standard error
#'
#' Function estimates standard error (sd / sqrt(n))
#' @param 
#' @keywords 
#' @export
#' @examples
#' uniques

se<-function(x) {sd(x) / sqrt(length(x))}
