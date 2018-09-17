#' Counting unique items
#'
#' Function counts unique values.
#' @param x = vector to be counted
#' @keywords 
#' @export
#' @examples
#' uniques

uniques<-function(x) {length (unique(x))}
