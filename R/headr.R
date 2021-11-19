#' Short cut to printing full head on a tibble
#'
#'
#' @param
#' @keywords
#' @export
#' @examples
#' headr(iris)

headr<-function(x){
	library(tidyverse)
 	x %>% head() %>% data.frame()
 }
