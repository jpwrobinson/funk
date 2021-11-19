#' Short cut to printing full head on a tibble
#'
#' Function estimates standard error (sd / sqrt(n))
#' @param 
#' @keywords 
#' @export
#' @examples
#' uniques

headr<-function(x) x %>% head() %>% data.frame()
