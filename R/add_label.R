#' Consistent labelling for base plots
#'
#' For adding easy labels to multipanel figures (a,b,c). Stolen from Sean Anderson (https://seananderson.ca/2013/10/21/panel-letters/)
#' @param xfrac is label position, in fraction from x axis
#' @param yfrac is label position, in fraction from y axis
#' @param label is text to add
#' @param pos is text alignment (see ?text)
#' @keywords plot
#' @export
#' @examples
#' add_label


add_label <- function(xfrac, yfrac, label, pos = 4, ...){
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}