#' Pairs plots for correlations
#'
#' panel.smooth2 creates scatterplot with LOESS for data diagnostics.
#' @param x is x covariate
#' @param y is y covariate
#' @keywords multimodel
#' @export
#' @examples
#' panel.smooth2

panel.smooth2<-function (x, y, col = par("col"), bg = alpha('grey',0.5), pch = 21, 
    cex = 1,  span = 2/3, iter = 3, ...) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if(abs(cor(x,y))>0.5) {col.smooth='red'} else {col.smooth='black'}
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, lwd=1.5, ...)
    }
 
