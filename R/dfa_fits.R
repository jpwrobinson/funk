#' Examine fits from a DFA object  
#'
#' For DFA fitted to time series datasets. This builds on DFA tutorial here: https://nwfsc-timeseries.github.io/AFTSLabbook/sec-dfa-plot-data.html
#' @param MLEobj is fitted MARSS DFA object
#' @param dd is explanatory covariates
#' @param alpha is significance level
#' @keywords DFA
#' @export
#' @examples
#' get_DFA_fits

get_DFA_fits <- function(MLEobj, dd = NULL, alpha = 0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type = "matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    H_inv <- varimax(ZZ)$rotmat
    ## check for covars
    if (!is.null(dd)) {
        DD <- coef(MLEobj, type = "matrix")$D
        ## model expectation
        fits$ex <- ZZ %*% H_inv %*% MLEobj$states + DD %*% dd
    } else {
        ## model expectation
        fits$ex <- ZZ %*% H_inv %*% MLEobj$states
    }
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for (tt in 1:TT) {
        RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, 
            , tt] %*% t(ZZ)
        SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% 
            t(MLEobj$states[, tt, drop = FALSE])
        VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1 - alpha/2) * SE + fits$ex
    fits$lo <- qnorm(alpha/2) * SE + fits$ex
    return(fits)
}

