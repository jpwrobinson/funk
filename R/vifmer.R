#' Multimodel averaging with t-values
#'
#' variance inflation factors for lme4 package: https://github.com/aufrank/R-hacks/blob/master/mer-utils.R. Adapted from rms::vif
#' @param fit = fitted model
#' @keywords
#' @export
#' @examples
#' vif.mer

vif.mer <- function (fit) {
    
    v <- vcov(fit)
    nam <- names(fixef(fit))

    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
        v <- v[-(1:ns), -(1:ns), drop = FALSE]
        nam <- nam[-(1:ns)]
    }
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
}
