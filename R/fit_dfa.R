#' Fit DFA to timeseries datasets
#'
#' This builds on DFA tutorial here: https://nwfsc-timeseries.github.io/AFTSLabbook/sec-dfa-plot-data.html
#' @param df is a matrix of time series
#' @param proc is the number of underlying processes, defaults to 3
#' @param id is the name of each time series
#' @param plot specifies whether fits and loadings should be plotted, defaults to TRUE
#' @param yr_first is the first year of the time series
#' @param yr_last is the last year of the time series
#' @param clr is a vector of colours applied to each time series if plot = TRUE
#' @keywords DFA
#' @export
#' @examples
#' fit_DFA




fit_DFA <- function(df, proc=3, id, plot=TRUE, yr_frst, yrlast, clr=clr){

    mm <- proc
    ## 'BB' is identity: 1's along the diagonal & 0's elsewhere
    BB <- "identity"  # diag(mm)
    ## 'uu' is a column vector of 0's
    uu <- "zero"  # matrix(0,mm,1)
    ## 'CC' and 'cc' are for covariates
    CC <- "zero"  # matrix(0,mm,1)
    cc <- "zero"  # matrix(0,1,wk_last)
    ## 'QQ' is identity
    QQ <- "identity"  # diag(mm)

    # get number of time series
    N_ts <- dim(df)[1]

    ## Fit the observation model
    ## 'ZZ' is loadings matrix
    ZZ <- matrix(0, nrow = N_ts, ncol = proc, byrow = TRUE)

        for(i in 1:N_ts){
            for (j in 1:proc){
            ZZ[i, j]<-paste('z', i, j, sep='')
        }}

    ZZ[upper.tri(ZZ)]<-0

    ## 'aa' is the offset/scaling
    aa <- "zero"
    ## 'DD' and 'd' are for covariates
    DD <- "zero"  # matrix(0,mm,1)
    dd <- "zero"  # matrix(0,1,wk_last)
    ## 'RR' is var-cov matrix for obs errors
    RR <- "diagonal and unequal"

    ## list with specifications for model vectors/matrices
    mod_list <- list(B = BB, U = uu, C = CC, c = cc, Q = QQ, Z = ZZ, 
        A = aa, D = DD, d = dd, R = RR)

    ## list with model inits
    init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
    ## list with model control parameters
    con_list <- list(maxit = 3000, allow.degen = TRUE)

    ## fit MARSS
    dfa <- MARSS(y = df, model = mod_list, inits = init_list, 
        control = con_list)

    # The first proc * n_ts parameter estimates Z.z## are the loadings of each observed time series on the proc hidden states. 
    #  The next n_ts estimates R.(,) are the variances of the observation errors  (vi,t)(vi,t) . 
    # The last proc values, x0.X#, are the estimates of the initial states at  t=0t=0 

    if(plot == TRUE){

        ## Rotate trends and loadings
        ## get the estimated ZZ
        Z_est <- coef(dfa, type = "matrix")$Z
        ## get the inverse of the rotation matrix
        H_inv <- varimax(Z_est)$rotmat

        ## rotate factor loadings
        Z_rot = Z_est %*% H_inv
        ## rotate processes
        proc_rot = solve(H_inv) %*% dfa$states

        ## now plot hidden trends and loadings 
        ylbl <- id
        w_ts <- seq(dim(df)[2])
        layout(matrix(c(1:(mm*2)), mm, 2), widths = c(2, 1))
        ## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
        par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
        ## plot the processes
        for (i in 1:mm) {
            ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
            ## set up plot area
            plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
                xlab = "", ylab = "", xaxt = "n")
            ## draw zero-line
            abline(h = 0, col = "gray")
            ## plot trend line
            lines(w_ts, proc_rot[i, ], lwd = 2)
            lines(w_ts, proc_rot[i, ], lwd = 2)
            ## add panel labels
            mtext(paste("State", i), side = 3, line = 0.5)
            axis(1, (0:dim(focal)[2]) + 1, yr_frst + 0:dim(focal)[2])
        }
        ## plot the loadings
        minZ <- 0
        ylm <- c(-1, 1) * max(abs(Z_rot))
        for (i in 1:mm) {
            plot(c(1:N_ts)[abs(Z_rot[, i]) > minZ], 
                as.vector(Z_rot[abs(Z_rot[, i]) > minZ, i]), 
                type = "h", lwd = 2, xlab = "", ylab = "", 
                xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
            for (j in 1:N_ts) {
                if (Z_rot[j, i] > minZ) {
                    text(j, -0.001, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
                        col = clr[j])
                }
                if (Z_rot[j, i] < -minZ) {
                    text(j, 0.001, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
                        col = clr[j])
                }
                abline(h = 0, lwd = 1.5, col = "gray")
            }
            mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
            }
        }
        print(paste('AIC = ', round(dfa$AIC, 2)))
        return(dfa)
}
