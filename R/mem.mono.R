#' Multiplicative effect modeling with monotonic treatment effect assumption
#'
#' @description \code{mem.mono} is used to estimate the multiplicative effect under monotonic treatment effect assumption. Treatment can be continuous or ordinal.
#'
#' @param y The response vector. Should only take values 0 or 1.
#' @param z The treatment vector. Should be numerical.
#' @param va The covariates matrix for the target model.
#' @param vb The covariates matrix for the nuisance model.
#' @param alpha.start Starting values for the parameters in the target model.
#' @param beta.start Starting values for the parameters in the nuisance model.
#' @param max.step The maximal number of iterations to be passed into the optim function. The default is 10,000.
#' @param thres Threshold for judging convergence. The default is 1e-6.
#'
#' @details The estimation is based on the models: \eqn{log(RR(z0,z;v,\alpha)) = \alpha*v(z-z0)},
#' \eqn{log(RR(z_min,z_max;v,\beta)) = \beta*v}, where \eqn{RR(z0,z) = P(Y=1|Z=z,z)/(Y=1|Z=z,z0)}.
#'
#' @return A list consisting of
#'  \itemize{
#'    \item point.est - the point estimates.
#'    \item se.est - he standard error estimates.
#'    \item cov - estimate of the covariance matrix for the estimates.
#'    \item conf.lower - the lower limit of the 95\% confidence interval.
#'    \item conf.upper - the upper limit of the 95\% confidence interval.
#'    \item p.value - the two sided p-value for testing zero coefficients.
#'    \item coefficients - the matrix summarizing key information: point estimate, 95\% confidence interval and p-value.
#'    \item param - the measure of association.
#'    \item converged - Did the maximization process converge?
#'    \item likelihood - The likelihood under the estimates.
#' }
#'
#' @export
mem.mono = function(y, z, va, vb,
                    alpha.start, beta.start,
                    max.step=10^4, thres=10^-6){

        if (is.numeric(z)) stop("'z' should be numerical!")

        pa = length(alpha.start)
        pb = length(beta.start)
        ny = length(y)

        z.uniq = sort(unique(z))
        nz = length(z.uniq)

        neg.log.value = function(alpha, beta){
                #Pzmin_Pzmax = t(mapply(getProbScalarRR, va %*% alpha * (nz-1), vb %*% beta))
                Pzmin_Pzmax = t(mapply(getProbScalarRR, va %*% alpha * (z.uniq[nz] - z.uniq[1]), vb %*% beta))
                P.mat = matrix(0, ncol = nz, nrow = ny)
                P.mat[, c(1, nz)] = Pzmin_Pzmax
                # P.mat[, -c(1, nz)] = Pzmin_Pzmax[,1] * exp(va %*% alpha %*% t(c(1: (nz-2))))
                P.mat[, -c(1, nz)] = Pzmin_Pzmax[,1] * exp(va %*% alpha %*% t(z.uniq)[-c(1,nz)])

                n_0.vec = apply(P.mat, 2, function(x) sum(x==0))

                if(identical(all.equal(n_0.vec, rep(0, nz)), TRUE)){
                        value = 0

                        for(i in 1: length(z.uniq)){
                                value = value - sum(y[z==z.uniq[i]] * log(P.mat[z==z.uniq[i],i]) + (1-y[z==z.uniq[i]]) * log(1-P.mat[z==z.uniq[i],i]))

                        }
                        #save(z, P.mat,alpha, beta,va,Pzmin_Pzmax,nz, ny,vb, file = "test.Rdata")
                        return(value)
                }
                else{
                        return(10000)
                }
        }

        neg.log.likelihood = function(pars){
                alpha = pars[1:pa]
                beta = pars[(pa + 1):(pa + pb)]

                return(neg.log.value(alpha, beta))
        }

        neg.log.likelihood.alpha = function(alpha){
                return(neg.log.value(alpha, beta))
        }

        neg.log.likelihood.beta = function(beta){

                return(neg.log.value(alpha, beta))
        }

        Diff = function(x, y) sum((x - y)^2)/sum(x^2 + thres)
        alpha = alpha.start
        beta = beta.start
        diff = thres + 1
        step = 0


        while(diff > thres & step < max.step){
                step = step + 1
                opt1 = stats::optim(alpha, neg.log.likelihood.alpha,
                                    control = list(maxit = max.step))
                diff1 = Diff(opt1$par, alpha)
                alpha = opt1$par

                opt2 = stats::optim(beta, neg.log.likelihood.beta, control = list(maxit = max.step))

                diff2 = Diff(opt2$par, beta)
                beta = opt2$par

                diff = max(diff1, diff2)
        }
        # have hessian matrix
        hessian.mat = try(optimHess(c(alpha,beta), neg.log.likelihood), silent = TRUE)

        if ("try-error" %in% class(hessian.mat)) {
                return()
                alpha.cov = matrix(NA, pa+pb, pa+pb)
        }

        cov = solve(hessian.mat)
        point.est = c(alpha, beta)

        likelihood = exp(-neg.log.likelihood(point.est))
        opt = orgEst(point.est = point.est,
                     cov = cov,
                     type = "monotone",
                     name = c(paste("alpha", 1:ncol(va), sep = ""), paste("beta", 1:ncol(vb), sep = "")),
                     va=va,
                     vb=vb,
                     coverged = step < max.step,
                     likelihood = likelihood)

        return(opt)
}
