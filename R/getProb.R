#' Get the probabilities
#'
#' @description By having logarithm of relative risks and logarithm of generalized odds product, the probability of event from each each treatment is calculated.
#'
#' @param log.rr A vetor of logarithm of relative risks.
#' @param log.gop A number oflogarithm of generalized odds product.
#'
#' @return A vector of probabilities
#'
#' @examples getProb(c(1, 10, 5), 2)
#'
#' @export
#'


getProb = function(log.rr, log.gop){

        if (lenght(log.gop) > 1) {
                stop("the length of 'log.gop' should be one!")
        }

        k = length(log.rr) # levels

        log.rr_max = max(c(0,log.rr))

        rr_update =exp(c(0, log.rr) - log.rr_max)

        f = function(x){
                return((k+1) * log(x) + sum(log(rr_update)) - sum(log(1 - x * rr_update)) - log.gop)
        }

        tol = 10
        if(f(1-0.1^tol)<0) p.max = 1 else{
                p.max = uniroot(f, lower=0, upper=1, tol=10^(-tol))$root
        }
        return(p.max * rr_update)
}

