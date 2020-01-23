#' Organize estimation result
#'
#' @noRd

orgEst = function(point.est, cov, type, name, va, vb, coverged, likelihood) {

        if (any(is.na(test))) {
                se.est = conf.lower = conf.upper = p.value = rep(NA, ncol(va) + ncol(vb))

        } else {
                se.est = sqrt(diag(cov))
                conf.lower = point.est + stats::qnorm(0.025) * se.est
                conf.upper = point.est + stats::qnorm(0.975) * se.est

                p.temp = stats::pnorm(point.est/se.est, 0, 1)
                p.value = 2 * pmin(p.temp, 1 - p.temp)
        }


        if (type == "monotone") {
                names(point.est) = names(se.est) = rownames(cov) = colnames(cov) = names(conf.lower) = names(conf.upper) = names(p.value) = name
                coefficients = cbind(point.est, se.est, conf.lower, conf.upper,
                                     p.value)
                linear.predictors = va %*% point.est[1:ncol(va)]
                param.est = exp(linear.predictors)

                output = list(type = type,
                              point.est = point.est,
                              se.est = se.est,
                              cov = cov,
                              conf.lower = conf.lower,
                              conf.upper = conf.upper,
                              p.value = p.value,
                              coefficients = coefficients,
                              param.est = param.est,
                              converged = converged,
                              likelihood = likelihood)

                class(output) = c("mem", type, "list")

        }


}
