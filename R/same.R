#' Decide if two values are equivalent
#'
#' @param x A numerical vector.
#' @param y A numerical vector, same dimension as \code{x}.
#' @param tolerance Difference tolerance, \code{tolerance = .Machine$double.eps^0.5}.
#' @return \code{TRUE} or \code{FALSE}.
#' @noRd

same = function (x, y, tolerance = .Machine$double.eps^0.5) {
        abs(x - y) < tolerance
}
