#' Title
#'
#' @param ff function. `ff(para, extra)` gives the function `f()` that
#'   \eqn{y=f(x)}
#' @param para numeric. see \code{ff}
#' @param extra numeric. see \code{ff}
ff_para_extra <- function() {}



#' Title
#'
#' @param data numeric matrix. `data[,1]` gives \eqn{x} and `data[,2]` gives
#'   \eqn{y}
#' @inheritParams ff_para_extra
#'
#' @return numeric scalar. \eqn{R^2} (coefficient of determination) showing how
#'   well the function fits `data`
#'
#' @export
#'
#' @examples
#' \donotrun {
#'     R_square(data1, f_protein, ideal.para11, extra11)
#'     R_square(data1, f_protein, ideal.para11, extra12)
#'     R_square(data2, f_AI2_out, ideal.para2, extra2)
#' }
R_square <- function(data, ff, para, extra) {
    f <- ff(para, extra);
    x <- data[,1];
    y <- data[,2];
    1 - sum((y - f(x)) ^ 2) / sum((y - mean(y)) ^ 2);
}



#' Title
#'
#' @param reality numeric. real parameters
#' @param ideal numeric. ideal parameters
#'
#' @return numeric scalar. the relative similarity between `reality` and `ideal`
#'
#' @export
#'
#' @examples
#' \donotrun {
#'     assess(ideal.para11, ideal.para11 * (1 + 0.1 * runif(3)))
#' }
assess <- function(reality, ideal) {
    result <- reality/ideal;
    result = ifelse(result > 1, result, 1 / result);
    prod(result) ^ (1 / length(result));
}





