#' @title get R square from data and fit function
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



#' @title assess how reality accord with the ideal
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



#' @title General-purpose fit
#'
#' @description you can almost fit any function you want
#'
#' @details even complex coupled differential equation that you can't get a analytic solution is ok.
#'
#' @param data numeric matrix. `data[i, 1]` gives \eqn{xi} and `data[i, 2]` gives
#'   \eqn{yi}
#' @param ff function. `ff(para, extra)` gives the function `f()` that
#'   \eqn{y=f(x)} where the range of para make a `space`.
#' @param extra numeric vector. see `ff`
#' @param space,partition,times,trim,enlarge see [bisec_optim]
#'
#' @return list of list.
#' 1. information about subrange1
#'     1. \eqn{R^2} using median of subrange1 as `para` for `ff`
#'     2. subrange1
#' 2. information about subrange2
#'     1. \eqn{R^2} using median of subrange2 as `para` for `ff`
#'     2. subrange1
#' 3. ...
#'
#' @export
#'
#' @seealso [biosec_optim]
#'
#' @examples
#' \donotrun {
#'     bisec_fit(data1, f_protein, extra11, space1, 2, 5, 1/3, 1);
#'     bisec_fit(data2, f_AI2_out, extra2, space2, 2, 5, 10, 1);
#' }
#'
#' @section to do:
#'     1. support multivariate function
bisec_fit <- function(data, ff, extra, space, partition, times, trim, enlarge) {
	fun <- function(para) {
		abs(1 - R_square(data, ff, para, extra));
	}

	bisec_optim(fun, space, partition, times, trim, enlarge)
}












