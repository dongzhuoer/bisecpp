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
#'     partition_fit_impl(data1, f_protein, space1, extra11, 2, 1, 1);
#'     partition_fit_impl(data2, f_AI2_out, space2, extra2, 2, 1, 1);
#' }
#'
#' @section to do:
#'     1. support multivariate function
bisec_fit <- function(data, ff, extra, space, partition, times, trim, enlarge) {
	df <- nrow(space[,,1]);
	n <- partition ^ df;
	space.len <- dim(space)[3];
	range <- space[,,1];
	sub.range <- range;
	para <- rowMeans(range);
	vec <- integer(df)
	middle <- matrix(NA, n * space.len, 3);
	reserve <- 1;
	if (trim > 0) {
		if (trim >= 1)
			reserve = trim
		else
			reserve = ceiling(n ^ trim);
	}
	result <- vector('list', reserve)
	i = 1;

	for (irange in seq(space.len)) {
		range = space[,,irange]
		for (ipartition in seq(0, n - 1)) {
			vec <- base(ipartition, partition, df);
			para = range[,1] + (vec + 1 / 2) / partition * (range[,2] - range[,1])
			middle[i, ] = c(R_square(data, ff, para, extra), ipartition, irange);
			i = i + 1;
		}
	}
	middle = middle[order(middle[,1], decreasing = T),,drop = FALSE];
	middle = middle[seq(reserve),, drop = FALSE];

	for (j in seq_along(middle[,1])) {
		range = space[,,middle[j,3]];
		vec = base(middle[j,2], partition, df);
		sub.range[,1] = range[,1] + vec * (range[,2] - range[,1]) / partition;
		sub.range[,2] = sub.range[,1] + (range[,2] - range[,1]) / partition;
		para = rowMeans(sub.range);
		sub.range[,1] = para - (para - sub.range[,1]) * enlarge;
		sub.range[,2] = para - (para - sub.range[,2]) * enlarge;
		result[[j]] = list(middle[j,1], sub.range)
	}
	result
}












