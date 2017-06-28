# to do: all sum((x-mean)^2) should be sum(x^2) - sum(x)^2 / n



#' Title
#'
#' @param data numeric matrix. `data[,1]` gives \eqn{x} and `data[,2]` gives
#'   \eqn{y}
#' @param ff function. `ff(para, extra)` gives the function `f()` that
#'   \eqn{y=f(x)}
#' @param space numeric matrix. \eqn{n*2*m} matirx, \eqn{n} is the number of
#'   parameters to fit, \eqn{m} is the number of ranges. `space[i,1,j]` and
#'   `space[i,2,j]` gives the min and max value of ith parameter in jth range
#'   respectively.
data_ff_space <- function() {}



#' Title
#'
#' @param extra usually numeric. see \code{ff}
#' @param partition integer scalar. how partitions are the original range
#'   divided per parameter.
extra_partition <- function() {}



#' Title
#'
#' @param trim integer scalar. determine what proportion of total subranges will
#'   be reserved for the next round of partition.
#' @param enlarge integer scalar. how many times is the selected subrange
#'   enlarged before the next round of partition.
trim_enlarge <- function() {}



#' Title
#'
#' @inheritParams data_ff_space
#' @inheritParams  extra_partition
#' @inheritParams trim_enlarge
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
#' @family partition fit functions
#'
#' @examples
#' \donotrun {
#'     partition_fit_impl(data1, f_protein, space1, extra11, 2, 1, 1);
#'     partition_fit_impl(data2, f_AI2_out, space2, extra2, 2, 1, 1);
#' }
partition_fit_impl <- function(data, ff, space, extra, partition, trim, enlarge) {
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



#' Title
#'
#' @inheritParams  extra_partition
#' @param times integer scalar. how many times should the initial range be
#'   partitioned.
#' @inheritParams trim_enlarge
extra_partition_times_trim_enlarge <- function() {}



#' Title
#'
#' @inheritParams data_ff_space
#' @inheritParams  extra_partition
#' @param times integer scalar. how many times should the initial range be partitioned.
#' @inheritParams trim_enlarge
#'
#' @return list.
#' 1. \eqn{R^2} using the final parameter `para` for `ff`
#' 2. final best fitted parameter
#'
#' @export
#'
#' @family partition fit functions
#'
#' @examples
#' \donotrun {
#'     partition_fit(data1, f_protein, space1, extra11, 2, 2, 1/2, 1);
#'     partition_fit(data2, f_AI2_out, space2, extra2, 3, 1, 1/3, 1);
#' }
partition_fit <- function(data, ff, space, extra, partition, times, trim, enlarge) {
	range <- space[,,1];
	result <- NULL;

	for (i in seq(times)) {
		result <- partition_fit_impl(data, ff, space, extra, partition, trim, enlarge);
		space <- vapply(result, function(x){x[[2]]},range);
		gc();
	}

	result <- result[[1]];
	result[[2]] = rowMeans(result[[2]]);
	result;
}



bisec_optim <- function(fun, space)