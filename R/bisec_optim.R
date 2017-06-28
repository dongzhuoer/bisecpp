#' @title decimal to base-n
#'
#' @description transform decimal integer to other number mode, such as binary
#'   (`radix = 2`)
#'
#' @param x integer scalar. The number to transform
#' @param radix integer scalar. radix of the target number, \eqn{2} for binary,
#'   \eqn{8} for octal, etc.
#' @param n.digit integer scalar. how many digits the result contains, unused
#'   uppper digits are filled by \eqn{0}
#'
#' @return integer. each scalar of the vector represent one digit of the
#'   transformed number
#'
#' @examples
#' base(443L, 3L, 6L)
base <- function(x, radix, n.digit) {
	result <- integer(n.digit);
	i = 1L;
	while (x > 0) {
		result[i] = x %% radix;
		x = x %/% radix;
		i = i + 1;
	}
	rev(result);
}


#' @param spaces numeric array. `spaces[ , , k]` should
#'   give a `space`, see `bisec_optim` below
biosec_optim_impl <- function(fun, spaces, partition, trim, enlarge) {
	df <- nrow(spaces[ , , 1]);		# how many parameters
	n <- partition ^ df;			# how many subspaces per space
	nspace <- dim(spaces)[3];		# how many spaces

	range <- spaces[ , , 1];
	sub.range <- range;
	para <- rowMeans(range);
	vec <- integer(df)
	middle <- matrix(NA, n * nspace, 3);
	reserve <- 1;
	if (trim > 0) {
		if (trim >= 1)
			reserve = trim
		else
			reserve = ceiling(n ^ trim);
	}
	result <- vector('list', reserve)
	i = 1;

	for (irange in seq(nspace)) {
		range = spaces[,,irange]
		for (ipartition in seq(0, n - 1)) {
			vec <- base(ipartition, partition, df);
			para = range[,1] + (vec + 1 / 2) / partition * (range[,2] - range[,1])
			middle[i, ] = c(fun(para), ipartition, irange);
			i = i + 1;
		}
	}
	middle = middle[order(middle[,1]),,drop = FALSE];
	middle = middle[seq(reserve),, drop = FALSE];

	for (j in seq_along(middle[,1])) {
		range = spaces[,,middle[j,3]];
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



#' @title General-purpose Optimization (Minimization)
#'
#' @description `bisec_optim()` found parameters which minimize target function in a given space
#'
#' @details
#'
#' @param fun function. the function to be optimized, i.e. find minimal value. `fun(c(para1, para2, ..., paran))` should return a numeric scalar where `para` is a numeric vector, such as
#' @param space numeric matrix. \eqn{n*2} matrix, where \eqn{n} is the number of
#'   parameters to fit. `space[i, ]` should give a the range of ith parameter, i.e. `c(para1.min, para2.max)`. `space` specifies a n-dimension space, this space is uniformly partitioned into subspaces, and those subspaces which have the highest score is reserved for further partitioning. Finally we get the best space, and the centre point of this space contains the parameters we want.
#' @param partition integer scalar. how partitions are the original range
#'   divided per parameter.
#' @param times integer scalar. how many times should the initial range be partitioned.
#' @param trim integer scalar. determine what proportion of total subranges will
#'   be reserved for the next round of partition.
#' @param enlarge integer scalar. how many times is the selected subrange
#'   enlarged before the next round of partition.
#'
#' @return numeric vector. The best parameters found.
#' @export
#'
#' @examples
#'
#' @section to do:
#'     1. partition can be a integer vector, i.e. different parameter can be partitioned for different times
bisec_optim <- function(fun, space, partition, times, trim, enlarge) {
	range <- space[,,1];
	result <- NULL;

	for (i in seq(times)) {
		result <- biosec_optim_impl(fun, space, partition, trim, enlarge);
		space <- vapply(result, function(x){x[[2]]},range);
		gc();
	}

	result <- result[[1]];
	result[[2]] = rowMeans(result[[2]]);
	result;
}