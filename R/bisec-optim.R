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
#' \dontrun{
#'     base(443L, 3L, 6L)
#' }
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



#' @title working horse for [biosec_optim]
#'
#' @param spaces numeric array. `spaces[ , , k]` should
#'   give a `space`, see `bisec_optim` below
#'
#' @return list. list of `list(fun.val, subsapce)` , first is best
#' 
#' @keywords internal
biosec_optim_impl <- function(fun, spaces, partition, trim, enlarge) {
	# perpare some variables which storage number
	df <- nrow(spaces[ , , 1]);		# how many parameters
	nsubspace <- partition ^ df;	# how many subspaces per space
	nspace <- dim(spaces)[3];		# how many spaces

	# adjust reserve according to trim
	if (trim >= 1)
		reserve <- as.integer(trim)
	else
		reserve <- ceiling(nsubspace^trim);

	# for each space, partition it into subspaces and eval `fun` with central point of each subspace
	central <- matrix(NA, nsubspace * nspace, 3);
	 		# storage information about the central point of all subsaces;
	for (ispace in seq(nspace)) {
		space <- spaces[ , , ispace];

		for (isubspace in seq(0, nsubspace - 1)) {
			vec <- base(isubspace, partition, df);
			para = space[ , 1] + (vec + 1 / 2) / partition * (space[ , 2] - space[ , 1]);
			central[isubspace + 1 + (ispace - 1)*nsubspace, ] = c(fun(para), isubspace, ispace);
		}
	}

	# sort and reserve best central point
	central = central[order(central[ , 1]),,drop = FALSE];
	central = central[seq(reserve), , drop = FALSE];

	# for earch central point, get the subspace it belongs and enlarge it (optional).
			# we only do this for reserved central points so we can save times.
	result <- vector('list', reserve);
	for (j in seq_along(central[ , 1])) {
		space <- spaces[ , , central[j, 3]];
		subspace <- space;

		vec <- base(central[j,2], partition, df);
		subspace[ , 1] = space[ , 1] + vec * (space[ ,2] - space[ ,1]) / partition;
		subspace[ , 2] = subspace[ , 1] + (space[ , 2] - space[ ,1]) / partition;

		para = rowMeans(subspace);
		subspace[ , 1] = para - (para - subspace[ ,1]) * enlarge;
		subspace[ , 2] = para - (para - subspace[ ,2]) * enlarge;

		result[[j]] = list(central[j,1], subspace);
	}

	# return
	result;
}



#' @title General-purpose Optimization (Minimization)
#'
#' @description `bisec_optim()` found parameters which minimize target function in a given space
#'
#' @details see paras.
#'
#' @param fun function. the function to be optimized, i.e. find minimal value. `fun(c(para1, para2, ..., paran))` should return a numeric scalar where `para` is a numeric vector, such as
#' @param space numeric matrix. \eqn{n*2} matrix, where \eqn{n} is the number of
#'   parameters to fit. `space[i, ]` should give a the range of ith parameter, i.e. `c(para1.min, para2.max)`. `space` specifies a n-dimension space, this space is uniformly partitioned into subspaces, and those subspaces which have the highest score is reserved for further partitioning. Finally we get the best space, and the centre point of this space contains the parameters we want.
#' @param partition integer scalar. how partitions are the original range
#'   divided per parameter.
#' @param times integer scalar. how many times should the initial range be partitioned.
#' @param trim numeric scalar. must be positive, determine what proportion of total subranges will
#'   be reserved for the next round of partition. Let `n` be the number of subspaces per space, if `trim` >= 1, `as.integer(trim)`; otherwise, `ceiling(n^trim)`.
#' @param enlarge integer scalar. how many times is the selected subrange
#'   enlarged before the next round of partition.
#'
#' @return list.
#' 1. value, numeric saclar. smallest function value found
#' 2. para, numeric vector. best parameters.
#' 3. space, numeric matrix. best solution space.
#' @export
#'
#' @examples
#' bisec_optim(sum, matrix(c(0,1,0,1,0,1), 3, byrow = T), 3, 10, 4, 1)
#' bisec_optim(sum, matrix(c(0,1,0,1,0,1), 3, byrow = T), 3, 10, 1/3, 1)
#'
#' @section to do:
#'     1. partition can be a integer vector, i.e. different parameter can be partitioned for different times
#'
bisec_optim <- function(fun, space, partition, times, trim, enlarge) {
	# transform `space` to `spaces` which contains only one space, for consistency in the following loop
	spaces <- array(space, dim = c(dim(space), 1));

	# loop `times` times
	for (i in seq(times)) {
		result <- biosec_optim_impl(fun, spaces, partition, trim, enlarge);
		spaces <- vapply(result, function(x) {x[[2]]}, space);
	}

	# make result
	best <- result[[1]];			# choose best `list(fun.val, subsapce)`
	list(value = best[[1]], para = rowMeans(best[[2]]), space = best[[2]])
}
