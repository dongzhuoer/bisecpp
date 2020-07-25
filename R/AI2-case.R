#' @title function for protein concentration
#'
#' @param para numeric. unknown parameters
#' @param extra numeric, given parameters
#'
#' @return function. returns \eqn{y} when giving \eqn{x} as argument
#'
#' @family AI-2 case functions
#'
#' @examples
#' NULL
#'
#' @export

f_protein <- function(para, extra) {
	function(x) {
		extra[[2]]*exp(-para[[3]]*x)*(exp((extra[[1]] + para[[3]])*x) - 1)*(para[[1]] + para[[2]])/(extra[[1]] + para[[3]])
	}
}


#' @title function for extracelluar AI2 concentration
#'
#' @inheritParams f_protein
#' @inherit f_protein return
#'
#' @family AI-2 case functions
#'
#' @examples
#' NULL
#'
#' @export
f_AI2_out <- function(para, extra) {
	k.in <- para[1];
	k.out <- para[2];
	k.p <- para[3];
	n = 270L;
	LsrACDB <- (if (extra[1]) protein_ACDB else if (extra[2]) protein_K else protein)(1:n);
	LsrK <- LsrACDB;
	AI2.in <- integer(n + 1);
	AI2.out <- integer(n + 1);
	AI2.out[1] = extra[3];
	for (i in seq(2, n + 1)) {
		AI2.in[i] = AI2.in[i - 1] + k.in * LsrACDB[i - 1] * AI2.out[i - 1] - k.out * AI2.in[i - 1] - k.p * LsrK[i - 1] * AI2.in[i - 1];
		AI2.out[i] = AI2.out[i - 1] - k.in * LsrACDB[i - 1] * AI2.out[i - 1] + k.out * AI2.in[i - 1];
		i = i + 1;
	}
	function(i) AI2.out[i + 1]
}