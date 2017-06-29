ideal.para11 <- c(K.i = 0.9, K.nat = 0.1, k.d = 0.02);        #protein (with IPTG)
ideal.para12 <- c(K.i = 0, K.nat = 0.1, k.d = 0.02);          #protein (without IPTG)
ideal.para2 <- c(k.in = 0.008, k.out = 0.045, k.p = 0.006);   #AI-2
devtools::use_data(ideal.para11, ideal.para12, ideal.para2, overwrite = TRUE)

extra11 <- c(mu = 0.0044, iOD600 = 0.5);                      #cell (with ACDB)
extra12 <- c(mu = 0.0056, iOD600 = 0.5);                      #cell (others)
extra2 <- c(IPTG.ACDB = TRUE, IPTG.K = TRUE, AI2.out.0 = 1); #sundries
devtools::use_data(extra11, extra12, extra2, overwrite = TRUE);

space1 <- cbind(integer(length(ideal.para11)), ideal.para11 * 3);
space2 <- cbind(integer(length(ideal.para2)), ideal.para2 * 3);
devtools::use_data(space1, space2, overwrite = TRUE)



#' Title
#'
#' @param para numeric. unknown parameters
#' @param extra numeric, given parameters
#'
#' @return function. returns \eqn{y} when giving \eqn{x} as argument
#'
#' @export
#'
#' @family AI-2 case functions
#'
#' @examples
#' donotrun {
#'     protein_ACDB <- f_protein(ideal.para11, extra11) #protein with ACDB gene and IPTG
#'     protein_K <- f_protein(ideal.para11, extra12) #protein with IPTG
#'     protein <- f_protein(ideal.para12, extra12) #protein with nothing
#' }
f_protein <- function(para, extra) {
	function(x) {
		extra[[2]]*exp(-para[[3]]*x)*(exp((extra[[1]] + para[[3]])*x) - 1)*(para[[1]] + para[[2]])/(extra[[1]] + para[[3]])
	}
}

protein_ACDB <- f_protein(ideal.para11, extra11) #protein with ACDB gene and IPTG
protein_K <- f_protein(ideal.para11, extra12) #protein with IPTG
protein <- f_protein(ideal.para12, extra12) #protein with nothing
devtools::use_data(protein_ACDB, protein_K, protein, overwrite = TRUE);



#' Title
#'
#' @inheritParams f_protein
#' @inherit f_protein return
#'
#' @export
#'
#' @family AI-2 case functions
#'
#' @examples
#' f_AI2_out(ideal.para2, extra2)(seq(0, 270, 30))
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


#' @include simulate-fit.R
data1 <- simulate_data(f_protein, ideal.para11, extra11, seq(60,240,60));
data2 <- simulate_data(f_AI2_out, ideal.para2, extra2, seq(0,270,30));
devtools::use_data(data1, data2, overwrite = TRUE);
