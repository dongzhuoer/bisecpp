#' @title generate simulated data from given model
#'
#' @param x numeric. used to generate `y` from `ff(para, extra)(x))`
#'
#' @return numeric matrix. `cbind(x, y)`
#'
#' @family simulate functions
#'
#' @examples
#' simulate_data(f_protein, ideal.para11, extra11, seq(60,240,60))
#' simulate_data(f_AI2_out, ideal.para2, extra2, seq(0,270,30))
#'
#' @export
simulate_data <- function(ff, para, extra, x) {
	cbind(x, y = ff(para, extra)(x));
}



#' @title use simulated data to evaluate fitting function
#'
#' @param x numeric. see [simulate_data]
#' @param ff function. see [simulate_data]
#' @param ideal.para numeric. see [simulate_data]. ideal value for `para`
#' @param times integer scalar. see [partition_fit]
#' @param scale integer scalar. how many times is `ideal.para` enlarged to
#'   generate the initial range
#' @param details logical scalar. whether result should include extra
#'   information.
#'
#' @return vector. see examples for details.
#'
#' @family simulate functions
#'
#' @examples
#' simulate_partition_fit(seq(60,240,60), f_protein, ideal.para11 , extra11, 10,3,1/3, 1, 3)
#' simulate_partition_fit(seq(0,270,30), f_AI2_out, ideal.para2 , extra2, 3,1,1/3, 1, 3)
#'
#' @export
simulate_partition_fit <- function(x, ff, ideal.para, extra, patition, times, trim, enlarge, scale, details = TRUE) {
	data <- simulate_data(ff, ideal.para, extra, x);
	range <- cbind(integer(length(ideal.para)),ideal.para*scale);
	start <- Sys.time();
	result <- bisec_fit(data, ff, extra, range, patition, times, trim, enlarge);
	result = c(time = as.numeric(Sys.time() - start),
			   assess = assess(result[[2]], ideal.para),
			   RSquare = 1 - result[[1]],
			   len = patition,
			   times = times,
			   trim = trim,
			   enlarge = enlarge,
			   scale = scale
	)
	if (details) c(result, ideal.para, extra) else result;
}
