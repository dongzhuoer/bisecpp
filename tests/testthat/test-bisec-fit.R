testthat::context('Testing bisec0-fit.R')
setwd(here::here(''))  # workspace is reset per file

testthat::test_that('R_square()', {
	testthat::expect_equal(R_square(data1, f_protein, ideal.para11, extra11), 1)
	testthat::expect_equal(R_square(data1, f_protein, ideal.para11, extra12), 0.577768039)
	testthat::expect_equal(R_square(data2, f_AI2_out, ideal.para2, extra2), 1)
})


testthat::test_that('assess()', {
	set.seed(0)
	testthat::expect_equal(assess(ideal.para11, ideal.para11 * (1 + 0.1 * runif(3))), 1.05078617)
})


testthat::test_that('bisec_fit()', {
	result <- list(
		value = 0.000144774101699885,
		para = c(K.i = 0.8859375, K.nat = 0.0984375, k.d = 0.0196875),
		space = matrix(
			c(0.84375, 0.09375, 0.01875, 0.928125, 0.103125, 0.020625),
			nrow = 3, dimnames = list(c("K.i", "K.nat", "k.d"), NULL)
		)
	)

	testthat::expect_equal(bisec_fit(data1, f_protein, extra11, space1, 2, 5, 1/3, 1), result)
	testthat::expect_equal(bisec_fit(data1, f_protein, extra11, space1, 2, 5, 4, 1), result)
})

