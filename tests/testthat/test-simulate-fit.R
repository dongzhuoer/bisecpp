testthat::context('Testing simulate-fit.R')
setwd(here::here(''))  # workspace is reset per file

testthat::test_that('simulate_data()', {
	testthat::expect_equal(
		simulate_data(f_protein, ideal.para11, extra11, seq(60,240,60)),
		matrix(
			c(60, 120, 180, 240, 20.510942303047, 32.8856534080412, 44.6820472621586, 58.7421889031268),
			ncol = 2, dimnames = list(NULL, c("x", "y"))
		)
	)

	testthat::expect_equal(
		simulate_data(f_AI2_out, ideal.para2, extra2, seq(0,270,30)),
		matrix(
			c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 1, 0.334298884487434,
			  0.038597082621433, 0.00157863290641664, 2.31624618400186e-05,
			  1.19488872914242e-07, 1.97233353895208e-10, 8.84708104129726e-14,
			  8.54174640118863e-18, 1.29719337541273e-22),
			ncol = 2, dimnames = list(NULL, c("x", "y"))
		)
	)
})


testthat::test_that('simulate_partition_fit()', {
	testthat::expect_equal(
		simulate_partition_fit(seq(60,240,60), f_protein, ideal.para11 , extra11, 10,3,1/3, 1, 3)[-1],
		c(assess = 1.0005, RSquare = 0.999999855696314, len = 10, times = 3,
		  trim = 0.333333333333333, enlarge = 1, scale = 3, K.i = 0.9,
		  K.nat = 0.1, k.d = 0.02, mu = 0.0044, iOD600 = 0.5)
	)
	testthat::expect_equal(
		simulate_partition_fit(seq(0,270,30), f_AI2_out, ideal.para2 , extra2, 3,1,1/3, 1, 3)[-1],
		c(assess = 1.77844665224503, RSquare = 0.999051495634497, len = 3,
		  times = 1, trim = 0.333333333333333, enlarge = 1, scale = 3,
		  k.in = 0.008, k.out = 0.045, k.p = 0.006, IPTG.ACDB = 1, IPTG.K = 1,
		  AI2.out.0 = 1)
	)
})

