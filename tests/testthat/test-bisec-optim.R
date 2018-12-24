testthat::context('Testing bisec-optim.R')
setwd(here::here(''))  # workspace is reset per file

testthat::test_that('base()', {
	testthat::expect_identical(base(443L, 3L, 6L), c(1L, 2L, 1L, 1L, 0L, 2L))
})


testthat::test_that('bisec_optim()', {
	result <- list(
		value = 2.54026317126454e-05,
		para = c(8.46754390421514e-06, 8.46754390421514e-06, 8.46754390421514e-06),
		space = matrix(
			c(0, 0, 0, 1.69350878084303e-05, 1.69350878084303e-05, 1.69350878084303e-05),
			nrow = 3
		)
	)

	testthat::expect_equal(
		bisec_optim(sum, matrix(c(0,1,0,1,0,1), 3, byrow = T), 3, 10, 4, 1),
		result
	)
	testthat::expect_equal(
		bisec_optim(sum, matrix(c(0,1,0,1,0,1), 3, byrow = T), 3, 10, 1/3, 1),
		result
	)
})

