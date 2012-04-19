library(ddgraph)
library(testthat)

ic.formula.mul.test = function(){	
	data = data.frame(gtools::permutations(2,3, repeats.allowed=T)-1)
	names(data) = c("a", "b", "c")
	
	target.formula = "a | (!a & b & !c)"
	nt = independent.contributions.formula.mul(data, c("a", "b", "c"), rep(1,3), rep(0,3), target.formula)
	
	test_that("IC multiple formula with 3 variables complex boolean function", {
		expect_equal(nt, c(0, 0, 1, 0, 1, 1, 1, 1))
	})
}

ic.formula.mul.test()

test_that("is.binary", {
	expect_true(is.binary(c(0, 0, 1, 1, 0, 1)))
	expect_true(is.binary(c(0, 0)))
	expect_true(is.binary(0))
	expect_false(is.binary(3))
	expect_false(is.binary(c(1,3,4,3)))
	expect_true(is.binary( data.frame(a=c(0,1,1), b=c(0,0,0)) ))
	expect_false(is.binary( data.frame(a=c(0,1,1), b=c(0,4,0)) ))
	expect_false(is.binary( data.frame(a=c(T,F,T), b=c(T,T,T)) ))
	expect_false(is.binary( cbind(c(0,1,1), c(0,4,0)) ))
	expect_false(is.binary( cbind(c(T,F,T), c(T,T,T)) ))
})

test_that("convertToFactor", {
	expect_true(is.factor(convertToFactor(c(0,1,1))))
	expect_true(is.factor(convertToFactor( matrix(c(0,1,1,0,0,1), nrow=2))[,1]))
	expect_true(is.factor(convertToFactor( data.frame(a=c(0,1,1), b=c(0,0,1)))[,1]))
	expect_true(is.null(convertToFactor(NULL)))
})


