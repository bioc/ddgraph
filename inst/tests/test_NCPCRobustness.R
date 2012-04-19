library(ddgraph)
library(testthat)

# test if the calculation of various robustness stats is correct

bin = makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "b"=c(1, 1, 0, 0, 0, 1), "c"=c(0, 0, 0, 1, 1, 1), 
	"d"=c(0, 1, 1, 0, 0, 0)), "binary example", c(0, 0, 0, 1, 1, 1))

raw = list(
	list(direct=c("a", "b"), joint=c(), indirect=c("c"), directAndJoint=c("a", "b"), jointIfNotDirect=c("a", "b"), positiveClassSize=3),
	list(direct=c("a"), joint=c("b"), indirect=c("c"), directAndJoint=c("a", "b"), jointIfNotDirect=c("a"), positiveClassSize=4),
	list(direct=c("a"), joint=c("b", "c"), indirect=c("d"), directAndJoint=c("a", "b", "c"), jointIfNotDirect=c("a"), positiveClassSize=5),
	list(direct=c("a"), joint=c("b"), indirect=c(), directAndJoint=c("a", "b"), jointIfNotDirect=c("a"), positiveClassSize=4)
)
r = ddgraph:::makeNCPCRobustness(bin, raw, list())
r = ddgraph:::calculateNCPCRobustnessStats(r)
freq.table = cbind(c(1, 0.25, 0), c(0, 0.75, 0.25), c(0, 0, 0.5))

test_that("calculateNCPCRobustnessStats", {
	expect_equal(r@runs, 4)
	
	expect_equal(as.matrix(r@tables$direct)[,1], c("a"=4, "b"=1))
	expect_equal(as.matrix(r@tables$joint)[,1], c("b"=3, "c"=1))
	expect_equal(as.matrix(r@tables$indirect)[,1], c("c"=2, "d"=1))	
	expect_equal(as.matrix(r@tables$directAndJoint)[,1], c("a"=4, "b"=4, "c"=1))
	expect_equal(as.matrix(r@tables$jointIfNotDirect)[,1], c("a"=4, "b"=1))
	expect_equal(as.matrix(r@tables$positiveClassSize)[,1], c("3"=1, "4"=2, "5"=1))
	
	rownames(freq.table) = c("a", "b", "c")
	colnames(freq.table) = c("direct", "joint", "indirect")
	expect_equal(as.matrix(r@enriched.pss[1:3, 1:3]), freq.table)
	
	expect_equal(as.character(r@enriched.pss[,4]), c("direct", "joint", "indirect"))
	expect_equal(r@enriched.pss[2,5], -0.25*log2(0.25)-0.75*log2(0.75))
	expect_equal(r@enriched.pss[,6], c(NA, NA, 3))
	expect_equal(r@enriched.pss[,7], c(0, 0, -0.25*log2(0.25)-0.75*log2(0.75)))
	
	expect_true(all(r@enriched.ps[,1] == freq.table[,1]+freq.table[,2]))
	expect_true(all(r@enriched.ps[,2] == freq.table[,3]))
	expect_equal(as.character(r@enriched.ps[,3]), c("direct and joint", "direct and joint", "indirect"))
	expect_true(all(r@enriched.ps[3,4] == c(2)))
	
	expect_true(all(r@not.enriched[1,] == c(0.25,  0.75, 3, -0.25*log2(0.25)-0.75*log2(0.75))))
	
	expect_equal(r@final.calls$name, c("a", "b", "c", "d"))
	expect_equal(r@final.calls$type, c("direct", "joint", "indirect", "no dependence"))
	expect_equal(r@final.calls$probability, c(1, 0.75, 0.5, 0.75))
})

