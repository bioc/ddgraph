library(ddgraph)
library(testthat)

bin = makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "b"=c(1, 1, 0, 0, 0, 1), 
                         "c"=c(0, 0, 0, 1, 1, 1)), "binary example", c(0, 0, 0, 1, 1, 1))
cont = makeDDDataSet(cbind("a"=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 
                          "b"=c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6)), 
                          "continuous example", c(0, 0, 0, 1, 1, 1))

test_that("DDDataSet constructor", {
	expect_equal(dataType(bin), "binary")
	expect_equal(bin$a, c(0, 0, 1, 0, 1, 0))
	expect_equal(bin$a, bin[["a"]])
	expect_equal(bin$a, bin[,1])
	expect_equal(bin$class, c(0, 0, 0, 1, 1, 1))
	
	expect_equal(dataType(cont), "continuous")
	expect_equal(cont$a, c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
})

bin.res = ciTest(bin, "a", "b", "c", test.type="x2")

test_that("DDDataSet::ciTest", {
	# verify if the right backend function is called, and error produced
	expect_equal(ciTest(bin, "a", "b", test.type="x2")$pValue, 
              bnlearn::ci.test(as.factor(bin$a), as.factor(bin$b), test="x2")$p.val)
	expect_equal(ciTest(bin, "a", "b", "c", test.type="x2")$pValue, 
              bnlearn::ci.test(as.factor(bin$a), as.factor(bin$b), as.factor(bin$c), test="x2")$p.val)
	expect_error(ciTest(bin, "a", "b", test.type="zf"))
	expect_equal(ciTest(bin, "a", "b", "class", test.type="x2")$pValue, 
              bnlearn::ci.test(as.factor(bin$a), as.factor(bin$b), as.factor(bin$class), test="x2")$p.val)
	
	expect_equal(ciTest(cont, "a", "b", test.type="zf")$pValue, 
              bnlearn::ci.test(cont$a, cont$b, test="zf")$p.val)
	expect_error(ciTest(cont, "a", "b", test.type="x2"))
	expect_equal(ciTest(cont, "a", "b", "class", test.type="zf")$pValue, 
              bnlearn::ci.test(cont$a, cont$b, cont$class, test="zf")$p.val)
	
	# verify that the return structure has all of stuff filled-out correctly
	expect_equal(bin.res$targetInx, 1)
	expect_equal(bin.res$sourceInx, 2)
	
	expect_equal(bin.res$targetName, "a")
	expect_equal(bin.res$sourceName, "b")
	
	expect_equal(bin.res$condSetInx, 3)
	expect_equal(bin.res$condSetName, "c")
	
	expect_equal(bin.res$testType, "x2")
	
	expect_equal(bin.res$pValue, bin.res[["pValue"]])
	
})

# generate some data to test extractCITestResult
res1 = ciTest(bin, "a", "b", "c", test.type="x2")
res2 = ciTest(bin, "a", "b", test.type="x2")

ciList = list()
ciList[[1]] = list(res1, res2)
ciList[[2]] = list(res1)

test_that("extractCITestResultProperty", {
	expect_equal(ddgraph:::extractCITestResultProperty(ciList, "pValue")[[1]], c(res1$pValue, res2$pValue))
	expect_equal(ddgraph:::extractCITestResultProperty(ciList, "pValue")[[2]], c(res1$pValue))
	expect_equal(length(ddgraph:::extractCITestResultProperty(ciList, "pValue")), 2)
	
	expect_equal(ddgraph:::extractCITestResultProperty(ciList, "condSetName")[[1]], 
              list("c", vector(mode="character")))
	expect_equal(ddgraph:::extractCITestResultProperty(ciList, "condSetInx")[[1]], 
              list(3, vector(mode="numeric")))
})

# the the pValueAfterMultipleTesting function
adjC.pvals.at.n = list(`1`=c(0.001, 0.003, 0.04, 0.5))

test_that("pValueAfterMultipleTesting", {
	expect_equal(ddgraph:::pValueAfterMultipleTesting(res1, 3, adjC.pvals.at.n, "none"), 
              res1$pValue)
	expect_equal(ddgraph:::pValueAfterMultipleTesting(res1, 3, adjC.pvals.at.n, "fdr"), 
              p.adjust(c(0.001, 0.003, res1$pValue, 0.5), "fdr")[3])
	expect_true(ddgraph:::pValueAfterMultipleTesting(res1, 3, adjC.pvals.at.n, "none") 
             != ddgraph:::pValueAfterMultipleTesting(res1, 3, adjC.pvals.at.n, "fdr"))
})
 
bin.var = makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "b"=c(1, 1, 1, 1, 1, 1), 
                       "c"=c(0, 0, 0, 1, 1, 1)), "binary example", c(0, 0, 0, 1, 1, 1),
                        removeZeroVar=T)

bin.names = makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "b"=c(1, 1, 0, 0, 0, 1), 
                        "c"=c(0, 0, 0, 1, 1, 1), "d"=c(0, 0, 0, 1, 1, 1)), "binary example",
                        classLabelsCol="d")

bin.names2 = makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "d"=c(0, 0, 0, 1, 1, 1),
                               "b"=c(1, 1, 0, 0, 0, 1), "c"=c(0, 0, 0, 1, 1, 1)), "binary example",
                        classLabelsCol="d")

bin.names3 = makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "d"=c(0, 0, 0, 1, 1, 1),
                               "b"=c(1, 1, 0, 0, 0, 1), "c"=c(0, 0, 0, 1, 1, 1)), "binary example",
                        classLabelsCol=2)

test_that("makeDDDataSet consistency checks", {
  expect_error(makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "b"=c(1, 1, 1, 1, 1, 1), 
                         "c"=c(0, 0, 0, 1, 1, 1)), "binary example", c(0, 0, 0, 1, 1, 1)))
                         
  expect_equal(names(bin.var), c("a", "c", "class"))
  expect_true(all(rawData(bin.names) == rawData(bin)))
  expect_true(all(rawData(bin.names) == rawData(bin.names2)))
  expect_true(all(rawData(bin.names) == rawData(bin.names3)))
  
  expect_error(makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "b"=c(0, 0, 1, 0, 1, 0), 
                         "c"=c(0, 0, 1, 0, 1, 0)), "binary example", c(0, 0, 0, 1, 1, 1)))
                         
  expect_error(makeDDDataSet(cbind("a"=c(0, 0, 1, 0, 1, 0), "a"=c(0, 1, 1, 1, 1, 1), 
                         "c"=c(0, 0, 0, 1, 1, 1)), "binary example", c(0, 0, 0, 1, 1, 1)))

})

