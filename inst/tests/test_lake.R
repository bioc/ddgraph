
# test case for the "Lake" example which illustrates the differences between DDGraph and DDGraph*

library(ddgraph)

library(testthat)

# fix the random seed
set.seed(311)

# the real graphs is as follows:
# Trees <- Lake
# Leaves <- Trees
# Fish <- Lake, Season
# People <- Season, Fish
# Lake, Season have no parents

data = matrix(0, ncol=6, nrow=5000)
colnames(data) = c("Fish", "Trees", "Leaves", "Season", "People", "Lake")
data = data.frame(data)

for(i in 1:nrow(data)){
	data$Season[i] = sample(c(0,1), 1)
	data$Lake[i] = sample(c(0,1), 1)
	if(data$Lake[i] == 1){
		data$Trees[i] = sample(c(0,1), 1, prob=c(0.2, 0.8))		
		if(data$Season[i] == 1)
			data$Fish[i] = sample(c(0,1), 1, prob=c(0.4, 0.6))
		else	
			data$Fish[i] = sample(c(0,1), 1, prob=c(0.1, 0.9))
	} else {
		data$Trees[i] = sample(c(0,1), 1, prob=c(0.6, 0.4))
		if(data$Season[i] == 1)
			data$Fish[i] = sample(c(0,1), 1, prob=c(0.6, 0.4))
		else
			data$Fish[i] = sample(c(0,1), 1, prob=c(0.1, 0.9))
	}
	
	if(data$Season[i] == 1 & data$Fish[i] == 1){
		data$People[i] = sample(c(0,1), 1, prob=c(0.1, 0.9))
	} else {
		data$People[i] = sample(c(0,1), 1, prob=c(0.9, 0.1))
	}
	
	if(data$Trees[i] == 1)
		data$Leaves[i] = sample(c(0,1), 1, prob=c(0.4, 0.6))
	else
		data$Leaves[i] = sample(c(0,1), 1, prob=c(0.9, 0.1))		
}

dd = data
names(dd)[ncol(dd)] = "class"

d = makeDDDataSet(dd, "Lake")

eg = ncpc(d, alpha=0.05, star=F, verbose=F, test.type="x2")
egs = ncpc(d, alpha=0.05, star=T, verbose=F, test.type="x2")

# only the DDGraph* can correctly identify only Fish and Trees are primary variables
# DDGraph (original) will have a false positive People even when there is infinite data

test_that("Lake example", {
	expect_equal(eg$direct, structure(c(1, 2, 5), names=c("Fish", "Trees", "People")))
	expect_equal(egs$direct, structure(c(1, 2), names=c("Fish", "Trees")))	
})


