library(ddgraph)
library(testthat)

# this test case calls system() so it will only work as intenden in linux
if(.Platform$OS.type == "unix"){
	fd = readFurlongData()

	test_that("readFurlongData data integrity", {
		expect_equal(ncol(signalMatrix(fd)), 15)
		expect_equal(nrow(signalMatrix(fd)), 310)
		expect_equal(names(fd), 
			c("tin_2.4", "tin_4.6", "tin_6.8", "bin_6.8", "bin_8.10",
			"bin_10.12", "twi_2.4", "twi_4.6", "twi_6.8", "bap_6.8",
			"mef2_2.4", "mef2_4.6", "mef2_6.8", "mef2_8.10", "mef2_10.12"))
	
		expect_true(is.factor(classLabels(fd)))
		expect_equal(levels(classLabels(fd)), c("neg", "Meso", "VM", "SM", "CM", "Meso_SM", "VM_SM"))	
		expect_equal(as.character(classLabels(fd)[1:10]), c("neg", "Meso", "neg", "VM", "neg", "SM", "CM", "VM_SM", "neg", "neg"))
	
	
	})

	meso = toDDDataSet(fd, "Meso")
	all.data = toDDDataSet(fd)
	tin_2.4.binary = signalMatrix(fd)[,"tin_2.4"]
	tin_2.4.binary[tin_2.4.binary>0] = 1

	test_that("DDDataSet creation from Furlong data", {
		expect_true(class(meso) == "DDDataSet")
		expect_true(class(all.data) == "list")
	
		expect_equal(dim(rawData(meso)), c(310,16))
		expect_equal(meso$class, as.numeric(classLabels(fd) == "Meso"))
		expect_equal(meso$tin_2.4, tin_2.4.binary)
	
		expect_equal(names(all.data), c("neg", "Meso", "VM", "SM", "CM", "Meso_SM", "VM_SM"))
		expect_equal(rawData(all.data$Meso), rawData(meso))
	
		expect_equal(dataType(meso), "binary")
		expect_equal(datasetName(meso), "Meso")
	})

	# run the inference tool using x2 as distribution so we always get the same results

	ncpc(all.data$Meso, alpha=0.05, test.type="x2", report.file="reports/meso_0.05.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$Meso, alpha=0.1, test.type="x2", report.file="reports/meso_0.1.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$Meso, alpha=0.2, test.type="x2", report.file="reports/meso_0.2.txt", verbose=FALSE, min.table.size=NULL)

	ncpc(all.data$VM, alpha=0.05, test.type="x2", report.file="reports/vm_0.05.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$VM, alpha=0.1, test.type="x2", report.file="reports/vm_0.1.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$VM, alpha=0.2, test.type="x2", report.file="reports/vm_0.2.txt", verbose=FALSE, min.table.size=NULL)

	ncpc(all.data$VM_SM, alpha=0.05, test.type="x2", report.file="reports/vm_sm_0.05.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$VM_SM, alpha=0.1, test.type="x2", report.file="reports/vm_sm_0.1.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$VM_SM, alpha=0.2, test.type="x2", report.file="reports/vm_sm_0.2.txt", verbose=FALSE, min.table.size=NULL)

	ncpc(all.data$Meso_SM, alpha=0.05, test.type="x2", report.file="reports/meso_sm_0.05.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$Meso_SM, alpha=0.1, test.type="x2", report.file="reports/meso_sm_0.1.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$Meso_SM, alpha=0.2, test.type="x2", report.file="reports/meso_sm_0.2.txt", verbose=FALSE, min.table.size=NULL)

	ncpc(all.data$SM, alpha=0.05, test.type="x2", report.file="reports/sm_0.05.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$SM, alpha=0.1, test.type="x2", report.file="reports/sm_0.1.txt", verbose=FALSE, min.table.size=NULL)
	ncpc(all.data$SM, alpha=0.2, test.type="x2", report.file="reports/sm_0.2.txt", verbose=FALSE, min.table.size=NULL)

	all.names = c("meso_0.05.txt", "meso_0.1.txt", "meso_0.2.txt", "vm_0.05.txt", "vm_0.1.txt", "vm_0.2.txt",
		"vm_sm_0.05.txt", "vm_sm_0.1.txt", "vm_sm_0.2.txt", "meso_sm_0.05.txt", "meso_sm_0.1.txt", "meso_sm_0.2.txt",
		"sm_0.05.txt", "sm_0.1.txt", "sm_0.2.txt")

	test_that("The report files correspond to gold standard reports", {
		for(name in all.names){
			cat("Checking", name, "\n")
			gold.standard.difference = system(paste("diff -u gold_standard/", name, " reports/", name, sep=""))
			expect_equal(gold.standard.difference, 0)		
		}
	})

	# run same with the star algorithm

	ncpc(all.data$Meso, alpha=0.05, test.type="x2", report.file="reports/star-meso_0.05.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$Meso, alpha=0.1, test.type="x2", report.file="reports/star-meso_0.1.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$Meso, alpha=0.2, test.type="x2", report.file="reports/star-meso_0.2.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)

	ncpc(all.data$VM, alpha=0.05, test.type="x2", report.file="reports/star-vm_0.05.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$VM, alpha=0.1, test.type="x2", report.file="reports/star-vm_0.1.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$VM, alpha=0.2, test.type="x2", report.file="reports/star-vm_0.2.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)

	ncpc(all.data$VM_SM, alpha=0.05, test.type="x2", report.file="reports/star-vm_sm_0.05.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$VM_SM, alpha=0.1, test.type="x2", report.file="reports/star-vm_sm_0.1.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$VM_SM, alpha=0.2, test.type="x2", report.file="reports/star-vm_sm_0.2.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)

	ncpc(all.data$Meso_SM, alpha=0.05, test.type="x2", report.file="reports/star-meso_sm_0.05.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$Meso_SM, alpha=0.1, test.type="x2", report.file="reports/star-meso_sm_0.1.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$Meso_SM, alpha=0.2, test.type="x2", report.file="reports/star-meso_sm_0.2.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)

	ncpc(all.data$SM, alpha=0.05, test.type="x2", report.file="reports/star-sm_0.05.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$SM, alpha=0.1, test.type="x2", report.file="reports/star-sm_0.1.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)
	ncpc(all.data$SM, alpha=0.2, test.type="x2", report.file="reports/star-sm_0.2.txt", verbose=FALSE, star=TRUE, min.table.size=NULL)

	all.files = paste("reports/", all.names, sep="")

	all.names = paste("star-", c("meso_0.05.txt", "meso_0.1.txt", "meso_0.2.txt", "vm_0.05.txt", "vm_0.1.txt", "vm_0.2.txt",
		"vm_sm_0.05.txt", "vm_sm_0.1.txt", "vm_sm_0.2.txt", "meso_sm_0.05.txt", "meso_sm_0.1.txt", "meso_sm_0.2.txt",
		"sm_0.05.txt", "sm_0.1.txt", "sm_0.2.txt"), sep="")

	test_that("The report files correspond to gold standard reports (star)", {
		for(name in all.names){
			cat("Checking", name, "\n")
			gold.standard.difference = system(paste("diff -u gold_standard/", name, " reports/", name, sep=""))
			expect_equal(gold.standard.difference, 0)		
		}
	})


	set.seed(100)
	rob = ncpcResampling(all.data$VM, "jackknife", 310/5, alpha=0.05, test.type="x2", verbose=F)

	test_that("Robustness estimates", {
		expect_equal(rob@raw[[1]]$direct, "bap_6.8")
		expect_equal(rob@raw[[1]]$joint, character(0))
		expect_equal(rob@raw[[1]]$indirect, c("tin_6.8", "bin_6.8", "bin_8.10", "twi_2.4"))
		expect_equal(rob@raw[[1]]$directAndJoint, "bap_6.8")
		expect_equal(rob@raw[[1]]$jointIfNotDirect, "bap_6.8")
	
		expect_equal(rob@raw[[2]]$direct, "bin_8.10")
		expect_equal(rob@raw[[2]]$joint, character(0))
		expect_equal(rob@raw[[2]]$indirect, c("bin_6.8", "bin_10.12"))
		expect_equal(rob@raw[[2]]$directAndJoint, "bin_8.10")
		expect_equal(rob@raw[[2]]$jointIfNotDirect, "bin_8.10")
	
		expect_equal(rob@raw[[3]]$direct, c("bin_6.8", "bin_8.10"))
		expect_equal(rob@raw[[3]]$joint, c("bap_6.8"))
		expect_equal(rob@raw[[3]]$indirect, c("tin_6.8", "bin_10.12", "twi_2.4", "twi_4.6"))
		expect_equal(rob@raw[[3]]$directAndJoint, c("bin_6.8", "bin_8.10", "bap_6.8"))
		expect_equal(rob@raw[[3]]$jointIfNotDirect, c("bin_6.8", "bin_8.10"))
	
		expect_equal(rob@raw[[4]]$direct, character(0))
		expect_equal(rob@raw[[4]]$joint, c("bin_8.10", "bin_10.12"))
		expect_equal(rob@raw[[4]]$indirect, c("bin_6.8"))
		expect_equal(rob@raw[[4]]$directAndJoint, c("bin_8.10", "bin_10.12"))
		expect_equal(rob@raw[[4]]$jointIfNotDirect, c("bin_8.10", "bin_10.12"))
	
		expect_equal(rob@raw[[5]]$direct, c("bin_6.8", "bin_8.10"))
		expect_equal(rob@raw[[5]]$joint, character(0))
		expect_equal(rob@raw[[5]]$indirect, c("tin_6.8", "bin_10.12", "twi_2.4", "twi_4.6", "bap_6.8"))
		expect_equal(rob@raw[[5]]$directAndJoint, c("bin_6.8", "bin_8.10"))
		expect_equal(rob@raw[[5]]$jointIfNotDirect, c("bin_6.8", "bin_8.10"))
	
	})

	# clean up
	all.files = c(all.files, paste("reports/", all.names, sep=""))
	for(file in all.files)
		unlink(file)
}
