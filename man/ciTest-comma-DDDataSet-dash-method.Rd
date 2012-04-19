\name{ciTest,DDDataSet-method}
\alias{ciTest,DDDataSet-method}
\alias{ciTest}
\title{Do conditional independence test on DDDataSet...}
\usage{\S4method{ciTest}{DDDataSet}(obj, var1, var2, cond, test.type="mc-x2-c", B, min.table.size, ...)
}
\description{Do conditional independence test on DDDataSet}
\details{This function does a conditional independence var1 indep var2 | cond. The following
test types are available (implemented by package \code{bnlearn}).

For binary data:
\itemize{
\item "fisher" - Fisher's exact test (only for unconditional independence)
\item "mi" - Mutual Information (discrete)
\item "mi-sh" - Mutual Information (discrete, shrinkage)
\item "mc-mi" - Mutual Information (discrete, Monte Carlo)
\item "aict" - AIC-like Test
\item "x2" - Pearson's X^2
\item "mc-x2" - Pearson's X^2 (Monte Carlo)
\item "mc-x2-c" - Pearson's X^2 (Monte Carlo) the corrected version
\item "g2" - G^2 test (requires pcalg package)
}
For continuous data:
\itemize{
\item "mi-g" - Mutual Information (Gaussian)
\item "mi-g-sh" - Mutual Information (Gaussian, shrinkage)
\item "mc-mi-g" - Mutual Information (Gaussian, Monte Carlo)
\item "cor" - Pearson's Linear Correlation
\item "mc-cor" - Pearson's Linear Correlation (Monte Carlo)
\item "zf" - Fisher's Z Test
\item "mc-zf" - Fisher's Z Test (Monte Carlo)
}}
\value{CITestResult object with the result of the test}
\arguments{\item{obj}{DDDataSet object on which (conditional) independence test needs to be done}
\item{var1}{the name or index of the first variable to be tested}
\item{var2}{the name or index of the second variable}
\item{cond}{the names or indexes of variables to condition on (defaults to NULL)}
\item{test.type}{the type of statistical test (defaults to mc-x2)}
\item{B}{the number of replicates for MC-based tests (default to NULL)}
\item{min.table.size}{the minimal number of samples in a contingency table per conditioning set 
(makes sense only for discrete data)}
\item{...}{unused}
}
\examples{data(mesoBin)
# test if tin_4.6 is independent of class labels
ciTest(mesoBin$Meso, "Tin 4-6h", "class")
# test if tin_4.6 is independent of class conditioned on twi_2.4
ciTest(mesoBin$Meso, "Tin 4-6h", "class", "Twi 2-4h")
# repeat the test using G2 asymptotic distribution
ciTest(mesoBin$Meso, "Tin 4-6h", "class", "Twi 2-4h", test.type="g2")}

