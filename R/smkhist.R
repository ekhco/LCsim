#' Sample SHG data
#'
#' This is a sample processed SHG dataset of 10,000 men born in 1950. The code for generating is shown below.
#' The cohort was generated using SHG 6.3.4
#'
#' @name smkhist
#' @docType data
#' @keywords data smoking history male 1950 SHG sample data set
#' @examples
#' library(LCRFsim)
#' library(splines)
#' runSHG("~/SHG3.6.4", 10000, "M", 1950, 1) # runs SHG to create output_1950.out in SHG directory
#' smkhist <- processSHG(file = "~/SHG3.6.4/output_1950.out", birth_cohort = 1950) 
#' OUT_m <- riskFactorSimulator(gender="M", birth_cohort = 1950, SHG_out = smkhist, ages = 45:90, seed = 1618)
#' summary(OUT_m$outputOnly[,1:5]) # look at results of simulation

NULL
