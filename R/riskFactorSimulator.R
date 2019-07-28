#' Risk Factor Simulator
#'
#' For a given SHG cohort with a specific birth year, for a given gender ("F" or "M"), for a range of ages (ie: 45:90), this function will return a simulated
#' set of risk factors for lung cancer (BMI, family history, personal history, COPD) for each individual in the SHG cohort, conditioned on their smoking
#' histories provided by runSHG().  Use this function after you have created a SHG cohort using runSHG().
#' You can also load a sample processed dataset of the SHG by typing data(smkhist). smkhist is a processed SHG data frame of 10000 men born in 1950
#' simulated using SHG3.6.4. This function requires a processed SHG object (SHG_out) that is created by processsSHG(). If you have your own simulated SHG
#' cohort (as an output_19X0.out file), you can process that using processSHG() to create a data.frame that can be inputted into this function.
#' @param gender         the gender (ie: "M" or "F")
#' @param birth_cohort   the birth year of the cohort (ie: 1950)
#' @param SHG_out        the name of the object that is the output from the SHG
#' @param ages           the range of longitudinal ages to be passed to reformatSHG (ie: 45:90)
#' @param seed           the randomization seed
#' @keywords             gender gender birth cohort seed risk factor simulator
#' @export
#' @examples
#' library(LCRFsim)
#' library(splines)
#' runSHG("~/SHG3.6.4", 10000, "M", 1950, 1) # runs SHG to create output_1950.out in SHG directory
#' smoking_history <- processSHG(file = "~/SHG3.6.4/output_1950.out", birth_cohort = 1950)
#' data(smkhist) # or you can load a sample processed SHG output of 10K men born in 1950
#' ANS_m <- riskFactorSimulator(gender="M", birth_cohort = 1950, SHG_out = smkhist, ages = 45:90, seed = 1618)
#' summary(ANS_m$outputOnly[,1:5]) # look at results of simulation

# -----------------------------------------------------------------------------

riskFactorSimulator <- function(gender, birth_cohort, SHG_out, ages, seed) {  # formatted_SHG=NULL

	Gender = gender

	# ---------------------
	# READ FROM SYSDATA.RDS
	# ---------------------
	if (gender=="M" & birth_cohort == 1950) { plco_race_edu <- PROB.edu.raceM.1950 }
	if (gender=="M" & birth_cohort == 1960) { plco_race_edu <- PROB.edu.raceM.1960 }
	if (gender=="F" & birth_cohort == 1950) { plco_race_edu <- PROB.edu.raceF.1950 }
	if (gender=="F" & birth_cohort == 1960) { plco_race_edu <- PROB.edu.raceF.1960 }

	# ---------------------
	# READ REGRESSION MODELS
	# FROM SYSDATA.RDS
	# ---------------------
	if (gender=="M" & birth_cohort == 1950) {
		OUT.bmi <- OUT.bmi.M.1950
		OUT.fh <- OUT.fh.M.1950
		OUT.ph <- OUT.ph.M.1950
		OUT.copd <- OUT.copd.M.1950
	}
	if (gender=="M" & birth_cohort == 1960) {
		OUT.bmi <- OUT.bmi.M.1960
		OUT.fh <- OUT.fh.M.1960
		OUT.ph <- OUT.ph.M.1960
		OUT.copd <- OUT.copd.M.1960
	}
	if (gender=="F" & birth_cohort == 1950) {
		OUT.bmi <- OUT.bmi.F.1950
		OUT.fh <- OUT.fh.F.1950
		OUT.ph <- OUT.ph.F.1950
		OUT.copd <- OUT.copd.F.1950
	}
	if (gender=="F" & birth_cohort == 1960) {
		OUT.bmi <- OUT.bmi.F.1960
		OUT.fh <- OUT.fh.F.1960
		OUT.ph <- OUT.ph.F.1960
		OUT.copd <- OUT.copd.F.1960
	}

			# createSmokingHIstories.R here
			# createSmokingHIstorries.R takes files in SHG and turns them into SHG_out (ie: m_smoke_big)
			if (gender == "F") {dat.f <- SHG_out }
			if (gender == "M") {dat.m <- SHG_out }


			# RUN REFORMATTING HERE
			# if (is.null(formatted_SHG)) { # they haven't reformatted yet
				# print("reformatting the SHG data ... ")
				if (gender == "F") { CIGDAT.f <- reformatSHG(ages, dat.f)}
				if (gender == "M") { CIGDAT.m <- reformatSHG(ages, dat.m)}
			# } # end if (formatted_SHG)
			# else { # they already reformatted - use these files
			# 	if (gender == "F") { CIGDAT.f <- formatted_SHG}
			# 	if (gender == "M") { CIGDAT.m <- formatted_SHG}
			# }


			if(gender=="F"){

				DAT=dat.f
				CIGDAT = CIGDAT.f

				PROB.edu.race = plco_race_edu
				PROB.edu.race

			}# end of gener==F


			if(gender=="M"){

					DAT=dat.m
					CIGDAT = CIGDAT.m


					PROB.edu.race = plco_race_edu
					PROB.edu.race


			} # if(gender=="M"){

	# -----------------------
	# RUN THE SIMULATION HERE


			ANS = riskFactorGenerator3_EC(DAT, CIGDAT, Gender, birth_cohort, PROB.edu.race, OUT.bmi, OUT.fh, OUT.ph, OUT.copd,     seed) # ref.race, ref.edu?

			# turn BMI vars into numeric from cat ---------------------------------
			# get the BMI variables
			has_bmi <- stringr::str_sub(names(ANS$outputOnly),1,3) == "BMI"
			BMI_vars <- names(ANS$outputOnly[has_bmi])
			# for each BMI var, replace it with a numeric version of itself
			for (var in BMI_vars) {
			  ANS$outputOnly[,var] <- as.numeric(as.character(ANS$outputOnly[,var]))
			}

			# turn BMI vars into numeric from cat ---------------------------------
			# get the BMI variables
			has_bmi <- stringr::str_sub(names(ANS$output_smk),1,3) == "BMI"
			BMI_vars <- names(ANS$output_smk[has_bmi])
			# for each BMI var, replace it with a numeric version of itself
			for (var in BMI_vars) {
				ANS$output_smk[,var] <- as.numeric(as.character(ANS$output_smk[,var]))
			}

	# return the simulated risk factors at the end
	return(ANS)
}

# ERIC SAYS: END riskFactorSimulator() here!
# TAG: END OF CODE COPIED TO R PACKAGE
