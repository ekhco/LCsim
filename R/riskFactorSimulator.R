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
#' runSHG("~/SHG3.6.4", 10000, "M", 1950, 1) # runs SHG to create output_1950.out in SHG directory
#' smoking_history <- processSHG(file = "~/SHG3.6.4/output_1950.out", birth_cohort = 1950)
#' data(smkhist) # or you can load a sample processed SHG output of 10K men born in 1950
#' ANS_m <- riskFactorSimulator(gender="M", birth_cohort = 1950, SHG_out = smkhist, ages = 45:90, seed = 1618)
#' summary(ANS_m$outputOnly[,1:5]) # look at results of simulation
#' diagnosticPlots(gender = "M", birth_cohort=1950, ANS_m$outputOnly) # diagnostic plots

# -----------------------------------------------------------------------------

riskFactorSimulator <- function(gender, birth_cohort, SHG_out, ages, seed) {  # formatted_SHG=NULL

	Gender = gender
	# TAG: START OF CODE COPIED TO R PACKAGE
	# gender="M"
	# gender="F"

	# work="~/QSU/lung/LCsim"
	work2="risk.lc"
	# setwd(work)

	####### get BMI, FH, and PH regression results first #############
	### run myBMI.cond2 to get the regression results for prediction ##
	# OUT.bmi=OUT.fh=NULL

	# plco_file, , impfile, nafile
	# if(gender=="M") { fit_PLCO <- fit_PLCO_BMI_male(plco_race_edu) }   # source("../LCsim/R/fit_PLCO.BMI.male.only.r")
	# if(gender=="F") { fit_PLCO <- fit_PLCO_BMI_female(plco_race_edu) } # source("../LCsim/R/fit_PLCO.BMI.female.only.r")

	# return(list(OUT.bmi, OUT.fh, OUT.ph, OUT.copd))
	# fit_PLCO[[1]] -> OUT.bmi
	# fit_PLCO[[2]] -> OUT.fh
	# fit_PLCO[[3]] -> OUT.ph
	# fit_PLCO[[4]] -> OUT.copd

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

	# setwd(work)
	# source("source.NLST.functions.R")
	# source("source.gwas.functions.R")

	Head="16_0406_1950.cohort"

	outdir = paste(work2,Head,sep="/")
	dir.create(outdir,showWarnings = F, recursive = T)

			# readThis=F
			# if(readThis==T){
			# 	###################### [1] Read inputfiles once ###################################
			# 	dat.m= read.csv(m_smoke_big_100000) # read.csv("Male_smoke100000_big.csv")
			# 	dat.f= read.csv(f_smoke_big_100000) # read.csv("Female_smoke100000_big.csv")
			#
			# 	CIGDAT.m = read.csv(m_smoke_big_CIGDAT) # read.csv("Male_smoke100000_big_CIGDAT2.csv")[,-1]  # from age 1 to 90
			# 	CIGDAT.f = read.csv(f_smoke_big_CIGDAT) # read.csv("Female_smoke100000_big_CIGDAT2.csv")[,-1]  # from age 1 to 90
			# 	### pack-years was generated in this document: 13_0221_US_extrap.summer2.r
			# } #end of readThis



			# ******************
			# ******* EDIT *****
			# 7/11/2019
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

			# # OLD - Eric 7/22/2019 for 1950 cohort
			# ### reference distribution of each variable obtained from literature (1950 cohort) or smoking history generator (smoking related variables)
			# ref.race=c(0.76, 0.10, 0.08, 0.05); names(ref.race)=c("W","B","H","A") ## info at 05_note_fit.models.PLCO2.ppt
			# ref.edu = c(0.78, 0.15, 0.07); names(ref.edu)=c("high","college","postcol")
			# ref.fh = 10.9  # 10.9% for age = 55 for 1950 birth cohort
			# ref.fh.ci=c(9.74, 11.90)
			# ref.ph = 6.14    ####### at age 55 from McClure ---> age = 55
			# ref.ph.ci = c(5.9, 6.4)
			# ref.copd = 7.4
			# ref.copd.ci = c(7.2, 7.6)


			if(gender=="F"){
				# # OLD FOR 1950 - Eric 7/22/2019
				# #### this is smoking status distributuion in the 1950 population (smoking history generator) ###
				# #### this is for plotting comparison plots #####
				# ref.smk = c(0.52135, 0.30737,0.17128 ); names(ref.smk)=c("Never","Former","Current")
				# ref.pky = c(0.66752 ,0.16092, 0.12864, 0.04292,0);names(ref.pky)=c("lst10", "10_30","30_60","60_100","gt100" )
				# ref.cigstop = c(0.52135, 0.17899, 0.15378, 0.13326, 0.01262);names(ref.cigstop)=c("never","cigstop.0", "cigstop.0_20", "cigstop.20_40",  "cigstop.gt40")
				# ref.bmi = c(28.3, 32.1, 30.4); names(ref.bmi)=c("W","B","H") # their average BMI at age = 50 "--Ogden et al. 2004"
				# ref.bmi.sd = c(0.4, 0.5, 0.5) ; names(ref.bmi.sd)=c("W","B","H") # ERic added

				DAT=dat.f
				CIGDAT = CIGDAT.f

			    ####### making CIGDAT.f: reformatting smoking history data ######

			    # doThis=F
			    # if(doThis==T) {
					#
					# ### reformated cigdat ####
					#
					# ages=45:90
					#
					# #### need all age from 1 to 100 ----> TAKES LONG TIME ###########
					# ages=1:90
					# DAT=dat.f
					# CIGDAT2 = reformatSHG(ages,DAT)
					# #setwd(work)
					# #write.csv(CIGDAT2, "Female_smoke100000_big_CIGDAT2.csv")
					#
					# 		# doCheckReference=F
					# 		# if(doCheckReference==T) {
					# 		# 		##################### (1) extract CPD at a given age 50 #########################
					# 		# 		####  "PY.50" is packyear at age 50
					# 		# 		####  "YSQ.50" is year since quit at age=50
					# 		# 		###   "age_smoke50"
					# 		# 		x=DAT[4,]
					# 		# 		AGE=50
					# 		# 		getCPD.small(x,AGE=50)
					# 		#
					# 		#
					# 		# 		#for(j in 1:nrow(DAT)){
					# 		# 		#
					# 		# 		 # tt=getCPD.small(DAT[j,],AGE=50)
					# 		# 		 # if(j==1) ans= tt
					# 		# 		 # if(j>1) ans=c(ans,tt)
					# 		#
					# 		# 		#}
					# 		#
					# 		# 		CPD50.m = ans=apply(DAT,1, getCPD.small, AGE=50)
					# 		# 		mean(CPD50.m)
					# 		#
					# 		# 		CPD60.m = ans=apply(DAT,1, getCPD.small, AGE=60)
					# 		# 		mean(CPD60.m)
					# 		# 		#[1] 4.031163
					# 		#
					# 		# tt2 = apply(DAT,1, getCPD.small, AGE=55); hist(tt2)
					# 		# ref.cpd.mean.55 = mean(tt2)
					# 		# ref.cpd.mean.55 =3.3804 # famale
					# 		# # ref.cpd.mean.55 = 5.088254	male
					# 		#
					# 		# 		##################### (1) extract cigyears at a given age 50 #########################
					# 		# 		#x=DAT[4,]
					# 		# 		#AGE=50
					# 		# 		getCPD.small(x,AGE=50)
					# 		# 		#   X race gend birth_yr age_init age_cess age_death_OC quintile age_smoke1     CPD1 age_smoke2     CPD2 age_smoke3     CPD3 age_smoke4     CPD4 age_smoke5     CPD5
					# 		# 		# 4 4    0    0     1950       16       82           89        3         16 5.744068         17 9.856438         18 11.81121         19 13.34136         20 14.62065
					# 		# 		#   age_smoke6     CPD6 age_smoke7     CPD7 age_smoke8     CPD8 age_smoke9     CPD9 age_smoke10    CPD10 age_smoke11    CPD11 age_smoke12    CPD12 age_smoke13    CPD13
					# 		# 		# 4         21 15.72104         22 16.68238         23 17.53028         24 18.28275          25 18.95324          26 19.55221          27 20.08808          28 20.56776
					# 		# 		j=1
					# 		# 		DAT[j,]
					# 		# 		getCigYears.sm(DAT[j,], AGE=50)
					# 		#
					# 		#
					# 		# 		tt = apply(DAT,1, getCigYears.sm, AGE=50)
					# 		# 		ref.cigyear.50.m = mean(tt)
					# 		# 		ref.cigyear.50.m = 11.6615
					# 		# 		#cigyears50.m.mean = 15.644 # male
					# 		#
					# 		#     tt = apply(DAT,1, getCigYears.sm, AGE=55)
					# 		# ref.cigyear.55.m = mean(tt)
					# 		# ref.cigyear.55.m = 12.83015 # female
					# 		# #ref.cigyear.55.m = 17.09977	 # male
					# 		# 		#  mean(cigyears50.m)
					# 		# 		# [1] 15.64421
					# 		#
					# 		# ####### pack-years ######
					# 		# 	summary(DAT[,"PY.55"])
					# 		# 	ref.pky.mean.55 = 11.42
					# 		# }#doThis
	        # } #end of dothis


				PROB.edu.race = plco_race_edu
				#PROB.edu.race = read.csv("~/Dropbox/work/research/projects/CISNET/cisnet.lung/risk.lc/14_0117_PLCO.fit.models_female/14_0117_PLCO.fit.models_female.probs.race.edu.csv")

				# row.names(PROB.edu.race)=as.character(PROB.edu.race[,1])
				# PROB.edu.race=PROB.edu.race[,-1]
				PROB.edu.race

					# doCheckReference=F
					# if(doCheckReference==T){
					#
					# 		##################### (1) extract CPD at a given age 50 #########################
					#
					# 		####  "PY.50" is packyear at age 50
					# 		####  "YSQ.50" is year since quit at age=50
					# 		###   "age_smoke50"
					# 		x=DAT[4,]
					# 		#AGE=50
					# 		getCPD.small(x,AGE=50)
					# 		#[1] 15.93271
					#
					#
					# 		CPD50.m = ans=apply(DAT,1, getCPD.small, AGE=50)
					# 		mean(CPD50.m)
					# 		# [1] [1] 4.192261
					# 		# > summary(CPD50.m)
					# 		#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
					# 		#   0.000   0.000   0.000   6.257  12.860  38.370
					#
					# 		CPD.50.mean = 4.192261
					#
					#
					# 		##################### (1) extract cigyears at a given age 50 #########################
					# 		#x=DAT[4,]
					# 		#AGE=50
					#
					# 		j=1
					# 		x=DAT[j,]
					# 		getCigYears.sm(x, AGE=50)
					#
					#
					# 		cigyears50.m = apply(DAT,1, getCigYears.sm, AGE=50)
					# 		mean(cigyears50.m)
					# 		#[1] 11.6615
					#
					# 		cigyears50.mean = 11.6615
					#
					# 	ref.pky.mean.55 = 17.61
					#
					# } #doThis

			}# end of gener==F


			if(gender=="M"){

					# # OLD - 1950 cohort
					# #### this is smoking status distributuion in the 1950 population (smoking history generator) ###
					# ref.smk = c(0.37, 0.40, 0.23); names(ref.smk)=c("Never","Former","Current")
					# #ref.pky = c(0.51, 0.21, 0.21, 0.07, 0); names(ref.pky)=c("lst10","10_30","30_60","60_100","gt100")# ref.pky  # from smoking history generator
					# ref.cigstop = c(0.37, 0.24, 0.19, 0.18, 0.02)  ; names(ref.cigstop)=c("never","cigstop.0", "cigstop.0_20", "cigstop.20_40",  "cigstop.gt40") #### reference (age=60 to be used for COPD calibration)
					# ref.bmi = c(28.7, 27.7, 28.9); names(ref.bmi)=c("W","B","H") # their average BMI at age = 50 "--Ogden et al. 2004"
					# ref.bmi.sd = c(0.3, 0.4, 0.3) ; names(ref.bmi.sd)=c("W","B","H") # Eric added
					# ref.pky.mean.55 = 17.61

					DAT=dat.m
					CIGDAT = CIGDAT.m

					##### [0] Process smk history generator data #######
					# make CPD by age (not by year)

					####  [1] Get the mean pky, year of smoking, and quit years at age = 50 --> don't have reference data, but need to use this for calculating mean BMI at age 50 for each race

					# doCheckReference=F
					# if(doCheckReference==T){
					#
					# 		##################### (1) extract CPD at a given age 50 #########################
					# 		getCPD.small(x,AGE=50)
					# 		CPD50.m = ans=apply(DAT,1, getCPD.small, AGE=50)
					# 		mean(CPD50.m)
					#
					# 		CPD60.m = ans=apply(DAT,1, getCPD.small, AGE=60)
					# 		mean(CPD60.m)
					# 		#[1] 4.031163
					#
					# 		tt2 = apply(DAT,1, getCPD.small, AGE=55); hist(tt2)
					# 		ref.cpd.mean.55 = mean(tt2)
					# 		ref.cpd.mean.55 = 5.088254
					#
					# 		##################### (1) extract cigyears at a given age 50 #########################
					# 		getCPD.small(x,AGE=50)
					#
					# 		j=1
					# 		DAT[j,]
					# 		getCigYears.sm(DAT[j,], AGE=50)
					#
					# 		cigyears50.m = apply(DAT,1, getCigYears.sm, AGE=50)
					#
					# 		tt = apply(DAT,1, getCigYears.sm, AGE=55)
					# 		ref.cigyear.55.m = mean(tt)
					# 		ref.cigyear.55.m = 17.09977
					#
					# 		cigyears50.m.mean = 15.644
					#
					# } # doCheckReference

					PROB.edu.race = plco_race_edu
					# row.names(PROB.edu.race)=as.character(PROB.edu.race[,1])
					# PROB.edu.race=PROB.edu.race[,-1]
					PROB.edu.race

					### packyear average ###
					# doThis=F
					# if(doThis==T){
					#
					# 	# > summary(DAT[,"PY.55"])
					# 	#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
					# 	#   0.000   0.000   9.125  17.610  29.820  84.850
					#
					# 	ref.pky.mean.55 = 17.61
					# }#end

			} # if(gender=="M"){




	Head = "16_0406_1950.cohort"

	#Head="14_0117_1950.cohort"

	outdir = paste(work2, Head,sep="/")

	Head2=paste(Head,gender,sep=".")
	plotfile1=paste("Plot_",Head2,".pdf",sep="")
	outfile=paste("riskfactorsData_",Head2,".csv",sep="")

	setwd(outdir)

	# -----------------------
	# RUN THE SIMULATION HERE
	pdf(plotfile1)

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

			# > dim(ANS)
			# [1] 100000    586
			write.csv(ANS$output_smk, outfile)
			## small sample ##
			#write.csv(ANS[1:200,], "outdata_1950.riskSim.1_200.csv")
	dev.off()

	# return the simulated risk factors at the end
	return(ANS)
}

# ERIC SAYS: END riskFactorSimulator() here!
# TAG: END OF CODE COPIED TO R PACKAGE
