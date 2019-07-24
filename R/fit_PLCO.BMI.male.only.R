
fit_PLCO_BMI_male <- function(plco_race_edu) {
## plco_file, impfile, nafile,
	OUT.bmi = OUT.fh = OUT.ph= NULL
	work2="risk.lc"
	Gender="M"

	#if(Gender=="M")  Head="14_0113_PLCO.fit.models_male"
	#if(Gender=="F")  Head="14_0113_PLCO.fit.models_female"

	if(Gender=="M")  Head="16_0303_PLCO.fit.models_male"
	if(Gender=="F")  Head="14_0303_PLCO.fit.models_female"

	#Head="14_0113_PLCO.fit.models_female"

	outdir = paste(work2,Head,sep="/")
	dir.create(outdir,showWarnings = F, recursive = T)

	## data is
	#~/Dropbox/work/research/projects/CISNET/cisnet.lung/Receive/PLCO.data.download_0909_2015/PLCO-146/

	doThis=FALSE
	if(doThis==T){

	###################### [1] Read inputfiles once ###################################

			#setwd(work)
			#mydat00=read.csv("PLCO.dat.mycopy_cisnetpop072312.sas7bdat_reduced_0703_2013.csv")

			mydat00=read.csv(plco_file) # Eric: used to be "plco146_jul15_090115.csv"
			mydat0 = as.matrix(mydat00); colnames(mydat0)= colnames(mydat00)
			dim(mydat0)
			#
			#[1] 148024    192
			# [1] 154900    191
					#mydat00=read.table("PLCO_popfile_072312_cleaned_Stanford2.txt",sep="\t",header=T)

	# age=as.numeric(as.character(mydat00[,"bq_age"]))
	# rr=mydat00[,"rndyear"]
	# birthyear = rr-age
	# hist(birthyear)

	NLST.elig=F

	CombineGender=F
	combineArm=T   ######## combine two arm! I will split into two arms after processing data togeter (data cleaning)

	Trial_arm= "CXR" #"UC" # doesn't matter now

			##################[2]Take the subset of data #########################################

			indic.NLST.elig = rep(F,nrow(mydat0))
			indic.NLST.elig[mydat0[,"nlst_flag"]=="1"]=T
			#indic.NLST.elig[mydat0[,"sct_flag_rnd"]=="1"]=T

			indic.subset = rep(T,nrow(mydat0))

			if(NLST.elig==T) indic.subset = indic.NLST.elig

			######## Selecting the subset of the data for corresponding population ##########
			var.rnd ="arm" #1: "CXR" 2: control(usual care)
			var.sex ="sex" #1: male, 2: female

			mydat=NULL

			if(CombineGender==F){

				if(combineArm==F){

					if(Gender=="M" & Trial_arm=="CXR") mydat = subset(mydat0, mydat0[,var.rnd]=="1" & mydat0[,var.sex]=="1" & indic.subset==T)
					if(Gender=="M" & Trial_arm=="UC") mydat = subset(mydat0, mydat0[,var.rnd]=="2" & mydat0[,var.sex]=="1" & indic.subset==T)
					if(Gender=="F" & Trial_arm=="CXR") mydat = subset(mydat0, mydat0[,var.rnd]=="1" & mydat0[,var.sex]=="2" & indic.subset==T)
					if(Gender=="F" & Trial_arm=="UC") mydat = subset(mydat0, mydat0[,var.rnd]=="2" & mydat0[,var.sex]=="2" & indic.subset==T)

				}#end of combineArm=F

				if(combineArm==T){

					if(Gender=="M") mydat = subset(mydat0, mydat0[,var.sex]=="1" & indic.subset==T)
							#if(Gender=="M" & Trial_arm=="UC") mydat = subset(mydat0, mydat0[,var.rnd]=="C" & mydat0[,var.sex]=="M" & indic.subset==T)
					if(Gender=="F") mydat = subset(mydat0, mydat0[,var.sex]=="2" & indic.subset==T)
							#if(Gender=="F" & Trial_arm=="UC") mydat = subset(mydat0, mydat0[,var.rnd]=="C" & mydat0[,var.sex]=="F" & indic.subset==T)

				}#end of combineArm=F
			}#if(CombineGender==F){

			if(CombineGender==T){

				if(combineArm==F){

					if(Trial_arm=="CXR") mydat = subset(mydat0, mydat0[,var.rnd]=="1" & indic.subset==T)
					if(Trial_arm=="UC") mydat = subset(mydat0, mydat0[,var.rnd]=="2" &  indic.subset==T)
					#if(Gender=="F" & Trial_arm=="CXR") mydat = subset(mydat0, mydat0[,var.rnd]=="I" & mydat0[,var.sex]=="F" & indic.subset==T)
					#if(Gender=="F" & Trial_arm=="UC") mydat = subset(mydat0, mydat0[,var.rnd]=="C" & mydat0[,var.sex]=="F" & indic.subset==T)
				}#end of combineArm=F

				if(combineArm==T){

					mydat = subset(mydat0, indic.subset==T)
					#if(Trial_arm=="UC") mydat = subset(mydat0, mydat0[,var.rnd]=="C" &  indic.subset==T)
					#if(Gender=="F" & Trial_arm=="CXR") mydat = subset(mydat0, mydat0[,var.rnd]=="I" & mydat0[,var.sex]=="F" & indic.subset==T)
					#if(Gender=="F" & Trial_arm=="UC") mydat = subset(mydat0, mydat0[,var.rnd]=="C" & mydat0[,var.sex]=="F" & indic.subset==T)
				}#end of combineArm=F
			}#if(CombineGender==T){




			dim(mydat);#table(mydat[,var.rnd]);table(mydat[,var.sex])
			# [1] 15770   109
			#[1] 76682 # male
			#[1]

			n.dat = nrow(mydat) # number of individuals
			n.dat

				#has_xry0-3   Had X-Ray in T0-3
				# ="Unknown/NA" 0="No"
				# 1="Yes"

	 	  mydata = data.frame(mydat)

			# > table(mydata[,"sex"])
			#
			#     1     2
			# 38340 39104
			# > table(mydata[,"arm"])
			#
			#     1
			# 77444

	}#end of doThis

	################### [2] Fit a series of models ###################################################
	# file is "dictionary_cisnet.mar12.051712.pdf"


	## risk factors described at: Tammemagi_2013_NEJM_NLST_riskModel.pdf
	# (1) Age (1-year), (2) Race (W/B/H/A/Inidian/NativeHawaii),(3)Education (level1-6), (4) BMI (continous)
	# (5) COPD (Y/N), (6) History of cancer (Y/N), (7) family history of lung cancer(Y/N),
	# (8) smoking staus (current/former), (9)smoking intensity (c/d), (10)smoking duration (year), (11)smoking quit time(years)

	#GENDER Gender of the participant. "F"='Female' "M"='Male'
	#age  Age at Randomization (updated) Age derived from rnddate and dob. Numeric
	#educat #Education
			# Question 3
			# .F="No Form" .M="Missing"
			# 1="Less than 8 years" 2="8-11 years"
			# 3="12 years or completed high school" 4="Post high school training other than college"
			# 5="Some college"
			# 6="College graduate" 7="Postgraduate"

			# 1: middle, 2: high dropout 3: high, 4: other than college, 5: some college, 6: collage graduate, 7: post graduate
	# race7 Race
			# Derived from questions 2, 2a (v3)
			# 1="White, Non-Hispanic" 2="Black, Non-Hispanic" 3="Hispanic"
			# 4="Asian"
			# 5="Pacific Islander" 6="American Indian" 7="Missing"
	# bmi_curr Current BMI
			# Derived from questions 22,23. People who are under 48 inches tall are set to .R. Females who are taller than 78 inches and males who are taller than 84 are also set to .R. Participants who report a weight less than 60 pounds are set to .R. BMI values
			# Numeric
			# .F="No Form" .M="Missing" .R="Out of Range"
	#bmi_curc
			#.F="No Form" .M="Missing" .R="Out of Range" 1="0 - <18.5" 2="18.5 - <25" 3="25 - <30" 4=">=30"
	#copd_f  Chronic Obstructive Pulmonary Disease (COPD)
			# Built from bronchit_f and emphys_f. If bronchit_f or emphys_f = 1 then copd_f = 1. If either is missing, then copd_f is missing. Else if both are no, then copd_f is no.
			# (format: ynf)
			# .F="No Form" .M="Missing" 0="No" 1="Yes"
	# lung_fh Family History of Lung Cancer?
			# Derived from question 21
			# .F="No Form"
			# .M="Missing"
			# 0="No"
			# 1="Yes, immediate family member" 9="Possibly - relative or cancer type not clear"

	# cig_stat Cigarette Smoking Status
			# Derived from questions 10 and 12. Set to .A if answers contradict or if one answer is missing but the other isn't
			# .A="Ambiguous"
			# .F="No Form"
			# .M="Missing"
			# 0="Never Smoked Cigarettes" 1="Current Cigarette Smoker" 2="Former Cigarette Smoker"
	#cig_per_day  of Cigarettes Smoked Per Day
			#The maximum number of cigarettes smoked per day, based on the range specified. The category 81+ is set to 100.

	# cig_years_rnd # of Years Smoked Cigarettes
			# Derived from questions 10-13. This differs from cig_years because it uses randomization age instead of BQ age.
			# Numeric .F="No Form" .M="Missing"
	#cig_stop_rnd # of Years Since Stopped Smoking Cigarettes
			# Derived from questions 12 and 13. This differs from cig_stop because it uses randomization age instead of BQ age.
			# Numeric
			# .F="No Form" .M="Missing" .N="Not Applicable"
	# ssmokea_f Age Stopped Smoking
			# Question 13, updated based on responses to other questions
			# Numeric
			# .F="No Form"
			# .M="Missing or Inconsistent Data" .N="Not Applicable"
			# .R="Out of Range"
			# smokea_f  Age Started Smoking
					# Question 11, updated based on responses to other questions
					# Numeric
					# .F="No Form"
					# .M="Missing or Inconsistent Data" .N="Not Applicable"
					# .R="Out of Range"
	# pack_years_rnd Pack Years
			# Maximum number of packs smoked per day times years smoked. This differs from pack_years because it uses randomization age instead of BQ age.
			# Numeric .F="No Form" .M="Missing"



	############## Prediction model building procedure ##############

	#####(1) Full model fitting ################################################
	#####(2) Variable selection (non-linear, interaction, etc) ###########################
	#### (3) Apparent validation--evaluating performance (various metrics, e.g. brier, calibration, c-index) ###
	#### (4) Validation (cross-validation and bootstrap validation  ###################
	#### (5) External validation



	##############[2.0] Reference distribution of each variable obtained from literature (1950 cohort) or smoking history generator (smoking related variables) ##

		ref.race=c(0.76, 0.10, 0.08, 0.05); names(ref.race)=c("W","B","H","A") ## info at 05_note_fit.models.PLCO2.ppt
				## from reference (1950 cohort)
				#  W    B     H    A
				# 0.76 0.10 0.08 0.05
		ref.edu = c(0.78, 0.15, 0.07); names(ref.edu)=c("high","college","postcol")
				### from reference (1950 cohort)
				#      high   college   postcol
				#      0.78     0.15    0.07
		ref.fh = 10.9  # 10.9% for age = 55 for 1950 birth cohort
		ref.fh.ci = c(9.74, 11.90)

		#ref.ph = 8.3
		#ref.ph.ci =  c(8, 8.5) # from McClure ---> age = 55

		### exclude LC and skin cancer ### LC is 13% of all cancer and skin is similar
			# > 8.3*(1-2*0.13)
			# [1] 6.142

	    ref.ph = 6.14    ####### at age 55 from McClure ---> age = 55
	    ref.ph.ci = c(5.9, 6.4)

		 ref.copd = 7.4
		 ref.copd.ci = c(7.2, 7.6)

		## calculation from Hewitt_2003_Prevalence_history.of.cancer.pdf  (Table 1)
		#0.0529=1474/(1473+26354)
		# this exclude skin cancer as well

		doThis=F
		if(doThis==T){
			p=0.053
			n=27827# total number
			se=sqrt(p*(1-p)/n)
			ci1 = p - 1.96*se
			ci2 = p + 1.96*se
			ci1;ci2
			# [1] 0.0503677
			# [1] 0.0556323
	     }#



	if(Gender=="M"){

		####### these are needed to calculate the "average" COPD based on the average of popultation in logisgic regression ##
		#### this is smoking status distributuion in the 1950 population (smoking history generator) ###
		ref.smk = c(0.37, 0.40, 0.23); names(ref.smk)=c("Never","Former","Current")
		#ref.pky = c(0.51, 0.21, 0.21, 0.07, 0); names(ref.pky)=c("lst10","10_30","30_60","60_100","gt100")# ref.pky  # from smoking history generator

		ref.bmi = c(28.7, 27.7, 28.9); names(ref.bmi)=c("W","B","H") # their average BMI at age = 50 "--Ogden et al. 2004"
		ref.bmi.sd = c(0.3, 0.4, 0.3) ; names(ref.bmi.sd)=c("W","B","H")

		CPD.50.mean = 6.257271  # from SHG
		cigyears50.mean=15.644


		ref.pky.mean.55 = 17.61
		ref.cigyear.mean.55 = 17.09
		ref.cpd.mean.55 = 5.088254

		ref.cigstop = c(0.37, 0.24, 0.19, 0.18, 0.02)  ; names(ref.cigstop)=c("never","cigstop.0", "cigstop.0_20", "cigstop.20_40",  "cigstop.gt40") #### reference (age=60 to be used for COPD calibration)


		ref.bmi.mean.55 = 27.85  # this is from 1950 birth cohort data (needed for calculating mean personal history given BMI ####

			# > mean(BMIdat[,"BMI.55"])
			# [1] 27.85902
			# > mean(BMIdat[,"BMI.50"])
			# [1] 28.07846



	}#end of if(Gender=="M"){

	if(Gender=="F"){

		#### this is smoking status distributuion in the 1950 population (smoking history generator) ###
		ref.smk = c(0.52135, 0.30737,0.17128 ); names(ref.smk)=c("Never","Former","Current")
		ref.pky = c(0.66752 ,0.16092, 0.12864, 0.04292,0);names(ref.pky)=c("lst10", "10_30","30_60","60_100","gt100" )
		ref.cigstop = c(0.52135, 0.17899, 0.15378, 0.13326, 0.01262);names(ref.cigstop)=c("never","cigstop.0", "cigstop.0_20", "cigstop.20_40",  "cigstop.gt40")

		ref.bmi = c(28.3, 32.1, 30.4); names(ref.bmi)=c("W","B","H") # their average BMI at age = 50 "--Ogden et al. 2004"
		ref.bmi.sd = c(0.4, 0.5, 0.5) ; names(ref.bmi.sd)=c("W","B","H")

	    ### these are needed for BMI prediction ##
		CPD.50.mean = 4.192261
		cigyears50.mean = 11.6615

		ref.pky.mean.55 = 11.42
		ref.cigyear.mean.55 = 12.83 #male: 17.09
		ref.cpd.mean.55 = 3.38  #male: 5.088254

		ref.bmi.mean.55= 28.92

	}#end of if(Gender=="F"){



	##############[2.1] Pr(Education & Race | smoking status) ###############################



	############[2.2] Pr(BMI | education, race, age and smoking status, stopage, smokingyears..) ###########################################
				#Estimate the distribution from PLCO: time varying BMI by smoking history & prevalence data for BMI for 1950 birth cohort
				#http://www.cdc.gov/nchs/data/ad/ad347.pdf



	##setwd(outdir)
	PROB.edu.race = read.csv(plco_race_edu) # paste(Head,".probs.race.edu.csv",sep=""))
			PROB.edu.race; sum(PROB.edu.race[,-1])
	 #[1] 3

	pdf("male.bmi.out.pdf")

	  age.ref = 50 # calculate average BMI at age=50
	  #OUT.bmi = myBMI.cond2(mydat,  age.ref, ref.bmi, Gender)
	  OUT.bmi = myBMI.cond2( PROB.edu.race, age.ref, ref.bmi, ref.bmi.sd, Gender,CPD.50.mean, cigyears50.mean)

	dev.off()

	############[2.3] Pr(Family history of lung cancer| gender, race, age, education, smoking) ##########################
				#logistic regression using PLCO data and offsets using prevalence data for FH



	OUT.fh = myFHL.cond2( ref.fh, ref.fh.ci, ref.pky.mean.55, ref.bmi.mean.55, ref.edu, ref.race, Gender)


	############ Run personal history ##################


	OUT.ph = myPH.cond( ref.ph, ref.ph.ci, ref.pky.mean.55, ref.bmi.mean.55, ref.edu, ref.race, ref.smk,Gender)


	############# run COPD ##############################

	OUT.copd = myCOPD.cond2( ref.copd, ref.copd.ci, ref.smk, ref.cigstop, ref.fh, ref.pky.mean.55, ref.cpd.mean.55, ref.cigyear.mean.55, ref.bmi.mean.55,  ref.edu, ref.race, Gender)

	print("fit PCLO BMI male successful!")

	# added for R package:
	return(list(OUT.bmi, OUT.fh, OUT.ph, OUT.copd))

}
