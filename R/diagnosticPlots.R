# Creates diagnostic plots
#
# This function creates diagnostic plots after the SHG, and LC simulation have been ran.
# @param riskFactorSimulator_out           The output after simulating the risk factors
# @param gender            either "F" for female or "M" for male
# @param birth_cohort      the year of birth of the cohort
# @keywords                diagnostic plot plots
# @examples
# SHG_out <- runSHG(100, gender="M", 1950, 1)
# ANS_f <- riskFactorSimulator(gender="F", birth_cohort = 1950, SHG_out = SHG_out, ages = 45:90, seed = 1)
# diagnosticPlots(gender="M", 1950, riskFactorSimulator_out = ANS_f$outputOnly)

diagnosticPlots <- function(gender, birth_cohort, riskFactorSimulator_out){


			########## reference ##################

			if(birth_cohort==1950){

				ref.race=c(0.76, 0.10, 0.08, 0.05); names(ref.race)=c("W","B","H","A") ## info at 05_note_fit.models.PLCO2.ppt
				ref.edu = c(0.78, 0.15, 0.07); names(ref.edu)=c("high","college","postcol")
				ref.fh = 10.9  # 10.9% for age = 55 for 1950 birth cohort
				ref.fh.ci=c(9.74, 11.90)
				ref.ph = 6.14    ####### at age 55 from McClure ---> age = 55
				ref.ph.ci = c(5.9, 6.4)
				ref.copd = 7.4
				ref.copd.ci = c(7.2, 7.6)

				if(gender=="M"){

						ref.smk = c(0.37, 0.40, 0.23); names(ref.smk)=c("Never","Former","Current")
						#ref.pky = c(0.51, 0.21, 0.21, 0.07, 0); names(ref.pky)=c("lst10","10_30","30_60","60_100","gt100")# ref.pky  # from smoking history generator
						ref.cigstop = c(0.37, 0.24, 0.19, 0.18, 0.02)  ; names(ref.cigstop)=c("never","cigstop.0", "cigstop.0_20", "cigstop.20_40",  "cigstop.gt40") #### reference (age=60 to be used for COPD calibration)
						ref.bmi = c(28.7, 27.7, 28.9); names(ref.bmi)=c("W","B","H") # their average BMI at age = 50 "--Ogden et al. 2004"
						ref.bmi.sd = c(0.3, 0.4, 0.3);names(ref.bmi.sd)=c("W","B","H") # #--Ogden et al. 2004"
						ref.pky.mean.55 = 17.61

				}# gender M

				if(gender=="F"){

						ref.smk = c(0.52135, 0.30737,0.17128 ); names(ref.smk)=c("Never","Former","Current")
						ref.pky = c(0.66752 ,0.16092, 0.12864, 0.04292,0);names(ref.pky)=c("lst10", "10_30","30_60","60_100","gt100" )
						ref.cigstop = c(0.52135, 0.17899, 0.15378, 0.13326, 0.01262);names(ref.cigstop)=c("never","cigstop.0", "cigstop.0_20", "cigstop.20_40",  "cigstop.gt40")
						ref.bmi = c(28.3, 32.1, 30.4); names(ref.bmi)=c("W","B","H") # their average BMI at age = 50 "--Ogden et al. 2004"
						ref.bmi.sd = c(0.4, 0.5, 0.5);names(ref.bmi.sd)=c("W","B","H") # #--Ogden et al. 2004"
				}# gender F

			}# end of if(birth_cohort==1950){



			#################### (1) Race-Education #########################################

			######### breaking into education and race ####################

			race=as.character(riskFactorSimulator_out[,"race"])
			edu=as.character(riskFactorSimulator_out[,"edu"])

			TB.race = table(race)/sum(table(race))
			TB.race
			# race
			#       A       B       H       W
			# 0.04015 0.05039 0.02170 0.88776
			TB.edu = table(edu)/sum(table(edu))
			TB.edu
			# educ
			#  <=high college postcol
			# 0.60296 0.18542 0.21162
			## remove "<="
			names(TB.edu) = gsub("<=","",names(TB.edu))
			#educ = gsub("<=","", educ)


			####### plot comparing references vs. simulated ###

			TAB=TB.race[names(ref.race)]
			TB.race2=TB.race[names(ref.race)]
			barplot(TAB,col=(heat.colors(4)),ylim=c(0,1),main="[Simulated] Race Distribution",cex.axis=1.5,cex.lab=2,cex.main=1.5)
			text(1:length(TAB),TAB+0.1,round(TAB,2),cex=2,col="blue")
			#axis(1, at=c(0.75,2), labels=c("Unadjusted", " Ref"),cex.axis=1.5)

			TAB=ref.race
			barplot(TAB,col=(heat.colors(4)),ylim=c(0,1),main="[Ref] Race Distribution",cex.axis=1.5,cex.lab=2,cex.main=1.5)
			text(1:length(TAB),TAB+0.1,round(TAB,2),cex=2,col="blue")

			TAB=TB.edu[names(ref.edu)]
			TB.edu2=TB.edu[names(ref.edu)]
			barplot(TAB,col=(heat.colors(3)),ylim=c(0,1),main="[Simulated] Education Distribution",cex.axis=1.5,cex.lab=2,cex.main=1.5)
			text(1:length(TAB),TAB+0.1,round(TAB,2),cex=2,col="blue")
			#axis(1, at=c(0.75,2), labels=c("Unadjusted", " Ref"),cex.axis=1.5)

			TAB=ref.edu
			barplot(TAB,col=(heat.colors(4)),ylim=c(0,1),main="[Ref] Education Distribution",cex.axis=1.5,cex.lab=2,cex.main=1.5)
			text(1:length(TAB),TAB+0.1,round(TAB,2),cex=2,col="blue")


			############# education comparing simulation vs. reference with CI from simulation ##########

		   par(mar=c(5,5,6,3))

		   mat=cbind(ref.edu, TB.edu2)
		   barplot(t(mat),col=c("red","yellow"),beside=T, cex.axis=1.5, cex.lab=1.5,ylab="Relative Frequency",cex.main=1.5, ylim=c(0,1),main="[Education] Observed vs. simulated distribution")
		  legend(x="topright", legend=c("Reference", "Simulated"), fill=c("red","yellow"), ncol=2,cex=1.5)


			############# Race comparing simulation vs. reference with CI from simulation ##########
	   	    mat=cbind(ref.race, TB.race2)
		   barplot(t(mat),col=c("red","yellow"),beside=T, cex.axis=1.5, cex.lab=1.5, main="[Race] Observed vs. simulated distribution",
		   ylab="Relative Frequency",cex.main=1.5, ylim=c(0,1))

 		  legend(x="topright", legend=c("Reference", "Simulated"), fill=c("red","yellow"), ncol=2,cex=1.5)



			######################## (2) BMI ###############################################################



			###### check to reference: mean BMI at age 50  #####
			# > ref.bmi
			#    W    B    H
			# 28.7 27.7 28.9

			BMI.50 =as.numeric( as.character(riskFactorSimulator_out[,"BMI.50"]))
			#tt=data.frame(cbind(BMI.50 = bmi50, race=race))
			bmi.ave = aggregate(BMI.50, by=list(race=race), mean)
			#   race        x
			# 1    A 26.76746
			# 2    B 29.16310
			# 3    H 28.84572
			# 4    W 28.88204

			#   race        x
			# 1    A 26.26849
			# 2    B 29.08601
			# 3    H 28.85768
			# 4    W 28.90377

			row.names(bmi.ave)=as.character(bmi.ave[,1])

				## make a table to plot
			tt1=t(bmi.ave)
			tt2=as.numeric(tt1[-1,])
			# > tt2
			# [1] "26.26849" "29.08601" "28.85768" "28.90377"
			names(tt2)=tt1[1,]
			#        A        B        H        W
			# 26.26849 29.08601 28.85768 28.90377

			#        A        B        H        W   --> female
			# 25.10848 32.41363 29.65636 28.60779

			#bmi.ave = tt2

			########### plot excluding asian (since reference doens't have asian) ####
			ci1=ref.bmi - 1.96*ref.bmi.sd
			ci2=ref.bmi + 1.96*ref.bmi.sd

			cbind(bmi.ave[c("W","B","H"),2],ref.bmi, ci1, ci2)
			#           ref.bmi    ci1    ci2
			# W 27.89465    28.7 28.112 29.288
			# B 27.96616    27.7 26.916 28.484
			# H 27.79131    28.9 28.312 29.488

		   tab = c(bmi.ave[bmi.ave[,1]=="W",2], bmi.ave[bmi.ave[,1]=="B",2], bmi.ave[bmi.ave[,1]=="H",2])
		   names(tab)=c("W","B","H")



		  	par(mar=c(5,5,5,3))

		   mat=cbind(ref.bmi, tab)
		   barplot(t(mat),col=c("red","yellow"),beside=T, cex.axis=1.5, cex.lab=1.5,ylab="BMI",
		   main="Reference vs. Simulated \nAverage BMI (age 50) by race", cex.main=1.5, ylim=c(0,35))

		   ##### confidence interval ###
		   xs = c(1.5,4.5,7.5)

		   for(u in 1:length(xs)){


		   lines(rep(xs[u],2), c(ci1[u],ci2[u]), lty="dotted", lwd=1.5)
		   lines(c(xs[u]-0.25, xs[u]+0.25), c(ci1[u],ci1[u]), lwd=1.5)
		   lines(c(xs[u]-0.25, xs[u]+0.25), c(ci2[u],ci2[u]), lwd=1.5)

		   }#end of

		  legend(x="topleft", legend=c("Reference", "Simulated"), fill=c("red","yellow"), ncol=2,cex=1)




	     ############### (3) FH #####################################################

	     FH.55=as.numeric(as.character(riskFactorSimulator_out[,"FH.55"]))

	     ### NA should in incorporated into the denominator
	     myprev = function(x){   sum(is.na(x)==F & (x==1))/ length(x) }
		 meanfh55=myprev(FH.55)


		############ (4) barplot compare to the reference data (age = 55) #############


		  par(mar=c(5,5,5,3))
			#ref.fh.ci = c(9.74, 11.90)
			ci1=ref.fh.ci[1]
			ci2=ref.fh.ci[2]

		   mat=cbind(ref.fh, meanfh55*100)
		   barplot(t(mat),col=c("red","yellow"),beside=T, cex.axis=1.5, cex.lab=1.5,ylab="FH prevalence (%)",
		   main="Prevalence (%) of FH (age 55)", cex.main=1.5, ylim=c(0,15))

		   ##### confidence interval ###
		   xs = c(1.5,2.5)#,7.5)

		   for(u in 1:length(xs)){


		   lines(rep(xs[u],2), c(ci1[u],ci2[u]), lty="dotted", lwd=1.5)
		   lines(c(xs[u]-0.25, xs[u]+0.25), c(ci1[u],ci1[u]), lwd=1.5)
		   lines(c(xs[u]-0.25, xs[u]+0.25), c(ci2[u],ci2[u]), lwd=1.5)

		   }#end of

		  legend(x="topleft", legend=c("Reference", "Simulated prob-based"), fill=c("red","yellow"), ncol=2,cex=1)


	    ######################  (5) PH ######################

	     PH.55=as.numeric(as.character(riskFactorSimulator_out[,"PH.55"]))
		 meanph55=mean(PH.55,na.rm=T)

		############ (4) barplot compare to the reference data (age = 55) #############

		  par(mar=c(5,5,5,3))
			#ref.fh.ci = c(9.74, 11.90)
			ci1=ref.ph.ci[1]
			ci2=ref.ph.ci[2]

		   mat=cbind(ref.ph, meanph55*100)
		   barplot(t(mat),col=c("red","yellow"),beside=T, cex.axis=1.5, cex.lab=1.5,ylab="PHC prevalence (%)",
		   main="Prevalence (%) of PHC (age 55)", cex.main=1.5, ylim=c(0,15))

		   ##### confidence interval ###
		   xs = c(1.5,2.5)#,7.5)

		   for(u in 1:length(xs)){


		   lines(rep(xs[u],2), c(ci1[u],ci2[u]), lty="dotted", lwd=1.5)
		   lines(c(xs[u]-0.25, xs[u]+0.25), c(ci1[u],ci1[u]), lwd=1.5)
		   lines(c(xs[u]-0.25, xs[u]+0.25), c(ci2[u],ci2[u]), lwd=1.5)

		   }#end of

		  legend(x="topleft", legend=c("Reference", "Simulated prob-based"), fill=c("red","yellow"), ncol=2,cex=1)

		  ################ (6)  COPD #########################

	     COPD.55=as.numeric(as.character(riskFactorSimulator_out[,"COPD.55"]))
		 meancopd55=mean(COPD.55,na.rm=T)



			 #simcopd55 = prev.est["COPD.55"]#-0.002
			#0.10779

			  par(mar=c(5,5,5,3))
				#ref.fh.ci = c(9.74, 11.90)
				ci1=ref.copd.ci[1]
				ci2=ref.copd.ci[2]

			   mat=cbind(ref.copd, meancopd55*100)
			   barplot(t(mat),col=c("red","yellow"),beside=T, cex.axis=1.5, cex.lab=1.5,ylab="COPD prevalence (%)",
			   main="Reference vs. Simulated \n Prevalence of COPD (age 55)", cex.main=1.5, ylim=c(0,15))

			   ##### confidence interval ###
			   xs = c(1.5,2.5)#,7.5)

			   for(u in 1:length(xs)){


			   lines(rep(xs[u],2), c(ci1[u],ci2[u]), lty="dotted", lwd=1.5)
			   lines(c(xs[u]-0.25, xs[u]+0.25), c(ci1[u],ci1[u]), lwd=1.5)
			   lines(c(xs[u]-0.25, xs[u]+0.25), c(ci2[u],ci2[u]), lwd=1.5)

			   }#end of

			  legend(x="topleft", legend=c("Reference", "Simulated-based"), fill=c("red","yellow"), ncol=2,cex=1)



}#end of
