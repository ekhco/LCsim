if (FALSE) {


# above code is in riskFactorSimulator.R
# run riskFactorSimulator.R before functsion below?
# ERIC SAYS: END riskFactorSimulator() here!
# TAG: END OF CODE COPIED TO R PACKAGE





















#################### [2] LC risk calculator 7/21/2017 using ##################################################################

work="~/Dropbox/work/research/projects/CISNET/cisnet.lung"
work2="~/Dropbox/work/research/projects/CISNET/cisnet.lung/risk.lc"
#work="/home/hans4/cisnet" # biowulf
#datfol1="~/Dropbox/work/research/paperDraft/Stanford/13_0617_PLCO_NLST_eligibles.from.Ayca"
#datfol2="~/Dropbox/work/research/paperDraft/Stanford/13_0617_PLCO_all_data.fromAyca"

setwd(work)
source("source.NLST.functions.R")
source("source.gwas.functions.R")


mydata=ANS$output_smk[1:100,]



		## reformat the education variable so it meets Tammemagi model implemented in PLCO.risk() ####
		# > levels(data$edu)
		# [1] "<=high"  "college" "postcol"

		xx=data$edu
		levels(xx)=c(2,5,6)
		table(xx)
		# xx
		#  2  5  6
		# 76 18  6


		mydata$edu2= as.numeric(xx)



		## relevel the education variable ##

 		ages = 45:90
 		#intensity=3 # average over the lifetime since the birth

		# intensity = 1: age-specific intensity (wrong)
		# intensity = 2: average intensity during which one smoked; also take into account 0 intensity during which he quit
		# intensity = 3: average intensity since age 1 to the given time: 0 intensity until the beginning of smoking is counted



### I was here
  		myrsk1 = LC.risk.by.age(mydata,ages,intensity=1)

		myrsk2 = myRiskPLCO.ages.cohort2(mydata,ages,intensity=2)
 		myrsk3 = myRiskPLCO.ages.cohort2(mydata,ages,intensity=3)





###################### [2] 10/11/2016: eligibility checking using Martin's model #########################################

work="~/Dropbox/work/research/projects/CISNET/cisnet.lung"
work2="~/Dropbox/work/research/projects/CISNET/cisnet.lung/risk.lc"
#work="/home/hans4/cisnet" # biowulf
#datfol1="~/Dropbox/work/research/paperDraft/Stanford/13_0617_PLCO_NLST_eligibles.from.Ayca"
#datfol2="~/Dropbox/work/research/paperDraft/Stanford/13_0617_PLCO_all_data.fromAyca"

setwd(work)
source("source.NLST.functions.R")
source("source.gwas.functions.R")


outdir = "~/Dropbox/work/research/projects/CISNET/cisnet.lung/risk.lc/16_0406_1950.cohort"

Gender="F"

if(Gender=="M") {

	outfile = "riskfactorsData_16_0406_1950.cohort.M.csv"
	setwd(work)
	CIGDAT = CIGDAT.m = read.csv("Male_smoke100000_big_CIGDAT2.csv")[,-1]  # from age 1 to 90

}#

if(Gender=="F") {
	outfile= "riskfactorsData_16_0406_1950.cohort.F.csv"
	setwd(work)
	CIGDAT = CIGDAT.f = read.csv("Female_smoke100000_big_CIGDAT2.csv")[,-1]  # from age 1 to 90

}#

setwd(outdir)
mydata0=read.csv(outfile)#[,-c(1:3)]
# > colnames(mydata)[1:5]
# [1] "X"        "gend"     "birth_yr" "age_init" "age_cess"


		mydata=NULL

		############ (0) making various CPD variables to be used in risk calculation ########################

		doThis1=T
		if(doThis1==T){

		########### (1) First option to use CPD is to use age-specific CPD as provided in SHG data ######

		########### (2) CPD2: Update CPD variable: average during the time smoked #########################


				#### modify intensity (CPD) so it is accmulative average intensity not intensity of the year ###
				### I was here
				## OK can I calculate CPD from pack-years?
				##  for current smokers: PY = (cig/day)/20 x duration years --> (overall cig/day) = OCPD =  (PY/dur.years)*20  for current smokers
				###### for former smokers,  (OCPD*dur.years + 0*ysq)/(dur.years + year.since.quit)

				ages=45:90

				### extract pack-years
				pynames=paste("PY",ages,sep=".")
				pydat=mydata0[,pynames]
				# > colnames(pydat)
				#  [1] "PY.45" "PY.46" "PY.47" "PY.48" "PY.49" "PY.50" "PY.51" "PY.52" "PY.53" "PY.54" "PY.55" "PY.56" "PY.57" "PY.58" "PY.59"
				# [16] "PY.60" "PY.61" "PY.62" "PY.63" "PY.64" "PY.65" "PY.66" "PY.67" "PY.68" "PY.69" "PY.70" "PY.71" "PY.72" "PY.73" "PY.74"
				# [31] "PY.75" "PY.76" "PY.77" "PY.78" "PY.79" "PY.80" "PY.81" "PY.82" "PY.83" "PY.84" "PY.85" "PY.86" "PY.87" "PY.88" "PY.89"
				# [46] "PY.90"

				cpddat = mydata0[,paste("CPD",ages,sep=".")]
				# > pydat[1:10,10:10]
				# > pydat[1:20,1:10]
				#        PY.45     PY.46     PY.47     PY.48     PY.49     PY.50     PY.51     PY.52     PY.53     PY.54
				# 1   0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				# 2   0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				# 3  28.053756 29.020226 29.981163 30.936569 31.886442 32.830783 33.769592 34.702868 35.630612 36.552824
				# 4  27.403327 28.369797 29.330735 30.286140 31.236013 32.180354 33.119163 34.052439 34.980184 35.902395
				# 5  24.382102 25.076264 25.757649 26.426256 27.082086 27.725139 28.355413 28.972911 29.577631 30.169573
				# 6   0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				# 7  43.645535 43.645535 43.645535 43.645535 43.645535 43.645535 43.645535 43.645535 43.645535 43.645535
				# 8   0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				# 9   0.107161  0.107161  0.107161  0.107161  0.107161  0.107161  0.107161  0.107161  0.107161  0.107161
				# 10 33.435355 34.653188 35.858235 37.050496 38.229970 39.396659 40.550562 41.691679 42.820009 43.935554
				# 11 19.501322 19.501322 19.501322 19.501322 19.501322 19.501322 19.501322 19.501322 19.501322 19.501322
				# 12 13.496734 14.190896 14.872281 15.540888 16.196718 16.839770 17.470045 18.087542 18.692262 19.284204
				# 13  3.391735  3.391735  3.391735  3.391735  3.391735  3.391735  3.391735  3.391735  3.391735  3.391735
				# 14 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054
				# 15 27.403327 27.403327 27.403327 27.403327 27.403327 27.403327 27.403327 27.403327 27.403327 27.403327
				# 16  8.991935  9.247441  9.497133  9.741013  9.979080 10.211334 10.437775 10.658403 10.873219 11.082221
				# 17  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				# 18  0.987791  0.987791  0.987791  0.987791  0.987791  0.987791  0.987791  0.987791  0.987791  0.987791
				# 19  7.870551  8.126057  8.375749  8.619629  8.857696  9.089950  9.316391  9.537019  9.751834  9.960837
				# 20  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000

				### extract smoking years ##
				durdat=mydata0[,paste("DUR",ages,sep=".")]
				# > durdat[1:20,1:10]
				#    DUR.45 DUR.46 DUR.47 DUR.48 DUR.49 DUR.50 DUR.51 DUR.52 DUR.53 DUR.54
				# 1       0      0      0      0      0      0      0      0      0      0
				# 2       0      0      0      0      0      0      0      0      0      0
				# 3      30     31     32     33     34     35     36     37     38     39
				# 4      29     30     31     32     33     34     35     36     37     38
				# 5      33     34     35     36     37     38     39     40     41     42
				# 6       0      0      0      0      0      0      0      0      0      0
				# 7      22     22     22     22     22     22     22     22     22     22  ---> wait
				# 8       0      0      0      0      0      0      0      0      0      0
				# 9       0      0      0      0      0      0      0      0      0      0
				# 10     26     27     28     29     30     31     32     33     34     35
				# 11     28     28     28     28     28     28     28     28     28     28
				# 12     16     17     18     19     20     21     22     23     24     25
				# 13      3      3      3      3      3      3      3      3      3      3
				# 14     20     20     20     20     20     20     20     20     20     20
				# 15     29     29     29     29     29     29     29     29     29     29
				# 16     31     32     33     34     35     36     37     38     39     40
				# 17      0      0      0      0      0      0      0      0      0      0
				# 18      1      1      1      1      1      1      1      1      1      1
				# 19     26     27     28     29     30     31     32     33     34     35
				# 20      0      0      0      0      0      0      0      0      0      0

				ix.smk = durdat!=0
				#  ix.smk[1:10,1:10]
				#       DUR.45 DUR.46 DUR.47 DUR.48 DUR.49 DUR.50 DUR.51 DUR.52 DUR.53 DUR.54
				#  [1,]  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE
				#  [2,]  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE
				#  [3,]   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE
				#  [4,]   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE
				#  [5,]   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE
				#  [6,]  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE
				#  [7,]   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE
				#  [8,]  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE
				#  [9,]  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE  FALSE
				# [10,]   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE

				ave.cpd.dat = matrix(0, ncol=ncol(durdat), nrow=nrow(durdat))
				ave.cpd.dat[ix.smk] = ((pydat/durdat)*20)[ix.smk]
				colnames(ave.cpd.dat)=paste("CPD2",ages,sep=".")

				# > ave.cpd.dat[1:20,1:10]
				#         CPD2.45   CPD2.46   CPD2.47   CPD2.48   CPD2.49   CPD2.50   CPD2.51   CPD2.52   CPD2.53   CPD2.54
				#  [1,]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				#  [2,]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				#  [3,] 18.702504 18.722726 18.738227 18.749436 18.756731 18.760447 18.760884 18.758307 18.752954 18.745038
				#  [4,] 18.898846 18.913198 18.923055 18.928838 18.930917 18.929620 18.925236 18.918022 18.908207 18.895998
				#  [5,] 14.777032 14.750744 14.718657 14.681254 14.638966 14.592178 14.541238 14.486455 14.428112 14.366463
				#  [6,]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				#  [7,] 39.677759 39.677759 39.677759 39.677759 39.677759 39.677759 39.677759 39.677759 39.677759 39.677759
				#  [8,]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				#  [9,]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				# [10,] 25.719504 25.669028 25.613025 25.552066 25.486647 25.417199 25.344101 25.267684 25.188241 25.106031
				# [11,] 13.929516 13.929516 13.929516 13.929516 13.929516 13.929516 13.929516 13.929516 13.929516 13.929516
				# [12,] 16.870917 16.695172 16.524756 16.358829 16.196718 16.037876 15.881859 15.728298 15.576885 15.427364
				# [13,] 22.611565 22.611565 22.611565 22.611565 22.611565 22.611565 22.611565 22.611565 22.611565 22.611565
				# [14,] 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054 16.694054
				# [15,] 18.898846 18.898846 18.898846 18.898846 18.898846 18.898846 18.898846 18.898846 18.898846 18.898846
				# [16,]  5.801248  5.779650  5.755838  5.730008  5.702331  5.672963  5.642041  5.609686  5.576009  5.541110
				# [17,]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
				# [18,] 19.755820 19.755820 19.755820 19.755820 19.755820 19.755820 19.755820 19.755820 19.755820 19.755820
				# [19,]  6.054270  6.019301  5.982678  5.944572  5.905131  5.864484  5.822744  5.780012  5.736373  5.691907
				# [20,]  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000

				##### people who quit smoking at given year should be modified--e.g. ID=10 quit smoking at age=59 and overall cpd should decrease overtime
				###### for former smokers,  (OCPD*dur.years + 0*ysq)/(dur.years + year.since.quit)

				ysqdat=mydata0[,paste("YSQ",ages,sep=".")]
				# > ysqdat[1:20,1:10]
				#    YSQ.45 YSQ.46 YSQ.47 YSQ.48 YSQ.49 YSQ.50 YSQ.51 YSQ.52 YSQ.53 YSQ.54
				# 1      NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
				# 2      NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
				# 3      NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
				# 4      NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
				# 5      NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
				# 6      NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
				# 7       6      7      8      9     10     11     12     13     14     15
				# 8      NA     NA     NA     NA     NA     NA     NA     NA     NA     NA
				# 9      28     29     30     31     32     33     34     35     36     37 <---- quit
				# 10     NA     NA     NA     NA     NA     NA     NA     NA     NA     NA

				ix.quit= is.na(ysqdat)==F
				###### for former smokers,  (OCPD*dur.years + 0*ysq)/(dur.years + year.since.quit)

				ave.cpd.dat2= ave.cpd.dat
				ave.cpd.dat2[ix.quit] = (ave.cpd.dat[ix.quit]*durdat[ix.quit] + 0 * ysqdat[ix.quit])/( durdat[ix.quit]+ ysqdat[ix.quit])
				# > ave.cpd.dat2[1:10,1:10]
				#        CPD2.45  CPD2.46  CPD2.47  CPD2.48  CPD2.49  CPD2.50  CPD2.51  CPD2.52  CPD2.53  CPD2.54
				#  [1,]  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
				#  [2,]  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
				#  [3,] 18.70250 18.72273 18.73823 18.74944 18.75673 18.76045 18.76088 18.75831 18.75295 18.74504
				#  [4,] 18.89885 18.91320 18.92305 18.92884 18.93092 18.92962 18.92524 18.91802 18.90821 18.89600
				#  [5,] 14.77703 14.75074 14.71866 14.68125 14.63897 14.59218 14.54124 14.48646 14.42811 14.36646
				#  [6,]  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
				#  [7,] 31.17538 30.10037 29.09702 28.15841 27.27846 26.45184 25.67384 24.94031 24.24752 23.59218 <---- goes down
				#  [8,]  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
				#  [9,]  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
				# [10,] 25.71950 25.66903 25.61302 25.55207 25.48665 25.41720 25.34410 25.26768 25.18824 25.10603

				#ave.cpd.dat2



			####### (3) CPD3: make average cpd over lifetime: method using all CPD from age 1 to 90 ######

			cpddat=NULL

			CIGDAT[1:10,1:10]
			# > CIGDAT[1:10,1:10]
			#    cigday.1 cigyears.1 cigstat.1 cigday.2 cigyears.2 cigstat.2 cigday.3 cigyears.3 cigstat.3 cigday.4
			# 1         0          0     Never        0          0     Never        0          0     Never        0
			# 2         0          0     Never        0          0     Never        0          0     Never        0
			# 3         0          0    Former        0          0    Former        0          0    Former        0
			# 4         0          0    Former        0          0    Former        0          0    Former        0
			# 5         0          0    Former        0          0    Former        0          0    Former        0
			# 6         0          0     Never        0          0     Never        0          0     Never        0
			# 7         0          0    Former        0          0    Former        0          0    Former        0
			# 8         0          0     Never        0          0     Never        0          0     Never        0
			# 9         0          0    Former        0          0    Former        0          0    Former        0
			# 10        0          0    Former        0          0    Former        0          0    Former        0



			## extract cigday
			AGES=1:90
			cpd = CIGDAT[,paste("cigday",AGES,sep=".")]
			# > cpd[1:10,1:15]
			#    cigday.1 cigday.2 cigday.3 cigday.4 cigday.5 cigday.6 cigday.7 cigday.8 cigday.9 cigday.10 cigday.11 cigday.12 cigday.13
			# 1         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 2         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 3         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 4         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 5         0        0        0        0        0        0        0        0        0         0         0  1.589957  5.278929
			# 6         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 7         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 8         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 9         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 10        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#    cigday.14 cigday.15
			# 1   0.000000  0.000000
			# 2   0.000000  0.000000
			# 3   0.000000  4.926697
			# 4   0.000000  0.000000
			# 5   7.205635  8.745075
			# 6   0.000000  0.000000
			# 7   0.000000  0.000000
			# 8   0.000000  0.000000
			# 9   0.000000  0.000000
			# 10  0.000000  0.000000

			# > cpd[1:10,1:15]
			# > cpd[1:10,1:15]
			#    cigday.1 cigday.2 cigday.3 cigday.4 cigday.5 cigday.6 cigday.7 cigday.8 cigday.9 cigday.10 cigday.11 cigday.12 cigday.13
			# 1         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 2         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 3         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 4         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 5         0        0        0        0        0        0        0        0        0         0         0  1.589957  5.278929
			# 6         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 7         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 8         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 9         0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# 10        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#    cigday.14 cigday.15
			# 1   0.000000  0.000000
			# 2   0.000000  0.000000
			# 3   0.000000  4.926697
			# 4   0.000000  0.000000
			# 5   7.205635  8.745075
			# 6   0.000000  0.000000
			# 7   0.000000  0.000000
			# 8   0.000000  0.000000
			# 9   0.000000  0.000000
			# 10  0.000000  0.000000
			cpdcum=t(apply(cpd, 1, cumsum))
			#cpdcum2[1:10,1:15]
			#      cigday.1 cigday.2 cigday.3 cigday.4 cigday.5 cigday.6 cigday.7 cigday.8 cigday.9 cigday.10 cigday.11 cigday.12 cigday.13
			#  [1,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#  [2,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#  [3,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#  [4,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#  [5,]        0        0        0        0        0        0        0        0        0         0         0  1.589957  6.868885
			#  [6,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#  [7,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#  [8,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#  [9,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			# [10,]        0        0        0        0        0        0        0        0        0         0         0  0.000000  0.000000
			#       cigday.14 cigday.15
			#  [1,]   0.00000  0.000000
			#  [2,]   0.00000  0.000000
			#  [3,]   0.00000  4.926697
			#  [4,]   0.00000  0.000000
			#  [5,]  14.07452 22.819596
			#  [6,]   0.00000  0.000000
			#  [7,]   0.00000  0.000000
			#  [8,]   0.00000  0.000000
			#  [9,]   0.00000  0.000000
			# [10,]   0.00000  0.000000

			ag=matrix(rep(AGES, nrow(cpd)),ncol=length(AGES),byrow=T)
			# > ag[1:5,1:5]
			#      [,1] [,2] [,3] [,4] [,5]
			# [1,]    1    2    3    4    5
			# [2,]    1    2    3    4    5
			# [3,]    1    2    3    4    5
			# [4,]    1    2    3    4    5
			# [5,]    1    2    3    4    5


			cpdmean = cpdcum/ag
			cpdmean[1:10,1:20]
			#       cigday.1 cigday.2 cigday.3 cigday.4 cigday.5 cigday.6 cigday.7 cigday.8 cigday.9 cigday.10 cigday.11 cigday.12 cigday.13
			#  [1,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			#  [2,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			#  [3,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			#  [4,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			#  [5,]        0        0        0        0        0        0        0        0        0         0         0 0.1324964 0.5283758
			#  [6,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			#  [7,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			#  [8,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			#  [9,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			# [10,]        0        0        0        0        0        0        0        0        0         0         0 0.0000000 0.0000000
			#       cigday.14 cigday.15 cigday.16 cigday.17 cigday.18 cigday.19 cigday.20
			#  [1,]  0.000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000  0.000000
			#  [2,]  0.000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000  0.000000
			#  [3,]  0.000000 0.3284465 0.8742671 1.4740443 2.0950310 2.7204271  3.340393
			#  [4,]  0.000000 0.0000000 0.3590042 0.9176768 1.5228730 2.1448989  2.768686
			#  [5,]  1.005323 1.5213064 2.0539422 2.5898833 3.1208173 3.6414896  4.148595
			#  [6,]  0.000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000  0.000000
			#  [7,]  0.000000 0.0000000 0.0000000 0.8034612 1.9980628 3.2788193  4.588181
			#  [8,]  0.000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000  0.000000
			#  [9,]  0.000000 0.0000000 0.0000000 0.1260718 0.1190678 0.1128010  0.107161
			# [10,]  0.000000 0.0000000 0.0000000 0.0000000 0.0000000 0.5538254  1.331813

			colnames(cpdmean)=paste("CPD3",AGES,sep=".")
		cpddat3 = cpdmean[,paste("CPD3",45:90,sep=".")]


		###### update mydata ######
		mydata=cbind(mydata0,ave.cpd.dat2, cpddat3)

		### there are three CPD:
		# 1) CPD: intensity for the given age
		# 2) CPD2: average intensity since the beginning of the smoking: if quit smoking 0 is used for during that times
		# 3) CPD3: average intensity whole lifetime

		#mydata=cbind(mydata,cpddat)

		}#end of dothis1



		############ (1) Evaluate individual risk by age using Martin's model ####################################


		ansUSPF = ansNLST = ansRisk = NULL

		# > dim(mydata)
		# [1] 100000    586
		# > dim(mydata)
		# [1] 100000    631
		# > colnames(mydata)
		#   [1]   "gend"         "birth_yr"     "age_init"     "age_cess"     "age_death_OC" "quintile"     "age_smoke1"   "CPD1"         "age_smoke2"
		#  [11] "CPD2"         "age_smoke3"   "CPD3"         "age_smoke4"   "CPD4"         "age_smoke5"   "CPD5"         "age_smoke6"   "CPD6"         "age_smoke7"
		# [241] "YSQ.90"       "YSQ.91"       "YSQ.92"       "YSQ.93"       "YSQ.94"       "YSQ.95"       "YSQ.96"       "YSQ.97"       "YSQ.98"       "YSQ.99"
		# [261] "PY.54"        "PY.55"        "PY.56"        "PY.57"        "PY.58"        "PY.59"        "PY.60"        "PY.61"        "PY.62"        "PY.63"
		# [341] "STAT.78"      "STAT.79"      "STAT.80"      "STAT.81"      "STAT.82"      "STAT.83"      "STAT.84"      "STAT.85"      "STAT.86"      "STAT.87"
		# [361] "CPD.52"       "CPD.53"       "CPD.54"       "CPD.55"       "CPD.56"       "CPD.57"       "CPD.58"       "CPD.59"       "CPD.60"       "CPD.61"
		# [421] "DUR.66"       "DUR.67"       "DUR.68"       "DUR.69"       "DUR.70"       "DUR.71"       "DUR.72"       "DUR.73"       "DUR.74"       "DUR.75"
		# [451] "BMI.48"       "BMI.49"       "BMI.50"       "BMI.51"       "BMI.52"       "BMI.53"       "BMI.54"       "BMI.55"       "BMI.56"       "BMI.57"
		# [571] "PH.76"        "PH.77"        "PH.78"        "PH.79"        "PH.80"        "PH.81"        "PH.82"        "PH.83"        "PH.84"        "PH.85"
		# [591] "COPD.50"      "COPD.51"      "COPD.52"      "COPD.53"      "COPD.54"      "COPD.55"      "COPD.56"      "COPD.57"      "COPD.58"      "COPD.59"

		# > x[-1]
		#         age         edu         bmi        copd hist.cancer          fh      race.w      race.b
		#          65           4          27           0           0           0           1           0
		#      race.h      race.a     race.ai     race.hw   smkstatus   intensity    duration         ysq
		#           0           0           0           0           1          20          40           0
		x=c(1, 65,4,27,0,0,0,
				1,0,0,0,0,0,
				1,20,40,0)
				 #names(x)=cnames

		PLCO.risk(x)
		#[1] 0.02941633

	   #######################(1.1) reformat some variables so it meet PLCO risk model ###########

	   ##### reformat education ###
		  # education, 1:less than high grad
		   #            2:highschool grad
		   #            3:post highschool
		   #            4:some college
		   #            5:college grad
		   #            6:postgrad

	    edu = as.character(mydata[,"edu"])
		# > table(mydata[,"edu"])
		#
		#  <=high college postcol
		#   78156   14916    6928
		edu2 = edu
		edu2[edu=="<=high"]= 2
		edu2[edu=="college"]=5
		edu2[edu=="postcol"]=6
		# > table(edu2)
		# edu2
		#     2     5     6
		# 78156 14916  6928
		edu2=as.numeric(edu2)

		# Martin's threshold >=0.0151


		########## also calculate overall CPD not yearly CPD variables ###################


		th = 0.0151

		################ (1.2) calculate risk using martin's model only for smokers and those didnt die from OCM at given age #######

 		ages = 45:90
 		#intensity=3 # average over the lifetime since the birth

		# intensity = 1: age-specific intensity (wrong)
		# intensity = 2: average intensity during which one smoked; also take into account 0 intensity during which he quit
		# intensity = 3: average intensity since age 1 to the given time: 0 intensity until the beginning of smoking is counted

  		myrsk1 = myRiskPLCO.ages.cohort2(mydata,ages,intensity=1)
		myrsk2 = myRiskPLCO.ages.cohort2(mydata,ages,intensity=2)
 		myrsk3 = myRiskPLCO.ages.cohort2(mydata,ages,intensity=3)


 		## write

 		setwd(outdir)

 		riskfile1=paste("1950.lifetime.LC.risk.CPD1",".",Gender,".csv",sep="")
 		riskfile2=paste("1950.lifetime.LC.risk.CPD2",".",Gender,".csv",sep="")
 		riskfile3=paste("1950.lifetime.LC.risk.CPD3",".",Gender,".csv",sep="")

 		write.csv(myrsk1, riskfile1)
 		write.csv(myrsk2, riskfile2)
 		write.csv(myrsk3, riskfile3)


 		#myrsk = myRiskPLCO.ages.cohort2(mydata,ages,doPlot=T, plotfile="LC.risk_trajectories.500_errorCorrected.pdf")
		myrsk[1:10,]
		#            RSK.45      RSK.46      RSK.47      RSK.48     RSK.49      RSK.50     RSK.51     RSK.52     RSK.53     RSK.54     RSK.55     RSK.56     RSK.57
		#  [1,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA         NA         NA         NA
		#  [2,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA         NA         NA         NA
		#  [3,] 0.005200578 0.008220257 0.009127020 0.010132629 0.01124758 0.012483456 0.01385326 0.01537142 0.01705389 0.03352574 0.03712629 0.04110169 0.04548876
		#  [4,] 0.018891911 0.020950215 0.023227204 0.025744809 0.02852687 0.031599349 0.03499096 0.03873305 0.04285960 0.04740739 0.05241608 0.05792825 0.06398943
		#  [5,] 0.001922351 0.003767177 0.004106542 0.004472341 0.00486594 0.005288673 0.00574189 0.00981835 0.01063186 0.01149824 0.01241869 0.01339401 0.01442447
		#  [6,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA         NA         NA         NA
		#  [7,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA         NA         NA         NA
		#  [8,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA         NA         NA         NA
		#  [9,] 0.000000000 0.000000000 0.000000000 0.000000000 0.00000000 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000         NA         NA
		# [10,] 0.005558774 0.009710826 0.010752403 0.011902977 0.01317341 0.014575621 0.01612285 0.01782967 0.01971205 0.02178748 0.03399669 0.03751601 0.04138361
		# > mean(myrsk,na.rm=T)
		# [1] 0.01298662

 		### female ##
		#  > mean(myrsk1,na.rm=T)
		# [1] 0.006729703
		# > mean(myrsk2,na.rm=T)
		# [1] 0.01606137
		# > mean(myrsk3,na.rm=T)
		# [1] 0.01150968

		mean(myrsk1,na.rm=T)
		#[1] 0.01298662
		mean(myrsk2,na.rm=T)
		#[1] 0.008595686
		mean(myrsk3,na.rm=T)
		#[1] 0.01727035

 		############## (1.2) Choose the right scale ###############

 		 #intensity=1
		 #intensity=2
		 intensity=3

		 if(intensity==1) mat.risk = myrsk1
		 if(intensity==2) mat.risk = myrsk2
		 if(intensity==3) mat.risk = myrsk3



 		############### (1.3) visualize to confirmation #########################################################



		  plotThis=T

		  ## use doPlot=T in the above
		  if(plotThis==T){

					setwd(outdir)
					plotfile=paste("LC.risk_trajectories.500_intensity",intensity,".", Gender,".pdf",sep="")
					pdf(plotfile)

						kk=500
						par(mar=c(5,6,5,5))

						for(u in 1:kk){

							rcc=as.character(mydata[u,"race"])
							ed=as.character(mydata[u,"edu"])
							#ag.quit=as.character(age.quit[u])
							ag.quit=as.character(mydata[u,"age_cess"])


								plot(ages,  mat.risk[u,], pch=16, col="red",cex.lab=1.5,cex.main=1.5, cex.axis=1.5, xlab="Age", ylab="6-year LC Risk", ylim=c(0,0.5),
										main=paste("Subject ID=",u,"\n(gender=",Gender,", race=",rcc,", edu=",ed,  ", & age.quit=", ag.quit,")",  sep=""),axes=F, cex=1.2)
								axis(1,cex.axis=1.5)
								axis(2,cex.axis=1.5)

								axis(4, at=seq(0,1.5, by=0.25), labels=seq(0,1.5, by=0.25)*100, cex.axis=1.5)
								lines(ages,  mat.risk[u,], lty="dotted", col="red")

								mtext("Cig/day", side=4, line=2.2, cex=1.5)

								### estimated probablity ###

								#lines(ages,COPD.prev[u,]*10,pch=16, col="dark grey")#, lty="dotted")

								 abline(v=mydata[u,"age_death_OC"]-1, col="blue", lty="dotted", lwd=1.5)
								 #points(0, ageOCM[u]

								### cig/day  ###

								### intensity choice ###
								if(intensity==2 | intensity==3) cpdnames =  paste(paste("CPD",intensity,sep=""),ages,sep=".")

								if(intensity==1) cpdnames =  paste("CPD",ages,sep=".")

								fhnames=paste("FH",ages,sep=".")
								phnames=paste("PH",ages,sep=".")
								copdnames=paste("COPD",ages,sep=".")


								lines(ages,mydata[u,cpdnames]/100,lty="dotted", col="dark green")#, lty="dotted")
								points(ages,mydata[u,cpdnames]/100,pch=18, col="dark green")#, lty="dotted")

# 								lines(ages,mydata[u,cpdnames]/100,lty="dotted", col="dark green")#, lty="dotted")
# 								points(ages,mydata[u,cpdnames]/100,pch=18, col="dark green")#, lty="dotted")

								#### family history yes ###
								points(ages, mydata[u,fhnames]*0.1, pch=19, col="purple",cex=0.5)
							    lines(ages,mydata[u,fhnames]*0.1,lty="dotted", col="purple")#, lty="dotted")


								#### COPD yes ###
								points(ages, mydata[u,copdnames]*0.09, pch=19, col="orange",cex=0.5)
							    lines(ages,mydata[u,copdnames]*0.09,lty="dotted", col="orange")#, lty="dotted")


								#### PH yes ###
								points(ages, mydata[u,phnames]*0.08, pch=19, col="pink",cex=0.5)
							    lines(ages,mydata[u,phnames]*0.08,lty="dotted", col="pink")#, lty="dotted")





							#mtext("Cig/Day", side=4, line=2.2,cex=1.5)

								legend(x="topleft",legend=c("LC risk","Cig/Day", "OCM age" , "FH", "COPD","PH"),pch=c(16,18, 16,19), col=c("red","dark green","blue", "purple","orange","pink"),cex=1,
									lty=c("solid","dotted", "dotted","dotted","dotted","dotted"))

							}

				  dev.off()

			}#end of if(doPlot==T){

		######### (1.4) apply different thresholds to estimate the percent eligible over 45 years #########

		# > myrsk[1:10,1:10]
		#            RSK.45      RSK.46      RSK.47      RSK.48     RSK.49      RSK.50     RSK.51     RSK.52     RSK.53     RSK.54
		#  [1,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA
		#  [2,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA
		#  [3,] 0.005200578 0.008220257 0.009127020 0.010132629 0.01124758 0.012483456 0.01385326 0.01537142 0.01705389 0.03352574
		#  [4,] 0.018891911 0.020950215 0.023227204 0.025744809 0.02852687 0.031599349 0.03499096 0.03873305 0.04285960 0.04740739
		#  [5,] 0.001922351 0.003767177 0.004106542 0.004472341 0.00486594 0.005288673 0.00574189 0.00981835 0.01063186 0.01149824
		#  [6,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA
		#  [7,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA
		#  [8,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA


		 if(intensity==1) mat.risk = myrsk1
		 if(intensity==2) mat.risk = myrsk2
		 if(intensity==3) mat.risk = myrsk3

		thrs = c(0.005, 0.0151, 0.02, 0.03)

		start.age=55
		stop.age=80

		ansRisk = myElig.thr2(mat.risk, thrs, start.age, stop.age)

		ansRisk
		#   threshold percentEl n.screenPerPerson n.screesPerElig
		# 1    0.0050   0.23544              4.18           17.74
		# 2    0.0151   0.16826              2.57           15.27
		# 3    0.0200   0.14878              2.11           14.21
		# 4    0.0300   0.11990              1.50           12.47



		ansRisk2 = myElig.thr2(mat.risk, thrs, start.age=55, stop.age=75)
		ansRisk2
		# 1    0.0050   0.27948              4.97           17.77
		# 2    0.0151   0.20167              3.26           16.16
		# 3    0.0200   0.18037              2.75           15.23
		# 4    0.0300   0.14747              2.01           13.60

		row.names(ansRisk)=paste(">=",ansRisk[,1],sep="")

		par(mar=c(5,6,6,5))

		barplot(ansRisk[,2]*100, beside=T, col=rev(terrain.colors(4)), ylim=c(0,27), cex.names=1.5,
					cex.lab=1.5, cex.axis=1.5,cex.main=1.5, ylab="Percent of population eligible (%)", xlab="Threshold",
			main="Percent of the U.S. population eligible for screening")
			xx=round(ansRisk[,2]*100,2)
			text(c(0.7,1.9,3.2, 4.2),xx+1, paste(xx,"%",sep=""), col="blue",cex=1.5)


		barplot(ansRisk[,3], beside=T, col=rev(terrain.colors(4)), ylim=c(0,5), cex.names=1.5, xlab="Threshold",
						cex.lab=1.5, cex.axis=1.5,cex.main=1.5, ylab="Number of screens per person",
				main="Number of screens per person in the cohort")
				xx=round(ansRisk[,3],2)
				text(c(0.7,1.9,3.2, 4.2),xx+0.5, paste(xx), col="blue",cex=1.5)


		######## (2) USPF eligibility  55-80-30-15 #############


		an1=myElig(mydata, op="USPSTF")
		# [1] "PercentEl = 21.804%"
		# [1] "Number of screens per person = 3.31"
		# [1] "Number of screens per person ever screened= 15.16"

		# > names(an1)
		# [1] "el"                "percentEl"         "n.screenPerPerson" "n.screesPerElig"

		######## (3) CMS eligibility  55-80-30-15 #############


		an2=myElig(mydata, op="CMS")
		# [1] "PercentEl = 21.798%"
		# [1] "Number of screens per person = 3.16"
		# [1] "Number of screens per person ever screened= 14.51"

		######## (4) NLST eligibility  55-80-30-15 #############

		an3=myElig(mydata, op="NLST")
		# [1] "PercentEl = 21.79%"
		# [1] "Number of screens per person = 3.04"
		# [1] "Number of screens per person ever screened= 13.93"


		####### (5)comparison plots  #############

		percentEl = c(an1$percentEl, an2$percentEl, an3$percentEl, ansRisk[2, "percentEl"])
		nscr = c(an1$n.screenPerPerson, an2$n.screenPerPerson, an3$n.screenPerPerson, ansRisk[2, "n.screenPerPerson"])
		nscr2 = c(an1$n.screesPerElig, an2$n.screesPerElig, an3$n.screesPerElig, ansRisk[2, "n.screesPerElig"])

		mat = cbind(percentEl, nscr, nscr2)
		row.names(mat)=c("USPSTF","CMS","NLST","Risk-based")
		#            percentEl nscr nscr2
		# USPSTF       0.21804 3.31 15.16
		# CMS          0.21798 3.16 14.51
		# NLST         0.21790 3.04 13.93
		# Risk-based   0.16728 2.19 13.10


		barplot(mat[,1]*100, beside=T, col=rev(heat.colors(4)), ylim=c(0,25), cex.names=1.5,
					cex.lab=1.5, cex.axis=1.5,cex.main=1.5, ylab="Percent of population eligible (%)",
			main="Percent of the U.S. population eligible for screening")
			xx=round(mat[,1]*100,2)
			text(c(0.7,1.9,3.2, 4.2),xx+1, paste(xx,"%",sep=""), col="blue",cex=1.5)


		barplot(mat[,2], beside=T, col=rev(heat.colors(4)), ylim=c(0,5), cex.names=1.5,
						cex.lab=1.5, cex.axis=1.5,cex.main=1.5, ylab="Number of screens per person",
				main="Number of screens per person in the cohort")
				xx=round(mat[,2],2)
				text(c(0.7,1.9,3.2, 4.2),xx+1, paste(xx), col="blue",cex=1.5)


		barplot(mat[,3], beside=T, col=rev(heat.colors(4)), ylim=c(0,20), cex.names=1.5,
						cex.lab=1.5, cex.axis=1.5,cex.main=1.5, ylab="Number of screens per person eligible",
				main="Number of screens per person ever screened")
				xx=round(mat[,3],2)
				text(c(0.7,1.9,3.2, 4.2),xx+1, paste(xx), col="blue",cex=1.5)



	 ##### (6) trajectory visualization ###########


    doThis=F
    if(doThis==T){

		el.usps = an1$el
		#       USPSTF.45 USPSTF.46 USPSTF.47 USPSTF.48 USPSTF.49 USPSTF.50 USPSTF.51 USPSTF.52 USPSTF.53 USPSTF.54
		#  [1,]         0         0         0         0         0         0         0         0         0         0
		#  [2,]         0         0         0         0         0         0         0         0         0         0
		#  [3,]         0         0         0         0         0         0         0         0         0         0
		#  [4,]         0         0         0         0         0         0         0         0         0         0
		#  [5,]         0         0         0         0         0         0         0         0         0         0
		#  [6,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
		#  [7,]        NA        NA        NA        NA        NA        NA        NA        NA        NA        NA

		# > myrsk[1:10,1:10]
		#            RSK.45      RSK.46      RSK.47      RSK.48     RSK.49      RSK.50     RSK.51     RSK.52     RSK.53     RSK.54
		#  [1,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA
		#  [2,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA
		#  [3,] 0.005200578 0.008220257 0.009127020 0.010132629 0.01124758 0.012483456 0.01385326 0.01537142 0.01705389 0.03352574
		#  [4,] 0.018891911 0.020950215 0.023227204 0.025744809 0.02852687 0.031599349 0.03499096 0.03873305 0.04285960 0.04740739
		#  [5,] 0.001922351 0.003767177 0.004106542 0.004472341 0.00486594 0.005288673 0.00574189 0.00981835 0.01063186 0.01149824
		#  [6,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA
		#  [7,]          NA          NA          NA          NA         NA          NA         NA         NA         NA         NA

		myth = 0.015

		ix.us = rowSums( el.usps > 0 & is.na(el.usps)==F )#eligible
		ix.rsk = rowSums(myrsk >= 0.015 & is.na(myrsk)==F) # eligible once



		doPlot=F
		if(doPlot==T){

				setwd(outdir)
				pdf("LC.risk_trajectories_risk.some_update_intensity3.pdf")

					ix = rep(T,nrow(mydata))
					sum(ix)

					kk=100
					myPlot.trj2(ix, kk, mydata, myrsk, el.usps, myth,intensity)

			  dev.off()



				setwd(outdir)
				pdf("LC.risk_trajectories_risk.yes.USPF.no_intensity3.pdf")

					ix = ix.rsk==T & ix.us==F
					sum(ix)

					kk=100
					myPlot.trj2(ix, kk, mydata, myrsk, el.usps, myth,intensity)

			  dev.off()


				setwd(outdir)
				pdf("LC.risk_trajectories_risk.no.USPF.yes.intensity3.pdf")

					ix =ix.rsk==F & ix.us==T
					sum(ix)

					kk=100
					myPlot.trj2(ix, kk, mydata, myrsk, el.usps, myth,intensity)

			  dev.off()


			setwd(outdir)
				pdf("LC.risk_trajectories_risk.yes.USPF.yes.intensity3.pdf")

					ix =ix.rsk==T & ix.us==T
					kk=15
					myPlot.trj2(ix, kk, mydata, myrsk, el.usps, myth,intensity)

			  dev.off()





		}#end of if(doPlot==T){

   }#end of




######################## 5/3/2017 comparing female vs. male ################

work="~/Dropbox/work/research/projects/CISNET/cisnet.lung"
work2="~/Dropbox/work/research/projects/CISNET/cisnet.lung/risk.lc"
#work="/home/hans4/cisnet" # biowulf
#datfol1="~/Dropbox/work/research/paperDraft/Stanford/13_0617_PLCO_NLST_eligibles.from.Ayca"
#datfol2="~/Dropbox/work/research/paperDraft/Stanford/13_0617_PLCO_all_data.fromAyca"

setwd(work)
source("source.NLST.functions.R")
source("source.gwas.functions.R")


outdir = "~/Dropbox/work/research/projects/CISNET/cisnet.lung/risk.lc/16_0406_1950.cohort"

#riskfile1=paste("1950.lifetime.LC.risk.CPD1",".",Gender,".csv",sep="")
#riskfile2=paste("1950.lifetime.LC.risk.CPD2",".",Gender,".csv",sep="")
#riskfile3=paste("1950.lifetime.LC.risk.CPD3",".",Gender,".csv",sep="") # lifetime average of CPD

setwd(outdir)

riskf = read.csv("1950.lifetime.LC.risk.CPD3.F.csv") # lifetime average of CPD
riskm = read.csv("1950.lifetime.LC.risk.CPD3.M.csv") # lifetime average of CPD


# > dim(riskf)
# [1] 100000     47
# > riskf[1:5,]
#   X      RSK.45      RSK.46      RSK.47      RSK.48      RSK.49      RSK.50      RSK.51      RSK.52      RSK.53     RSK.54
# 1 1          NA          NA          NA          NA          NA          NA          NA          NA          NA         NA
# 2 2          NA          NA          NA          NA          NA          NA          NA          NA          NA         NA
# 3 3 0.004211647 0.006764605 0.007618956 0.008570128 0.009628503 0.010805488 0.012113599 0.013566547 0.015179475 0.03011859
# 4 4 0.010101519 0.011399842 0.012844863 0.014452050 0.016238316 0.018222127 0.013901537 0.014160084 0.014423503 0.01469257
# 5 5 0.001268804 0.002592268 0.002940951 0.003330195 0.003764390 0.004248377 0.004787483 0.008499081 0.009548486 0.01071466

# > tt=apply(riskf,1,median,na.rm=T)
# > hist(tt)
# > range(tt)
# [1] 6.351113e-219  1.000000e+05

hist(riskf[,"RSK.65"])
hist(riskm[,"RSK.65"])




ag=80
rname=paste("RSK",ag,sep=".")

mydat=data.frame(rbind(cbind(Gender=rep("F",nrow(riskf)), risk65=riskf[,rname]),
				cbind(Gender=rep("M",nrow(riskf)), risk65=riskm[,rname])))
mydat[1:5,]
#  > mydat[1:5,]
#   Gender             risk65
# 1      F               <NA>
# 2      F               <NA>
# 3      F  0.143189125545885
# 4      F 0.0182710492674228
# 5      F 0.0120904973005482


mydat[,2]=as.numeric(as.character(mydat[,2]))
mydat[1:5,1]
mydat[1:5,2]

with(mydat,boxplot(risk65~Gender))

library(ggplot2)
ggplot(mydat, aes(x=risk65, fill=Gender)) +    geom_density(alpha=0.25) +ggtitle("Gender comparison: LC risk at age 65") +
		theme(axis.text=element_text(size=15),axis.title=element_text(size=14,face="bold"),title=element_text(size=14,face="bold"))




########### other comparison #########

### smoking status at age 55

	ref.smk.f = c(0.52135, 0.30737,0.17128 ); names(ref.smk.f)=c("Never","Former","Current")

 	ref.smk.m = c(0.37, 0.40, 0.23); names(ref.smk.m)=c("Never","Former","Current")


 mat=cbind(rev(ref.smk.f), rev(ref.smk.m))
 colnames(mat)=c("Female","Male")

 barplot(mat,beside=T, col=(heat.colors(3)),main="Smoking Status in 1950 birth cohort (age=55)")
 legend(x="topright",legend=rev(c("Never","Former","Current")),fill=heat.colors(3))







} # FALSE wrapper
