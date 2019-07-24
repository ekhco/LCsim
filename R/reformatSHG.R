#' Reformats Smoking History Generator Output
#'
#' This function takes a SHG output from runSHG() and reformats it by longitudinal age. It can take some time to run. It returns these columns: cigday, cigyears, cigstat, for each age in the range of ages. It is called internally by riskFactorSimulator(), and was not intended for direct use by the user. see riskFactorSimulator().
#' @param ages              the range of ages for each person in the SHG cohort (ie: 45:90)
#' @param DAT               The SHG output data, see runSHG()
#' @keywords                reformat smoking history generator cigday cigyears citstat ages DAT
#' @export
#' @examples
#' SHG_output <- runSHG(100, gender="M", 1950, 1)
#' ages <- 45:90
#' reformatted_SHG_output <- reformatSHG(ages, SHG_output)

reformatSHG=function(ages,DAT){
        cat("reformatting SHG data ...\n")
        for(u in 1:length(ages)){

            ag = ages[u]

            #### (1) age ###
            age=rep(ag, nrow(DAT))

            ### (2) Extract cigday at age 45 ###
            CPD.at.this.age = apply(DAT,1, getCPD.small,AGE=ag)

            ###(3) Extract cigyears at this age ##
            cigyears.at.this.age = apply(DAT,1,getCigYears.sm,AGE=ag)

            #### (4) cigstat at this age ############
            cigstat.at.this.age = rep(NA, nrow(DAT))
            cigstat.at.this.age[CPD.at.this.age > 0] = "Current"

            ## never smoker ####################
            cigstat.at.this.age[DAT[,"age_init"]=="-999"]= "Never"
            cigstat.at.this.age[CPD.at.this.age == 0 & DAT[,"age_init"]!="-999"]="Former"

            #### combine #######
            cigdat = data.frame(cigday=CPD.at.this.age, cigyears= cigyears.at.this.age, cigstat = cigstat.at.this.age)

            colnames(cigdat)=paste(c("cigday","cigyears","cigstat"),rep(ag,3),sep=".")

            if(u==1) ans = cigdat
            if(u>1)  ans=data.frame(ans,cigdat)
        } # end of for loop
        CIGDAT = ans
        return(CIGDAT)
}
