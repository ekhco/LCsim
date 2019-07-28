# Creates Smoking histories
#
# This function post-processes the SHG output. It is called by processSHG.R, after decompose_output_files.py has processed the SHG outputs. It reads the people_ages_mat.txt, people_cpd_mat.txt, and people_mat.txt files created by the decompose[...].py.  This function is not intended for direct use by the user, see processSHG().
# @param year               the year of birth for the cohort
# @keywords                 smoking history generator
# @examples
# setwd("~/SHG/") # navigate to the SHG directory first
# createSmokingHistories(year = 1950) # Do not run this code, see processSHG()

createSmokingHistories <- function(year) {
#for (year in seq(1890, 2100,by=10)){
# updated version from Jacob (Dec 11, 2017)

    cat("creating smoking histories... ", year,"\n", sep="")
    # Eric edit: no _female in txt string
    ages.list <- lapply(strsplit(readLines(paste("people_ages_mat_",toString(year),".txt",sep=""))," "),as.integer)
    cpd.list <- lapply(strsplit(readLines(paste("people_cpd_mat_",toString(year),".txt",sep=""))," "),as.integer)
    people.list <- lapply(strsplit(readLines(paste("people_mat_",toString(year),".txt",sep=""))," "),as.integer)

    num_subjects <- length(ages.list)
    smoking.histories <- matrix(0,num_subjects,6+4*100)

    # "quintile" is not included as a column
    colnames(smoking.histories) <- c("race","gend","birth_yr","age_init","age_cess","age_death_OC","age_smoke1","CPD1","age_smoke2","CPD2","age_smoke3","CPD3","age_smoke4","CPD4","age_smoke5","CPD5","age_smoke6","CPD6","age_smoke7","CPD7","age_smoke8","CPD8","age_smoke9","CPD9","age_smoke10","CPD10","age_smoke11","CPD11","age_smoke12","CPD12","age_smoke13","CPD13","age_smoke14","CPD14","age_smoke15","CPD15","age_smoke16","CPD16","age_smoke17","CPD17","age_smoke18","CPD18","age_smoke19","CPD19","age_smoke20","CPD20","age_smoke21","CPD21","age_smoke22","CPD22","age_smoke23","CPD23","age_smoke24","CPD24","age_smoke25","CPD25","age_smoke26","CPD26","age_smoke27","CPD27","age_smoke28","CPD28","age_smoke29","CPD29","age_smoke30","CPD30","age_smoke31","CPD31","age_smoke32","CPD32","age_smoke33","CPD33","age_smoke34","CPD34","age_smoke35","CPD35","age_smoke36","CPD36","age_smoke37","CPD37","age_smoke38","CPD38","age_smoke39","CPD39","age_smoke40","CPD40","age_smoke41","CPD41","age_smoke42","CPD42","age_smoke43","CPD43","age_smoke44","CPD44","age_smoke45","CPD45","age_smoke46","CPD46","age_smoke47","CPD47","age_smoke48","CPD48","age_smoke49","CPD49","age_smoke50","CPD50","age_smoke51","CPD51","age_smoke52","CPD52","age_smoke53","CPD53","age_smoke54","CPD54","age_smoke55","CPD55","age_smoke56","CPD56","age_smoke57","CPD57","age_smoke58","CPD58","age_smoke59","CPD59","age_smoke60","CPD60","age_smoke61","CPD61","age_smoke62","CPD62","age_smoke63","CPD63","age_smoke64","CPD64","age_smoke65","CPD65","age_smoke66","CPD66","age_smoke67","CPD67","age_smoke68","CPD68","age_smoke69","CPD69","age_smoke70","CPD70","age_smoke71","CPD71","age_smoke72","CPD72","age_smoke73","CPD73","age_smoke74","CPD74","age_smoke75","CPD75","age_smoke76","CPD76","age_smoke77","CPD77","age_smoke78","CPD78","age_smoke79","CPD79","age_smoke80","CPD80","age_smoke81","CPD81","age_smoke82","CPD82","age_smoke83","CPD83","age_smoke84","CPD84","age_smoke85","CPD85","age_smoke86","CPD86","age_smoke87","CPD87","age_smoke88","CPD88","age_smoke89","CPD89","age_smoke90","CPD90","age_smoke91","CPD91","age_smoke92","CPD92","age_smoke93","CPD93","age_smoke94","CPD94","age_smoke95","CPD95","age_smoke96","CPD96","age_smoke97","CPD97","age_smoke98","CPD98","age_smoke99","CPD99","age_smoke100","CPD100","YSQ.1", "YSQ.2","YSQ.3","YSQ.4","YSQ.5","YSQ.6","YSQ.7","YSQ.8","YSQ.9","YSQ.10","YSQ.11","YSQ.12","YSQ.13","YSQ.14","YSQ.15","YSQ.16","YSQ.17","YSQ.18","YSQ.19","YSQ.20","YSQ.21","YSQ.22","YSQ.23","YSQ.24","YSQ.25","YSQ.26","YSQ.27","YSQ.28","YSQ.29","YSQ.30","YSQ.31","YSQ.32","YSQ.33","YSQ.34","YSQ.35","YSQ.36","YSQ.37","YSQ.38","YSQ.39","YSQ.40","YSQ.41","YSQ.42","YSQ.43","YSQ.44", "YSQ.45","YSQ.46","YSQ.47","YSQ.48","YSQ.49","YSQ.50","YSQ.51","YSQ.52","YSQ.53","YSQ.54","YSQ.55","YSQ.56","YSQ.57","YSQ.58","YSQ.59","YSQ.60","YSQ.61","YSQ.62","YSQ.63","YSQ.64","YSQ.65","YSQ.66","YSQ.67","YSQ.68","YSQ.69","YSQ.70","YSQ.71","YSQ.72","YSQ.73","YSQ.74","YSQ.75","YSQ.76","YSQ.77","YSQ.78","YSQ.79","YSQ.80","YSQ.81","YSQ.82","YSQ.83","YSQ.84","YSQ.85","YSQ.86","YSQ.87","YSQ.88","YSQ.89","YSQ.90","YSQ.91","YSQ.92","YSQ.93","YSQ.94","YSQ.95","YSQ.96","YSQ.97","YSQ.98","YSQ.99","YSQ.100","PY.1", "PY.2","PY.3","PY.4","PY.5","PY.6","PY.7","PY.8","PY.9","PY.10","PY.11","PY.12","PY.13","PY.14","PY.15","PY.16","PY.17","PY.18","PY.19","PY.20","PY.21","PY.22","PY.23","PY.24","PY.25","PY.26","PY.27","PY.28","PY.29","PY.30","PY.31","PY.32","PY.33","PY.34","PY.35","PY.36","PY.37","PY.38","PY.39","PY.40","PY.41","PY.42","PY.43","PY.44","PY.45","PY.46","PY.47","PY.48","PY.49","PY.50","PY.51","PY.52","PY.53","PY.54","PY.55","PY.56","PY.57","PY.58","PY.59","PY.60","PY.61","PY.62","PY.63","PY.64","PY.65","PY.66","PY.67","PY.68","PY.69","PY.70","PY.71","PY.72","PY.73","PY.74","PY.75","PY.76","PY.77","PY.78","PY.79","PY.80","PY.81","PY.82","PY.83","PY.84","PY.85","PY.86","PY.87","PY.88","PY.89","PY.90","PY.91","PY.92","PY.93","PY.94","PY.95","PY.96","PY.97","PY.98","PY.99","PY.100")
    CPD_final = matrix(0,num_subjects,100)
    for (i in 1:num_subjects){
        # print(i)

        age_and_cpd <- NA*matrix(1,1,2*100)
        CPD <- NA*matrix(1,1,100)
        for (k in 1:100){
            if (people.list[i][[1]][4] != -999){
              if (k == people.list[i][[1]][4]){
                if (people.list[[i]][5] != -999){
                  age_and_cpd[seq(2*k-1,2*(k+people.list[i][[1]][5]-people.list[i][[1]][4])-1, by = 2)] <- ages.list[i][[1]]
                  age_and_cpd[seq(2*k, 2*(k+people.list[i][[1]][5]-people.list[i][[1]][4]), by = 2)] <- cpd.list[i][[1]]
                  CPD[seq(k,(k+people.list[i][[1]][5]-people.list[i][[1]][4]), by = 1)] <- cpd.list[i][[1]]
                } else {
                  age_and_cpd[seq(2*k-1,2*(k+99-people.list[i][[1]][4])-1, by = 2)] <- ages.list[i][[1]]
                  age_and_cpd[seq(2*k, 2*(k+99-people.list[i][[1]][4]), by = 2)] <- cpd.list[i][[1]]
                  CPD[seq(k,(k+99-people.list[i][[1]][4]), by = 1)] <- cpd.list[i][[1]]
                }
                break
              }
            }
        }
        CPD_final[i,] = CPD


        quit_age <- people.list[i][[1]][5]
        if (quit_age == -999){
            quit_age <- 200 # any number 100 or greater will work here
        }

        years_since_quitting <- seq(1,100,by=1) - quit_age*matrix(1,1,100)
        for (k in 1:length(years_since_quitting)){
            if (years_since_quitting[k] < 1){
                years_since_quitting[k] <- NA
            }
        }

        pack_years_final <- matrix(1:100,1,100)
        # if (!is.na(match(1,ages.list[i][[1]]))){
        if (people.list[i][[1]][4] != -999){
            pack_years <- cumsum(cpd.list[i][[1]])/20
            pack_years_temp1 <- pack_years[match(pack_years_final,ages.list[i][[1]])]
            pack_years_temp1[which(!is.na(pack_years_temp1))[length(which(!is.na(pack_years_temp1)))]:length(pack_years_temp1)] <- pack_years_temp1[which(!is.na(pack_years_temp1))[length(which(!is.na(pack_years_temp1)))]]

            pack_years_final[1:length(pack_years_temp1)] <- pack_years_temp1
        } else {
          pack_years_final <- NA*matrix(1,1,100)
        }

        smoking.histories[i,] <- c(people.list[i][[1]],as.vector(age_and_cpd),as.vector(years_since_quitting),as.vector(pack_years_final))
    }

    # Eric edit here: no female in the .csv string
    write.csv(smoking.histories,file=paste("smoking_histories_",toString(year),".csv",sep=""))


    # Eric added for function
    return(smoking.histories)
}
