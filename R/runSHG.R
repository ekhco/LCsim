#' Custom SHG Call
#'
#' The SHG is written in Python and is maintained seperately from this package by CISNET.
#' This function writes a custom "run_tests.py" into the current director for the SHG so
#' that only selected birth cohort years can be run instead of all of the 1890 - 2110 range
#' that is hardcoded into the current run_tests.py distributed with SHG 6.3.4. This function
#' is called by runSHG() and does not need to be called directly. It creates a python script
#' in the SHG directory, called 'run_custom_tests.py'
#' @param birth_cohort       the year of birth for the cohort (ie: 1950)
#' @keywords                 smoking history generator run_tests
#' @export
#' @examples
#' setwd("~/SHG/") # navigate to the SHG directory first
#' SHG_output <- makeCustomRun(birth_cohort = 1950)
makeCustomRun <- function(birth_cohort) {
  # working directory should have been set as SHG directory before this

  # write the following code into a custom python
  run_custom_tests <- c(
    stringr::str_c("from __future__ import division"),
    stringr::str_c("from ipdb import set_trace"),
    stringr::str_c("import numpy as np"),
    stringr::str_c("import os"),
    stringr::str_c("from shutil import move"),
    stringr::str_c("import socket"),
    stringr::str_c("from sys import platform"),
    stringr::str_c(""),
    stringr::str_c(""),
    stringr::str_c("def recompile_source():"),
    stringr::str_c("    if platform == 'linux' or platform == 'linux2':"),
    stringr::str_c("        os.system('sh install.sh')"),
    stringr::str_c("        try:"),
    stringr::str_c("            move('lbc_smokehist.exe', 'lbc_smokehist_linux.exe')"),
    stringr::str_c("        except:"),
    stringr::str_c("            pass"),
    stringr::str_c("    elif platform == 'darwin':"),
    stringr::str_c("        os.system('sh install.sh')"),
    stringr::str_c("        try:"),
    stringr::str_c("            move('lbc_smokehist.exe', 'lbc_smokehist_osx.exe')"),
    stringr::str_c("        except:"),
    stringr::str_c("            pass"),
    stringr::str_c("    elif platform == 'win32' or platform == 'cygwin':"),
    stringr::str_c("        os.system('install.sh')"),
    stringr::str_c("        try:"),
    stringr::str_c("            move('lbc_smokehist.exe', 'lbc_smokehist_win.exe')"),
    stringr::str_c("        except:"),
    stringr::str_c("            pass"),
    stringr::str_c(""),
    stringr::str_c("def delete_old_io_files():"),
    stringr::str_c("    [os.remove(test_file) for test_file in os.listdir('.') if test_file.startswith('input') or test_file.startswith('output') or test_file.startswith('errors')]"),
    stringr::str_c(""),
    stringr::str_c("def create_input_file(cohort):"),
    stringr::str_c("    seed1 = np.int(np.random.random() * 2147483647); seed2 = np.int(np.random.random() * 2147483647)"),
    stringr::str_c("    seed3 = np.int(np.random.random() * 2147483647); seed4 = np.int(np.random.random() * 2147483647)"),
    stringr::str_c("    open('input_{cohort}.txt'.format(**locals()), 'w').write(open('template_input_05032017.txt').read().format(**locals()))"),
    stringr::str_c(""),
    stringr::str_c("def run_shg(cohort):"),
    stringr::str_c("    if platform == 'linux' or platform == 'linux2':"),
    stringr::str_c("        command = './lbc_smokehist_linux.exe input_{cohort}.txt'.format(**locals())"),
    stringr::str_c("    elif platform == 'darwin':"),
    stringr::str_c("        command = './lbc_smokehist_osx.exe input_{cohort}.txt'.format(**locals())"),
    stringr::str_c("        print command"),
    stringr::str_c("    elif platform == 'win32' or platform == 'cygwin':"),
    stringr::str_c("        command = 'lbc_smokehist_win.exe input_{cohort}.txt'.format(**locals())"),
    stringr::str_c("    os.system(command)"),
    stringr::str_c(""),
    stringr::str_c("def parse_results(year):"),
    stringr::str_c("    global sims"),
    stringr::str_c("    people = []"),
    stringr::str_c("    sims = open('output_' + str(year) + '.out', 'r').read()"),
    stringr::str_c("    if sims.find('<RUN>') != -1:"),
    stringr::str_c("        sims = sims[sims.find('<RUN>') + 5:sims.find('</RUN>')].split('", "\\", "n')"),
    stringr::str_c("    else:"),
    stringr::str_c("        sims = sims.split('","\\","n')"),
    stringr::str_c("    sims = [sim.split(';')[0:-1] for sim in sims if len(sim) != 0]"),
    stringr::str_c("    for sim in sims:"),
    stringr::str_c("        person = {}"),
    stringr::str_c("        person['race'] = sim[0]"),
    stringr::str_c("        person['sex'] = sim[1]"),
    stringr::str_c("        person['yob'] = sim[2]"),
    stringr::str_c("        person['init_age'] = sim[3]"),
    stringr::str_c("        person['cess_age'] = sim[4]"),
    stringr::str_c("        person['ocd_age'] = sim[5]"),
    stringr::str_c("        if person['init_age'] != '-999':"),
    stringr::str_c("            smoking_history = np.array(sim[6:])"),
    stringr::str_c("            years = int(smoking_history.shape[0] / 2)"),
    stringr::str_c("            person['smoking_history'] = smoking_history.reshape(years, 2)"),
    stringr::str_c("        people.append(person)"),
    stringr::str_c("    return people"),
    stringr::str_c(""),
    stringr::str_c("def comment_on_prevalence(people, year):"),
    stringr::str_c("    smokers = 0"),
    stringr::str_c("    for person in people:"),
    stringr::str_c("        if 'smoking_history' in person.keys():"),
    stringr::str_c("            if person['ocd_age'] > person['smoking_history'][0][0] or person['ocd_age'] == '-999':"),
    stringr::str_c("                smokers += 1"),
    stringr::str_c("    print 'year = ' + str(year)"),
    stringr::str_c("    print 'prevalence = ' + str(smokers / len(people))"),
    stringr::str_c(""),
    stringr::str_c("def comment_on_ages(people, key, year):"),
    stringr::str_c("    init_ages = [int(person[key]) for person in people]"),
    stringr::str_c("    init_ages = [d for d in init_ages if d != -999]"),
    stringr::str_c("    local_min = min(init_ages)"),
    stringr::str_c("    local_max = max(init_ages)"),
    stringr::str_c("    local_median = np.median(init_ages)"),
    stringr::str_c(""),
    stringr::str_c("    if key is 'smoking_history':"),
    stringr::str_c("        local_min = min(init_ages)"),
    stringr::str_c("        local_max = max(init_ages)"),
    stringr::str_c("        local_median = np.median(init_ages)"),
    stringr::str_c(""),
    stringr::str_c("    print key + ' = ' + str(local_min) + ', ' + str(local_max) + ', ' + str(local_median) + ' (min, max, median)'"),
    stringr::str_c(""),
    stringr::str_c("    # Addressing minimums"),
    stringr::str_c("    # Minimum init age is normally within some delta of 8 years old, using 6 years here"),
    stringr::str_c("    if key is 'init_age' and abs(local_min - 8) > 7:"),
    stringr::str_c("        print 'WARNING: weird min init age, ' + str(local_min)"),
    stringr::str_c("    # Minimum cess age is normally within some delta of 15 years old, using 6 years here"),
    stringr::str_c("    if key is 'cess_age' and abs(local_min - 15) > 7:"),
    stringr::str_c("        print 'WARNING: weird min cess age, ' + str(local_min)"),
    stringr::str_c("    # Infant mortality in the first year is normally unavoidable"),
    stringr::str_c("    if key is 'ocd_age' and local_min is not 0:"),
    stringr::str_c("        print 'WARNING: weird min ocd age, ' + str(local_min)"),
    stringr::str_c(""),
    stringr::str_c("    # Addressing maximums"),
    stringr::str_c("    # Maximum init age is normally within some delta of 8 years old, using 6 years here"),
    stringr::str_c("    # Maximum cess age is normally within some delta of 15 years old, using 6 years here"),
    stringr::str_c("    if key is 'cess_age' and abs(local_max - 99) > 9 and local_max + year != 2050:"),
    stringr::str_c("        print 'WARNING: weird max cess age, ' + str(local_max)"),
    stringr::str_c("    # Infant mortality in the first year is normally unavoidable"),
    stringr::str_c("    if key is 'ocd_age' and abs(local_max - 99) > 9 and local_max + year != 2050:"),
    stringr::str_c("        print 'WARNING: weird max ocd age, ' + str(local_max)"),
    stringr::str_c(""),
    stringr::str_c("def flatten(list_of_lists):"),
    stringr::str_c("    return [item for sublist in list_of_lists for item in sublist]"),
    stringr::str_c(""),
    stringr::str_c("def comment_on_smoking_histories(people, key, year):"),
    stringr::str_c("    smoking_histories = [person['smoking_history'] for person in people if 'smoking_history' in person.keys()]"),
    stringr::str_c("    smoking_ages = [list(np.array(smoking_history[:, 0], dtype='int')) for smoking_history in smoking_histories]"),
    stringr::str_c("    smoking_amounts = [list(np.array(smoking_history[:, 1], dtype='float')) for smoking_history in smoking_histories]"),
    stringr::str_c("    smoking_averages = [np.mean(np.array(smoking_amount)) for smoking_amount in smoking_amounts]"),
    stringr::str_c("    num_switchers = len([d for d in smoking_averages if float(d) != int(d)])"),
    stringr::str_c("    min_smoking_age = min(flatten(smoking_ages))"),
    stringr::str_c("    max_smoking_age = max(flatten(smoking_ages))"),
    stringr::str_c("    min_smoking_amount = min(flatten(smoking_amounts))"),
    stringr::str_c("    average_smoking_amount = np.mean(np.array(flatten(smoking_amounts)))"),
    stringr::str_c("    max_smoking_amount = max(flatten(smoking_amounts))"),
    stringr::str_c("    fraction_of_switchers = num_switchers / len(smoking_histories)"),
    stringr::str_c("    print 'Smoking ages: min = {min_smoking_age}, max = {max_smoking_age}'.format(**locals())"),
    stringr::str_c("    print 'Smoking amounts: min = {min_smoking_amount}, max = {max_smoking_amount}, mean = {average_smoking_amount}'.format(**locals())"),
    stringr::str_c("    print 'Fraction of switchers = {fraction_of_switchers}'.format(**locals())"),
    stringr::str_c(""),
    stringr::str_c("if __name__ == '__main__':"),
    stringr::str_c("    delete_old_io_files()"),
    stringr::str_c("    if socket.gethostname() == 'BenRacine-PC':          ## Add name of your pc in order to rebuild"),
    stringr::str_c("        recompile_source()"),
    stringr::str_c("    # for year in range(1890, 2110, 10):"),
    stringr::str_c("    year = ",birth_cohort),
    stringr::str_c("    create_input_file(year)"),
    stringr::str_c("    print"),
    stringr::str_c("    print"),
    stringr::str_c("    run_shg(year)"),
    stringr::str_c("    people = parse_results(year)"),
    stringr::str_c("    comment_on_prevalence(people, year)"),
    stringr::str_c("    comment_on_ages(people, 'init_age', year)"),
    stringr::str_c("    comment_on_ages(people, 'cess_age', year)"),
    stringr::str_c("    comment_on_ages(people, 'ocd_age', year)"),
    stringr::str_c("    comment_on_smoking_histories(people, 'ocd_age', year)"),
    stringr::str_c("    print"),
    stringr::str_c("")
            )

  # write to file
  fileConn <- file("run_custom_tests.py")
      writeLines(run_custom_tests, fileConn)
  close(fileConn)
}





#' Runs Smoking History Generator
#'
#' The SHG is written in Python and is maintained seperately from this package by CISNET. This function calls the SHG from within R,
#' creates the SHG output_19X0.out file, and also returns the SHG cohort. The SHG needs to be installed properly prior to use.
#' To use runSHG(), you pass it the directory where the SHG is installed.  Tested using SHG6.3.4 and Python v2.7.10 on MacOSX.
#' SHG6.3.4.zip was obtained from cisnet.flexkb.net, see their readme.md for setup instructions.  In python, you'll likely want these
#' dependencies: six, ipdb, enum, ipython_genutils, ipython.  To install the SHG, run install.sh.  To run the SHG from R, run:
#' (1) python run_tests.py directly, or use this R wrapper function to create a custom 'run_custom_tests.py'. No object is returned
#' by this wrapper function. You will be prompted to look for the appropriate output_19X0.out, input_19X0.txt, and errors_19X0.txt
#' in the SHG directory.
#' If you do not wish to use the SHG, you can load a processed smoking history dataset by typing: data(smkhist) that is included
#' with this package. The smkhist dataset is for a 1950 birth cohort of 10000 men from SHG3.6.4.

#' @param path               a string indicating the directory that the SHG is installed in (ie: "~/SHG6.3.4")
#' @param n                  the size of the cohort (ie: 10000)
#' @param gender             the gender (ie: "M" == 0 or "F" == 1)
#' @param birth_cohort       the year of birth for the cohort (ie: 1950 or 1960)
#' @param seed               the seed for random number generator
#' @keywords                 smoking history generator
#' @export
#' @examples
#' library(LCRFsim)
#' runSHG("~/SHG3.6.4", 10000, "F", 1950, 1) # runs SHG to create output_1950.out in SHG directory
#' smoking_history <- processSHG(file = "~/SHG3.6.4/output_1950.out", birth_cohort = 1950)
#' data(smkhist) # or you can load a sample processed SHG output

runSHG <- function(path, n, gender, birth_cohort, seed){

    oldwd <- getwd()
    # the path of the SHG
    setwd(path)

    if (gender %in% c("F","f","female")) {SEX = 1}
    if (gender %in% c("M","m","male")) {SEX = 0}

    # fix weird 5x10e6 problem
    n <- as.integer(n)

    cat(stringr::str_c("\nrunning SHG for n = ",n,", gender = ",SEX," (",gender,")",", cohort birth year = ",birth_cohort, ", seed = ",seed,"... "))

    # runs twice, once for M, once for F
    # ---------------------------------------------------------
    # step 0: produce the template input file that the SHG uses
    # filename: template_input_05032017.txt
    # ---------------------------------------------------------
    CONFIG <- c(stringr::str_c("SEED_INIT=",seed),
                stringr::str_c("SEED_CESS=",seed+1),
                stringr::str_c("SEED_OCD=",seed+2),
                stringr::str_c("SEED_MISC=",seed+3),
                stringr::str_c("RACE=0"),
                stringr::str_c("SEX=",SEX),
                stringr::str_c("YOB=", birth_cohort),
                stringr::str_c("CESSATION_YEAR=0"),
                stringr::str_c("REPEAT=", n),
                stringr::str_c("INIT_PROB=data/05_03_2017/lbc_shg_initiation05032017.txt"),
                stringr::str_c("CESS_PROB=data/05_03_2017/lbc_shg_cessation05032017.txt"),
                stringr::str_c("OCD_PROB=data/05_03_2017/lbc_smokehist_oc_mortality_02212016.txt"),
                stringr::str_c("CPD_DATA=data/05_03_2017/lbc_shg_cpd05032017.txt"),
                stringr::str_c("OUTPUTFILE=output_{cohort}.out"),
                stringr::str_c("ERRORFILE=errors_{cohort}.txt"))

    # write to file
    fileConn <- file("template_input_05032017.txt")
        writeLines(CONFIG, fileConn)
    close(fileConn)

    # ---------------------------------------------------------
    # step 1 - 3: run_tests.py, decompose.py, create_SH.R
    # ---------------------------------------------------------

    # make a "run_custom_tests.py" in SHG directory for the user selected brith cohort
    makeCustomRun(birth_cohort)

    # call_SHG <- "run_tests.py"
    call_SHG <- "run_custom_tests.py"   # ASSUMED TO BE IN SHG DIR

    # --------------------------------
    # step 1: run_tests.py (ships with SHG), actually run_custom_tests.py for birth_cohort year
    cmd <- "python"
    system2(cmd, call_SHG, stdout=TRUE) # for one birth_cohort year

    cat("DONE!\n")
    cat("check the SHG directory for the appopriate SHG output files: \noutput_", birth_cohort,".out\ninput_", birth_cohort,".txt\nerrors_", birth_cohort,".txt\n", sep="")
    setwd(oldwd)

}










#' Process raw SHG output
#'
#' The SHG creates an out_1950.out file which is the raw output. The following
#' function runs python code that decomposes that output into a few other files
#' that will create an R-readable table. A SHG file is saved to "smoking_histories.csv"
#' in the working directory. The function also returns the smoking history as a data frame object in R.
#' If you have a output_XXXX.out file from your own SHG simulations, you can process it using
#' this function to load an R object that can be used for simulating the risk factors.
#'
#' @param file               the location of the output_1950.out file (ie: "shg/output_1950.out")
#' @param birth_cohort       the year of birth for the cohort (ie: 1950)
#' @keywords                 smoking history generator run_tests process processSHG
#' @export
#' @examples
#' smoking_history <- processSHG(file = "~/SHG3.6.4/output_1950.out", birth_cohort = 1950)

processSHG <- function(file, birth_cohort){
  # decompose output output_19X0.out into a file using Iakovos's python code
  cat("processing SHG output file ... \n")

  # --------------------------------
  # step 2: decompose_output_files in R package /inst folder
  decompose_out <- paste(system.file(package="LCRFsim"), "decompose_output_files.py", sep="/")
  # run the python code that decomposes the SHG output
  cmd <- "python"
  system2(cmd, c(decompose_out, file, birth_cohort), stdout=TRUE)
  # creates 3 files: people_ages_mat_1950.txt, people_cpd_mat_1950.txt, people_mat1950.txt

  # ---------------------------------
  # step 3: run create_smoking_histories.R file (genderless version) for birth_cohort year
  result <- createSmokingHistories(birth_cohort)
  cat("smoking history saved as 'smoking_histories_", birth_cohort,".csv'\n", sep="")

  # if never started smoking, NA, so set to 0
  result[is.na(result)] <- 0

  # turn it into a data.frame
  result <- data.frame(result)

  # return the smoking history
  return(result)
}
