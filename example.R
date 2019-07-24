# -----------------------------------
# Eric Chow, 2017
# this code tests the LCRFsim PACKAGE
# -----------------------------------

rm(list=ls()); gc()
library(devtools) # for loading library


# ------------------------------------------------------
# STEP 0: LOAD LCRFsim Package
# ------------------------------------------------------

    # load LCRFsim
    setwd("~/QSU/LCsim/LCRFsim/R")
    start <- Sys.time(); load_all(); message(round((Sys.time()-start),2), " seconds to load LCRFsim package ...")
    document() # creates the Rd files
    # in terminal: build the pdf:
    # cd ~/QSU/LCsim
    # R CMD Rd2pdf LCRFsim

    # load LCRFsim another way
    if (FALSE) {
      setwd("~/QSU/LCsim")
      start <- Sys.time(); install("LCRFsim"); message(round((Sys.time()-start),2), " seconds to load LCRFsim package ...")
      library(LCRFsim)
    }
# ------------------------------------------------------
# STEP 1: RUN SMOKING HISTORY GENERATOR
# these files are produced by the SHG
# write wrapper code that will produce these!
# ------------------------------------------------------

        setwd("~/QSU/LCsim/test_LCRFsim")

        # call the SHG -------------------------------------------------------------
        # MALE
        runSHG(path="~/QSU/LCsim/SHG/SHG6.3.4", n=10000, gender="M", birth_cohort=1950, seed=1234)
        # process the SHG output (saves to smoking_histories_1950.csv)
        m_smoke_big <- processSHG(file="~/QSU/LCsim/SHG/SHG6.3.4/output_1950.out", birth_cohort=1950)
          # smkhist <- m_smoke_big
          # # for package!
          # save(smkhist, file="~/QSU/LCsim/LCRFsim/data/smkhist.RData")
        # --------------------------------------------------------------------------

        # call the SHG -------------------------------------------------------------
        # FEMALE
        runSHG(path="~/QSU/LCsim/SHG/SHG6.3.4", n=1000, gender="F", birth_cohort=1950, seed=1234)
        # process the SHG output (saves to smoking_histories_1950.csv)
        f_smoke_big <- processSHG(file="~/QSU/LCsim/SHG/SHG6.3.4/output_1950.out", birth_cohort=1950)
        # --------------------------------------------------------------------------




        # --------------------------------------------------------------------------
        # use the 10K M 1950 cohort that is packaged with LCRFsim (already processed)
        if (FALSE) {
          data(smkhist); m_smoke_big <- smkhist
          head(m_smoke_big); dim(m_smoke_big)
        }
        # --------------------------------------------------------------------------





        # --------------------------------------------------------------------------
        # uses Summer's sim results (100,000) from previously
        if (FALSE) {
            m_smoke_big = read.csv("~/QSU/LCsim/test_LCRFsim/Male_smoke100000_big.csv")
            f_smoke_big = read.csv("~/QSU/LCsim/test_LCRFsim/Female_smoke100000_big.csv")
        }
        # --------------------------------------------------------------------------


# ------------------------------------------------------
# STEP 2: RUN RISK FACTOR SIMULATOR
# This uses the reformatted SHG files to run the sim
# ------------------------------------------------------

    setwd("~/QSU/LCsim/test_LCRFsim")

    # set parameters

    #######################################
    # DEBUGGING:
    # riskFactorSimulator calls reformatSHG
    #  - take it out and put it here
    #######################################

    # gender="M"
    # birth_cohort = 1950
    # SHG_out = m_smoke_big
    # ages = 45:90
    # seed = 1618

    # Run LC sim here -------------------------------------------------------------
    ANS_m <- riskFactorSimulator(gender="M", birth_cohort = 1950, SHG_out = m_smoke_big, ages = 45:90, seed = 1618)

    ANS_f <- riskFactorSimulator(gender="F", birth_cohort = 1950, SHG_out = f_smoke_big, ages = 45:90, seed = 1618)
    # -----------------------------------------------------------------------------



    summary(ANS_m$outputOnly[,1:5]) # look at results of simulation
    summary(ANS_f$outputOnly[,1:5]) # look at results of simulation

# ------------------------------------------------------
# STEP 3: RUN DIAGNOSTIC PLOTS
# Use output from ariskFactorSimulator
# ------------------------------------------------------

    # test the diagnostic plot functions
    setwd("~/QSU/LCsim/test_LCRFsim")

    pdf("figs/diagnostic_plots_M.pdf", width=6, height=6)
      diagnosticPlots(gender = "M", birth_cohort=1950, ANS_m$outputOnly)
    dev.off()

    pdf("figs/diagnostic_plots_F.pdf", width=6, height=6)
      diagnosticPlots(gender = "F", birth_cohort=1950, ANS_f$outputOnly)
    dev.off()


#      ~ fin ~
