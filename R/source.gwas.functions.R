# source("wga.R") ###### data management #######
# source("wga_read.R") 
# source("wga_util.R") 
# source("wga_stat.R")
# source("R_myExtractGenotype.gwas.r") 

# source("R_myReduce.r")
# source("R_myOR.CI4.r")
# source("R_myMAF.r")
# source("R_runAssoc.small.r")
# source("R_runAssoc.r")
# source("R_myStratAnalysis3.r")



# History Feb 14 2008  Initial coding
#         Feb 19 2008  Add function convert.GenABEL
#                      Change recode.vec to recode.geno
#         Feb 22 2008  Add the name "file" to snp.list and remove
#                      the names "prefix" and "extension". 
#                      Change "in.dir" to "dir" in snp.list
#         Mar 04 2008  Make convert.snpMatrix more efficient
#         Mar 07 2008  Allow snps from different chromosomes to be part
#                      of the same snpMatrix data object.
#                      Move chrms, start.vec, stop.vec into snp.list
#         Mar 12 2008  Add seperate functions to check each list
#                      Add snpNames field to snp.list
#         Mar 13 2008  Change the purpose of the chrms vector in snp.list
#                      Allow reading other flat files
#         Mar 14 2008  Change dir and file fields in locusMap.list
#                      Add snp.list$row1.omitN
#                      Add funtion readPhenoData
#                      Change the function getPhenoData
#                      Add readTable function.
#         Mar 17 2008  Change recode.geno
#                      Rename readTable to table2Format
#                      Rename readPhenoData to readTable
#                      Allow phenotype data and locus map data to be type 1
#         Mar 18 2008  Check for numeric vectors when there is no header.
#                      Do not always check the phenotype list for snpMatrix
#                      Add function to transform a sas file (SAS2Format).
#         Mar 19 2008  Add getLocusMap to read the locus map data set
#                      In convert.GenABEL, delete temporary files as they
#                      are no longer needed instead of at the end.
#         Mar 20 2008  Make scanFile more general.
#                      Add file.type = 4 (SAS data sets)
#         Mar 21 2008  Change in convert.GenABEL when adding the snp data.
#                      Read in the temporary files as lines instead of each
#                      element at a time.
#                      Allow readTable not to have an id variable
#         Mar 24 2008  Use loadData function to load the phenotype and
#                      locusMap data.
#                      Allow an option to pass in the snp data file into
#                      convert.snpMatrix
#                      Change in convert.GenABEL to subset by subject ids.
#                      Change in convert.snpMatrix to subset by subject ids.
#         Mar 25 2008  Add createDummy function and createDummy option in the
#                      phenotype list.
#         Mar 26 2008  Add option to save the snp output files
#                      Re-write sas2Format
#         Mar 27 2008  Get the order of the snps in the snp files
#                      Break functions up ino seperate files
#         Mar 31 2008  Change the snpNames option
#         Apr 01 2008  Add error checks.
#         Apr 03 2008  Change in getPhenoData to check if the subject ids
#                      are unique. If not, then create a variable with
#                      unique ids. Return a list of the data and new var.
#                      Change is complete for GenABEL.
#         Apr 04 2008  Make changes in convert.snpMatrix for non-unique
#                      subject ids.
#                      Change is complete for snpMatrix
#         Apr 05 2008  Change in getDataObject for which = 1.
#                      Change the subject ids when the ids are not unique
#         Apr 08 2008  Change in loadData.type1 to return the subject ids
#                      in the order they appear in the data
#         Apr 09 2008  Add getSubjIds function
#                      Set pdata.ids to be the subject ids even when the
#                      ids are unique in getPheno.info
#         Apr 11 2008  Add getDelimiter
#         Apr 15 2008  Change to sas.list
#         Apr 17 2008  Put the names delete, id, and dir in sas.list
#                      Change in getPhenoData
#                      Put function addToPhenoData in wga_unused.R
#         Apr 19 2008  Change update.snpNames function.
#         Apr 30 2008  Add option to remove rows with missing values
#                      in getPhenoData.
#         May 02 2008  Check for same number of genotypes for each SNP.
#         May 08 2008  Add option so that loadData.type1 returns the
#                      phenotype data.
#                      Remove snp.list$snpNames after the snpNames
#                      vector is defined.
#         May 09 2008  Add inheritance option
#         May 12 2008  Remove remove.miss from getPhenoData
#         Jun 10 2008  Change inheritance to genetic.model
#         Jun 27 2008  Have getLocusMap call getColumns
#         Jun 27 2008  Call update.snp.list to get all snp ids to use
#         Jun 28 2008  Extract snpNames after loadData is called
#                      Add loadData.stream
#         Jul 08 2008  When subjects are not found, save them in a 
#                      global object bad.subject.ids.global
#         Jul 30 2008  Add function to get the intersection of the
#                      subject ids on the phenotype and genotype files.
#         Aug 14 2008  Move intersectSubIds call into getData.1
#         Sep 17 2008  Fix bug in intersectSubIds for streamed input
#         Oct 02 2008  Set snpNames and snpNames.list to NULL in 
#                      intersectSubIds function
#         Oct 10 2008  Add code not to recode genotypes
#         Oct 29 2008  Recode before subsetting the subjects (changes in
#                      loadData.type1 and orderSNP
#         Dec 11 2008  Add code for formulas in pheno.list
#         Dec 15 2008  Move snpMatrix, GenABEL code to wga_unused.R
#         Feb 06 2009  Change in getPheno.info to return the subject ids
#                      which are controls
#         Feb 19 2009  Add MAF calculation in loadData.type1
#         Mar 05 2009  Do not produce an error if cc.var is not specified in
#                      loadData.type1
#                      Compute MAF in orderSNP
#         Mar 20 2009  Check cc.var is 0-1 in getPheno.info
#         Apr 07 2009  Leave ids as variable in the phenotype data
#         May 01 2009  Update getLocusMap function to return snps from a 
#                      given chromosome in a selected range.
#         Jul 21 2009  Call checkDelimiter function in getData.1
#         Aug 10 2009  Do not stop if snpNames were not found. Add error option.
#                      Use options out.miss and out.delimiter in snp.list
#         Sep 17 2009  Add option in getDelimiter for the output snp data
#         Oct 16 2009  Set up changes:
#                      getPhenoData
#                      loadData.type1
#                      getData.1
#                      getPheno.info
#                      intersectSubIds
#         Oct 16 2009  Change in recode.geno
#         Dec 30 2009  Add code for getting variables from formulas
#         Feb 25 2010  Add genotype frequencies and number of missing genotypes
#                      to loadData.type1
#         Mar 25 2010  Add gene.var to getLocusMap

# Function to read the locus map data set
# This function return a list with the names "snp", "chrm", and "loc".
# "snp" and "chrm" are character vectors, and "loc" is a numeric vector.
getLocusMap <- function(file, locusMap.list, temp.list=NULL, op=NULL) {

  # file           File to read
  #                No default
  # locusMap.list  List of options (See locusMap.list.wordpad)
  #                No default
  # op             List with names chrm, start, stop

  # Set up the input arguments to getColumns()
  vars <- c(locusMap.list$snp.var, locusMap.list$chrm.var,
            locusMap.list$loc.var, locusMap.list$alleles.var,
            locusMap.list$gene.var)
  file.list      <- locusMap.list
  file.list$file <- file

  # Check for errors
  chrmFlag <- !is.null(locusMap.list[["chrm.var", exact=TRUE]])
  locFlag  <- !is.null(locusMap.list[["loc.var", exact=TRUE]])
  snpFlag  <- !is.null(locusMap.list[["snp.var", exact=TRUE]])
  allFlag  <- !is.null(locusMap.list[["alleles.var", exact=TRUE]])
  geneFlag <- !is.null(locusMap.list[["gene.var", exact=TRUE]])

  opFlag   <- !is.null(op)

  data <- getColumns(file.list, vars, temp.list=temp.list)

  # Return a list
  ret <- list() 
  
  if (snpFlag) ret$snp     <- data[[locusMap.list$snp.var]]
  if (chrmFlag) ret$chrm   <- data[[locusMap.list$chrm.var]]
  if (locFlag) ret$loc     <- as.numeric(data[[locusMap.list$loc.var]])
  if (allFlag) ret$alleles <- data[[locusMap.list$alleles.var]]
  if (geneFlag) ret$gene   <- data[[locusMap.list$gene.var]]


  # Determine if only certain snps in a range is desired
  if (opFlag) {
    temp <- rep(TRUE, times=length(data[[1]]))
    chr <- getListName(op, "chrm")
    if ((!is.null(chr)) & (chrmFlag)) {
      temp <- temp & (ret$chrm %in% chr)
    }
    start <- getListName(op, "start")
    if ((!is.null(start)) & (locFlag)) {
      temp <- temp & (ret$loc >= as.numeric(start))
    }
    stop <- getListName(op, "stop")
    if ((!is.null(stop)) & (locFlag)) {
      temp <- temp & (ret$loc <= as.numeric(stop))
    }
    temp[is.na(temp)] <- FALSE
    if (snpFlag) ret$snp     <- ret$snp[temp]
    if (chrmFlag) ret$chrm   <- ret$chrm[temp]
    if (locFlag) ret$loc     <- ret$loc[temp]
    if (allFlag) ret$alleles <- ret$alleles[temp]
    if (geneFlag) ret$gene   <- ret$gene[temp]
  }

  ret
 
} # END: getLocusMap

# Function to read and modify the phenotype data based on which
#  package is being used. 
getPhenoData <- function(p, temp.list=NULL) {

  # The options sex.var, sex.female, and sex.male are used
  #  with which = 2 and 3 (see below).
  
  # p is a list with the folowing fields:
  #############################################################
  # file           Phenotype data file. This file must have an
  #                id variable.
  #                No default
  # file.type      1, 3, 4  1 is for an R object file created with the
  #                save() function. 3 is for a table that will be read in
  #                with read.table(). 4 is for a SAS data set.
  #                The default is 3
  # header         0 or 1 . Set to 0 if the file does not contain a header.
  #                The default is 1.
  # delimiter      The delimiter in the table.
  #                The default is ""
  # id.var         Name or column number of the id variable. 
  #                No default.
  # keep.vars      Vector of variable names or column numbers to keep. 
  #                The default is NULL, so that all variables will be kept.
  # remove.vars    Vector of variable names or column numbers to remove. 
  #                The default is NULL, so that all variables will be kept.
  #                Both use.vars and remove.vars cannot be specified. 
  # factor.vars    Vector of variable names or column numbers to convert
  #                into factors.
  #                The default is NULL.
  # make.dummy     0 or 1 to make dummy variables for factor.vars
  #                The default is 0.
  # keep.ids       Vector of ids or row numbers to keep.
  #                The default is NULL.
  # remove.ids     Vector of ids or row numbers to remove.
  #                Both keep.ids and remove.ids cannot be specified
  #                The default is NULL.
  # in.miss        Vector of character strings to define the missing values.
  # remove.miss    0 or 1 to remove rows with missing values
  #                The default is 0.
  # sas.list       See below
  ###############################################################

  # Check the list
  p <- check.pheno.list(p)

  id.var <- p$id.var

  if (!is.null(p$keep.vars)) {
    if (!(p$id.var %in% p$keep.vars)) stop("ERROR: keep.vars must contain id.var")
  } 

  # Set options 
  p$transpose    <- 0
  p$start.row    <- 1
  p$stop.row     <- -1
  p$stream       <- 0
  p$include.row1 <- p$header
  p$snpNames     <- NULL
  p$method       <- 1
  if (p$file.type == 4) {
    p$sas.list$delimiter <- "|"
    temp.list <- check.temp.list(temp.list)
    p$sas.list$temp.list <- temp.list
  }

  # Check if missing values are to be removed.
  # NOTE: a formula can create missing values when the variable
  #  does not have any: log(x)
  removeFlag <- p$remove.miss

  # Check if any formulas are to be applied
  formulas <- getFormulas(p)
  formFlag <- length(formulas) && removeFlag

  if (formFlag) p$remove.miss <- 0

  # The data has been loaded, call readTable
  x <- readTable(p$file, p) 
  p$data <- NULL

  if (formFlag) {
    # Update the data for any formulas
    x <- applyFormulas(x, formulas)

    # Remove missing values
    vars <- getAllVars(p)
    x    <- removeMiss.vars(x, vars=vars, miss=p$in.miss)
    
  } else {
    if (removeFlag) x <- removeMiss(x, miss=p$in.miss)
  }

  nr <- nrow(x)

  # Get the id var
  id <- x[, p$id.var]

  # Check if the ids are unique
  orig.id <- NULL
  uniq.id <- unique(id)
  n.uid   <- length(uniq.id)
  if (n.uid != nr) {
    # Print a warning message
    warning("The subject ids are not all unique")
    
    cnames  <- colnames(x)
    orig.id <- "orig_id_57839"
    
    x[, orig.id] <- x[, id.var]
 
    # Get the new ids
    count        <- rep(0, n.uid)
    names(count) <- uniq.id
    temp         <- x[, orig.id]
    if (is.factor(temp)) temp <- as.character(levels(temp))[temp]
    for (i in 2:nr) {
      if (temp[i] %in% temp[-i]) {
        name <- temp[i]
        if (count[name]) {
          temp[i]     <- paste(name, "_", count[name], sep="")
          count[name] <- count[name] + 1
        } else {
          count[name] <- 2
        }
      }
    }

    x[, id.var] <- temp
  } # END: if (n.uid != nr)

  # Set the row names
  rownames(x) <- as.character(x[, id.var])
  #x[, id.var] <- NULL
 
  # Return a list
  list(data=x, orig.id=orig.id)

} # END: getPhenoData

# Function to return the snp data for which = 1
loadData.type1 <- function(snp.list, pheno.list, temp.list, op=NULL) {

  # snp.list
  # pheno.list
  #####################################################
  # op               A list with the following names
  #  include.row1    0 or 1 to include the header
  #                  The default is 1
  #  include.snps    0 or 1 to include the snp names as the
  #                  first element of field.
  #                  The default is 0.
  #  return.type     1 or 2   1 is a vector of characters
  #                  2 is a matrix. The returned matrix will have
  #                  the column names as subject ids and row names
  #                  as snpnames.
  #                  The default is 1
  #  missing         0 or 1 to return a logical vector to determine
  #                  which snps had missing values
  #                  The default is 1.
  #  snpNames        0 or 1 to return a vector of snp names
  #                  The default is 1.
  #  subjIds         0 or 1 to return a vector of subject ids corresponding
  #                  to the order they appear. (Subjects are columns)
  #                  The default is 0
  #  orderByPheno    0 or 1 to order the columns of the snps according
  #                  to the order of the subject ids in the phenotype
  #                  data. 
  #                  The default is 1.
  #  return.pheno    0 or 1 to return the phenotype data
  #                  The default is 1.
  #  useControls     0 or 1 to use only the subset of controls to determine
  #                  the minor allele. If set to 1, then there must be a name
  #                  "cc.var" in pheno.list which gives the name of the 0-1
  #                  case-control variable.
  #                  The default is 1.
  #  MAF             0 or 1 to compute the MAF for each SNP. If useControls = 1,
  #                  then only the controls will be used.
  #                  The default is 0.
  #  alleles         0 or 1 to return the major/minor alleles
  #                  The default is 1
  #  genoFreqs       0 or 1 to return genotype frequencies as a string 0-1-2
  #                  The default is 0
  #  n.miss          0 or 1 to return number of missing genotypes
  #                  The default is 0.

  # Check the lists. pheno.list is checked in getPheno.info
  snp.list  <- check.snp.list(snp.list)
  temp.list <- check.temp.list(temp.list)

  # Check the options list
  if (is.null(op)) op <- list()
  op <- default.list(op, 
         c("include.row1", "include.snps", "return.type", "missing",
            "snpNames", "subjIds", "orderByPheno", "return.pheno", 
            "useControls", "MAF", "stopOnError", "alleles",
            "genoFreqs", "n.miss"),
         list(1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0))
  if (op$return.type == 2) op$include.row1 <- 1

  # Check for error condition
  if (op$useControls) {
    ccvar <- getListName(pheno.list, "cc.var")
    if (is.null(ccvar)) {
      #print("##########################################################")
      #print("NOTE: pheno.list$cc.var not specified")
      #print("##########################################################")
      op$useControls <- 0
    }
  }

  # Character vector containing the snp data. First element contains the
  # subject ids. Remaining elements contain the snp name and genotypes.
  # This character vector must be delimited.
  i              <- 1
  data.obj       <- NULL
  delimiter      <- getDelimiter(snp.list)
  in.miss        <- get.in.miss(snp.list)
  heter.codes    <- snp.list$heter.codes 
  mode           <- snp.list$genetic.model
  codes          <- getInheritanceVec(mode, recode=snp.list$recode)
  missFlag       <- op$missing
  phenoData.list <- NULL

  # Get phenotype information
  temp       <- getPheno.info(pheno.list, snp.list, temp.list=temp.list)
  pheno.id   <- temp$pheno.id
  pheno.uid  <- temp$pheno.uid
  pdata.flag <- temp$pdata.flag
  pdata.ids  <- temp$pdata.ids
  nsubj      <- temp$nsubj
  nsubj.u    <- temp$nsubj.u
  controls   <- getListName(temp, "controls")

  if (op$return.pheno) phenoData.list <- temp$phenoData.list
  rm(pheno.list, temp)
  temp <- gc(verbose=FALSE)

  # Update snp.list
  snp.list <- update.snp.list(snp.list)

  # Define a list to call loadData
  tList <- getSnpLoadList(snp.list, temp.list)

  snpNames   <- getListName(snp.list, "snpNames")
  total.nsnp <- length(snpNames)
  snpsFlag   <- 1 - is.null(snpNames)
  snp.list$snpNames <- NULL
  count      <- 0
  index      <- 0
  rows       <- NULL
  firstFlag  <- 0
  stopFlag   <- 0
  missing    <- NULL
  snps       <- NULL
  sFlag      <- op$snpNames
  inc.snps   <- op$include.snps
  ret.type   <- op$return.type
  temp.matrix <- NULL
  save.subs  <- NULL
  ret.order  <- NULL
  temp.snp   <- NULL
  MAF.flag   <- op$MAF
  MAF        <- NULL
  allFlag    <- op$alleles
  alleles    <- NULL
  if (!is.null(snp.list$heter.codes)) allFlag <- 0
  gfreqFlag  <- op$genoFreqs
  nmissFlag  <- op$n.miss
  genoFreqs  <- NULL
  n.miss     <- NULL

  # If the snp names were specified, then set sFlag to 1 so that
  #  the snp names will be kept
  if (snpsFlag) sFlag <- 1
  if (ret.type == 2) sFlag <- 1

  # Function to convert character to numeric
  if (mode != 3) {
    convert.fun <- as.integer
  } else {
    convert.fun <- as.numeric
  }

  out.miss <- snp.list$out.miss
  if (snp.list$file.type == 4) {
    out.sep <- "|"
  } else {
    out.sep <- snp.list$out.delimiter
  }

  # Loop over each file and combine the objects
  for (file in snp.list$file) {
    index <- index + 1

    # Get the input file name
    temp <- paste(snp.list$dir, file, sep="")

    # Update the list
    tList$start.row <- snp.list$start.vec[index]
    tList$stop.row  <- snp.list$stop.vec[index]
    tList$snpNames  <- snpNames

    # Get the snp data
    snpData <- loadData(temp, tList)
    nsnps <- length(snpData)
    if (!firstFlag) {
      if (nsnps <= 1) next
    } else {
      if (!nsnps) next
    }

    # Get the subjects and check for errors
    temp     <- getSubjIds(snpData[1], snpData[2], delimiter, pheno.uid, 
                           controls=controls)
    subjects <- temp$subjects
    subj.ids <- temp$subj.ids
    subjFlag <- temp$subjFlag
    total.nsubjects <- temp$total.nsubjects
    control.ids     <- getListName(temp, "control.ids")

    # Check for controls
    if (!any(control.ids)) {
      #print("##########################################################")
      #print("NOTE: No controls found to determine the minor allele")
      #print("##########################################################")
      control.ids <- NULL
    }

    # Define a matrix for return type = 2
    if (ret.type == 2) temp.matrix <- matrix(data=NA, nrow=nsnps-1, ncol=nsubj)

    # Save the subjects to preserve the order
    if (!firstFlag) {
    
      if (pdata.flag) {
        # The order will be as in the phenotype data
        subj.order  <- pheno.id
        subj.order2 <- getOrder(pheno.id, subjects)

        # The new subject ids will be written out
        subjects    <- pdata.ids
      } else {
        if (op$orderByPheno) {
          subj.order <- pheno.id
        } else {
          subj.order <- subjects
        }
        subj.order2  <- getOrder(subj.order, subjects) 
        subjects     <- subjects[subj.order2]
      }
 
      if (ret.type == 2) colnames(temp.matrix) <- subjects
      if (op$subjIds) save.subs <- subjects
        
      rm(pdata.ids)
    } else {
      subj.order2 <- getOrder(subj.order, subjects)
    }

    # Set the first row
    snpData[1] <- paste(subjects, sep="", collapse=out.sep)

    # For missing vector
    if (missFlag) temp.miss <- rep(FALSE, times=nsnps-1)

    # snp names vector
    if (sFlag) temp.snp <- character(nsnps-1)

    # MAF
    if (MAF.flag) {
      temp.MAF <- numeric(nsnps-1)
    } else {
      temp.MAF <- NULL
    }

    # alleles
    if (allFlag) {
      temp.all <- character(nsnps-1)
    } else {
      temp.all <- NULL
    }

    # Geno freqs
    if (gfreqFlag) {
      temp.gfreq <- character(nsnps-1)
    } else {
      temp.gfreq <- NULL
    }

    # n.miss
    if (gfreqFlag) {
      temp.nmiss <- numeric(nsnps-1)
    } else {
      temp.nmiss <- NULL
    }

    for (i in 2:nsnps) {

      temp <- getVecFromStr(snpData[i], delimiter=delimiter)

      # Remove the snp name
      snp.name <- temp[1]
      temp     <- temp[-1]
      if (sFlag) temp.snp[i-1] <- snp.name

      # Check for error
      if (length(temp) != total.nsubjects) stop(paste("ERROR with SNP", snp.name))

      # Use the 0, 1, 2 codes
      temp2 <- recode.geno(temp, in.miss=in.miss, out.miss=out.miss,
              out.genotypes=codes, heter.codes=heter.codes, subset=control.ids)

      # Vector of recoded genotypes
      temp <- temp2$vec

      # Get the MAF
      if (MAF.flag) temp.MAF[i-1] <- getMAF(temp, sub.vec=control.ids, controls=TRUE) 

      if (allFlag) temp.all[i-1] <- temp2$alleles

      # Geno Frequencies
      if (gfreqFlag) temp.gfreq[i-1] <- paste(table(temp, exclude=NA), collapse="|", sep="")

      # n.miss
      if (nmissFlag) temp.nmiss[i-1] <- sum(is.na(temp))

      # Get the correct subjects
      if (subjFlag) temp <- temp[subj.ids]
  
      # Get the correct order
      temp <- temp[subj.order2]

      # Determine if the snp has missing values
      if (missFlag) {
        if (any(is.na(temp))) temp.miss[i-1] <- TRUE
      }

      if (ret.type == 2) {
        temp.matrix[i-1, ] <- convert.fun(temp)
      } else {
        snpData[i] <- paste(temp, sep="", collapse=out.sep)
        if (inc.snps) snpData[i] <- paste(snp.name, out.sep, snpData[i], sep="")
      }
    } # END: for (i in 2:nsnps)

    # Get the rows by searching for the snp names
    if (snpsFlag) {
      snpNames <- update.snpNames(snpNames, temp.snp)
      count    <- count + length(snpData) - 1

      if (!length(snpNames)) stopFlag <- 1
    } # END: if (snpsFlag)

    # Update
    if (ret.type == 2) {
      rownames(temp.matrix) <- temp.snp
      data.obj <- convert.fun(rbind(data.obj, convert.fun(temp.matrix)))
    } else {
      if (firstFlag) snpData <- snpData[-1]
      data.obj  <- c(data.obj, snpData)
    }
    if (missFlag) missing <- c(missing, temp.miss)
    if (sFlag) snps <- c(snps, temp.snp)
    if (MAF.flag) MAF <- c(MAF, temp.MAF)
    if (allFlag) alleles <- c(alleles, temp.all)
    if (gfreqFlag) genoFreqs <- c(genoFreqs, temp.gfreq)
    if (nmissFlag) n.miss <- c(n.miss, temp.nmiss)
    firstFlag <- 1

    if (stopFlag) break
  } # END: for (file in snp.list$file) 

  # Check for error
  if (snpsFlag) {
    if ((total.nsnp != count) && (length(snpNames))) {
      print("Some SNPs were not found")
      if (op$stopOnError) {
        print(snpNames)
        stop("ERROR: The above SNPs were not found")
      }
    }
  }

  rm(subjects, snpNames, tList, subj.ids, pheno.id, pheno.uid, 
    temp.matrix, subj.order, total.nsubjects, temp.MAF, temp.all, temp2,
    temp.gfreq, temp.nmiss)
  temp <- gc(verbose=FALSE)

  if (!op$include.row1) data.obj <- data.obj[-1]

  list(data=data.obj, missing=missing, snpNames=snps, nsubjects=nsubj,
       subjects=save.subs, order=ret.order, phenoData.list=phenoData.list,
       MAF=MAF, alleles=alleles, n.miss=n.miss, genoFreqs=genoFreqs)

} # END: loadData.type1

# Function to return phenotype data and file id
loadData.stream <- function(snp.list, pheno.list, temp.list, op=NULL) {

  snp.list$stream <- 1
  op <- default.list(op, c("file.index", "useControls"), 
                     list(1, 1), error=c(0, 0))
  index <- op$file.index

  if (index == 1) {
    # Check the lists. pheno.list is checked in getPheno.info
    snp.list  <- check.snp.list(snp.list)
    temp.list <- check.temp.list(temp.list)

    # Update snp.list for snpNames
    snp.list <- update.snp.list(snp.list)

    # Get phenotype information
    pinfo <- getPheno.info(pheno.list, snp.list, temp.list=temp.list)
  } else {
    pinfo <- NULL
  } 
 
  ccvar <- getListName(pheno.list, "cc.var")
  if (is.null(ccvar)) {
    #print("##########################################################")
    #print("WARNING in loadData.stream: pheno.list$cc.var not specified")
    #print("##########################################################")
    op$useControls <- 0
  }
  if (op$useControls == 0) pinfo$controls <- NULL

  # Define a list to call loadData
  tList <- getSnpLoadList(snp.list, temp.list, op=op)

  temp <- paste(snp.list$dir, snp.list$file[index], sep="")

  # Open the data and get the subject ids
  fid <- loadData(temp, tList)

  list(phenoData.list=pinfo$phenoData.list, pheno.id=pinfo$pheno.id,
      pheno.uid=pinfo$pheno.uid, fid=fid, controls=pinfo$controls)

} # END: loadData.stream

# Function to get the data or the file id
getData.1 <- function(snp.list, pheno.list, temp.list, op=NULL) {

  # Check the lists
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)

  # Check the delimiters
  temp <- snp.list
  temp$file <- temp$file[1]
  if (!checkDelimiter(temp)) {
    stop("ERROR: Check the delimiter in snp.list")
  }
  rm(temp)
  if (!checkDelimiter(pheno.list)) {
    stop("ERROR: Check the delimiter in pheno.list")
  }

  if (snp.list$stream) {
    ret <- loadData.stream(snp.list, pheno.list, temp.list, op=op)
  } else {
    ret <- loadData.type1(snp.list, pheno.list, temp.list, op=op)
  }

  ret

} # END: getData.1

# Function to get the next observation 
scanNextObs <- function(fid, fid.row, get.row, snpFlag=0, sep="|",
                      snpNames=NULL) {
 
  if (!snpFlag) {
    row <- scan(file=fid, what="character", nlines=1, sep=sep,
                      skip=get.row-fid.row, quiet=TRUE)
  } else {
    # Loop and check each snp name
    cont  <- 1
    while (cont) {
      row <- scan(file=fid, what="character", nlines=1, sep=sep, quiet=TRUE)
      if (!length(row)) break
      if (row[1] %in% snpNames) break 
    }
  }

  row

} # END: scanNextObs

# Function to get the next observation and close the file if no more
#   reading is to be done
getNextObs <- function(i, snpFid, snpFlag, snpNames, tempfile, delete,
                  start.stream, stop.stream, delimiter) {

  # Check the snpNames vector
  if ( ((snpFlag) && (!length(snpNames))) || (i > stop.stream) ) {
    closeFile(snpFid, file=tempfile, delete=delete)
    return(NULL)
  }

  # Read in the next snp
  snp <- scanNextObs(snpFid, start.stream, start.stream, 
           sep=delimiter, snpFlag=snpFlag, snpNames=snpNames)
  if (!length(snp)) {
    closeFile(snpFid, file=tempfile, delete=delete)
    return(NULL)
  }
     
  snp

} # END: getNextObs

# Function to order and subset a SNP
orderSNP <- function(snp, snp.name, subj.order2, total.nsubjects, 
                     in.miss, heter.codes, subjFlag=0, subj.ids=NULL, 
                     out.genotypes=0:2, control.ids=NULL) {

  # Check for error
  if (length(snp) != total.nsubjects) stop(paste("ERROR with SNP", snp.name))
 
  # Use the 0, 1, 2 codes
  temp <- recode.geno(snp, in.miss=in.miss, out.miss=NA,
                 out.genotypes=out.genotypes, heter.codes=heter.codes,
                 subset=control.ids)
  snp <- temp$vec
  alleles <- temp$alleles

  # Compute the MAF
  maf <- getMAF(snp, sub.vec=control.ids, controls=TRUE)

  # Get the correct subjects
  if (subjFlag) snp <- snp[subj.ids]
  
  # Get the correct order
  snp <- snp[subj.order2]

  list(SNP=snp, MAF=maf, alleles=alleles)

} # END: orderSNP

# Function to call for stream input after getData.1 was called
setUp.stream <- function(obj, snp.list, tempfile, delete, controls=NULL) {
  
  # obj        Return object from getData.1
  
  pheno.uid <- obj$pheno.uid
  pheno.id  <- obj$pheno.id
  snpFid    <- obj$fid
  start.vec <- snp.list[["start.vec", exact=TRUE]]
  stop.vec  <- snp.list[["stop.vec", exact=TRUE]]
  start     <- max(2, start.vec[1])
  # Get stop in terms of i = 1, 2, ....
  stop      <- stop.vec[1] - start + 1
  if (stop < 0) stop <- Inf
  delimiter <- getDelimiter(snp.list)
  snames    <- snp.list[["snpNames", exact=TRUE]]
  snpFlag   <- 1 - is.null(snames)

  # Read in the subjects and first snp
  temp <- scanNextObs(snpFid, 1, 1, sep=delimiter)
  snp  <- scanNextObs(snpFid, 2, start, sep=delimiter, 
                     snpFlag=snpFlag, snpNames=snames)
  if (!length(snp)) {
     closeFile(snpFid, file=tempfile, delete=delete)
     stop("ERROR: Check snp.list options")
  }
 
  # Get the subjects id vector and order
  temp            <- getSubjIds(temp, snp, delimiter, pheno.uid, controls=controls)
  subjects        <- temp$subjects
  subjFlag        <- temp$subjFlag
  subj.ids        <- temp$subj.ids
  total.nsubjects <- temp$total.nsubjects
  control.ids     <- getListName(temp, "control.ids")
  if (!any(control.ids)) {
    #print("##########################################################")
    #print("NOTE: No controls found to determine the minor allele")
    #print("##########################################################")
    control.ids <- NULL
  }
  subj.order2     <- getOrder(pheno.id, subjects)

  list(subjFlag=subjFlag, subj.ids=subj.ids, start.stream=start,
       total.nsubjects=total.nsubjects, subj.order2=subj.order2,
       snp=snp, stop.stream=stop, control.ids=control.ids)

} # END: setUp.stream

# Function to get info from the phenotype data
getPheno.info <- function(pheno.list, snp.list, temp.list=NULL) {

  id.var <- pheno.list$id.var

  # Load the phenotype data
  temp <- list(file.type=pheno.list$file.type, header=pheno.list$header,
               delimiter=pheno.list$delimiter, id.var=id.var,
               remove.miss=0, in.miss=pheno.list$in.miss)
  data <- loadData(pheno.list$file, temp)

  # Get all of the control ids
  ccvar <- pheno.list[["cc.var", exact=TRUE]]
  if (!is.null(ccvar)) {
    # Check that it is 0-1
    temp <- (data[, ccvar] %in% c(0, 1))
    temp[is.na(temp)] <- TRUE
    if (!all(temp)) stop("ERROR in getPheno.info: pheno.list$cc.var is not coded as 0-1")

    temp     <- (data[, ccvar] == 0)
    temp[is.na(temp)] <- FALSE
    controls <- unique(makeVector(data[temp, id.var]))
  } else {
    controls <- NULL
  }

  # Update pheno.list
  pheno.list$data        <- data
  pheno.list$is.the.data <- 1

  # Get the intersecting subject ids with the genotype data
  ids <- intersectSubIds(snp.list, pheno.list, temp.list=temp.list)

  # Keep these ids
  temp <- data[, id.var] %in% ids
  temp[is.na(temp)] <- FALSE
  data <- removeOrKeepRows(data, temp)
  pheno.list$data <- data
  rm(data, ids)
  gc()

  # Call getPhenoData for the other pheno.list options to get the 
  # data to be used in the analysis
  getPheno.list <- getPhenoData(pheno.list, temp.list=temp.list)
  phenoData     <- getPheno.list$data
  pheno.list$data <- NULL
  pdata.flag    <- !is.null(getPheno.list$orig.id)
  pheno.list$orig.id <- getPheno.list$orig.id

  # pdata.flag is the flag for non-unique subject ids
  if (!pdata.flag) {
    # Get the subject ids
    pheno.id  <- rownames(phenoData) 
    pdata.ids <- pheno.id
  } else {
    pheno.id   <- as.character(phenoData[, pheno.list$orig.id])
    pdata.ids  <- rownames(phenoData) 
  }

  # Get the total number of subjects
  nsubj <- length(pheno.id)

  # Get the unique ids
  pheno.uid <- unique(pheno.id)
   
  # Get the number of unique subjects
  nsubj.u <- length(pheno.uid)

  # Return a list of info
  list(pheno.id=pheno.id, pheno.uid=pheno.uid, pdata.ids=pdata.ids,
       pdata.flag=pdata.flag, nsubj=nsubj, nsubj.u=nsubj.u,
       phenoData.list=getPheno.list, controls=controls)

} # END: getPheno.info

# Function to check for errors is the subject ids
getSubjIds <- function(snpData1, snpData2, delimiter, pheno.uid, controls=NULL) {

  # snpData1      Row 1      
  # snpData2      Row 2
  # delimiter     delimiter
  # pheno.uid     Unique (original) phenotype ids
  # controls      Unique (original) control ids

  # Get the subject ids from the snp data
  subs <- getSubject.vec(snpData1, snpData2, delimiter)
  total.nsubjects <- length(subs)

  # Get the number of unique original subject ids
  nsubj.u <- length(pheno.uid)

  # Check for an error
  if (sum(pheno.uid %in% subs) != nsubj.u) {
    temp <- !(pheno.uid %in% subs)
    temp <- pheno.uid[temp]
    #bad.subject.ids.global <<- temp
    print(temp)
    #print("bad.subject.ids.global")
    stop("The above subject ids were not found in the genotype data")
  }

  # Determine if any subjects should be removed
  subj.ids <- subs %in% pheno.uid

  if (sum(subj.ids) != nsubj.u) {
    print("The subject ids in the genotype data may not all be unique")
    stop("ERROR with subject ids")
  }
  
  # Get the logical vector for controls. Recall: MAF is determined by the 
  #  entire sample of subjects.
  if (!is.null(controls)) {
    control.ids <- subs %in% controls
  } else {
    control.ids <- NULL
  }

  if (total.nsubjects == nsubj.u) {
    # All subjects are here, so subj.ids is no longer needed
    subj.ids <- NULL
    subjFlag <- 0
  } else {
    subjFlag <- 1
    subs     <- subs[subj.ids]
  }

  # Return list
  list(subjects=subs, subjFlag=subjFlag, subj.ids=subj.ids,
       total.nsubjects=total.nsubjects, control.ids=control.ids)

} # END: getSubjIds

# Function to return the vector of subject ids from a snp file
getSubject.vec <- function(snpData1, snpData2, delimiter) {

  # snpData1    Row 1 (header) 
  # snpData2    Row 2
  # delimiter

  # Get the subject ids from the snp data
  if (length(snpData1) == 1) {
    subs <- getVecFromStr(snpData1, delimiter=delimiter)
  } else {
    subs <- snpData1
  }

  # Get the number of fields to remove
  n.omit <- getRow1.omitN(subs, snpData2, delimiter=delimiter) 
  if (n.omit > 0) subs <- subs[-(1:n.omit)] 

  subs

} # END: getSubject.vec

# Function to determine the number of fields to remove in row1 of the snp
#  data 
getRow1.omitN <- function(row1, row2, delimiter="|") {

  # row1  
  # row2  
  # delimiter

  n1 <- length(row1)
  if (length(row2) == 1) row2 <- getVecFromStr(row2, delimiter=delimiter)
  n2 <- length(row2)

  n.omit <- n1 - n2 + 1
  if (n.omit < 0) stop("ERROR: in the snp files")
  if (n.omit > 1) warning("Possible error in the SNP data")

  n.omit
 
} # END: getRow1.omitN

# Function to return the delimiter used in the file read in
getDelimiter <- function(snp.list, output=0) {

  delimiter <- snp.list$delimiter
  if ((output == 1) && (!snp.list$stream)) delimiter <- snp.list$out.delimiter
  if (snp.list$file.type %in% c(4)) {
    delimiter <- "|"
  } 
  delimiter

} # END: getDelimiter

# Function to get the vector of missing values
get.in.miss <- function(snp.list) {

  ret <- snp.list$in.miss
  if (snp.list$file.type == 4) ret <- c(ret, "-9")
  ret

} # END: get.in.miss

# Function to update the snpNames vector when snp names are specified
update.snpNames <- function(snpNames, temp.snp) {

  # snpNames    Vector of all snp names desired
  # temp.snp    Vector of snp names on current data

  # Remove the snp names from the snpNames vector
  rows     <- as.logical(1 - (snpNames %in% temp.snp))
  snpNames <- snpNames[rows]
  
  snpNames

} # END: update.snpNames

# Function to return a list of parameter options for loading the 
#  snp data.
getSnpLoadList <- function(snp.list, temp.list, op=NULL) {

  op <- default.list(op, c("file.index"), list(1))

  ret <- list(file.type=snp.list$file.type, 
         delimiter=snp.list$delimiter, read.n=snp.list$read.n,
         sas.list=snp.list$sas.list, transpose=1, include.row1=1,
         id.var=snp.list$id.var, 
         start.row=snp.list$start.vec[op$file.index],
         stop.row=snp.list$stop.vec[op$file.index],
         snpNames.keep=snp.list$snpNames.keep)
  ret$snpNames <- getListName(snp.list, "snpNames")
  stream <- snp.list$stream
  if (stream) {
    ret$stream  <- 1
    ret$outfile <- op$outfile
  }
  if (ret$file.type == 4) ret$sas.list$temp.list <- temp.list
  
  ret

} # END: getSnpLoadList

# Function to get the common subject ids.
# Returns the updated pheno.list
intersectSubIds <- function(snp.list, pheno.list, temp.list=NULL) {

  snp.list$snpNames      <- NULL
  snp.list$snpNames.list <- NULL
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)
  temp.list  <- check.temp.list(temp.list)
  delimiter  <- getDelimiter(snp.list)
  snp.list   <- update.snp.list(snp.list)
  stream     <- getListName(snp.list, "stream")

  # Get the unique files
  files <- unique(snp.list$file)

  # Define a list to call loadData
  tList           <- getSnpLoadList(snp.list, temp.list)
  tList$start.row <- 1
  tList$stop.row  <- 2
  isub            <- NULL
  index           <- 1

  # Loop over each file and combine the objects
  for (file in files) {

    # Get the data
    temp <- loadData(paste(snp.list$dir, file, sep=""), tList)

    # Get the first 2 rows
    if (!stream) {
      row1 <- temp[1]
      row2 <- temp[2]
    } else {
      row1 <- scanNextObs(temp, 1, 1, sep=delimiter, snpFlag=0, 
                          snpNames=NULL)
      row2 <- scanNextObs(temp, 1, 1, sep=delimiter, snpFlag=0, 
                          snpNames=NULL)
      close(temp)
    }

    # Get the subject ids
    temp <- getSubject.vec(row1, row2, delimiter)

    if (index == 1) {
      isub <- temp
    } else {
      isub <- intersect(isub, temp)
    }
  }

  if (!length(isub)) {
    stop("ERROR: No intersecting subject ids in the genotype files")
  }

  # Get the phenotype data
  dflag <- !is.null(getListName(pheno.list, "is.the.data"))
  if (dflag) {
    temp <- unique(makeVector(pheno.list$data[, pheno.list$id.var]))
  } else {
    temp <- NULL # ADD code here
  }

  isub <- intersect(isub, temp)
  if (!length(isub)) {
    stop("ERROR: No intersecting subject ids in the genotype and phenotype files")
  }

  if (length(temp) != length(isub)) {
    temp <- temp[!(temp %in% isub)]
    print(temp)
    #bad.subject.ids.global <<- temp
    print("The above subject ids were not found in the genotype data")
    warning("Some subject ids were not found in the genotype data")
  }

  isub

} # END: intersectSubIds




# History Mar 28 2008  Allow infile to be a connection in scanFile
#                      In readTable, call createDummy using the factors
#                      variables.
#         Mar 31 2008  Change snpNames option    
#         Apr 01 2008  Add option in loadData for file.type = 3 to call
#                      scanFile instead of readTable.  
#         Apr 03 2008  Change the option creatDummy to make.dummy in readTable.
#         Apr 04 2008  Add in.miss to readTable
#                      Remove makeIdsNames from readTable.
#         Apr 11 2008  Add getNcols function.
#                      Use scanFile in table2Format
#         Apr 15 2008  Change in sas.list
#         May 12 2008  Add remove.miss to readTable
#                      Add baseline option for readTable
#         Jun 27 2008  Add function getColumns
#         Jun 27 2008  Use the tempfile function for generating temporary
#                      file names
#         Jun 28 2008  Do not extract snp names in loadData
#         Jun 30 2008  Add code for streamed input
#         Jul 02 2008  Add code for type 5 (.zip) files
#         Jul 03 2008  Add code for type 6, 8
#                      Add getNrows function
#         Jul 08 2008  Add changeLevels and subsetData options to readTable
#                      Unfactor all variables in the beginning of readTable
#         Jul 09 2008  Remove baseline option from readTable and
#                      keep all factors
#         Sep 17 2008  In loadData, check for outfile when it is needed.
#         Sep 24 2008  In getNcols, add option to return the columns
#         Mar 05 2009  Generalize getColumns function
#         Apr 20 2009  Change getNrows function
#         Jun 12 2009  Change default delimiter in loadData to "\t" 
#         Jun 15 2009  In loadData, method=2 file.type=3,6,8 change
#                      default returnMatrix to 1    
#         Jun 15 2009  For snpNames option in loadData, use the substring
#                      option in extractByStr   
#         Jul 20 2009  Fix bug in getNrows. Set delimiter to "\n".   
#         Aug 03 2009  Change in readTable: set as.is=FALSE in readTable
#                      Remove changeLevels option    
#         Aug 10 2009  Pass snp.list option snpNames.keep into extractByStr function  
#         Aug 19 2009  Add checkVars function   
#         Aug 27 2009  Set as.is=TRUE in readTable
#         Oct 16 2009  Changes for set up:
#                      readTable
#         Dec 24 2009  In readTable, call unfactor.all for file.type=1
#         Dec 28 2009  Add function getFileType
#                      Use getFileType in loadData 
#         Dec 29 2009  Add getFileDelim function
#         Jan 06 2010  Call getFileType in getFID
#         Jan 19 2010  Add function loadData.table
#         Mar 23 2010  Update getFileDelim

# Function to load and return the data object (or subset)
loadData <- function(infile, p) {

  # infile         The R object file containing the data. Only 1 object should
  #                be in the file. The load() function is used to read in
  #                the data.
  #                No default
  #######################################################################
  # p         A list with the folowing names:
  #
  # file.type      1-8  Type 1 is for a file created with the save()
  #                command. Type 2 is for other flat files that have row 1
  #                as subject ids and then the snps as rows.
  #                Type 3 is for tables where column 1 is the subject id
  #                and the remining columns are snps. See delimiter,
  #                transpose and method below.
  #                Type 4 is a sas data set. See sas.list below.
  #                No default
  # start.row      The starting row. The default value is 1
  # stop.row       The end row or -1 to return all the data from row
  #                start.row onwards.
  #                The default is -1
  # include.row1   0 or 1 to include the first row of the data. The
  #                first row may be a header containing variable names.
  #                This option overrides start.row
  # snpNames       Character vector of the snp names to use. If NULL, then
  #                start.row and stop.row will be used.
  #                The default is NULL.
  # snpNames.keep  0 or 1 to remove or keep the snps in snpNames
  #                The default is 1 (keep).
  # read.n         Integer   This option is only used with file.type=3.
  #                It specifies the maximum number of lines to read in
  #                each iteration.
  #                The default value is -1, so that the entire table will
  #                be read in all at once.
  # delimiter      This option is only for file.type = 3. 
  #                The default is "|".
  # transpose      0 or 1 to transpose the data. Set to 1 for the snp data.
  #                This option is only for file.type = 3 and 4.
  #                The default is 0.
  # method         1 or 2  This option is only for file.type = 3 and 
  #                transpose = 0.
  #                1: call readTable
  #                2: call scanFile
  #                The default is 1.
  # stream         0 or 1 for streamed input
  #                The default is 0
  # outfile        Name of the output (temporary) file if one needs it to
  #                be created with stream = 1.
  #                The default is NULL. 
  # id.var         For type 3 and 4 only
  #                No default
  # For file.type = 4
  # temp.list
  # sas.list       Only for file.type = 4. A list with the following names:
  #         sas.exe
  #         sas.file
  #         shell
  #########################################################################

  # Set defaults
  p <- default.list(p, c("start.row", "stop.row", 
         "include.row1", "read.n", "delimiter", "transpose", 
         "stream", "snpNames.keep"), 
       list(1, -1, 1, -1, "\t", 0, 0, 1), 
       error=c(0, 0, 0, 0, 0, 0, 0, 0))

  type <- p[["file.type", exact=TRUE]]
  if (is.null(type)) type <- getFileType(infile, default=-1)
  if (type == -1) stop("ERROR in loadData: file.type must be specified")
  p$file.type <- type

  # Check for errors with start.row and stop.row
  if ((!p$include.row1) && (p$start.row == 1)) p$start.row <- 2

  # Check for errors
  if (p$start.row < 1) stop("ERROR: with start.row")
  if ((p$stop.row > 0) && (p$start.row > p$stop.row)) {
    stop("ERROR: with start.row and/or stop.row")
  }

  type      <- p$file.type
  stream    <- p$stream
  outfile   <- p[["outfile", exact=TRUE]]
  p$outfile <- NULL

  # Set default for some file type
  temp <- p[["method", exact=TRUE]]
  if (is.null(temp)) {
    if (type %in% c(6, 8)) {
      p$method <- 2
    } else {
      p$method <- 1
    }
  } 
  if ((type %in% c(3, 6, 8)) & (p$method == 2)) {
    p <- default.list(p, c("returnMatrix", "what"), list(1, "character"))
  }

  # Load the data
  if (type == 1) {
    # R object file (variables are rows)
    dat <- loadFile(infile, p)

  } else if (type == 2) {
    # Flat file (variables are rows)
    # If stream, then return file id
    if (stream) return(file(infile, "r"))

    p$delimiter <- "\n"
    dat <- scanFile(infile, p)

  } else if (type == 3) {
    # Flat file (variables are columns)
    if (p$transpose) {
      p$include.header <- p$include.row1
      p$vars           <- p[["snpNames", exact=TRUE]]
      p$start.col      <- p$start.row
      p$stop.col       <- p$stop.row
      dat              <- table2Format(infile, p)
    } else {
      if (p$method == 2) {
        dat <- scanFile(infile, p)
      } else {
        dat <- readTable(infile, p)
      }
    }
  } else if (type == 4) {
    # SAS data set
    p$vars            <- p$snpNames
    p$start.col       <- p$start.row
    p$stop.col        <- p$stop.row
    p$sas.list$id.var <- p$id.var

    dat <- getSASData(infile, p, transpose=p$transpose) 

  } else if (type == 5) {
    # Type 2 compressed zip file
    dat <- readFile.zip2(infile, p)
    if (stream) return(dat)

  } else if (type == 7) {
    # Type 2 compressed gz file
    dat <- readFile.gz2(infile, p)
    if (stream) return(dat)

  } else if (type %in% c(6, 8)) {
    # Type 3 compressed file
    dat <- readFile.gzip3(infile, p)
    if (stream) return(dat)

  } else {
    # GLU format
    dat <- glu.transform(infile, inFormat=type)
    return(dat)
  } 

  # Extract snps
  temp <- p[["snpNames", exact=TRUE]]
  if (!is.null(temp)) {
    op <- list(include.row1=p$include.row1, substr.vec=c(1, 15), 
               keep=p$snpNames.keep)
    dat <- extractByStr(dat, temp, op=op)
  }

  # For streamed input
  if (stream) {
    if (is.null(outfile)) stop("ERROR: outfile must be specified")
    writeLines(dat, con=outfile)
    dat <- file(outfile, "r")
  }

  dat

} # END: loadData

# Function to load an object from a file that was created with the 
#  save() function
loadFile <- function(infile, p) {

  # infile       Path to the file to load
  ############################################################
  # p      A list of options with the names:
  # name         The name of the object in the file
  #              For now, the file should only contain 1 object
  #              The default value is "dat"
  # start.row    The default is 1
  # stop.row     The default is -1
  # include.row1 The default is 1

  p <- default.list(p, c("start.row", "stop.row", "include.row1", 
                        "name", "stream"), 
       list(1, -1, 1, "dat", 0))

  # Initialize 
  dat <- NULL

  # Load the object
  temp <- load(infile)
 
  if (length(temp) > 1) stop("ERROR: infile contains more than 1 object")

  # Rename the object
  if (temp != p$name) {
    dat <- eval(parse(text=temp))
    eval(parse(text=paste("rm(", temp, ")", sep="")))
    temp <- gc(verbose=FALSE)
  }

  # Get the number of rows
  n <- length(dat)

  rows <- NULL
  if ((p$stop.row < 0) || (p$stop.row > n)) p$stop.row <- n
  if ((p$include.row1) && (p$start.row <= 2)) {
    p$start.row    <- 1
    p$include.row1 <- 0
  } 

  if ((p$stop.row < n) || (p$start.row > 1)) {
    rows <- (p$start.row):(p$stop.row)
  }
  if (p$include.row1) rows <- c(1, rows)
  temp <- length(rows)
  if ((temp) && (temp != n)) dat <- dat[rows]

  dat

} # END: loadFile

# Function to read a data set using scan()
scanFile <- function(infile, p) {

  # infile        File to read. Infile can also be a connection.
  #               If it is a connection, then include.row1 gets set to
  #               0 if start.row > 2 and start.row gets set to 1.
  #               No default
  #################################################################
  # p      A list with the following names
  # start.row     Starting row to read.
  #               The default is 1.
  # stop.row      The stopping row to read.
  #               The default is -1, so that up to the end of the file
  #               will be read.
  # include.row1  0 or 1 to include row1 
  #               The default is 1.
  # what          Type of data to read
  #               The default is "character"
  # delimiter     The delimiter
  #               The default is "\n"
  # returnMatrix  0 or 1 to return a matrix. If set to 1 and include.row1=1,
  #               then it is assumed that row 1 is a header and the column
  #               names will be assigned.
  #               The default is 0 
  # ncol          If returnMatrix = 1 and infile is a connection, then
  #               ncol must be set to the number of columns in the returned
  #               matrix. ncol does not need to be specified otherwise.
  #               The default is NULL
  ##################################################################

  p <- default.list(p, c("start.row", "stop.row", "what", 
         "include.row1", "delimiter", "returnMatrix"), 
       list(1, -1, "character", 1, "\n", 0), 
       error=c(0, 0, 0, 0, 0, 0))
    
  start.row    <- floor(p$start.row)
  stop.row     <- floor(p$stop.row)
  include.row1 <- p$include.row1
  connFlag     <- 0
 
  if (stop.row > 0) {
    if ((start.row < 0) || (start.row > stop.row)) start.row <- 1
    if ((stop.row == 1) && (include.row1)) include.row1 <- 0
  }

  row1Flag <- 0
  temp <- (p$what == "character")
  if (!length(temp)) temp <- 0
  if (temp) {
    if (include.row1) {
      row1Flag <- 1
      if (start.row == 1) start.row <- 2
    }
    #if ((include.row1) && (start.row == 2)) start.row <- 1
    #if ((include.row1) && (start.row > 1)) row1Flag <- 1
  } else {
    if ((include.row1) && (p$returnMatrix)) row1Flag <- 1
    if ((include.row1) && (p$returnMatrix) && (start.row == 1)) start.row <- 2
  }

  # See if infile is a connection
  if ("connection" %in% class(infile)) {
    connFlag <- 1
    if ((p$returnMatrix) && (is.null(p$ncol))) stop("ERROR: ncol must be set")
    if (include.row1) {
      if (start.row <= 2) {
        row1Flag  <- 1
        start.row <- 1
      } else {
        row1Flag     <- 0
        include.row1 <- 0
      }
    }
    skip   <- start.row - 1
    nlines <- stop.row - start.row + 1
    if ((row1Flag) && (start.row == 1)) nlines <- nlines - 1
    if (nlines == 0) nlines <- 1
  } else {
    skip   <- start.row - 1
    nlines <- stop.row - start.row + 1 
  } 

  # Get row1
  if (row1Flag) {
    row1 <- scan(file=infile, what="character", nlines=1, 
                 sep=p$delimiter, quiet=TRUE)
  } else {
    row1 <- NULL
  }

  # Get the rest
  dat <- scan(file=infile, what=p$what, skip=skip, quiet=TRUE,
                 nlines=nlines, sep=p$delimiter)

  if (p$returnMatrix) {
    # Get the number of columns
    if (!is.null(p$nc)) {
      nc <- p$nc
    } else {
      if (!include.row1) {
        nc <- length(scan(file=infile, what=p$what, nlines=1, 
                     sep=p$delimiter, quiet=TRUE))
      } else {
        nc <- length(row1)
      }
    }

    # Convert to a matrix
    dat <- matrix(data=dat, ncol=nc, byrow=TRUE)

    # Assign column names
    if (include.row1) colnames(dat) <- row1
  } else {
    if (!is.null(row1)) dat <- c(row1, dat)
  }

  dat

} # END: scanFile

# Function to read in a table and convert data to the proper format
table2Format <- function(infile, p) {

  # infile         Path to the file. The file must have a header.
  #                infile can also be a data frame or matrix
  #                No default
  ##################################################################
  # p        A list with the following names:
  # delimiter      Delimiter in infile (infile is a file)
  #                The default is " "
  # out.delimiter  Delimiter to use in transposed data
  #                The default is "|"
  # id.var         Name or column number of the id variable
  #                No default
  # include.header 0 or 1 to include the header as the first column
  #                in the transposed data.
  #                The default is 1.
  # read.n         Number of rows to read each iteration.
  #                The default is that all rows will be read.
  # vars           Vector of variable names or column numbers to read.
  #                If column numbers, then the vector must be numeric.
  #                The default is NULL
  # start.col      Starting column to read
  #                The default is 1
  # stop.col       Ending column to read.
  #                The default is -1.
  # outfile        Output file
  #                The default is NULL
  # out.type       1 = compressed R object file (using save() function)
  #                2 = flat file  
  #                The default is 1.  
  ###################################################################

  # Local function
  f2 <- function(col) {
    paste(col, sep="", collapse=delimiter)
  }

  p <- default.list(p, c("id.var", "start.col", "stop.col", "read.n", 
         "include.header", "out.delimiter", "out.type", "delimiter"), 
       list("ERROR", 1, -1, -1, 1, "\t", 1, " "),
        error=c(1, 0, 0, 0, 0, 0, 0, 0))

  # Initialize 
  include.header <- p$include.header
  delimiter      <- p$out.delimiter
  start.col      <- p$start.col
  stop.col       <- p$stop.col
  vars           <- p$vars
  stop           <- 0
  index          <- 1
  ids            <- NULL
  rflag          <- 0
  readFlag       <- 0
  idFlag         <- 0
  fid            <- NULL
  tlist          <- list()

  if ((!is.data.frame(infile)) && (!is.matrix(infile))) {
    # Get the number of columns
    temp <- getNcols(infile, p)

    # Open a connection
    fid <- file(infile, "r")

    # Define a list for scanFile
    tlist <- list(what="character", returnMatrix=1, delimiter=p$delimiter,
                  include.row1=1, start.row=1, stop.row=p$read.n,
                  ncol=temp)
  }

  while (!stop) {
    if (!is.null(fid)) {
      # Read in the data
      infile <- scanFile(fid, tlist)
    
      # Check for error
      if (class(infile) == "try-error") {
        if (readFlag) {
          break
        } else {
          # Set read.n and continue
          if (read.n < 1) {
            read.n <- 1000
          } else {
            read.n <- ceiling(read.n/2)
          }
          next
        }
      } else {
        if (!length(infile)) break
      }
    } # END: if (!is.null(fid))
    readFlag <- 1
    tlist$include.row1 <- 0

    if (index == 1) {
      nc <- ncol(infile)
      cnames <- colnames(infile)
      if (nc < 2) stop("ERROR: the table has too few columns")

      # Get the name of the id variable
      if (is.numeric(p$id.var)) {
        if ((p$id.var < 1) || (p$id.var > nc)) stop("ERROR with id.var")
        idname <- cnames[p$id.var]
      } else {
        idname   <- p$id.var
        p$id.var <- match(idname, cnames)
        if (is.na(p$id.var)) stop("ERROR: with id.var") 
      }

      # Get the column numbers for the vars
      if (!is.null(vars)) {
        rflag <- 1
        if (!is.numeric(vars)) {
          vars <- match(vars, cnames)

          # Remove NAs
          vars <- vars[!is.na(vars)]
          
          if (!length(vars)) {
            if (!is.null(fid)) close(fid)
            return("")
          }
 
        } else {
          if ((any(vars < 1)) || (any(vars > nc))) 
            stop("ERROR: column numbers are not correct")
        }

        # Set start row
        start.col <- min(vars)
      } else {
        if (start.col < 1) start.col <- 1
        if ((start.col == 2) && (include.header)) start.col <- 1
        if (stop.col < 1) stop.col   <- nc
        if (stop.col > nc) stop.col  <- nc
        #if ((start.col==1) && (stop.col==1)) include.header <- 1

        if ((stop.col < nc) || (start.col > 1)) {
          rflag <- 1
          vars  <- start.col:stop.col 
        }
      }
      # See if the id variable is the only variable
      if ((length(vars) == 1) && (vars == p$id.var)) include.header <- 1

    } # END if (index == 1)

    # Get the ids  
    ids <- c(ids, as.character(infile[, p$id.var]))
  
    # Keep columns
    if (rflag) {
      infile <- removeOrKeepCols(infile, vars, which=1)
      if (index == 1) {
        cnames <- colnames(infile)
        nc     <- ncol(infile)
        if (idname %in% cnames) {
          idFlag <- 1
          idPos  <- match(idname, cnames)
        }
      }
    } else {
      idFlag <- 1
      idPos  <- p$id.var
    }
 
    # Remove the id column if it was not removed before
    if (idFlag) {
      infile <- removeOrKeepCols(infile, idPos, which=-1)
      if (index == 1) {
        cnames <- colnames(infile)
        nc     <- length(cnames)
      }
    }

    # Get rows 
    temp <- try(apply(infile, 2, f2), silent=TRUE)

    # Check for error
    if (class(temp) == "try-error") {
      # It failed, so loop over each column
      temp <- character(nc)

      for (i in 1:nc) temp[i] <- f2(infile[, i])
    } 

    if (index == 1) {
      data <- temp
    } else {
      data <- paste(data, temp, sep=delimiter)
    }

    if (is.null(fid)) stop <- 1
    index <- index + 1

  } # END: while (!stop)

  # Close the file
  if (!is.null(fid)) close(fid)

  rm(infile, vars, tlist, fid)
  temp <- gc(verbose=FALSE)

  # First row are the subject ids
  if (include.header) {
    data   <- c(f2(ids), data)
    cnames <- c(idname, cnames) 
  }

  # Add the snp ids
  data <- paste(cnames, data, sep=delimiter)

  # Save data
  temp <- getListName(p, "outfile")
  if (!is.null(temp)) {
    if (p$out.type == 1) {
      save(data, file=temp)
    } else {
      writeLines(data, con=temp)
    }
  } 

  data

} # END: table2Format

# Function to return a text file from a sas data set. 
# Missing values will be represented as -9. The delimiter will be "|".
getSASData <- function(sas.data, p, transpose=0) {

  # sas.data     The complete path to the SAS data set.
  #              No default
  # p            A list with the names "sas.list", "temp.list" and
  #              other names as in scanFile (for transpose = 1) or
  #              readTable (for transpose = 0)
  #              No default
  # transpose    0 or 1  1 is for retrieving the snp data.
  #              0 is for reading a sas table with the read.table()
  #              function.
  #              The default is 0.

  # Check the lists
  if (is.null(p$sas.list)) stop("ERROR: SAS options list is NULL")
  p$sas.list <- default.list(p$sas.list, 
   c("sas.exe", "sas.file", "id.var", "shell"),
   list("ERROR", "ERROR", "ERROR", "bash"),
    error=c(1, 1, 1, 0))
  p$sas.list$delimiter <- "|"

  p$sas.list$temp.list <- check.temp.list(p$sas.list$temp.list)

  # Initialize
  dir    <- p$sas.list$temp.list$dir
  id     <- p$sas.list$temp.list$id
  delete <- p$sas.list$temp.list$delete
  str    <- paste("_sas", id, "_", sep="")

  # Define the file names
  p$sas.list$outfile     <- getTempfile(dir, prefix=c(str, "out"), ext=".txt")
  p$sas.list$temp.file   <- getTempfile(dir, prefix=c(str, "tmp"), ext=".sas") 
  p$sas.list$temp.script <- getTempfile(dir, prefix=c(str, "bat"), ext=".bat") 
  p$sas.list$idFile      <- getTempfile(dir, prefix=c(str, "ids"), ext=".txt") 
  p$sas.list$delete      <- delete 

  if (transpose) {
    # Create the text files
    p$sas.list$out.miss <- -9
    ret <- SAS2Format(sas.data, p)

    # Read in the text file
    p$what         <- "character"
    p$delimiter    <- "\n"
    p$returnMatrix <- 0
    p$start.row    <- 1
    p$stop.row     <- -1
    p$snpNames     <- NULL
    dat            <- scanFile(p$sas.list$outfile, p)

    # Read in the file of subject ids
    ids <- scanFile(p$sas.list$idFile, p)
    ids <- paste(ids, sep="", collapse="|")

    # Combine
    dat <- c(ids, dat)

    # Delete file
    if (delete) file.remove(p$sas.list$idFile)

  } else {
    # Create the text file
    ret <- exportSAS(sas.data, p$sas.list)

    # Set the list
    p$file.type <- 3
    p$header    <- 1
    p$delimiter <- "|"
    p$in.miss   <- c(p$in.miss, "", " ", "  ")

    # Read the table
    dat <- readTable(p$sas.list$outfile, p) 
  }

  # Delete file
  if (delete) file.remove(p$sas.list$outfile)

  dat

} # END: getSASData

# Function to export a SAS data set
exportSAS <- function(sas.data, p) {

  # See SAS2Format for the options  

  p <- default.list(p, 
   c("sas.exe", "sas.file", "outfile", "temp.file", "delimiter", 
     "temp.script", "delete", "shell"),
   list("ERROR", "ERROR", "outfile.txt", "tempfile.sas", " ", 
        "temp.bat", 1, "bash"),
    error=c(1, 1, 0, 0, 0, 0, 0, 0))

  # Check the existence of the sas data set
  if (check.files(sas.data)) stop()

  # Check the existence of the sas source file
  if (check.files(p$sas.file)) stop()

  # Get the path to the library
  libname <- dirname(sas.data)

  # Get the name of the sas data set
  temp <- strsplit(basename(sas.data), ".", fixed=TRUE)[[1]]
  data <- temp[length(temp)-1] 

  # Open a connection to the temporary file
  fid <- try(file(p$temp.file, "w"), silent=TRUE)
  if (class(fid)[1] == "try-error") stop("ERROR: temp.file is invalid")
  
  temp <- paste('libname _tmp4125 "', libname, '"; \n', sep="")
  cat(temp, file=fid)

  # % include file
  temp <- paste('%include "', p$sas.file, '"; \n', sep="")
  cat(temp, file=fid)

  # Define macro variables
  temp <- paste("%let data = _tmp4125.", data, "; \n", sep="")
  cat(temp, file=fid)

  temp <- paste('%let outfile = "', p$outfile, '"; \n', sep="")
  cat(temp, file=fid)

  temp <- paste('%let delimiter = "', p$delimiter, '"; \n', sep="")
  cat(temp, file=fid)

  # Define the macro call
  temp <- "%export(data=&data, outfile=&outfile, \n"
  cat(temp, file=fid)
  temp <- "delimiter=&delimiter); \n"
  cat(temp, file=fid)

  # Close the file
  close(fid)

  # Create the batch file and call sas
  ret <- createScript(p$temp.script, p$sas.exe, p$temp.file, shell=p$shell,
                 delete=p$delete) 

  # Delete file
  if (p$delete) file.remove(p$temp.file)

  # Check the return code
  if (ret) stop("ERROR: calling exportSAS")

  ret

} # END: exportSAS

# Function to transform a sas data set to a text file, where the columns
# of the sas data set will be the rows.
SAS2Format <- function(sas.data, p) {

  # sas.data         The complete path to the SAS data set
  #                  No default
  ###################################################################
  # p    A list with the following names
  # start.col
  # stop.col
  # vars
  ################################################################
  # sas.list  A list with the following names
  # sas.exe          The complete path to the executable file to 
  #                  start SAS. If there is a command (eg sas) to start
  #                  SAS from any directory, then use the command instead.
  #                  No default
  # sas.file         The file containing the SAS macros needed to
  #                  convert the SAS data set.
  #                  No default
  # outfile          The output file that will contain the data in
  #                  the new format.
  #                  The default is "outfile.txt"
  # temp.file        The temporary file that will contain SAS commands
  #                  to run the SAS macro %table2Format.
  #                  The default is "tempfile.sas"
  # id.var           The id variable on the SAS data set
  #                  The default is NULL
  # delimiter        The delimiter to be used in the outfile.
  #                  The default is "|".
  # temp.script      The temporary script file that will call SAS.
  #                  The default is "temp.bat"
  # delete           0 or 1 to delete the temporary files temp.file and
  #                  temp.script.
  #                  The default is 1.
  # shell            The default is "bash"
  # out.miss         The numeric value to denote missing values.
  #                  Use "." for a regular SAS missing value
  #                  The default is -9999
  # idFile           Output file containing the SNP ids.
  #                  The default is "id.txt"
  ##################################################################

  p$sas.list <- default.list(p$sas.list, 
   c("sas.exe", "sas.file", "outfile", "temp.file", "delimiter", 
     "temp.script", "delete", "shell", "out.miss", "idFile"),
   list("ERROR", "ERROR", "outfile.txt", "tempfile.sas", "|", 
        "temp.bat", 1, "bash", -9999, "id.txt"),
    error=c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0))
 
  p <- default.list(p, c("start.col", "stop.col"), list(1, -1))
  s <- p$sas.list

  # Check the existence of the sas data set
  if (check.files(sas.data)) stop()

  # Check the existence of the sas source file
  if (check.files(s$sas.file)) stop()

  # Convert vars to a string
  if (!is.null(p$vars)) p$vars <- paste(p$vars, collapse=" ")

  # Get the path to the library
  libname <- dirname(sas.data)

  # Get the name of the sas data set
  temp <- strsplit(basename(sas.data), ".", fixed=TRUE)[[1]]
  data <- temp[length(temp)-1] 

  # Open a connection to the temporary file
  fid <- try(file(s$temp.file, "w"), silent=TRUE)
  if (class(fid)[1] == "try-error") stop("ERROR: temp.file is invalid")
  
  temp <- paste('libname _tmp4125 "', libname, '"; \n', sep="")
  cat(temp, file=fid)

  # % include file
  temp <- paste('%include "', s$sas.file, '"; \n', sep="")
  cat(temp, file=fid)

  # Define macro variables
  temp <- paste("%let data = _tmp4125.", data, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste('%let outfile = "', s$outfile, '"; \n', sep="")
  cat(temp, file=fid)
  temp <- paste("%let outMiss = ", as.numeric(s$out.miss), "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste("%let idVar = ", s$id.var, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste('%let delimiter = "', s$delimiter, '"; \n', sep="")
  cat(temp, file=fid)
  temp <- paste("%let startCol = ", p$start.col, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste("%let stopCol = ", p$stop.col, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste("%let vars = ", p$vars, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste('%let idFile = "', s$idFile, '"; \n', sep="")
  cat(temp, file=fid)

  # Define the macro call
  temp <- "%table2Format(data=&data, idVar=&idVar, outfile=&outfile, \n"
  cat(temp, file=fid)
  temp <- "outMiss=&outMiss, delimiter=&delimiter, idFile=&idFile,\n"
  cat(temp, file=fid)
  temp <- "startCol=&startCol, stopCol=&stopCol, vars=&vars); \n"
  cat(temp, file=fid)

  # Close the file
  close(fid)

  # Create the batch file and call sas
  ret <- createScript(s$temp.script, s$sas.exe, s$temp.file, shell=s$shell,
                 delete=s$delete) 

  # Delete file
  if (s$delete) file.remove(s$temp.file)

  # Check the return code
  if (ret) stop("ERROR: calling SAS2Format")

  ret
  
} # END: SAS2Format 

# Function to create and execute a script file
createScript <- function(script.file, exe, run.file, shell="bash",
                 delete=1) {

  # script.file   The batch file to create.
  #               No default
  # exe           The executable file
  #               No default
  # run.file      The file that the executable program will run
  #               No default
  # shell         Shell to use (UNIX)
  #               The default is "bash"
  # delete        0 or 1 to delete script.file
  #               The default is 1

  # Open the script file
  fid <- try(file(script.file, "w"), silent=TRUE)
  if (class(fid)[1] == "try-error") stop("ERROR: script file is invalid")

  # Get the path seperator
  sep <- .Platform$file.sep

  # Get the directory names of the executable program
  dirs <- strsplit(dirname(exe), sep, fixed=TRUE)[[1]]
  file <- basename(exe)

  # Determine the platform
  os      <- .Platform$OS.type
  winFlag <- (os == "windows")

  if (winFlag) {
    temp <- paste(dirs[1], " \n", sep="")
    cat(temp, file=fid)
    cat("cd \\ \n", file=fid)
    dirs <- dirs[-1]
  }
  for (dir in dirs) {
    temp <- paste("cd ", dir, " \n", sep="") 
    cat(temp, file=fid)
  }
  
  # Define the command to call sas
  temp <- paste(file, ' "', run.file, '"', sep='')
  cat(temp, file=fid)

  # Close the file
  close(fid)

  # Call the script file
  if (winFlag) {
    ret <- shell(script.file)
  } else {
    ret <- system(paste(shell, script.file, sep=" "))
  }

  # Delete the file
  if (delete) file.remove(script.file)
 
  # Check the return code
  if (ret > 1) stop("ERROR: with system call")

  0

} # END: createScript

# Function to read in a table (data frame)
readTable <- function(file, p) {

  # file           Data file. 
  #                No default.
  #################################################################
  # p is a list with the folowing fields:
  # file.type      1 or 3  1 is for an R object file created with the
  #                save() function. 3 is for a table that will be read in
  #                with read.table()
  #                The default is 3
  # header         0 or 1 . Set to 0 if the file does not contain a header.
  #                The default is 1.
  # delimiter      The delimiter used in the file.
  #                The default is "".
  # id.var         Name or column number of the id variable. 
  #                Use NULL if there is no id variable.
  #                No default.
  # keep.vars      Vector of variable names or column numbers to keep. 
  #                The default is NULL, so that all variables will be kept.
  # remove.vars    Vector of variable names or column numbers to remove. 
  #                The default is NULL, so that all variables will be kept.
  #                Both use.vars and remove.vars cannot be specified. 
  # factor.vars    Vector of variable names or column numbers to convert
  #                into factors.
  #                The default is NULL.
  # make.dummy     0 or 1 to create dummy variables for the factor variables.
  #                The default is 0.
  # keep.ids       Vector of ids or row numbers to keep.
  #                The default is NULL.
  # remove.ids     Vector of ids or row numbers to remove.
  #                Both keep.ids and remove.ids cannot be specified
  #                The default is NULL.
  # in.miss        Character vector of strings to define missing values
  #                in read.table (na.strings)
  #                The default is "NA"
  # remove.miss    0 or 1 to remove rows with missing values.
  #                Rows are removed after keep.vars or remove.vars has
  #                been applied.
  #                The default is 0.
  # subsetData     List of sublists. Each sublist must contain the names
  #                "var", "operator", and "value". 
  #                See the function subsetData.list in the file wga_util.R
  #                Ex: list(list(var="GENDER", operator"==", value="MALE"))
  #                The default is NULL
  ##################################################################

  # Set defaults
  p <- default.list(p, 
      c("header", "file.type", "delimiter", "make.dummy", "in.miss",
        "remove.miss", "is.the.data"),
      list(1, 3, "", 0, "NA", 0, 0))

  # Check for errors
  if ((!is.null(p$keep.vars)) && (!is.null(p$remove.vars))) {
    stop("ERROR: keep.vars and remove.vars cannot both be specified.")
  }
  if ((!is.null(p$keep.ids)) && (!is.null(p$remove.ids))) {
    stop("ERROR: keep.ids and remove.ids cannot both be specified.")
  }

  # Read in the data
  if (p$is.the.data) {
    x <- p$data
    p$data <- NULL
  } else if (p$file.type == 1) {
    x <- loadFile(file)
  } else {
    x <- read.table(file, header=p$header, sep=p$delimiter, 
          na.strings=p$in.miss, as.is=TRUE)
  }

  # See if there is an id variable
  idFlag <- 1 - is.null(p$id.var)

  # Check the variables
  nc     <- ncol(x)
  cnames <- colnames(x)
  temp   <- list(maxValue=nc, minValue=1, checkList=cnames)
  check.vec(p$keep.vars, "keep.vars", temp)
  check.vec(p$remove.vars, "remove.vars", temp)
  check.vec(p$factor.vars, "factor.vars", temp)
  temp$len <- 1

  if (idFlag) check.vec(p$id.var, "id.var", temp)

  # Get the rows to keep
  if ((!is.null(p$remove.ids)) || (!is.null(p$keep.ids))) {
  
    if (!is.null(p$remove.ids)) {
      temp.ids <- p$remove.ids
      flag     <- 1
    } else {
      temp.ids <- p$keep.ids
      flag     <- 0
    }

    temp.ids <- as.character(unique(temp.ids))
    nrows    <- nrow(x)
    if (idFlag) {
      id <- as.character(x[, p$id.var])
    } else {
      id <- as.character(1:nrows)
    }
  
    temp <- id %in% temp.ids
    temp[is.na(temp)] <- FALSE
    x <- removeOrKeepRows(x, temp, which=flag)
  }

  # Subset data by variables
  sub <- getListName(p, "subsetData")
  if (!is.null(sub)) {

    # Get var names and unfactor
    for (i in 1:length(sub)) {
      temp <- getListName(sub[[i]], "var")
      x[, temp] <- unfactor(x[, temp], fun=NULL)
    } 

    x <- subsetData.list(x, sub)
  }

  # Remove variables
  temp <- getListName(p, "remove.vars")
  if (!is.null(temp)) {
    temp <- unique(temp)
    x    <- removeOrKeepCols(x, temp, which=-1) 
  }

  # Keep variables
  temp <- getListName(p, "keep.vars")
  if (!is.null(temp)) {
    temp <- unique(temp)
    x    <- removeOrKeepCols(x, temp, which=1) 
  }

  # Un-factor all factors
  if (p$file.type == 1) x <- unfactor.all(x)

  # Factor variables
  if (!is.null(p$factor.vars)) {
    factorFlag <- 1
    p$factor.vars <- unique(p$factor.vars)
    for (var in p$factor.vars) x[, var] <- factor(x[, var])

    if ((p$make.dummy) && (is.numeric(p$factor.vars))) {
      # If factor.vars is numeric, then get the variable names
      factor.char <- cnames[p$factor.vars]
    } else {
      factor.char <- p$factor.vars
    }
  } else {
    factorFlag <- 0
  }

  # Remove missing values
  if (p$remove.miss) x <- removeMiss(x)

  # Create dummy variables
  if ((p$make.dummy) && (factorFlag)) {
    n <- length(factor.char)
    b <- rep("x_!!_@#^&#*", times=n)
    d <- rep(1, times=n)
    x <- createDummy(x, vars=factor.char, baseline=b, keep.factor=d)$data
  }

  x

} # END: readTable

# Function to get the number of columns in a file or data object
getNcols <- function(infile, p) {

  p <- default.list(p, c("delimiter", "return.cols"), list("\t", 0))

  # See if infile is a data frame or matrix
  temp <- dim(infile)
  if (!is.null(temp)) return(temp[2])

  p$open <- "r"
  fid    <- getFID(infile, p)
  temp   <- scan(file=fid, what="character", sep=p$delimiter, 
               nlines=1, quiet=TRUE)
  close(fid)

  if (p$return.cols) {
    return(temp)
  } else {
    return(length(temp))
  }

} # END: getNcols

# Function to read in a file and return certain columns
getColumns <- function(file.list, vars, temp.list=NULL, op=NULL) {

  # file.list   (See file.list.wordpad)
  # vars        Vector of variable names or column numbers
  #############################################################
  # op         List with names
  #   return.type   0, 1, 2 0=list, 1=matrix, 2=data frame

  # Returns a list with each element is a name form vars, 
  # or if vars is numeric the position in the list.

  file.list <- default.list(file.list, 
                            c("file", "file.type", "delimiter", "header"), 
                            list("ERROR", 3, "\t", 1), error=c(1, 0, 0, 0))
  op <- default.list(op, c("return.type"), list(0))

  type <- getListName(file.list, "file.type")
  file <- getListName(file.list, "file")
  vars <- unique(vars)

  if (type == 1) {
    data <- loadFile(file)
  } else if (type == 3) {
    file.list$include.row1 <- file.list$header
    file.list$what         <- "character"
    file.list$returnMatrix <- 1

    data <- scanFile(file, file.list)

  } else if (type == 4) {
    temp.list <- check.temp.list(temp.list)
    file.list$sas.list$delimiter <- "|"
    file.list$sas.list$temp.list <- temp.list
   
    # Define an id variable to prevent an error
    file.list$sas.list$id.var <- "id"
    data <- getSASData(file, file.list, transpose=0)

  } else if (type %in% c(6, 8)) {
    file.list$method       <- 2
    file.list$include.row1 <- file.list$header
    file.list$what         <- "character"
    file.list$returnMatrix <- 1

    data <- readFile.gzip3(file, file.list)
  }

  # Check for errors
  oplist <- list(checkList=colnames(data), minValue=1, maxValue=ncol(data))
  check.vec(vars, "vars", oplist)

  return.type <- op$return.type
  if (return.type == 0) {
    ret <- list()
    for (var in vars) ret[[var]] <- as.vector(data[, var])
  } else {
    ret <- removeOrKeepCols(data, vars, which=1)
    if (return.type == 1) {
      ret <- as.matrix(ret)
    } else {
      ret <- data.frame(ret)
    }
  }

  ret
  
} # END: getColumns

# Function for type 2 zip files
readFile.zip2 <- function(infile, p) {

  p <- default.list(p, c("zipFile", "start.row", "stop.row", "stream"),
       list("ERROR", 1, -1, 0), error=c(1, 0, 0, 0))

  fid <- unz(infile, p$zipFile, open="r")
  if (p$stream) {
    return(fid)
  } else {
    return(readFile.type2(fid, p))
  }

} # END: readFile.zip2

# Function for type 2 zip files
readFile.gz2 <- function(infile, p) {

  p <- default.list(p, c("start.row", "stop.row", "stream"),
       list(1, -1, 0), error=c(0, 0, 0))

  fid <- gzfile(infile, open="r")
  if (p$stream) {
    return(fid)
  } else {
    return(readFile.type2(fid, p))
  }

} # END: readFile.gz2

# Function for readFile.zip2 and readFile.gz2
readFile.type2 <- function(fid, p) {

  # Read row 1
  row1 <- scan(file=fid, what="character", nlines=1, sep="\n", quiet=TRUE)

  # Update start.row and stop.row 
  p$start.row    <- max(p$start.row - 1, 1)
  if (p$stop.row > 0) p$stop.row <- max(p$stop.row - 1, 1)
  p$what         <- "character"
  p$delimiter    <- "\n"
  p$include.row1 <- 0
  p$returnMatrix <- 0
  dat            <- scanFile(fid, p)
  close(fid)
  dat            <- c(row1, dat)
  dat

} # END: readFile.type2

# Function for opening a file
getFID <- function(infile, p) {

  p <- default.list(p, c("open"), list("r"))
  if (is.null(p[["file.type", exact=TRUE]])) p$file.type <- getFileType(infile)

  type <- p$file.type
  if (type %in% c(2, 3)) {
    fid <- file(infile, open=p$open)
  } else if (type %in% c(5, 6)) {
    fid <- unz(infile, p$zipFile, open=p$open)
  } else if (type %in% c(7, 8)) {
    fid <- gzfile(infile, open=p$open)
  }
  fid

} # END: getFID

# Function for type 3 zip or gz files
readFile.gzip3 <- function(infile, p) {

  p <- default.list(p, c("file.type", "method", "delimiter"),
       list("ERROR", 1, "\t"), error=c(1, 0, 0),
       checkList=list(c(6, 8), 1:2, NA))

  # Call scanFile or readTable
  # NOTE: read.table automatically closes the connection
  if (p$method == 2) {
    # Set up list for scanFile
    p$open         <- "r"
    p$returnMatrix <- 1
    p$ncol         <- getNcols(infile, p)
    fid            <- getFID(infile, p) 
    dat            <- scanFile(fid, p)
    close(fid)
  } else {
    p$open         <- ""
    fid            <- getFID(infile, p) 
    dat            <- readTable(fid, p)
  }
  
  dat

} # END: readFile.gzip3

# Function to get the number of rows in a data set
getNrows <- function(infile, file.type=3, delimiter="\n") {

  fid <- getFID(infile, list(file.type=file.type, open="r"))
  i   <- 0
  while (1) {
    temp <- scan(fid, what="character", sep="\n", nlines=1, 
                 n=1, quiet=TRUE)
    if (!length(temp)) break
    i <- i + 1
  }
  close(fid)
  return(i)

} # END: getNrows

# Function to check variables in a data set
checkVars <- function(flist, vars) {

  flist <- default.list(flist, c("file", "file.type", "delimiter", "header"), 
                        list("ERROR", 3, "\t", 1), error=c(1, 0, 0, 0))
  flist$return.cols <- 1
  f <- flist$file
  cols <- getNcols(f, flist) 
  if (is.numeric(vars)) {
    op <- list(minValue=1, maxValue=length(cols))
  } else {
    op <- list(checkList=cols)
  }
  
  ret <- check.vec(vars, paste("Variables for ", f, sep=""), op) 
  if (ret) stop("ERROR in checkVars")

  0

} # END: checkVars

# Function to return the file type
getFileType <- function(f, default=3) {

  ret  <- default
  vec  <- getVecFromStr(f, delimiter=".")
  n    <- length(vec)
  ext1 <- ""
  if (n) {
    ext0 <- vec[n]
    if (n > 1) ext1 <- vec[n-1]
    if (ext0 == "gz") {
      if (ext1 %in% c("txt", "xls", "csv", "dat", "data", "sdat", "def")) {
        ret <- 8
      } else if (ext1 %in% c("ldat")) {
        ret <- 7
      }
    } else if (ext0 %in% c("txt", "xls", "csv", "dat", "data", "sdat", "def")) {
      ret <- 3
    } else if (ext0 %in% c("rda")) {
      ret <- 1
    } else if (ext0 %in% c("ldat")) {
      ret <- 2
    } else if (ext0 %in% c("zip")) {
      if (ext1 %in% c("txt")) {
        ret <- 6
      } else if (ext1 %in% c("ldat")) {
        ret <- 5
      }
    } else if (ext0 %in% c("lbat", "sbat")) {
      ret <- ext0
    }
  }

  ret

} # END: getFileType

# Function for guessing what the delimiter in a file is
getFileDelim <- function(f, type=NULL, default="\t") {

  delim <- default
  if (is.null(type)) type <- getFileType(f)
  if (type %in% c("lbat", "sbat")) return(delim)
  fid <- getFID(f, list(file.type=type))
  x <- scan(fid, what="character", nlines=2, sep="\n")
  close(fid)
  if (length(x) < 2) return(delim)

  x1 <- getVecFromStr(x[1], delimiter="")
  x1 <- sort(table(x1, exclude=NULL), decreasing=TRUE)
  x2 <- getVecFromStr(x[2], delimiter="")
  x2 <- sort(table(x2, exclude=NULL), decreasing=TRUE)
  n1 <- names(x1)
  n2 <- names(x2)
  temp <- (n1 %in% n2) & (x1 %in% x2)
  if (!any(temp)) return("\n")

  x1    <- x1[temp]
  n1    <- n1[temp]
  check <- c("\t", " ", "\n", "|", ",")
  len   <- length(n1)
  delim <- n1[1]
  if (len > 1) {
    for (i in 2:len) {
      temp <- n1[i]
      if (temp %in% check) {
        delim <- temp
        break
      }
    }
  }
  if (!(delim %in% check)) delim <- "\n"
  
  delim 

} # END: getFileDelim

# Function to load a table
loadData.table <- function(flist) {

  if (!is.list(flist)) {
    if (!is.character(flist)) stop("ERROR in loadData.table: incorrect type of input argument")
    temp <- flist
    flist <- list(file=temp)
  }
  flist <- check.file.list(flist)
  flag  <- 0
  type  <- flist$file.type
  if (type == 1) {
    # Load the object
    temp <- load(flist$file)
    name <- flist[["name", exact=TRUE]]
    if (is.null(name)) name <- "data"
    if ((length(temp) > 1) || (name != "data")) {
    
      # Rename the object
      data <- eval(parse(text=name))
      eval(parse(text=paste("rm(", temp, ")", sep="")))
      temp <- gc(verbose=FALSE)
    }
  } else if (type == 3) {
    method <- flist[["method", exact=TRUE]]
    if (is.null(method)) method <- 1
    if (method == 2) {
      flist$returnMatrix <- 1
      data <- scanFile(flist$file, flist)
    } else {
      data <- readTable(flist$file, flist)
      flag <- 1
    }
  } else if (type == 4) {
    # SAS data set
    flist$sas.list$id.var <- flist$id.var
    data <- getSASData(flist$file, flist, transpose=0) 
  } else if (type %in% c(6, 8)) {
    # Type 3 compressed file
    method <- flist[["method", exact=TRUE]]
    if (is.null(method)) flist$method <- 2
    data <- readFile.gzip3(flist$file, flist)
  } else {
    stop("ERROR in loadData.table: file.type must be 1, 3, 4, 6, or 8")
  } 

  # Apply other arguments
  if (!flag) {
    flist$is.the.data <- 1
    flist$data <- data
    data <- readTable(flist$file, flist)
  }

  data

} # END: loadData.table



# History  Mar 28 2008 Add option to createDummy
#                      Add getOrder function
#          Mar 31 2008 Change snpNames option. snpNames are added
#                      in check.snp.list
#          Apr 07 2008 Add getNames function.
#                      Check for file seperator in check list functions.
#          Apr 11 2008 Add function removeOrKeepCols
#          Apr 15 2008 Change to sas.list
#          Apr 16 2008 Add removeOrKeepRows function
#          Apr 17 2008 Change in default.list
#          Apr 25 2008 Add unfactor function
#          Apr 30 2008 Add getVarNames function
#          May 01 2008 Add factorVars function
#                      Allow logical vector in removeOrKeepRows
#          May 05 2008 Change getSnpNames to getIdsFromFile
#          May 06 2008 Change snp.col, chrm.col, and loc.col to
#                      snp.var, chrm.var, and loc.var
#          May 09 2008 Add getInheritanceVec function
#                      Add addInterVars function
#          May 13 2008 Add option to createDummy to keep the original
#                      factor variable in the returned data frame.
#          May 31 2008 Add function matchNames
#          Jun 07 2008 Unfactor factors in changeStrata
#          Jun 10 2008 Change inheritance to genetic.model
#          Jun 11 2008 Add checkList option to default.list function
#          Jun 13 2008 Add getCommandArg function
#          Jun 14 2008 Remove nvars in getVarNames
#          Jun 20 2008 Add function to define lists for genfile
#                      Make snp.list$in.miss and snp.list$heter.codes
#                      to be character vectors in sheck.snp.list
#                      Check if input list is NULL in default.list
#                      Fix bug in recode.geno
#          Jun 21 2008 Change createDummy to allow for input matrices
#                      Return a list with added variable names
#          Jun 25 2008 Make createDummy more efficient
#                      Add start.n to changeStrata
#                      Make removeMiss more general
#                      Add sort2D function
#          Jun 27 2008 Change createDummy function
#          Jun 27 2008 Use getColumns function
#          Jun 27 2008 Add update.snp.list function
#                      Add getTempfile
#          Jun 28 2008 Add function extractByStr
#          Jun 30 2008 Add function closeFile
#          Jul 02 2008 Add code for type 5 files in check.snp.list
#                      Allow for GLU formats
#                      Add makeVector function
#          Jul 06 2008 Add initDataFrame function
#          Jul 08 2008 Add functions for subsetting data
#          Jul 11 2008 Add renameVar function
#          Jul 15 2008 Make default delimiter in snp.list a tab
#          Jul 17 2008 Change addInterVars
#                      Add getVarNames.int
#          Jul 24 2008 Add getSNP.type2
#          Jul 30 2008 Change to checkForSep, check.snp.list, 
#                       check.locusMap.list
#          Aug 07 2008 Add checkForConstantVar function
#          Aug 11 2008 Add unfactor.all function
#          Aug 15 2008 Redo getSNPdata
#          Aug 18 2008 Add strataMatrix
#          Aug 19 2008 Add rep.rows, rep.cols, rep.mat
#          Oct 10 2008 Add recode.data function
#          Oct 11 2008 Fix bug in subsetData.var for vector value
#          Oct 21 2008 Fix bug in addInterVars
#          Oct 22 2008 Extend recode.data function
#                      Add mergeData function
#          Oct 29 2008 Make changeLevels.var more general
#          Oct 31 2008 Allow for user defined genotypes in recode.data
#          Nov 04 2008 Extend removeOrKeepRows for character rows
#          Nov 05 2008 Extend subsetData.var to return rows
#          Nov 06 2008 Change the functionality subsetData.list
#                      Extend subsetData.var for character vector values
#          Nov 12 2008 Add callOS function
#          Nov 26 2008 Add GetColsFromCharVec function
#          Nov 26 2008 Generalize removeMiss function
#          Dec 08 2008 Add orderGenotypes function
#          Dec 10 2008 Keep rownames in removeOrKeepCols/removeOrKeepRows
#                      Add changeStr.names function
#          Dec 11 2008 Add applyFormulas function
#                      Add removeMiss.vars function
#          Jan 08 2009 Add changeAlleles function
#          Jan 21 2009 Add addColumn function
#          Jan 30 2009 Update changeAlleles function.
#                      Add writeTable function
#          Jan 31 2009 Update subsetData.var for missing values
#          Feb 02 2009 Add getPartition function
#          Feb 04 2009 Add crossTab function
#          Feb 05 2009 Change in recode.geno to allow for determining the major
#                      and minor alleles from a subset. 
#          Feb 07 2009 Add mergePhenoGeno function
#          Feb 23 2009 Fix bug in subsetData.var with NAs
#          Feb 26 2009 Extend crossTab function
#          Mar 05 2009 Make addColumn more efficient
#          Mar 16 2009 Update crossTab for 1 variable
#                      Add parse.vec function
#          Mar 19 2009 Add error checks in subsetData.list
#          Mar 20 2009 Change to subsetData.list to return the rows
#          Mar 23 2009 Update sort2D for 1 NA row
#          Mar 24 2009 Add replaceStr.var
#          Mar 25 2009 Add warnings to addColumn function
#          Mar 27 2009 Add option returnRows to subsetData.list
#          Apr 01 2009 Add option for NAs in subsetData
#          Apr 07 2009 Fix bug in recode.geno
#          Apr 08 2009 Let file be a data frame in addColumn
#          Apr 14 2009 Generalize subsetData.list
#          May 07 2009 Add option to mergePhenoGeno to include original genotypes
#          May 07 2009 Change replaceStr.var function
#                      Combine addColumn and addCol2
#                      Remove addCol2 function
#          May 12 2009 Add function orderVars 
#          May 13 2009 Add option to addColumn to check for leading zeros
#                      in the id variables.
#          May 29 2009 Update removeOrKeepCols to print out cols not found
#          Jun 01 2009 Update addColumn
#          Jun 05 2009 Add getTopHits function
#          Jun 05 2009 Add option to addColumn for not replacing variables.
#                      Add function mergeTopHits
#          Jun 08 2009 Add getRank.file function
#          Jun 09 2009 Change in subsetData.var for a numeric variable
#          Jun 15 2009 Add option to extractByStr for efficiency
#          Jun 25 2009 Force input argument data to be a data frame in addColumn
#                      Create function addRow
#          Jul 02 2009 Return a matrix in addColumn if the input object is a matrix
#          Jul 02 2009 Add option to sort2D
#          Jul 10 2009 Fix bug in getVarNames with NULL input object
#          Jul 21 2009 Add option to removeOrKeepCols to not throw error
#                      when columns are not found
#                      Move out older code to wga_util2.R
#                      Add function to check the delimiter in a file
#          Jul 22 2009 Add flag to snp.list, pheno.list, etc to check if
#                      the check.xxx.list was called.
#                      Add name mergeDirFile to snp.list
#          Jul 24 2009 Add writeVec function
#          Aug 03 2009 Changes in check.pheno.list
#          Aug 04 2009 Add removeWhiteSpace function
#                      Move less used code to wga_util2.R
#          Aug 10 2009 Check snpNames.list in check.snp.list
#          Aug 10 2009 Add options to snp.list in check.snp.list
#                      Add option in extractByStr to keep or remove
#                      matched values
#          Aug 18 2009 Add debug.time function
#          Aug 19 2009 Update check.vec to return flag
#          Aug 27 2009 Add exclude option in crossTab
#          Sep 24 2009 Add function to compute all interactions
#                      between 2 sets of variables in a data frame
#          Oct 01 2009 Change var to "var" in checkForConstantVar function
#          Oct 16 2009 Set up changes for efficiency:
#                      check.snp.list
#                      check.pheno.list
#          Oct 16 2009 Return major/minor alleles in recode.geno
#          Oct 16 2009 Change in closeFile to check that the file exists
#                      before deleting it.
#          Oct 20 2009 Add fun option to unfactor.all
#          Oct 23 2009 Add function checkTryError
#          Oct 27 2009 Change else statement in getInheritanceVec
#          Oct 28 2009 Add option removeVars from checkForConstantVar
#          Oct 28 2009 Fix bug in checkForConstantVar
#          Nov 03 2009 Add column names in getColsFromCharVec
#          Nov 10 2009 Update genfile.list for functions
#          Nov 10 2009 Update crossTab so that snp.list can be NULL
#          Nov 12 2009 Fix bug in recode.geno when out.genotypes is NULL
#          Nov 12 2009 Check vars in addColumn
#          Nov 20 2009 NO CHANGE
#          Dec 10 2009 Fix bug in determining the major/minor alleles in recode.geno
#          Dec 17 2009 Add option to getVarNames.int for seperating 
#                       interaction terms      
#          Dec 22 2009 Initialize a1, a2 to NULL in recode.geno       
#          Dec 29 2009 Use getFileDelim in check.xxx.list 
#          Dec 30 2009 Add code for getting variables from formulas    
#                      Add miss option to removeMiss.vars      
#                      Set default value for in.miss in check.pheno.list 
#          Dec 31 2009 Add option to callOS
#          Jan 04 2010 Update applyFormulas
#          Jan 13 2010 Add function check.file.list
#          Jan 15 2010 Add check.subsetData.list, update check.file.list.
#                      Add getSubsetDataVars, check.subsetData.list
#          Jan 19 2010 Update check.pheno.list with check.file.list
#          Feb 01 2010 Set extended = FALSE in replaceStr.var
#          Mar 04 2010 Add normVarNames function
#          Mar 05 2010 Update getIdsFromFile to call loadData.table
#          Mar 14 2010 Check vector length in recode.geno
#          Mar 18 2010 Fix bug in genfile.list with an empty list

# Function to pull apart each string from the data object
getVecFromStr <- function(string, delimiter="|") {

  # string       String to break apart. No default
  # delimiter    Delimiter used in string. The default is "|". 

  strsplit(string, delimiter, fixed=TRUE)[[1]]

} # END: getVecFromStr

# Function to recode a vector of genotypes and missing values
recode.geno <- function(vec, in.miss=c("  "), out.miss=NA,
               out.genotypes=c(0,1,2), heter.codes=NULL, subset=NULL) {

  # vec           Input vector 
  #               No default.
  # in.miss       Vector of categories which denote missing values. 
  #               These categories will be changed to out.miss,
  #               if out.miss is not NULL.
  #               The default is "  "  (2 spaces)
  # out.miss      New category for missing values. Use NULL if you want to
  #               keep the original missing categories.
  #               The default is NA
  # out.genotypes Vector of length 3 containing the new genotypes,
  #               the first element is for the major homozygous
  #               genotype and the second element is for the 
  #               heterozygous genotype. 
  #               Use NULL for no recoding of the genotypes
  # heter.codes   Vector of codes used for the heterozygous genotype.
  #               If NULL, then it is assumed that the heterozygous 
  #               genotype is of the form "AB", "Aa", "CT", etc (a 2-character
  #               string with different characters... case sensitive!!!)
  #               The default is NULL.
  # subset        NULL or a logical vector of length length(vec) to be used
  #               in determining the major and minor homozygous genotypes.
  #               The default is NULL so that all (non-missing) elements
  #               of vec will be used. 

  # Get the rows that contain missing values
  # Do not change these values now, because we need to know the other codes.
  tempMiss <- vec %in% in.miss
  if (!is.null(out.miss)) {
    mFlag    <- 1
  } else {
    mFlag    <- 0
  }

  # Determine if a subset is to be used
  subFlag <- !is.null(subset)
  if (subFlag) {
    if (!is.logical(subset)) stop("ERROR in recode.geno: subset is not a logical vector")
    if (length(subset) != length(vec)) stop("ERROR in recode.geno: length(subset) != length(vec)")
    if (!any(subset)) subFlag <- 0
    vec0 <- vec
  }

  # Initialize
  index   <- Inf
  subSum  <- Inf
  alleles <- "  "
  hflag   <- 1
  a1      <- NULL
  a2      <- NULL

  if (!is.null(out.genotypes)) {

    # Get the frequency counts
    if (subFlag) {
      tab <- sort(table(vec[subset], exclude=in.miss), decreasing=TRUE)
    } else {
      tab <- sort(table(vec, exclude=in.miss), decreasing=TRUE)
    }

    # Remove the missing values
    if (!is.null(in.miss)) {
      genos <- names(tab)
      tab   <- tab[!(genos %in% in.miss)]
    }

    # Check for error
    if (length(tab) > 3) {
      print(tab)
      stop("ERROR: in recode.geno")
    }

    # Get the genotypes
    genos   <- names(tab)
    flag    <- 0
    flag2   <- 0
    minFlag <- 0
    hflag   <- !is.null(heter.codes)  

    # Get the ids for each genotype
    ids    <- list()
    index  <- 1
    if (subFlag) subSum <- sum(tempMiss)
    for (geno in genos) {
      ids[[index]] <- vec == geno
      if (subFlag) subSum <- subSum + sum(ids[[index]])
      index        <- index + 1   
    } 

    # tab is sorted in descending order
    index <- 1
    for (geno in genos) {

      # Determine if this genotype is heterozygous
      if (hflag) {
        heterFlag <- geno %in% heter.codes
      } else {
        a1        <- substr(geno, 1, 1)
        a2        <- substr(geno, 2, 2)
        heterFlag <- a1 != a2
      }

      if (!heterFlag) {
        # Homozygous, but check if major homozygous has been assigned
        if (!flag) {
          vec[ids[[index]]] <- out.genotypes[1]
          flag <- 1
          hom1 <- a1
        } else {
          # Minor homozygous
          vec[ids[[index]]] <- out.genotypes[3]
          minFlag <- 1
          hom2    <- a1
        }
      } else {
        # Heterozygous
        # Check for error 
        if (flag2) {
          print(tab)
          stop("ERROR: in recode.geno, flag2=1")
        }
        vec[ids[[index]]] <- out.genotypes[2]
        flag2 <- 1
        hetg  <- geno
      }
      index <- index + 1

    } # END: for (geno in genos)

  } # END: if (!is.null(out.genotypes))

  # Change missing values
  if (mFlag) vec[tempMiss] <- out.miss

  # Get the major/minor alleles
  if (!hflag) {
    if (flag) {
      # Major
      if (minFlag) {
        # There was major and minor genotypes
        alleles <- paste(hom1, hom2, sep="")
      } else if (flag2) {
        # Use heterozygous genotype
        a2 <- substr(hetg, 2, 2)
        if (hom1 != a2) {
          alleles <- paste(hom1, a2, sep="")
        } else {
          a1 <- substr(hetg, 1, 1)
          alleles <- paste(hom1, a1, sep="")
        }
      } else {
        # Only major homozygous
        alleles <- paste(hom1, hom1, sep="") 
      }
    } else {
      # Only heterozygous genotype
      if (flag2) alleles <- hetg
    }
  }

  # Final error check if a subset was used
  if ((subFlag) && (index < 4) && (subSum < length(vec))) {
    temp <- as.numeric(!flag) + as.numeric(!flag2) + as.numeric(!minFlag)
    if (temp > 1) {
      # Use all elements
      temp <- recode.geno(vec0, in.miss=in.miss, out.miss=out.miss,
               out.genotypes=out.genotypes, heter.codes=heter.codes, subset=NULL)
      return(temp)
    }

    # Get the ids that were not mapped
    temp <- tempMiss
    for (i in 1:length(ids)) temp <- temp | ids[[i]]
      
    if (!minFlag) vec[!temp] <- out.genotypes[3]
    if (!flag2)   vec[!temp] <- out.genotypes[2]
    if (!flag)    vec[!temp] <- out.genotypes[1]
  }

  list(vec=vec, alleles=alleles)

} # END: recode.geno

# Function to assign a default value to an element in a list
default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} # END: default.list

# Function to get a name from a list (without partial matching)
getListName <- function(inList, name) {

  if (name %in% names(inList)) {
    return(inList[[name]])
  } else {
    return(NULL)
  }

} # END: getListName

# Function to check snp.list
check.snp.list <- function(snp.list) {

  if (is.null(snp.list)) stop("ERROR: snp.list must be specified")

  # Check the names in the list
  snp.list <- default.list(snp.list, 
            c("file", "in.miss", 
        "read.n", "genetic.model", "stream", "recode",
        "alreadyChecked", "out.miss", "out.delimiter", "snpNames.keep"), 
        list("ERROR", "  ", -1, 0, 0, 1, 0, NA, "\t", 1), 
            error=c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
         checkList=list(NA, NA, NA, 0:3, 0:1, 0:1, NA, NA, NA, NA))

  if (snp.list$alreadyChecked == 1) return(snp.list)

  # Check file vector
  nfiles <- length(snp.list$file)
  snp.list$dir <- checkForSep(snp.list[["dir", exact=TRUE]])
  temp <- paste(snp.list$dir, snp.list$file, sep="")
  if (check.files(temp)) stop("ERROR: in check.snp.list") 
  snp.list$file <- temp
  snp.list$dir <- ""

  if (is.null(snp.list[["file.type", exact=TRUE]])) {
    snp.list$file.type <- getFileType(snp.list$file[1], default=7)
  }
  if (is.null(snp.list[["delimiter", exact=TRUE]])) {
    snp.list$delimiter <- getFileDelim(snp.list$file[1], type=snp.list$file.type, default="\t")
  }

  # Check start.vec and stop.vec
  snp.list <- default.list(snp.list, c("start.vec", "stop.vec"), 
               list(rep(1, nfiles), rep(-1, nfiles)), error=c(0, 0))
  temp <- snp.list[["snpNames", exact=TRUE]]
  if (is.null(temp) && is.null(snp.list$snpNames.list)) {
    a1 <- check.vec.num(snp.list$start.vec, "snp.list$start.vec", len=nfiles)
    a2 <- check.vec.num(snp.list$stop.vec, "snp.list$stop.vec", len=nfiles)
    if (sum(a1+a2)) stop() 
  } else {
    snp.list$start.vec <- rep(1, times=nfiles)
    snp.list$stop.vec  <- rep(-1, times=nfiles)
  }
  temp <- (snp.list$stop.vec == 1)
  if (any(temp)) snp.list$stop.vec[temp] <- 2
  
  # Check for id variable for type 3 and 4
  if (snp.list$file.type %in% c(3, 4)) {
    if (is.null(snp.list$id.var)) stop("snp.list$id.var is not specified")
  }

  # Check for sas.list
  if (snp.list$file.type == 4) {
    if (is.null(snp.list$sas.list)) stop("snp.list$sas.list is not specified")
  }

  # Check for zip file with type 5
  if (snp.list$file.type == 5) {
    if (is.null(snp.list$zipFile)) stop("snp.list$zipFile is not specified")
  }

  # Check in.miss
  vec  <- as.character(snp.list$in.miss)
  temp <- is.na(vec)
  if (any(temp)) vec[temp] <- "NA"
  snp.list$in.miss <- vec

  # Check heter.codes
  vec <- snp.list$heter.codes
  if (!is.null(vec)) {
    vec <- as.character(vec)
    temp <- is.na(vec)
    if (any(temp)) vec[temp] <- "NA"
    snp.list$heter.codes <- vec
  }

  # Check file.type
  temp <- snp.list$file.type
  if (is.numeric(temp)) {
    if (!(temp %in% 1:8)) {
      temp <- paste("ERROR:", temp, "is not a valid value for snp.list$file.type")
      stop(temp) 
    }
  } else {
    # GLU format
    if (snp.list$stream == 0) {
      warning("Assuming snp.list$file.type is a GLU format. Setting snp.list$stream to 1")
      snp.list$stream <- 1
    }
  }

  # Check stream and nfiles
  if ((snp.list$stream) && (nfiles > 1)) {
    snp.list$stream <- 0
    stop("ERROR: snp.list$stream = 1 and length(snp.list$file) > 1")
  } 

  # For file type 1, stream must be 0
  if (snp.list$file.type == 1) {
    print("snp.list$stream is set to 0 for file.type = 1")
    snp.list$stream <- 0
  }
   
  # Check snpNames.list
  temp <- snp.list[["snpNames.list", exact=TRUE]]
  if (!is.null(temp)) {
    temp <- default.list(temp, c("file", "file.type", "delimiter", "header", "method"),
                         list("ERROR", 3, "\n", 0, 2), error=c(1, 0, 0, 0, 0))
    snp.list$snpNames.list <- temp
  }

  snp.list$alreadyChecked <- 1

  snp.list

} # END: check.snp.list

# Function to check pheno.list
check.pheno.list <- function(pheno.list) {

  if (is.null(pheno.list)) stop("ERROR: pheno.list must be specified")

  pheno.list <- default.list(pheno.list, 
   c("file", "id.var", "header", "remove.miss", "alreadyChecked", 
     "is.the.data", "in.miss"),
        list("ERROR", "ERROR", 1, 0, 0, 0, c(NA, "NA", NaN, "NaN", ".")), 
               error=c(1, 1, 0, 0, 0, 0, 0))

  if (pheno.list$alreadyChecked == 1) return(pheno.list)

  temp  <- getAllVars(pheno.list)
  pheno.list <- check.file.list(pheno.list, op=list(exist=1, vars=temp))

  if (!(pheno.list$file.type %in% c(1, 3, 4, 6, 8)))
    stop("ERROR: pheno.list$file.type must be 1, 3, 4, 6 or 8")

  # Check for sas.list
  if (pheno.list$file.type %in% c(4)) {
    if (is.null(pheno.list$sas.list)) stop("pheno.list$sas.list is not specified")
  }

  # Check for zip file with type 5
  if (pheno.list$file.type == 5) {
    if (is.null(pheno.list$zipFile)) stop("pheno.list$zipFile is not specified")
  }

  # Check id names list
  temp <- pheno.list[["keep.ids.list", exact=TRUE]]
  if (!is.null(temp)) {
    pheno.list$keep.ids <- getIdsFromFile(temp, 
                           id.vec=getListName(pheno.list, "keep.ids"))
    pheno.list$keep.ids.list <- NULL
  }

  # Check id names list
  temp <- pheno.list[["remove.ids.list", exact=TRUE]]
  if (!is.null(temp)) {
    pheno.list$remove.ids <- getIdsFromFile(temp, 
                           id.vec=getListName(pheno.list, "remove.ids"))
    pheno.list$remove.ids.list <- NULL
  }

  temp <- pheno.list[["keep.vars", exact=TRUE]]
  if (!is.null(temp)) {
    temp2 <- pheno.list[["factor.vars", exact=TRUE]]
    if (!is.null(temp2)) {
      if (!all(temp2 %in% temp)) stop("ERROR in check.pheno.list: all factor.vars not in keep.vars")
    }
  }

  pheno.list$alreadyChecked <- 1

  pheno.list

} # END: check.pheno.list

# Function to check locusMap.list
check.locusMap.list <- function(locusMap.list) {

  if (is.null(locusMap.list)) stop("ERROR: locusMap.list must be specified")

  locusMap.list <- default.list(locusMap.list, 
     c("file", "header", "snp.var","chrm.var", "loc.var", 
       "alreadyChecked"),
              list("ERROR", 0, "ERROR", "ERROR", "ERROR", 0), 
              error=c(1, 0, 1, 1, 1, 0))

  if (locusMap.list$alreadyChecked == 1) return(locusMap.list)

  locusMap.list$dir <- checkForSep(getListName(locusMap.list, "dir"))
  temp <- paste(locusMap.list$dir, locusMap.list$file, sep="")
  if (check.files(temp)) stop("ERROR: in check.locusMap.list") 

  if (is.null(locusMap.list[["file.type", exact=TRUE]])) {
    locusMap.list$file.type <- getFileType(temp[1], default=3)
  }
  if (is.null(locusMap.list[["delimiter", exact=TRUE]])) {
    locusMap.list$delimiter <- getFileDelim(temp[1], type=locusMap.list$file.type, default="")
  }


  if (!(locusMap.list$file.type %in% c(1, 3, 4, 6, 8)))
    stop("ERROR: locusMap.list$file.type must be 1, 3, 4, 6 or 8")

  if (!locusMap.list$header) {
    a1 <- check.vec.num(locusMap.list$snp.var, "locusMap.list$snp.var", len=1)
    a2 <- check.vec.num(locusMap.list$chrm.var, "locusMap.list$chrm.var", len=1)
    a3 <- check.vec.num(locusMap.list$loc.var, "locusMap.list$loc.var", len=1)
    if (sum(a1+a2+a3)) stop() 
  } else {
    locusMap.list$header <- 1
  }

  # Check for sas.list
  if (locusMap.list$file.type %in% c(4)) {
    if (is.null(locusMap.list$sas.list)) stop("locusMap.list$sas.list is not specified")
  }

  locusMap.list$alreadyChecked <- 1

  locusMap.list

} # END: check.locusMap.list

# Function to check temp.list
check.temp.list <- function(temp.list) {

  if (is.null(temp.list)) temp.list <- list()

  temp.list <- default.list(temp.list, c("dir", "delete", "id", "alreadyChecked"),
                               list(" ", 1, 1, 0))

  if (temp.list$alreadyChecked == 1) return(temp.list)

  # Check the directory
  temp.list$dir <- checkForSep(temp.list$dir) 
  if (temp.list$dir != "") {
    warn <- getOption("warn")
    # Change the warn option so a warning turns to an error
    options(warn=2)
    temp <- try(dir(temp.list$dir), silent=TRUE)
    if (class(temp) == "try-error") {
      temp <- paste("ERROR: the directory ", temp.list$dir, " does not exist", sep="")
      options(warn=warn)
      stop(temp)
    }
    options(warn=warn)
  }

  temp.list$alreadyChecked <- 1

  temp.list

} # END: check.temp.list

# Function to check if a file exists
check.files <- function(files) {
  # files   Character vector of files
  
  temp <- !(file.exists(files))
  if (any(temp)) {
    ret <- 1
    files <- files[temp]
    for (file in files) print(paste("The file ", file, " does not exist", sep=""))
  } else {
    ret <- 0
  }
  
  ret

} # END: check.files

# Function to check a numeric vector
check.vec.num <- function(vec, name, len=NULL, maxValue=NULL, 
                 minValue=NULL) {

  if (is.null(vec)) return(0)
  ret <- 0
  if (!is.numeric(vec)) {
    print(paste("ERROR: ", name, " must be of type numeric", sep=""))
    ret <- 1
  }

  if (!is.null(len)) {
    if (length(vec) != len) {
      print(paste("ERROR: ", name, " must be of length ", len, sep=""))
      ret <- 1
    }
  }

  if (!is.null(maxValue)) {
    if (max(vec) > maxValue) {
      print(paste("ERROR: each element of ", name, " must <= ", maxValue, sep=""))
      ret <- 1
    }
  }

  if (!is.null(minValue)) {
    if (min(vec) < minValue) {
      print(paste("ERROR: each element of ", name, " must >= ", minValue, sep=""))
      ret <- 1
    }
  }

  ret

} # END: check.vec.num

# Function to check a character vector
check.vec.char <- function(vec, name, len=NULL, checkList=NULL) {

  if (is.null(vec)) return(0)
  ret <- 0
  if (!is.character(vec)) {
    print(paste("ERROR: ", name, " must be of type character", sep=""))
    ret <- 1
  }

  if (!is.null(len)) {
    if (length(vec) != len) {
      print(paste("ERROR: ", name, " must be of length ", len, sep=""))
      ret <- 1
    }
  }

  if (!is.null(checkList)) {
    temp <- is.na(match(vec, checkList))
    if (sum(temp)) {
      ret  <- 1
      temp <- vec[temp]
      for (x in temp) {
        print(paste("ERROR: ", name, " contains the invalid element ", x, sep=""))
      }
    }
  }

  ret

} # END: check.vec.char

# Function to check for character vectors is a list
check.list.vec.char <- function(inList, names=NULL) {

  # inList    List to check
  # names     Character vector of names to check. If NULL, then
  #           all names are checked
  #           The default is NULL
  
  if (is.null(names)) names <- names(inList)
  sum <- 0
  for (name in names) {
    if (!is.null(inList[[name]])) {
      sum <- sum + check.vec.char(inList[[name]], name)
    }
  }
  if (sum) stop()
  0

} # END check.list.vec.char

# Function to check a vector
check.vec <- function(vec, name, opList) {

  opList <- default.list(opList, c("stopOnError"), list(1))

  if (is.numeric(vec)) {
    ret <- check.vec.num(vec, name, len=opList$len, 
            maxValue=opList$maxValue, minValue=opList$minValue)
  } else {
    ret <- check.vec.char(vec, name, len=opList$len, 
            checkList=opList$checkList)
  }
  if (ret) {
    if (opList$stopOnError) stop()
  }
  ret

} # END: check.vec

# Function to get the (approximate) number of rows to read
getRead.n <- function(file, what="character") {

  fid  <- file(file, "r")
  temp <- scan(file=fid, what="character", sep="\n", nlines=1, quiet=TRUE)
  p1   <- seek(fid, where=0, origin="start")  
  temp <- scan(file=fid, what="character", sep="\n", nlines=2, quiet=TRUE)
  p2   <- seek(fid, where=0, origin="start")
  p3   <- seek(fid, where=0, origin="end")
  p3   <- seek(fid, where=0, origin="end")
  close(fid)
  nrow <- p3/(p2-p1)
 
  

} # END: getRead.n

# Function to create dummy variables. The returned data object will
# have the original factor variables removed from it (see keep.factor)
# and new dummy variables added. For example if data[, "x1"] is a factor 
# with levels 0, 1, and 2, then the new variables will be "x1_1" and
# "x1_2", provided 0 is the baseline category.
# If you do not want to keep the baseline variable too, then set
#  baseline to a value that is not a level of the factor.
createDummy <- function(data, vars=NULL, baseline=NULL, 
                         keep.factor=NULL) {
  # data           Data frame, matrix, or vector
  # vars           Character vector of variable names or numeric vector
  #                of column numbers for  the factor variables that will be
  #                turned into dummy variables. If data is a vector, then
  #                vars does not have to be specified. 
  #                If NULL, then all the 
  #                factors will turn into dummy variables.
  #                The default is NULL.
  # baseline       Vector of baseline categories. The length of
  #                this vector must be equal to the length of vars or 
  #                the number of factors in data.
  #                The default is NULL.
  # keep.factor    Logical vector to keep the original factors in the 
  #                returned object. If NULL, then all factors are
  #                removed. The length of this vector must equal the
  #                length of vars or the number of factors. 
  #                The default is NULL.

  # Check for vector and column names
  if (is.null(dim(data))) {
    dim(data) <- c(length(data), 1)
    vars      <- "VAR1"
  }
  cnames <- colnames(data)
  nc     <- ncol(data)
  if (is.null(cnames)) {
    cnames <- paste("VAR", 1:nc, sep="")
    colnames(data) <- cnames
  }

  # Check for numeric vector vars
  if (is.numeric(vars)) vars <- cnames[vars]

  # Get the variables to factor
  if (is.null(vars)) {
    for (i in nc) {
      if (is.factor(data[, i])) vars <- c(vars, cnames[i])
    }
  } else {
    # Check the variable names
    temp <- check.vec.char(vars, "vars", len=NULL, checkList=colnames(data))
    if (temp) stop("ERROR in vars")
  }
  nvar <- length(vars)
  if (!nvar) return(list(data=data))

  if (!is.null(baseline)) {
    if (length(baseline) != nvar) stop("ERROR: baseline is incorrect")
    baseFlag <- 1
  } else {
    baseFlag <- 0
  }
  if (!is.null(keep.factor)) {
    if (length(keep.factor) != nvar) stop("ERROR: keep.factor is incorrect")
  } else {
    keep.factor <- rep(FALSE, times=nvar)
  }

  # Initialize a list
  newVars   <- list()
  ret       <- data
  nr        <- nrow(data)
  addedCols <- NULL

  for (i in 1:nvar) {
    var <- vars[i]

    # Get the levels
    temp <- data[, var]
    if (!is.factor(temp)) temp <- factor(temp)
    levels  <- levels(temp)
    if (length(levels) == 1) {
      keep.factor[i] <- TRUE
      next
    }

    # Get the baseline category
    if (!baseFlag) {
      base <- levels[1]
    } else {
      base <- baseline[i]
    }
     
    # Remove the baseline
    temp    <- !(levels == base) 
    levels  <- levels[temp]
    nlevels <- length(levels)
 
    # Initialize a dummy matrix
    temp           <- matrix(data=NA, nrow=nr, ncol=nlevels)
    addedVars      <- paste(var, "_", levels, sep="")
    colnames(temp) <- addedVars

    # Add the new variables names to the return list
    newVars[[var]] <- addedVars

    # Get the binary values
    for (j in 1:nlevels) {
      temp[, j] <- as.numeric(data[, var] == levels[j])
    }

    # Add the new variables
    addedCols <- cbind(addedCols, temp)  
  }

  # Add the new variables
  ret <- cbind(ret, addedCols)

  # Free memory
  rm(data, addedCols)
  temp <- gc(verbose=FALSE)

  # Remove original factors
  if (sum(keep.factor) != nvar) {
    temp <- vars[!keep.factor]
    ret  <- removeOrKeepCols(ret, temp, which=-1)
  }

  list(data=ret, newVars=newVars)

} # END: createDummy

# Function to match elements of 2 vectors
getOrder <- function(baseVec, newVec, removeMiss=0, errorIfMiss=1) {

  # baseVec        Baseline vector
  # newVec
  # removeMiss
  # errorIfMiss

  ret <- match(baseVec, newVec)
  if (errorIfMiss) {
    if (any(is.na(ret))) stop("ERROR: in matching vectors")
  }
  if (removeMiss) ret <- ret[!is.na(ret)]

  ret

} # END: getOrder

# Function to read in files and save as an R object file
dat2rda <- function(infile, outfile) {
                   
  dat <- readLines(infile)
  save(dat, file=outfile)

  0
} # END: dat2rda

# Function to read in file (numeric matrix) and save as an R object file
matrix2rda <- function(infile, outfile, delimiter="\t", stopOnError=1) {
                   
  tlist <- list(returnMatrix=1, include.row1=1, what=double(0),
                delimiter=delimiter, start.row=1, stop.row=-1)

  options(warn=2)
  dat <- try(scanFile(infile, tlist), silent=TRUE)
  options(warn=1)
  if (class(dat) == "try-error") {
    temp <- paste("ERROR with file ", infile, sep="")
    if (stopOnError) stop(temp)
    return(0)
  }
  save(dat, file=outfile)

  0
} # END: matrix2rda

# Function to get all ids requested
getIdsFromFile <- function(file.list, id.vec=NULL) {

  if (!is.null(file.list)) {
    file.list <- check.file.list(file.list)
    file.list <- default.list(file.list, c("id.var"), list(1))

    if (file.list$id.var == -1) {
      fid  <- getFID(file.list$file, file.list)
      temp <- scan(file=file.list$file, what="character", 
                   sep=file.list$delimiter)
      close(fid)
      id.vec <- c(id.vec, temp) 
    } else {
      var    <- file.list$id.var
      temp   <- loadData.table(file.list)
      id.vec <- c(id.vec, makeVector(temp[, var]))
    }
  }
  id.vec <- unique(id.vec)
  if (!length(id.vec)) stop(paste("No ids in file", file.list$file))

  id.vec

} # END: getIdsFromFile

# Function to return the names of a vector or array
getNames <- function(obj) {

  if (is.null(dim(obj))) {
    ret <- names(obj)
  } else {
    ret <- dimnames(obj)
    ret <- ret[[1]]
  }
  ret

} # END: getNames

# Function to return or create variable names
getVarNames <- function(obj, prefix="VAR") {

  if (is.null(obj)) return(NULL)
  if (is.data.frame(obj)) return(colnames(obj))
  if (is.matrix(obj)) {
    ret <- colnames(obj)
    if (is.null(ret)) ret <- paste(prefix, 1:ncol(obj), sep="")
    return(ret)
  }
  ret <- names(obj)
  if (is.null(ret)) ret <- paste(prefix, 1:length(obj), sep="")
  ret

} # END: getVarNames

# Function to check if a directory ends with a slash.
# If not, then it will be added.
checkForSep <- function(dir) {

  if (is.null(dir)) return("")

  if (dir %in% c("", " ", "  ")) {
    dir <- paste(".", .Platform$file.sep, sep="")
    return("")
  }

  n <- nchar(dir)
  if (substring(dir, n, n) %in% c("\\", "/")) {
    return(dir)
  } else {
    dir <- paste(dir, .Platform$file.sep, sep="")
    return(dir)
  }

} # END: checkForSep

# Function to remove columns or variables from a matrix or data frame
removeOrKeepCols <- function(x, vars, which=1, stopOnError=1) {

  # x      Matrix or data frame 
  # vars   Vector of variable names or column numbers.
  # which  1 or -1 to keep or remove cols
  #        The default is 1 
  # stopOnError  0 or 1 to call stop() if columns are not found

  n      <- dim(x)
  dfFlag <- is.data.frame(x)
  if (is.null(n)) stop("ERROR: x should be 2 dimensional")
  if (which != 1) which <- -1

  if (!is.numeric(vars)) {
    cols <- match(vars, colnames(x))
    temp <- is.na(cols)
    if (any(temp)) {
      if (stopOnError) {
        print(vars[temp])
        print("The above columns were not found in the data")
        stop("ERROR in removeOrKeepCols")
      }
      cols <- cols[!temp]
      if (!length(cols)) {
        if (which == 1) {
          return(NULL)
        } else {
          return(x)
        }
      }
    }
  } else {
    cols <- vars
  }
  cols <- which*cols

  # Keep or remove columns
  ret <- x[, cols]
  
  # Check for NULL dimension
  if (is.null(dim(ret))) {
    dim(ret) <- c(n[1], length(ret)/n[1])
    if (dfFlag) ret <- data.frame(ret)
    cnames <- colnames(x)
    if (!is.null(cnames)) colnames(ret) <- cnames[cols]
    rownames(ret) <- rownames(x)
  }

  ret

} # END: removeOrKeepCols

# Function to remove columns or variables from a matrix or data frame
removeOrKeepRows <- function(x, rows, which=1) {

  # x      Matrix or data frame 
  # rows   Vector of row numbers, character vector of row names,
  #        or logical vector.
  # which  1 or -1 to keep or remove rows
  #        The default is 1 

  n      <- dim(x)
  dfFlag <- is.data.frame(x)
  if (is.null(n)) stop("ERROR: x should be 2 dimensional")
  if (which != 1) which <- -1
  rnames <- rownames(x)

  if (is.logical(rows)) {
    if (length(rows) != n[1]) stop("ERROR with logical vector rows")
    if (which != 1) rows <- !rows
  } else if (is.character(rows)) {
    rows <- match(rows, rnames)
    if (any(is.na(rows))) stop("ERROR in removeOrKeepRows")
    rows <- which*rows
  } else {
    rows <- which*rows
  }
  cnames <- colnames(x)

  # Keep or remove rows
  x <- x[rows, ]
  
  # Check for NULL dimension
  if (is.null(dim(x))) {
    dim(x) <- c(length(x)/n[2], n[2])
    if (dfFlag) x <- data.frame(x)
    if (!is.null(cnames)) colnames(x) <- cnames
    rownames(x) <- rnames[rows]
  }

  x

} # END: removeOrKeepRows

# Function to write snp rows to an open file
writeSnpLines <- function(snpNames, fid, snpData, nsub, orderFlag=0,
                    order=NULL, delimiter="|", sep="|") {

  # snpNames       Character vector of snps (no header)
  # fid            File connection
  # snpData        Character vector of the snp data
  # nsub           Number of subjects
  # orderFlag      0 or 1 if the subjects are already ordered
  #                0 = not ordered
  # order          Vector of integers for orderFlag = 1
  # delimiter      Input delimiter
  # sep            Output delimiter

  # Define a local function to search for character strings
  f1 <- function(str) {

    grep(str, snpData, value=FALSE)

  } # END: f1

  # Get the rows by searching for the snp names
  temp <- unlist(lapply(snpNames, f1))

  if (length(temp)) {
    snpData <- snpData[temp]
    nsnps   <- length(snpData)
  } else {
    return(0)    
  }

  nsubP1 <- nsub + 1
  ret    <- 0
  for (i in 1:nsnps) {
    temp  <- getVecFromStr(snpData[i], delimiter=delimiter)
    snp   <- temp[1]
    temp  <- temp[-1]
    if (!orderFlag) temp  <- temp[order]
    if (snp %in% snpNames) {
      write(c(snp, temp), file=fid, ncolumns=nsubP1, sep=sep)
      ret <- ret + 1
    }
  }

  ret

} # END: writeSnpLines

# Function to create factors in a data frame
factorVars <- function(data, vars) {

  # data    Data frame
  # vars    Vector of variables names or column numbers

  if (!is.data.frame(data)) stop("ERROR: data must be a data frame")
  if (is.numeric(vars)) vars <- colnames(data)[vars]

  for (var in vars) {
    data[, var] <- factor(data[, var])
  }

  data

} # END: factorVars

# Function to unfactor all columns ij a matrix or data frame
unfactor.all <- function(data, fun=NULL) {

  nc <- ncol(data)
  if (is.null(nc)) {
    nc <- 1
    dim(data) <- c(length(data), 1)
  }
  for (i in 1:nc) data[, i] <- unfactor(data[, i], fun=fun)
  data

} # END: unfactor.all

# Function to un-factor a factor
unfactor <- function(fac, fun=NULL) {

  # fac   Factor
  # fun   Function like as.character or as.numeric, etc

  if (is.factor(fac)) {
    ret <- levels(fac)[fac]
  } else {
    ret <- fac
  }

  if (!is.null(fun)) ret <- fun(ret)

  ret

} # END: unfactor

# Function to multiply each row or column of a matrix by a vector
matrixMultVec <- function(mat, vec, by=2) {

  # by    1 or 2  1 = rows, 2 = columns
  
  d <- dim(mat)
  if (by == 1) {
    vec <- rep(vec, each=d[1])
  } else {
    vec <- rep(vec, times=d[2])
  }

  dim(vec) <- d
  ret <- mat*vec
  ret

} # END: matrixMultVec

# Function to divide each row or column of a matrix by a vector
matrixDivideVec <- function(mat, vec, by=2) {

  # by    1 or 2  1 = rows, 2 = columns

  d <- dim(mat)
  if (by == 1) {
    vec <- rep(vec, each=d[1])
  } else {
    vec <- rep(vec, times=d[2])
  }

  dim(vec) <- d

  ret <- mat/vec
  ret

} # END: matrixDivideVec

# Function to get the column number of a 1 for each row in a matrix
#  of dummy variables. The matrix must contain only one 1 in each row.
# If the matrix has a 1 in the (i,j)th element, then in the returned
# vector ret, ret[i] = j.logistic.dsgnMat
getColNumber <- function(mat) {

  nc  <- ncol(mat)
  nr  <- nrow(mat)
  ret <- rep(NA, times=nr)
  for (i in 1:nc) {
    temp      <- mat[, i] == 1
    ret[temp] <- i
  }
  ret

} # END: getColNumber

# Function to return a block diagonal matrix
blockDiag <- function(mat.list) {
  
  # mat.list   List of matrices to form the block diagonal matrix

  # Get the dimensions of each matrix
  d <- matrix(unlist(lapply(mat.list, dim)), byrow=TRUE, ncol=2)

  # Initialize 
  ret  <- matrix(data=0, nrow=sum(d[,1]), ncol=sum(d[,2]))
  row0 <- 1
  row1 <- 0
  col0 <- 1
  col1 <- 0

  # Set each block
  for (i in 1:length(mat.list)) {
    row1 <- row1 + d[i,1]
    col1 <- col1 + d[i,2]
    ret[row0:row1, col0:col1] <- mat.list[[i]]
    row0 <- row1 + 1
    col0 <- col1 + 1
  }

  ret

} # END: blockDiag

# Function to add intercept column to a matrix or vector
addIntercept <- function(x, nrow=NULL) {

  if (is.null(x)) {
    if (is.null(nrow)) stop("Cannot add intercept")
    ret <- matrix(1, nrow=nrow, ncol=1)
  } else {
    if (is.null(nrow)) nrow <- nrow(x)
    if (is.null(nrow)) nrow <- length(x)
    nc <- ncol(x)
    if (is.null(nc)) {
      nc     <- 1
      dim(x) <- c(nrow, nc)
    } 
    ret <- cbind(rep(1, times=nrow), x)  
  }

  ret

} # END: addIntercept

# Function to use integers 1, ..., n for a categorical vector
changeStrata <- function(vec, start.n=1) {

  if (is.factor(vec)) vec <- unfactor(vec)
  uvec <- sort(unique(vec))
  ret  <- rep.int(0, times=length(vec))
  for (u in uvec) {
    ret[vec == u] <- start.n
    start.n <- start.n + 1
  }

  ret

} # END: changeStrata

# Function to remove rows that contain at least 1 missing value from
#  a data frame or matrix or vector
removeMiss <- function(x, miss=NA) {
  
  # x
  # miss    Vector of missing values

  d <- dim(x)

  # For a vector
  if (is.null(d)) {
    flag   <- 1
    cnames <- names(x)
    cFlag  <- !is.null(cnames)
    d      <- c(length(x), 1)
    dim(x) <- d
  } else {
    flag  <- 0
    cFlag <- 0
  }

  # Be careful for a data frame
  if (is.data.frame(x)) {
    temp <- matrix(data=TRUE, nrow=d[1], ncol=d[2])
    for (i in 1:d[2]) temp[, i] <- !(x[, i] %in% miss)
  } else {
    temp <- !(x %in% miss)
  }
  if (cFlag) cnames <- cnames[temp]

  dim(temp) <- d
  temp <- rowSums(temp)
  temp <- temp == d[2]
  x    <- removeOrKeepRows(x, temp, which=1)
  
  if (flag) {
    dim(x) <- NULL
    if (cFlag) names(x) <- cnames
  }
  
  x

} # END: removeMiss

# Function to remove missing values of a data frame from
#  certain variables
removeMiss.vars <- function(x, vars=NULL, miss=NA) {

  if (is.null(vars)) {
    x <- removeMiss(x, miss=miss)
    return(x)
  }

  temp <- rep(TRUE, times=nrow(x)) 
  for (var in vars) temp <- temp & !(x[, var] %in% miss)
  
  x <- removeOrKeepRows(x, temp, which=1)
  x

} # END: removeMiss.vars

# Function to write a (named) vector to a file. For type 3 and close = 0,
#  the file conection is returned.
writeVecToFile <- function(vec, file, colnames=NULL, type=3, close=1,
                           sep=" ") {

  # vec
  # file
  # colnames
  # type        1 or 3  
  # close       0 or 1 (for type 3 only)

  if (type == 1) {
    if (!is.null(colnames)) names(vec) <- colnames
    save(vec, file=file)
    return(0)
  } else {
    fid <- file(file, "w")
    n   <- length(vec)
    if (!is.null(colnames)) write(colnames, file=fid, ncolumns=n, sep=sep)
    write(vec, file=fid, ncolumns=n, sep=sep)
    if (close) {
      close(fid)
      fid <- 0
    }
    return(fid)
  }
  0

} # END: writeVecToFile 

# Function for writing vectors to files
writeVec <- function(vec, fileOrFID, colnames=NULL, isFID=0, sep="\t", 
                     close=1, type=3) {

  if (type == 1) {
    if (!is.null(colnames)) names(vec) <- colnames
    save(vec, file=fileOrFID)
    return(0)
  } 

  if (!isFID) fileOrFID <- file(fileOrFID, "w")
  if (!is.null(colnames)) {
    temp <- paste(colnames, collapse=sep, sep="")
    write(temp, file=fileOrFID, ncolumns=1)
  }
  temp <- paste(vec, collapse=sep, sep="")
  write(temp, file=fileOrFID, ncolumns=1)
  if (close) {
    close(fileOrFID)
    fileOrFID <- 0
  }

  fileOrFID

} # END: writeVec

# Function to return the vector of genotypes for the different modes
#  of inheritance
getInheritanceVec <- function(which, recode=1) {

  # which   NULL, 0, 1, 2, 3
  #         0 = trend
  #         1 = dominant
  #         2 = recessive
  #         3 = factor
  
  if ((!is.null(recode)) && (!recode)) return(NULL)

  if ((is.null(which)) || (!which)) {
    return(c(0, 1, 2))
  } else if (which == 1) {
    return(c(0, 1, 1))
  } else if (which == 2) {
    return(c(0, 0, 1))
  } else {
    return(c(0, 1, 2))
  }

} # END: getInheritanceVec

# Function to add interaction variables to a matrix or data frame
# Returns a list with the names data and newVars. newVars is a 
# character vector containing the variables added.
addInterVars <- function(data, vec, inter.vars, prefix="SNP") {

  # data           Data frame or matrix
  # vec            Numeric vector or factor for interactions
  #                If a matrix with more than 1 column, then
  #                it is assumed the columns are dummy variables
  # inter.vars     Matrix or data frame of variables that 
  #                will interact with var. These variables cannot
  #                be factors.
  # prefix         Variable prefix for interaction variables added
  #                The default is "SNP".

  cnames2 <- colnames(inter.vars) 
  if (is.null(cnames2)) cnames2 <- paste("V", 1:ncol(inter.vars), sep="")
 
  facFlag <- is.factor(vec)
  ncVec   <- ncol(vec)
  if (is.null(ncVec)) ncVec <- 0
  mFlag   <- ncVec > 1

  if (facFlag || mFlag) {
    newVars <- NULL
    if (facFlag) {
      vec <- data.frame(vec)
      colnames(vec) <- prefix
      vec <- createDummy(vec)$data
    }
    cnames <- colnames(vec)
    #temp2  <- NULL

    for (i in 1:ncol(vec)) {
      temp   <- matrixMultVec(inter.vars, vec[, i], by=2)
      tnames <- paste(cnames[i], "_", cnames2, sep="")
      newVars <- c(newVars, tnames) 
      colnames(temp) <- tnames
      if (i == 1) {
        temp2 <- temp
      } else {
        temp2 <- cbind(temp2, temp)
      }
    }
    data <- cbind(data, temp2)
  } else {
    temp <- matrixMultVec(inter.vars, vec, by=2)
    newVars <- paste(prefix, "_", cnames2, sep="")
    colnames(temp) <- newVars
    data <- cbind(data, temp)
  }
 
  list(data=data, newVars=newVars)

} # END: addInterVars

# Function to return variable names and positions from a vector
matchNames <- function(vec, name, exact=0) {

  # vec       Character vector to search from.
  # name      Character vector of names to search for.
  # exact     0 or 1 to perform exact matching
  #           The default is 0.

  names <- NULL
  pos   <- NULL

  if (exact) {
    temp <- match(name, vec)
    temp <- temp[!is.na(temp)]
    if (length(temp)) {
      names <- vec[temp]
      pos   <- temp
    } 
  } else {
    for (n in name) {
      temp <- grep(n, vec)
      if (length(temp)) {
        pos   <- c(pos, temp)
        names <- c(names, vec[temp])
      }
    }
    temp  <- !duplicated(pos)
    pos   <- pos[temp]
    names <- names[temp]
  }

  list(names=names, pos=pos)  

} # END: matchNames

# Function to get an option from the command arguments
getCommandArg <- function(optName, fun=NULL) {
  
  # optName      Unique option name in the command arguements
  # fun          Function to apply

  # Get the vector of command arguments
  comArgs <- commandArgs()

  # Get the one that we need
  temp <- grep(optName, comArgs)

  if (length(temp) == 1) { 
    ret <- sub(optName, "", comArgs[temp])
    if (!is.null(fun)) ret <- fun(ret)
  } else {
    temp <- paste("ERROR with ", optName, " in commandArgs()", sep="")
    stop(temp)
  }

  ret

} # END: getCommandArg

# Function to define list in swarm generator files
genfile.list <- function(inList, listName, fid) {

  # If the field is inList is a function or family,
  # put it in quotes and set the comment to "FUNCTION".

  names <- names(inList)
  llen  <- length(inList)
  if (is.null(names)) {
    flag  <- 1   
    names <- 1:llen
    str1  <- '[['
    str2  <- ']] <- '
  } else {
    flag <- 0
    str1 <- '[["'
    str2 <- '"]] <- '
  }
  temp  <- paste("\n ", listName, " <- list() \n", sep="")
  cat(temp, file=fid)

  if (!llen) return(NULL)
  for (name in names) {
    temp <- inList[[name]]
    cmm  <- comment(temp)
    if (is.null(cmm)) cmm <- "NULL"

    # For a list
    if (is.list(temp)) {
      if (flag) {
        name2 <- paste("list", name, sep="")
      } else {
        name2 <- name
      }
      genfile.list(temp, name2, fid)
      temp <- paste(listName, str1, name, str2, name2, ' \n', sep='')
      cat(temp, file=fid)
      next
    }
    if (length(temp) == 1) {
      if (is.character(temp)) { 
        if (cmm == "FUNCTION") {
          temp <- paste(listName, str1, name, str2, temp, ' \n', sep='')
        } else {
          temp <- paste(listName, str1, name, str2, '"', temp, '" \n', sep='')
        }
      } else {
        temp <- paste(listName, str1, name, str2, temp, ' \n', sep='')
      }
      cat(temp, file=fid) 
    } else {
      if (flag) {
        field <- paste(listName, "$", "'", name, "'", sep="")
      } else {
        field <- paste(listName, "$", name, sep="")
      }
      genfile.vec(temp, field, fid)
    }  
  }

  NULL

} # END: genfile.list

# Function to define a vector in swarm generator files
genfile.vec <- function(vec, name, fid) {

  n       <- length(vec)
  vecType <- "character"
  cflag   <- 1
  if (is.numeric(vec)) {
    vecType <- "numeric"
    cflag   <- 0
  }
  temp <- paste(name, " <- ", vecType, "(", n, ") \n", sep="")
  cat(temp, file=fid)
  k <- 1
  for (temp in vec) { 
    if (cflag) {
      temp <- paste(name, '[', k, '] <- "', temp, '" \n', sep='')
    } else {
      temp <- paste(name, '[', k, '] <- ', temp, ' \n', sep='')
    }
    cat(temp, file=fid)
    k <- k + 1
  }
  
  NULL

} # END: genfile.vec

# Function to sort a matrix or data frame by a column
sort2D <- function(data, col, dec=FALSE, fun=NULL) {

  # data   Data frame or matrix
  # col    Column number or variable name to sort on
  # dec    0 or 1 for decreasing
  #        The default is 0
  # fun    Function to apply to col before sorting
  #        The default is NULL

  vec  <- makeVector(unfactor(data[, col]))
  if (!is.null(fun)) vec <- fun(vec)
  temp <- is.na(vec)
  keep <- removeOrKeepRows(data, temp, which=1)
  data <- removeOrKeepRows(data, temp, which=-1)
  temp <- sort(vec, decreasing=dec, index.return=TRUE)$ix
  data <- rbind(data[temp, ], keep)
  return(data)

} # END: sort2D

# Function to update snp.list
update.snp.list <- function(snp.list, where=0) {

  # where    Integer specifying where in the program

  if (where == 1) {
    stream <- getListName(snp.list, "stream")
    type   <- getListName(snp.list, "file.type")
    if (stream) {
      if (type != 2) {
        n <- length(getListName(snp.list, "file"))
        snp.list$start.vec <- rep(1, times=n)
        snp.list$stop.vec  <- rep(-1, times=n)
      }
    }
  }

  temp <- getListName(snp.list, "snpNames.list")
  if (!is.null(temp)) {
    snp.list$snpNames <- getIdsFromFile(temp, 
                           id.vec=getListName(snp.list, "snpNames"))
    snp.list$snpNames <- unique(snp.list$snpNames)
    snp.list$snpNames.list <- NULL
  }

  temp <- getListName(snp.list, "snpNames")
  if (!is.null(temp)) {
    #n <- length(snp.list$file)
    #snp.list$start.vec <- rep(1, times=n)
    #snp.list$stop.vec  <- rep(-1, times=n)
  }

  snp.list

} # END: update.snp.list

# Function to return a temporary file name
getTempfile <- function(dir, prefix=NULL, ext=NULL) {

  pattern <- paste(prefix, collapse="", sep="")
  ret     <- tempfile(pattern=pattern, tmpdir="")
  ret     <- paste(dir, ret, ext, sep="") 

  ret

} # END: getTempfile

# Function to extract elements from a character vector
extractByStr <- function(dat, search, op=NULL) {

  # dat           Character vector to search in
  # search        Character vector for strings to search for in dat
  ################################################################
  # op            List with names:
  #  include.row1 0 or 1
  #               The default is 1
  #  substr.vec   Vector of max length 2 for start and stop options
  #               in the substring function
  #               The default is NULL.
  #  keep         0 or 1 to remove or keep matched strings
  #               The default is 1

  # Define a local function to search for character strings
  f1 <- function(str) {

    grep(str, dat, value=FALSE)

  } # END: f1

  op <- default.list(op, c("include.row1", "keep"), list(1, 1))

  sub.vec <- getListName(op, "substr.vec")
  if (!is.null(sub.vec)) {
    subFlag <- 1
    a       <- sub.vec[1]
    b       <- sub.vec[2]
    save    <- dat
    dat     <- substr(dat, a, b)
    if (is.na(b)) b <- Inf
  } else {
    subFlag <- 0
  }

  # Get the rows by searching for the snp names
  rows <- unlist(lapply(search, f1))
  # Get a logical vector of rows to keep or drop
  temp <- (1:length(dat)) %in% rows
  keep <- op$keep
  if (keep == 0) temp <- !temp 

  # Add the first row if needed
  if (op$include.row1) {
    temp[1] <- TRUE
  } else {
    temp[1] <- FALSE
  }
 
  # Get the subset
  if (subFlag) {
    return(save[temp])
  } else {
    return(dat[temp])
  }

} # END: extractByStr

# Function to close and delete a file
closeFile <- function(fid, file=NULL, delete=0) {

  close(fid)
  ret <- 0
  if ((delete) && (!is.null(file))) {
    if (file.exists(file)) ret <- file.remove(file)
  }
  ret

} # END: closeFile

# Function to create a vector from a matrix
makeVector <- function(x) {

  d <- dim(x)
  if (is.null(d)) return(x)

  nn <- NULL
  if (d[1] == 1) {
    nn <- colnames(x)
  } else if (d[2] == 1) {
    nn <- rownames(x)
  }
  dim(x) <- NULL
  if (!is.null(nn)) names(x) <- nn
  if ((!is.vector(x)) && (is.list(x))) x <- unlist(x)
  x
 
} # END: makeVector

# Function to initialize a data frame
initDataFrame <- function(nrow, columns, rownames=NULL, colnames=NULL,
                  initChar="", initNum=NA) {
  
  cc   <- character(nrow)
  nn   <- double(nrow)
  nn[] <- initNum
  cc[] <- initChar
  cc   <- data.frame(cc)
  nn   <- data.frame(nn)
  nc   <- length(columns)

  columns <- toupper(columns)
  if (columns[1] == "C") {
    data <- cc
  } else {
    data <- nn
  }
  
  if (nc > 1) { 
    for (x in columns[-1]) {
      if (x == "C") {
        data <- cbind(data, cc)
      } else {
        data <- cbind(data, nn)
      }
    }
  }  
  data <- data.frame(data)
  if (!is.null(rownames)) rownames(data) <- rownames
  if (!is.null(colnames)) colnames(data) <- colnames

  for (i in 1:nc) {
    if (is.factor(data[, i])) data[, i] <- unfactor(data[, i])
  }
  data

} # END: initDataFrame

# Function to subset data by a single variable
subsetData.var <- function(data, var, operator, value, which=1, 
                           returnRows=0, na.value=FALSE) {

  # data
  # var
  # operator
  # value
  # which        1 or -1 to keep or drop
  #              The default is 1
  # returnRows   0 or 1
  #              The default is 0
  # na.value     TRUE or FALSE on what to do with NAs
  #              The default is FALSE

  lenv <- length(value)
  if (any(is.na(value))) {
    if (lenv > 1) stop("ERROR: IN subsetData.var: with NA in value")
    naFlag  <- 1
    numFlag <- 0
  } else {
    numFlag <- is.numeric(value)
    naFlag  <- 0
  }

  # Be careful with a vector for value
  if (lenv > 1) {
    if (numFlag) {
      value <- paste(value, collapse=",", sep="")
      value <- paste("c(", value, ")", sep="")
    } else {
      # Change the operator
      if (operator == "%in%") {
        operator <- "=="
      } else { 
        stop("ERROR: value cannot be a character vector")
      }
    }
  }

  if (naFlag) {
    temp <- ""
    if (operator == "!=") temp <- "!"  
    callStr <- paste("(", temp, "is.na(data[, var]))", sep=" ") 
  }
  else if (numFlag) {
    callStr <- paste("(as.numeric(data[, var])", operator, value, ")", sep=" ") 
  } else {
    callStr <- ""
    for (i in 1:lenv) {
      temp <- paste('(data[, var] ', operator, ' "', value[i], '")', sep='')
      if (i < lenv) temp <- paste(temp, " | ", sep="")
      callStr <- paste(callStr, temp, sep="")
    }
  }

  rows <- eval(parse(text=callStr))
  rows[is.na(rows)] <- na.value
  if (returnRows) return(rows)

  data <- removeOrKeepRows(data, rows, which=which)
  data

} # END: subsetData.var

# Function to subset data. Returns the subsetted data, or if op$which = 0,
#   returns the vector of TRUE/FALSE to subset
subsetData.list <- function(data, slist, returnRows=0) {

  # data
  ###################################################################
  # slist       List of sublists with names:
  #  var
  #  operator
  #  value
  #  which      1 or -1 to include or NOT within each list
  #             The default is 1
  #  logic.op   Logical operator 
  #             Only used starting from the second sublist
  #             The default is "&"
  #  na.value   TRUE or FALSE
  #             The default is FALSE
  #  last.which 1 or -1 If -1, the rows defined by all the sublists will
  #             be negated. The value for last.which in the last sublist
  #             is the only one applied.
  #             The default is 1.
  ####################################################################
  # returnRows  Set to 1 to return the logical vector of rows instead of
  #             the data.
  #             The default is 0

  n      <- length(slist)
  wvec   <- rep(9999, times=n)
  cnames <- colnames(data)
  cflag  <- !is.null(cnames)

  for (i in 1:n) {

    tlist    <- slist[[i]]
    if (!is.list(tlist)) stop("ERROR in subsetData.list: Input list is incorrect")
    tlist    <- default.list(tlist, 
                  c("var", "operator", "value", "which", "logic.op", "na.value", "last.which"), 
                  list("ERROR", "ERROR", "ERROR", 1, "&", FALSE, 1), error=c(1,1,1,0,0,0,0))
    var      <- getListName(tlist, "var")
    if ((cflag) && (is.character(var))) {
      if (!(var %in% cnames)) {
        temp <- paste("ERROR in subsetData.list: ", var, " not in data", sep="")
        print(temp)
        stop()
      }
    }
    operator <- getListName(tlist, "operator")
    value    <- getListName(tlist, "value")
    which    <- getListName(tlist, "which")
    last     <- getListName(tlist, "last.which")
    wvec[i]  <- which
    na.value <- getListName(tlist, "na.value")
    temp     <- subsetData.var(data, var, operator, value, returnRows=1, na.value=na.value)
    if (which == -1) temp <- !temp
    if (i == 1) {
      rows <- temp
    } else {
      logic.op <- getListName(tlist, "logic.op")
      callStr  <- paste("rows ", logic.op, " temp", sep="") 
      rows     <- eval(parse(text=callStr))      
    }
  }

  if (last == -1) rows <- !rows
  if (returnRows) return(rows)

  data <- removeOrKeepRows(data, rows, which=1)

  data

} # END: subsetData.list

# Function to rename a variable on a matrix or data frame
renameVar <- function(data, old.var, new.var) {

  cnames <- colnames(data)
  i      <- match(old.var, cnames)
  if (is.na(i)) return(data)
  cnames[i] <- new.var
  colnames(data) <- cnames
  data

} # renameVar

# Function to return a list of variable names for a SNP
getVarNames.snp <- function(prefix="SNP_", genetic.model=0) {

  # genetic.model = 3 is for factors, assuming 0 is the baseline
  if (genetic.model != 3) {
    return(prefix)
  } else {
    return(paste(prefix, 1:2, sep=""))
  }

} # END: getVarNames.snp

# Function to return a list of variable names for interactions
getVarNames.int <- function(V, prefix="SNP_", genetic.model=0, sep="_") {

  # V   Matrix of interactions

  vnames <- getVarNames(V, prefix="V")

  # genetic.model = 3 is for factors, assuming 0 is the baseline
  if (genetic.model != 3) return(paste(prefix, sep, vnames, sep="")) 
  
  ret <- c(paste(prefix, 1, sep, vnames, sep=""),
           paste(prefix, 2, sep, vnames, sep=""))

  ret

} # END: getVarNames.int

# Function to check for a constant variable in a matrix or data frame
# Returns a list of the data with the constant variables removed, and
#  a vector of variable names removed
checkForConstantVar <- function(data, msg=1, removeVars=1) {

  # data
  # msg          0, 1, 2  0 = no message, 1 = warning, 2 = error
  # removeVars   0 or 1 to remove the constant variables from the
  #              data.
  #              The default is 1

  temp   <- apply(data, 2, "var", na.rm=TRUE)
  remove <- NULL
  temp2  <- temp == 0
  temp2[is.na(temp2)] <- FALSE
  ret    <- (1:length(temp))[temp2]

  if (length(ret)) {
    temp <- colnames(data)
    if (!is.null(temp)) {
      temp <- temp[ret]
    } else {
      temp <- ret
    }
    remove <- temp
    if (removeVars) data <- removeOrKeepCols(data, remove, which=-1)    

    # Message
    temp <- paste(temp, collapse=",")
    if (msg == 1) {
      temp <- paste("WARNING: Variables ", temp, " are constant and have been removed from the data", sep="") 
      warning(temp)
    } else if (msg == 2) {
      temp <- paste("ERROR: Variables ", temp, " are constant", sep="") 
      stop(temp)
    }
  } 

  list(data=data, remove=remove) 

} # END: checkForConstantVar

# Function to return a logical matrix from a stratification variable
strataMatrix <- function(strata) {

  # strata    Stratification vector

  strata <- unfactor(strata)
  nr     <- length(strata)
  us     <- unique(strata)
  nc     <- length(us)
  ret    <- matrix(data=FALSE, nrow=nr, ncol=nc)
  colnames(ret) <- us
  
  for (i in 1:nc) ret[, i] <- (strata == us[i])
  
  ret

} # END: strataMat

rep.cols <- function(vec, times) {
  # Function to replicate columns
  len      <- length(vec)
  dim(vec) <- c(len, 1)
  mat      <- rep(vec, times)
  dim(mat) <- c(len, times)
  mat
} # END: rep.cols

rep.mat <- function(scalar, nrow=NULL, ncol=NULL)
{
  mat      <- rep(scalar, each=nrow, times=ncol)
  dim(mat) <- c(nrow, ncol)
  mat
}

rep.rows <- function(vec, times) {
  # Function to replicate rows
  len      <- length(vec)
  dim(vec) <- c(1, len)
  mat      <- rep(vec, each=times)
  dim(mat) <- c(times, len)
  mat
} # END: rep.rows

# Function to merge 2 data frames 
mergeData <- function(base.data, new.data, new.vars, 
                      base.id="id", new.id="id") {

  # base.data   Data frame
  # new.data    Data frame
  # new.vars
  # base.id
  # new.id

  new.vars <- unique(new.vars)
  base.id  <- unique(base.id)
  new.id   <- unique(new.id)
  n.base   <- length(base.id)
  n.new    <- length(new.id)
  flag     <- 0

  if (n.base != n.new) stop("Number of id variables do not match")
  #if (nrow(new.data) != length(unique(new.data[, new.id]))) {
  #  stop("Ids not unique in new.data")
  #}
  if (nrow(base.data) != length(unique(base.data[, base.id]))) {
    flag <- 1
  }
  if (n.base > 1) flag <- 1

  temp <- match(new.id, new.vars)
  if (!is.na(temp)) new.vars <- new.vars[-temp]
  temp <- new.vars %in% colnames(base.data)
  if (any(temp)) {
     new.vars2 <- c(new.vars[!temp], paste(new.vars[temp], "_2", sep=""))
  } else {
     new.vars2 <- new.vars
  }

  # Add variables
  for (var in new.vars2) base.data[, var] <- NA
  nvars <- length(new.vars2)

  if (!flag) {
    rows <- base.data[, base.id] %in% new.data[, new.id]
    temp <- match(new.data[, new.id], base.data[, base.id])
    temp <- temp[!is.na(temp)]
    for (i in 1:nvars) {
      base.data[rows, new.vars2[i]] <- new.data[temp, new.vars[i]]
    } 
  } else {
    n <- nrow(base.data)
    for (i in 1:nrow(new.data)) {
      temp <- rep.int(TRUE, times=n)
      for (j in 1:n.base) {
        temp <- temp & (base.data[, base.id[j]] == new.data[i, new.id[j]])
      }
      if (any(temp)) base.data[temp, new.vars2] <- new.data[i, new.vars]
    }
  }

  base.data

} # END: mergeData

# Function to call the operating system
callOS <- function(command, intern=FALSE) {

  # Determine the platform
  os      <- .Platform$OS.type
  winFlag <- (os == "windows")

  if (winFlag) {
    ret <- shell(command, intern=intern)
  } else {
    ret <- system(command, intern=intern)
  }
  ret

} # END: callOS

# Function to get columns from a character vector
getColsFromCharVec <- function(vec, cols, delimiter="\t", colNames=NULL,
                        fun=as.numeric) {

  # vec
  # cols        Numeric or character vector of columns 
  # colNames    NULL or ordered vector of column names in vec
  #             if cols is character
  #             The default is NULL

  cnamesFlag <- !is.null(colNames)

  if (is.character(cols)) {
    cols <- match(cols, colNames)
  }

  n <- length(vec)
  ret <- matrix(data=NA, nrow=n, ncol=length(cols))
  if (cnamesFlag) colnames(ret) <- colNames[cols] 
  for (i in 1:n) {
    temp <- getVecFromStr(vec[i], delimiter=delimiter)
    ret[i, ] <- fun(temp[cols])
  }

  ret

} # END: getColsFromCharVec 

# Function to search and replace strings within the names of an oject
changeStr.names <- function(obj, search, replace="") {

  # obj
  # search
  # replace

  d <- dim(obj)
  if (is.null(d)) {
    n <- names(obj)
    if (!is.null(n)) names(obj) <- gsub(search, replace, n, fixed=TRUE)
  } else {
    n <- colnames(obj)
    if (!is.null(n)) colnames(obj) <- gsub(search, replace, n, fixed=TRUE)
    n <- rownames(obj)
    if (!is.null(n)) rownames(obj) <- gsub(search, replace, n, fixed=TRUE)
  }
  obj

} # END: changeStr.names

# Function to get the correct phenotype data when there are formulas
#  to be applied. The original phenotype data will be returned with
applyFormulas <- function(data, formulas, remove=NULL) {

  # data       Data frame
  # formulas   List of formulas with variables in the data frame

  flag   <- 0
  rflag  <- !is.null(remove)
  data2  <- data
  ids    <- 1:nrow(data2)
  rownames(data2) <- ids 
  for (i in 1:length(formulas)) {
    f <- formulas[[i]]
    if ("formula" %in% class(f)) {
      # Get the design matrix
      temp <- model.matrix(f, data=data2)
      ids  <- intersect(ids, rownames(temp))
      flag <- 1
      if (rflag) {
        temp <- removeMiss(data.frame(temp), miss=remove)
        ids  <- intersect(ids, rownames(temp))
      }
    }
  } 

  if (flag) data <- removeOrKeepRows(data, as.numeric(ids), which=1)
  data

} # END: applyFormulas

# Function to return the formulas from a list
getFormulas <- function(inlist) {

  # inlist    List

  ret   <- list()
  index <- 1
  for (i in 1:length(inlist)) {
    if ("formula" %in% class(inlist[[i]])) {
      ret[[index]] <- inlist[[i]]
      index        <- index + 1
    }
  }
  ret

} # END: getFormulas

# Function to return the variables from a particular object
getAllVars <- function(obj, names=c("response.var", "main.vars", "int.vars",
                                    "strata.var", "start.var", "stop.var",
                                    "partition.var", "group.var", "id.var",
                                    "keep.vars", "remove.vars", "factor.vars")) {

  if (is.null(obj)) return(NULL)
  clss <- class(obj)

  # Character vector
  if ((is.vector(obj)) && ("character" %in% clss)) return(unique(obj))

  # Formula
  if ("formula" %in% clss) return(unique(all.vars(obj)))

  # List
  if ("list" %in% clss) {
    ret <- NULL
    if (is.null(names)) names <- 1:length(obj)
    for (nn in names) {
      obj.n <- obj[[nn, exact=TRUE]]
      if (is.null(obj.n)) next

      clss <- class(obj.n)      
      if ((is.vector(obj)) && ("character" %in% clss)) {
        ret <- c(ret, obj.n)
      } else if ("formula" %in% clss) {
        ret <- c(ret, all.vars(obj.n))
      }
    }
    slist <- obj[["subsetData", exact=TRUE]]
    if (!is.null(slist)) ret <- c(ret, getSubsetDataVars(slist))

    return(unique(ret))
  }

  stop("ERROR in formulaVars: obj is of wrong type")

} # END: getAllVars

# Function to change the alleles in a data frame
changeAlleles <- function(data, alleleVars, rows=NULL, newVars=NULL) {

  # data
  # alleleVars
  # rows         Logical vector of length = nrow(data) or NULL
  #              The default is NULL, so that all rows will be changed
  # newVars      New variable names. If NULL, then allelVars will be used
  #              The order of newVars must match alleleVars.
  #              The default is NULL.

  if (is.null(rows)) rows <- rep(TRUE, times=nrow(data))
  if (is.null(newVars)) newVars <- alleleVars
  i <- 1
  for (var in alleleVars) {
    tempA <- (data[, var] == "A") & rows
    tempT <- (data[, var] == "T") & rows
    tempG <- (data[, var] == "G") & rows
    tempC <- (data[, var] == "C") & rows
    nv    <- newVars[i]
    data[tempA, nv] <- "T"
    data[tempT, nv] <- "A"
    data[tempG, nv] <- "C"
    data[tempC, nv] <- "G"
    i <- i + 1
  }
  data

} # END: changeAlleles

# Function to write a table
writeTable <- function(x, outfile) {

  write.table(x, file=outfile, sep="\t", row.names=FALSE, quote=FALSE)

} # END: writeTable

# Function to compute cross tab for ids data in the analysis
crossTab <- function(snp.list, pheno.list, varlist, op=NULL) {

  # snp.list
  # pheno.list
  # varlist       List of character vectors of length 2 
  ###########################################################################
  # op            List with names
  #  outfile
  #  by.var       By variable in pheno.list data
  #  exclude      Character string of values to exclude
  #  temp.list

  snpFlag <- !is.null(snp.list)
  if ((!snpFlag) && (is.null(getListName(pheno.list, "id.var")))) {
    pheno.list$id.var <- 1
  }

  # Check the input lists
  pheno.list <- check.pheno.list(pheno.list)

  outfile <- getListName(op, "outfile")
  outflag <- !is.null(outfile)
  by.var  <- getListName(op, "by.var")
  byflag  <- !is.null(by.var)

  if (snpFlag) {
    snp.list  <- check.snp.list(snp.list)
    temp.list <- getListName(op, "temp.list")
  
    pheno.list$remove.miss <- 0
    pheno.list$make.dummy  <- 0
    temp <- getListName(pheno.list, "keep.vars")
    if (!is.null(temp)) {
      pheno.list$keep.vars <- c(pheno.list$id.var, temp, by.var)
    }
    snp.list$stop.vec <- 2

    # Get the data vector of snps
    tlist <- list(include.row1=0, include.snps=0, return.type=1,
                missing=1, snpNames=1, orderByPheno=1, return.pheno=1)

    temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
    if (class(temp) == "try-error") {
      print(temp)
      stop("ERROR loading data")
    }

    # Get the phenotype data
    phenoData.list <- temp$phenoData.list
    x              <- phenoData.list$data

    rm(temp, phenoData.list, tlist)
    temp <- gc()
  } else {
    x <- loadData(pheno.list$file, pheno.list)
    pheno.list$data        <- data
    pheno.list$is.the.data <- 1
    x <- getPhenoData(pheno.list)$data 
    rm(pheno.list)
    gc()
  }

  if (outflag) sink(outfile)
  if (byflag) {
    byLevels <- unique(unfactor(x[, by.var]))
  } else {
    byLevels <- 1
  }

  len <- length(varlist)
  n   <- nrow(x)
  exclude <- getListName(op, "exclude")
  if (is.null(exclude)) exclude <- "NULL"
  suffix <- paste(", exclude=", exclude, ")", sep="") 

  for (level in byLevels) {
    if (byflag) {
      temp <- (x[, by.var] == level)
      print(level)
    } else {
      temp <- rep(TRUE, times=n)
    }
    
    for (i in 1:len) {
      vars <- varlist[[i]]
      str  <- parse.vec(length(vars), vec.prefix="x[temp, vars[", vec.suffix="]]",
                   prefix="table(", suffix=suffix, delimiter=", ")
      tab  <- eval(parse(text=str))
      print(vars)
      print(tab) 
    }
  }
  
  if (outflag) sink()

  0

} # END: crossTab

# Function to merge the phenotype and genotype data
mergePhenoGeno <- function(snp.list, pheno.list, temp.list=NULL, op=NULL) {

  # snp.list
  # pheno.list
  # temp.list
  # op
  #   outfile
  #   which    0-2  0=0-1-2 coding, 1=original coding, 2=both
  #            Then default is 0

  op <- default.list(op, c("which"), list(0))

  # Check the input lists
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)

  which <- op$which
  if (which) snp.list$recode <- 0

  # Get the data vector of snps
  tlist <- list(include.row1=0, include.snps=0, return.type=1,
                missing=1, snpNames=1, orderByPheno=1, return.pheno=1)

  temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
  if (class(temp) == "try-error") {
    print(temp)
    stop("ERROR loading data")
  }

  snpData   <- temp$data
  snpNames  <- temp$snpNames
  delimiter <- getDelimiter(snp.list)
  nsnps     <- length(snpData)

  # Get the phenotype data
  phenoData.list <- temp$phenoData.list
  phenoData0     <- phenoData.list$data

  # Determine if cc.var was specified
  ccVar <- getListName(pheno.list, "cc.var")
  if (!is.null(ccVar)) {
    subset <- phenoData0[, ccVar] == 0
    subset[is.na(subset)] <- FALSE
  } else {
    subset <- rep(TRUE, times=nrow(phenoData0))
  }

  for (i in 1:nsnps) {
    new <- snpNames[i]
    if (which == 2) new <- paste(new, "_GENO", sep="")
    phenoData0[, new] <- getVecFromStr(snpData[i], delimiter=delimiter)

    if (which == 2) {
      phenoData0[, snpNames[i]] <- recode.geno(phenoData0[, new], in.miss="NA", subset=subset)$vec
    }
  }

  out <- getListName(op, "outfile")
  if (!is.null(out)) writeTable(phenoData0, out)

  phenoData0

} # END: mergePhenoGeno

# Function to parse a vector. Returns the string to evaluate.
parse.vec <- function(vec.len, vec.prefix="", vec.suffix="",
                      prefix="", suffix="", delimiter=", ") {

  # vec.len
  # vec.prefix    Constant string (for now)
  # vec.suffix    Constant string (for now)

  str <- prefix 
  for (j in 1:vec.len) {
    str <- paste(str, vec.prefix, j, vec.suffix, sep="")
    if (j < vec.len) str <- paste(str, delimiter, sep="")
  }
  str <- paste(str, suffix, sep="")
  str

} # END: parse.vec

# Function to replace strings in levels of a categorical variable
replaceStr.var <- function(data, var, str=" ", newStr="") {

  data[, var] <- gsub(str, newStr, data[, var], extended=FALSE)

  data

} # END: replaceStr.var

# Function to replace strings in a data frame
replaceStr.list <- function(data, varlist) {

  for (i in 1:length(varlist)) {
    tlist <- varlist[[i]]
    tlist <- default.list(tlist, c("var", "str", "newStr"),
                          list("ERROR", " ", ""), error=c(1, 0, 0))
    data[, tlist$var] <- gsub(tlist$str, tlist$newStr, data[, tlist$var], extended=FALSE)
  }
  data

} # END: replaceStr.var

# Function to add columns to a data frame from another file or data frame
# A data frame will be returned
addColumn <- function(data, id.var, file.list, op=NULL) {

  # data        Data frame or matrix. 
  # id.var      ID variables on data (either 1 or 2 variables)
  #             No default
  # file.list   List of type file.list with the additional names
  #   id.var    Id variables to match against id.var
  #             The order of this vector must match id.var.
  #             No default
  #   vars      Vars to copy to data.
  #             No default
  #   names     New variable names of vars to copy to data
  #             The default is vars
  #   type      Vector of "C" or "N" for character/numeric variables
  #             The default is NULL
  #######################################################################
  # op          List with names:
  #  leading0   0 or 1 to check for leading zeros in the id variable
  #             The default is 0
  #  initValue  Initial value for new columns
  #             The default is NA.
  #  replace    0 or 1 to replace the new variables
  #             The default is 1

  op <- default.list(op, c("leading0", "initValue", "replace"), list(0, NA, 1))

  file.list <- default.list(file.list,
               c("file", "file.type", "delimiter", "header",
                 "id.var", "vars"), 
               list("ERROR", 3, "\t", 1, "ERROR", "ERROR"), 
               error=c(1, 0, 0, 0, 1, 1))

  check.vec(id.var, "id.var", list(checkList=colnames(data), maxValue=ncol(data), minValue=1))

  type  <- getListName(file.list, "type")
  vars  <- file.list$vars
  id    <- file.list$id.var
  names <- file.list$names
  if (is.null(names)) names <- vars
  nvars  <- length(vars)
  if (!nvars) return(data)
  if (nvars != length(names)) {
    stop("ERROR with vars/names field in file.list")
  }
  typeFlag <- !is.null(type)
  if (typeFlag) {
    type <- toupper(type)
    if (length(type) == 1) type <- rep(type, times=nvars)
    if (nvars != length(type)) {
      stop("ERROR with vars/type field in file.list")
    }
  }
  nid <- length(id.var)
  if (length(id) != nid) {
    stop("ERROR in addColumn: length(id.var) != length(file.list$id.var)")
  }

  if (!is.data.frame(data)) {
    dfFlag <- 0
    if (typeFlag) {
      if (length(unique(type)) > 1) typeFlag <- 0
    }
  } else {
    dfFlag <- 1
  }

  # Read in the file
  if (file.list$file.type %in% c(3, 6, 8)) {
    file.list$method       <- 2
    file.list$what         <- "character"
    file.list$returnMatrix <- 1
    file.list$include.row1 <- file.list$header
  }
  if (is.data.frame(file.list$file)) {
    x <- file.list$file
    file.list$file <- NULL
    temp <- gc()
  } else {
    temp <- checkVars(file.list, file.list$vars)
    x <- loadData(file.list$file, file.list)
  }

  # Check id
  if (nid == 1) {
    if (length(unique(x[, id])) != nrow(x)) {
      print("WARNING: ids are not unique in file")
    }
    if (length(unique(data[, id.var])) != nrow(data)) {
      print("WARNING: ids are not unique in data")
    }
  } else {
    # Compare the id variables
    if (length(unique(data[, id.var[2]])) > length(unique(data[, id.var[1]]))) {
      # Switch
      temp      <- id.var[1]
      id.var[1] <- id.var[2]
      id.var[2] <- temp
      temp      <- id[1]
      id[1]     <- id[2]
      id[2]     <- temp
    } 
  }

  # Check variable names
  dimx <- dim(x)
  temp <- colnames(x)
  check.vec(vars, "file.list$vars", list(checkList=temp, maxValue=dimx[2], minValue=1))

  for (v in c(id, vars)) x[, v] <- unfactor(x[, v])
  for (v in id.var) data[, v] <- unfactor(data[, v])

  # Remove leading zeros
  if (op$leading0) {
    for (v in id.var) {
      data[, v] <- removeLeading0(data[, v])
    }
    for (v in id) {
      x[, v] <- removeLeading0(x[, v])
    }
  }

  # Get unique id.var[2]
  if (nid > 1) {
    levels <- unique(data[, id.var[2]])
  } else {
    levels <- 1
  }
  ret     <- NULL
  replace <- op$replace
  for (lev in levels) {
    if (nid > 1) {
      temp <- data[, id.var[2]] == lev
      d2   <- removeOrKeepRows(data, temp, which=1)
      temp <- x[, id[2]] == lev
      x2   <- removeOrKeepRows(x, temp, which=1)
    } else {
      d2 <- data
      x2 <- x
    } 

    # Match ids
    temp <- x2[, id[1]] %in% d2[, id.var[1]]
    x2   <- removeOrKeepRows(x2, temp, which=1)
    rows <- match(d2[, id.var[1]], x2[, id[1]])
    temp <- !is.na(rows)
    rows <- rows[temp]
    nr   <- length(rows)

    # For a matrix, define a matrix to combine
    if (!dfFlag) {
      i <- matrix(data=op$initValue, nrow=nrow(d2), ncol=nvars)
      colnames(i) <- names
      d2 <- cbind(d2, i)
    }

    # Add the columns
    for (i in 1:nvars) {
      if ((replace) && (dfFlag)) d2[, names[i]] <- op$initValue
      if (nr) {
        vec <- x2[rows, vars[i]] 
        if (typeFlag) {
          if (type[i] == "N") {
            vec <- as.numeric(vec)
          } else if (type[i] == "C") {
            vec <- as.character(vec)
          }
        }
        d2[temp, names[i]] <- vec
      }
    }
    ret <- rbind(ret, d2)

  } # END: for (lev in levels)

  ret

} # END: addColumn

# Function to order columns in a matrix or data frame
orderVars <- function(data, order) {

  # data     matrix or data frame with column names
  # order    Character vector

  if (ncol(data) == 1) return(data)
  cnames <- colnames(data)
  temp   <- order %in% cnames
  order  <- order[temp]
  temp   <- !(cnames %in% order)
  if (any(temp)) order <- c(order, cnames[temp])
    
  data <- data[, order]
  data

} # END: orderVars

# Function to remove leading zeros in character vectors
removeLeading0 <- function(vec) {

  # vec    Character vector

  if (is.numeric(vec)) return(vec)

  ret    <- vec
  numvec <- as.numeric(vec)
  temp   <- !is.na(numvec)
  if (any(temp)) ret[temp] <- as.character(numvec[temp])
  ret

} # END: removeLeading0


# Function to check that the specified delimiter at leasts exists in the first
#   row of the file.
checkDelimiter <- function(file.list) {

  file.list <- default.list(file.list, c("file", "file.type", "delimiter"),
                 list("ERROR", "ERROR", "ERROR"), error=c(1, 1, 1))

  if (file.list$file.type == 1) return(1)

  # Open file
  fid <- getFID(file.list$file, file.list)

  # Read 1 row
  x <- scan(fid, what="character", nlines=1, sep="\n")
  close(fid)

  ret <- grep(file.list$delimiter, x)
  if (!length(ret)) ret <- 0

  ret

} # END: checkDelimiter

# Function to remove leading/trailing white space
removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

} # END: removeWhiteSpace

# Function for changing levels in a matrix or data frame
changeLevels.var <- function(data, var, old.levels, new.level,
                             new.var=NULL) {

  # data
  # var           Variable name or column number.
  #               Use NULL if data is a vector
  # old.levels    Current level(s) to change to new.level
  # new.level
  # new.var       New variable name or the old one if NULL
 
  # NA works with %in%

  if (is.null(dim(data))) {
    temp <- data %in% old.levels
    temp[is.na(temp)] <- FALSE
    data[temp] <- new.level
  } else {
    if (is.null(new.var)) new.var <- var
    temp <- data[, var] %in% old.levels
    temp[is.na(temp)] <- FALSE
    data[temp, new.var] <- new.level
  }

  data  

} # END: changeLevels.var

# Function for debugging 
debug.time <- function(time0, str=NULL) {

  if (!is.null(time0)) print(proc.time()-time0)
  if (!is.null(str)) print(str)
  return(proc.time())

} # END: debug.time

# Function to compute all interactions between 2 sets of vars
getInteractions <- function(data, vars1, vars2) {

  vars1 <- unique(vars1)
  vars2 <- unique(vars2)
  nv1   <- length(vars1)
  nv2   <- length(vars2)
  nv    <- nv1*nv2
  flip  <- 0
  if ((nv2 == 1) && (nv1 > 1)) {
    # Flip vars for efficiency
    temp  <- vars1
    vars1 <- vars2
    vars2 <- temp
    temp  <- nv1
    nv1   <- nv2
    nv2   <- temp
    flip  <- 1
  }
 
  new     <- matrix(data=NA, nrow=nrow(data), ncol=nv)
  newVars <- character(nv)
  if ((nv1 == 1) && (nv2 == 1)) {
    new[, 1] <- makeVector(data[, vars1])*makeVector(data[, vars2])
    newVars  <- paste(vars1, ".", vars2, sep="") 
  } else {
    start <- 1
    stop  <- nv2
    for (i in 1:nv1) {
       new[, start:stop] <- matrixMultVec(as.matrix(data[, vars2]), 
                              makeVector(data[, vars1[i]]), by=2)
       if (flip) {
         newVars[start:stop] <- paste(vars2, ".", vars1[i], sep="")
       } else {
         newVars[start:stop] <- paste(vars1[i], ".", vars2, sep="")
       }
       start <- stop + 1
       stop  <- stop + nv2
    }
  }

  colnames(new) <- newVars 
  data <- cbind(data, new) 
  
  list(data=data, newVars=newVars)

} # END: getInteractions

# Function to check for an error with try function
checkTryError <- function(obj, conv=1) {

  classObj <- class(obj)
  if ("try-error" %in% classObj) return(1)
  ret <- 0  

  # Check for convergence
  if (conv) {
    if (("glm" %in% classObj) || ("lm" %in% classObj)) {
      ret <- 1 - obj$converged
    } else if ("vglm" %in% classObj) {
      temp <- try(obj@criterion$loglikelihood, silent=TRUE)
      if ("try-error" %in% classObj) return(1)
      if ((!is.finite(temp)) || (is.null(temp))) ret <- 1
    } else if ("snp.logistic" %in% classObj) {
      if (is.null(obj$UML)) return(1)
    }
  }
  
  ret

} # END: checkTryError

# Function to check a list of type file.list
check.file.list <- function(flist, op=NULL) {

  # op          List with names 
  #  exist
  #  vars
 
  flist <- default.list(flist, c("file", "header"), list("ERROR", 1), 
                        error=c(1, 0))
  if (is.null(flist[["file.type", exact=TRUE]])) {
    flist$file.type <- getFileType(flist$file)
  }
  if (is.null(flist[["delimiter", exact=TRUE]])) {
    flist$delimiter <- getFileDelim(flist$file, type=flist$file.type)
  }

  op <- default.list(op, c("exist"), list(1)) 
  if (op$exist) {
    if (check.files(flist$file)) stop()
  }
  slist <- op[["subsetData", exact=TRUE]]
  if (!is.null(slist)) {
    slist <- check.subsetData.list(slist) 
    svars <- getSubsetDataVars(slist) 
  } else {
    svars <- NULL
  }
  vars <- op[["vars", exact=TRUE]]
  vars <- unique(c(vars, svars))
  if (!is.null(vars)) checkVars(flist, vars) 

  flist

} # END: check.file.list

# Function to return the variables in a list of type subsetData
getSubsetDataVars <- function(slist) {

  n   <- length(slist)
  ret <- character(n)
  for (i in 1:n) {
    ret[i] <- slist[[i]]$var
  }

  ret

} # END: getSubsetDataVars

# Function to check a list of type subsetData
check.subsetData.list <- function(slist) {

  n <- length(slist)
  for (i in 1:n) {
    temp <- default.list(slist[[i]], c("var", "operator", "value"), 
                         list("ERROR", "ERROR", "ERROR"), error=c(1, 1, 1))
  }

  slist

} # END: check.subsetData.list

# Function to normalize variable names
normVarNames <- function(cvec, op=NULL) {
 
  ret <- cvec
  ret <- gsub(">=", "_GTEQ_", ret, perl=TRUE)
  ret <- gsub("<=", "_LTEQ_", ret, perl=TRUE)
  ret <- gsub("==", "_EQ_", ret, perl=TRUE)
  ret <- gsub(">", "_GT_", ret, perl=TRUE)
  ret <- gsub("<", "_LT_", ret, perl=TRUE)
  ret <- gsub("%in%", "_IN_", ret, perl=TRUE)
  ret <- gsub("!=", "_NEQ_", ret, perl=TRUE)
  ret <- gsub(" ", ".", ret, perl=TRUE)

  str <- "[~`'!@#$%^&*()-+={}|\ ;<>?//]"
  ret <- gsub(str, "", ret, perl=TRUE)
  ret

} # END: normVarNames



# History: May 02 2008 Add getWaldTest
#          May 14 2008 Make getWald test more general
#          May 30 2008 Add getMAF function
#          Jun 03 2008 Use matchNames function in getWaldTest.
#          Jun 04 2008 Add getSummary function
#          Jun 10 2008 Add setUpSummary function
#          Jun 16 2008 Rename getLoglike to getLoglike.glm
#          Jun 18 2008 Redo getWaldTest
#          Jun 30 2008 Let getWaldTest call setUpSummary
#          Jul 11 2008 Add getPermutation function
#          Jul 15 2008 In getWaldTest, use the try function
#          Jul 17 2008 Add getGenoCounts
#          Jul 18 2008 Remove missing names from getWaldTest
#                      Remove inc.Beta.se option
#          Jul 21 2008 Generalize likelihood ratio to lists
#                      Add getXBeta function
#          Jul 25 2008 Add getDesignMatrix function
#          Aug 05 2008 Add effects functions
#          Aug 08 2008 Compute standard errors for effects
#                      Add getLinearComb.var
#          Aug 12 2008 Change getDesignMatrix
#          Aug 13 2008 Add getModelData
#          Aug 13 2008 Use getModelData in callGLM
#          Aug 14 2008 Change to getPermutation function
#          Aug 21 2008 Change to effects.init function
#          Aug 25 2008 Add effects function
#          Aug 29 2008 Add waldTest.main
#          Aug 29 2008 Add getSummary.main
#          Sep 29 2008 Add getORfromLOR and getCI
#          Oct 01 2008 Fix bug in computing LR p-value
#          Oct 22 2008 Fix bug in getDesignMatrix
#          Oct 23 2008 Add pvalue.normal function
#          Nov 18 2008 Add unadjustedGLM.counts
#          Dec 05 2008 Catch errors in unadjustedGLM.counts
#          Jan 15 2009 Add GC.adj.pvalues and inflationFactor functions
#          Feb 03 2009 Change in getMAF
#          Mar 26 2009 Add swap2cols.cov function
#          Apr 04 2009 Fix bug in getDesignMatrix
#          Apr 27 2009 Modify effects.init return list
#          Apr 28 2009 Change in effects.init for the first row of a
#                      stratified effects table.
#          Jun 04 2009 Add heterTest function
#          Jun 12 2009 Add freqCounts.var function
#                      Add standardize.z function
#          Jun 19 2009 Fix bug in getDesignMatrix for interaction
#                      matrices of 1 column, and setting the colnames
#          Jul 15 2009 Add inflation factor argument to GC.Adj.Pvalues function
#          Jul 20 2009 Update freqCounts.var for left endpoint = right endpoint
#          Jul 28 2009 Update freqCounts.var for left endpoint = right endpoint
#          Jul 30 2009 Update freqCounts.var to include data frame
#          Aug 03 2009 Change output of effects.init
#          Oct 06 2009 Change levels in freqCounts.var for leftEndClosed
#          Oct 09 2009 Update freqCounts.var to pass in labels
#          Oct 23 2009 Add function to check the convergence of an object
#          Oct 26 2009 Add functions for likelihood ratio test, Wald test
#          Oct 28 2009 Add function dsgnMat
#          Dec 28 2009 Let a colon be the seperator in snp.effects for the interaction
#          Dec 29 2009 Add option removeInt to dsgnMat
#          Mar 03 2010 Update getWaldTest, getEstCov for snp.matched class
#          Mar 03 2010 Update snp.effects for the snp.matched class
#          Mar 13 2010 Add method option to getSummary, getWaldTest
#                      Compute stratified effects in snp.effects
#                      Add function to print snp.effects object
#          Mar 15 2010 Add method option to snp.effects
#          Mar 18 2010 Add function getGenoStats

# Function to return point estimates and covariance matrix from an object
getEstCov <- function(fit) {

  clss <- class(fit)
  if (any(clss == "snp.logistic")) {
    methods <- c("UML", "CML", "EB")
    ret <- list(methods=methods)
    for (method in methods) {
      temp <- fit[[method, exact=TRUE]]
      if (!is.null(temp)) {
        ret[[method]] <- list(estimates=temp$parms, cov=temp$cov)
      } 
    }
    return(ret)
  } else if (any(clss == "glm")) {
    parms <- fit$coefficients
    fit   <- summary(fit)
    cov   <- fit$cov.scaled
  } else if (any(clss == "coxph")) {
    parms         <- fit$coefficients
    cov           <- fit$var
    vnames        <- names(parms)
    rownames(cov) <- vnames
    colnames(cov) <- vnames
  } else if (any(clss == "vglm")) {
    parms <- fit@coefficients
    fit   <- summary(fit)
    cov   <- fit@cov.unscaled
  } else if (any(clss == "snp.matched")) {
    methods <- c("CLR", "CCL", "HCL")
    ret <- list(methods=methods)
    for (method in methods) {
      temp <- fit[[method, exact=TRUE]]
      if (!is.null(temp)) {
        ret[[method]] <- list(estimates=temp$parms, cov=temp$cov)
      } 
    }
    return(ret)
  } else {
    parms <- fit$parms
    if (is.null(parms)) parms <- fit$coefficients 
    cov <- fit$cov
    if (is.null(cov)) cov <- fit$cov.scaled 
  }

  list(estimates=parms, cov=cov)

} # END: getEstCov

# Function to get the log-likelihood from a glm object
getLoglike.glm <- function(model) {
  # model

  # AIC = -2*loglike + 2*(# of parms)
  (2*model$rank - model$aic)/2 

} # END: getLoglike.glm

# Main function for likelihood ratio test
likelihoodRatio.main <- function(ll1, rank1, ll2, rank2) {

  # ll1, ll2   log-likelihood values

  df     <- abs(rank1 - rank2)
  test   <- 2*abs(ll1 - ll2)
  if (!df) {
    pvalue <- 1
  } else {
    pvalue <- pchisq(test, df, lower.tail=FALSE)
  }
  list(test=test, df=df, pvalue=pvalue)

} # END: likelihoodRatio.main

# Function to return the log-likelihood and rank of an object
loglikeAndRank <- function(fit) {

  clss <- class(fit)
  if (any(clss == "glm")) {
    rank <- fit$rank
    ll   <- (2*rank - fit$aic)/2
  } else if (any(clss == "coxph")) {
    rank <- sum(!is.na(fit$coefficients))
    ll   <- max(fit$loglik)
  } else if (any(clss == "vglm")) {
    rank <- fit@rank
    ll   <- fit@criterion$loglikelihood
  } else {
    rank <- fit$rank
    ll   <- fit$loglike
  } 

  list(loglike=ll, rank=rank)

} # END: loglikeAndRank

# Function to do a likelihood ratio test
likelihoodRatio <- function(model1, model2) {
  # model1    Return object from glm, lm, snp.logistic, coxph, vglm or list
  #           with names "loglike" and "rank"
  # model2

  l1 <- loglikeAndRank(model1)
  l2 <- loglikeAndRank(model2)

  ret <- likelihoodRatio.main(l1$loglike, l1$rank, l2$loglike, l2$rank) 

  ret

} # END: likelihoodRatio

# Function to compute the Wald test (2 - sided) 
waldTest.main <- function(parms, cov, parmNames) {

  # parms      Parameter vector
  # cov        Covariance matrix
  # parmNames  Character or numeric vector of parameters to test

  df     <- length(parmNames)
  nrcov  <- nrow(cov)
  vnames <- names(parms)

  if (is.numeric(parmNames)) {
    vpos <- parmNames

    # Update parms and cov
    parms <- parms[vpos]
    cov   <- cov[vpos, vpos]
    np    <- length(parms)

  } else {
    # Remove missing names in parmNames
    vnames <- vnames[vnames %in% parmNames]
     
    # Update the parameter vector (the name is kept if length(vnames) = 1)
    parms <- parms[vnames]

    # Remove missing values
    temp   <- !is.na(parms)
    parms  <- parms[temp]
    vnames <- vnames[temp]

    # Check for error
    np <- length(parms)
    if (!np) return(list(test=NA, df=0, pvalue=NA))
    
    # Update cov
    cov <- cov[vnames, vnames]
  }

  # See if matrix is invertible
  temp <- try(solve(cov), silent=TRUE)
  if (class(temp) == "try-error") {
    return(list(test=NA, df=np, pvalue=NA))
  } 

  # Get the test statistic
  dim(parms) <- c(np, 1)
  test       <- t(parms) %*% temp %*% parms
  dim(test)  <- NULL

  list(test=test, df=np, pvalue=pchisq(test, df=np, lower.tail=FALSE))

} # END: waldTest.main

# Function to compute the Wald test (2 - sided)
getWaldTest <- function(fit, parmNames, method=NULL) {

  # fit        Return object from glm, list with names "coefficients"
  #            and "cov.scaled", return object from snp.logistic or snp.matched.
  # parmNames  Character or numeric vector of parameters to test
  # method     

  estcov <- getEstCov(fit)
  methods <- estcov[["methods", exact=TRUE]]
  if (!is.null(methods)) {
    ret <- list()
    for (m in methods) {
      temp <- estcov[[m, exact=TRUE]]
      if (!is.null(temp)) ret[[m]] <- waldTest.main(temp$estimates, temp$cov, parmNames)
    } 
    if (!is.null(method)) ret <- ret[[method, exact=TRUE]]
    return(ret)
  } 

  return(waldTest.main(estcov$estimates, estcov$cov, parmNames))    
  
} # END: getWaldTest

# Function to compute the threshold p-value
threshold.trunc.prod <- function(pvals, threshold=0.05) {

  # pvals      Vector or matrix of p-values   
  # threshold  The default is 0.05

  dim(pvals) <- NULL

  # Get the pvalues less than threshold
  pvals <- pvals[((pvals < threshold) & (!is.na(pvals)))]

  ret <- exp(sum(log(pvals)))

  ret

} # END: rankTruncate.pvals

# Function to compute rank truncated p-value
rank.trunc.prod <- function(pvals, k=10) {

  # pvals   Vector or matrix of p-values   
  # k       Maximum number of p-values to use

  # Get the sorted p-values
  pvals <- sort(pvals)
  n     <- min(length(pvals), k)
  ret   <- exp(sum(log(pvals[1:n])))

  ret

} # END: rankTruncate.pvals

# Function to compute score test for logistic reg
score.logReg <- function(fit, mat) {

  # fit      Return object from glm with x=TRUE and y=TRUE in the call
  # mat      A single factor, matrix, or data frame of variables to test

  # See if mat is a factor
  if (is.factor(mat)) mat <- data.frame(mat)
  if (is.data.frame(mat)) {
    mat <- as.matrix(createDummy(mat)$data)
  }

  # Get the number of columns of mat
  df <- ncol(mat)
  if (is.null(df)) df <- 1

  temp <- exp(fit$linear.predictors)
  p    <- temp/(1 + temp)  

  # Add the new vector to x
  x  <- cbind(fit$x, mat)
  
  # Get the gradient
  U      <- colSums(matrixMultVec(x, fit$y-p, by=2)) 
  n      <- length(U)
  dim(U) <- c(n, 1)

  # Free memory
  rm(fit, mat)
  temp <- gc(verbose=FALSE)

  # Let p = p*(1-p)
  p <- p*(1 - p)

  # Get the negative Hessian
  hess <- matrix(data=NA, nrow=n, ncol=n)
  for (i in 1:n) {
    temp      <- p*x[, i]
    hess[i, ] <- colSums(matrixMultVec(x, temp, by=2)) 
    hess[, i] <- hess[i, ]
  }

  # Invert the hessian
  hess <- try(solve(hess), silent=TRUE)
  if (class(hess) == "try-error") {
    warning("Singular hessian matrix")
    return(list(test=NA, df=df, pvalue=NA))
  } 

  # Compute the chi-squared test statistic
  test <- t(U) %*% hess %*% U
  dim(test) <- NULL

  pvalue <- pchisq(test, df=df, lower.tail=FALSE)

  list(test=test, df=df, pvalue=pvalue)

} # score.logReg

# Function to compute minor allele frequency
getMAF <- function(genotype, sub.vec=NULL, controls=0) {

  # genotype       Vector of genotypes coded as 0, 1, 2, NA
  #                No default
  # sub.vec        NULL or case/control vector for only using
  #                a subset of the genotypes.
  #                The default is NULL
  # controls       Vector of values describing the controls in 
  #                sub.vec.
  #                The default is 0.

  temp <- !is.na(genotype)
  if (!is.null(sub.vec)) temp <- temp & (sub.vec %in% controls)
  genotype <- genotype[temp]
  ng <- length(genotype)
  if (!ng) return(NA)

  freq <- (sum(genotype==1) + 2*sum(genotype==2))/(2*ng)
  freq

} # END: getMAF 

# Function to return summary information for parameters
getSummary.main <- function(parms, cov, sided=2) {

  # parms   Vector of parameters
  # cov     Covariance matrix 
  # sided   1 or 2 

  if (sided != 1) sided <- 2
  cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")

  n    <- length(parms)
  ret  <- matrix(data=NA, nrow=n, ncol=4)
  pnames <- names(parms)
  rownames(ret) <- pnames
  colnames(ret) <- cols
  ret[, 1] <- parms
  
  cols <- colnames(cov)
  cov  <- sqrt(diag(cov))
  names(cov) <- cols
  
  # Get the correct order
  if (is.null(pnames)) pnames <- 1:n
  cov <- cov[pnames]
  ret[, 2] <- cov

  ret[, 3] <- parms/cov 
  ret[, 4] <- sided*pnorm(abs(ret[, 3]), lower.tail=FALSE)
  ret

} # END: getSummary.main

# Function to return summary information for parameters
getSummary <- function(fit, sided=2, method=NULL) {

  # fit     Return object from glm, snp.logistic, or a list
  #         with names "parms" and "cov".
  # sided   1 or 2 

  clss <- class(fit)

  # snp.logistic
  if (any(clss %in% "snp.logistic")) {
    if (is.null(method)) {
      methods <- c("UML", "CML", "EB")
    } else {
      methods <- toupper(method)
    }
    ret <- list()
    for (m in methods) {
      temp <- fit[[m, exact=TRUE]]
      if (!is.null(temp)) {
        ret[[m]] <- getSummary.main(temp$parms, temp$cov, sided=sided)
      } 
    }
    return(ret)
  } 

  # snp.matched
  if (any(clss %in% "snp.matched")) {
    if (is.null(method)) {
      methods <- c("CLR", "CCL", "HCL")
    } else {
      methods <- toupper(method)
    }
    ret <- list()
    for (m in methods) {
      temp <- fit[[m, exact=TRUE]]
      if (!is.null(temp)) {
        ret[[m]] <- getSummary.main(temp$parms, temp$cov, sided=sided)
      } 
    }
    return(ret)
  } 

  # GLM
  if ("glm" %in% clss) fit <- summary(fit)
  if (class(fit) == "summary.glm") {
    cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")
    if (sided != 1) sided <- 2

    fit$coefficients[, 4] <- 
       sided*pnorm(abs(fit$coefficients[, 3]), lower.tail=FALSE)
    colnames(fit$coefficients) <- cols
    return(fit$coefficients)
  }

  # List
  parms <- fit$parms
  if (is.null(parms)) {
    parms <- fit$coefficients
    if (is.null(parms)) return(NULL)
  }
  cov <- fit$cov
  if (is.null(cov)) {
    cov <- fit$cov.scaled
    if (is.null(cov)) return(NULL)
  }
  return(getSummary.main(parms, cov, sided=sided))      

} # END: getSummary

# Function to take a parameter vector and covariance matrix and
#  output an nx2 matrix that can be used with getWaldTest
setUpSummary <- function(parms, cov) {

  cnames         <- colnames(cov)
  nc             <- ncol(cov)
  if (is.null(nc)) nc <- 1
  temp           <- matrix(data=NA, nrow=nc, ncol=2)
  colnames(temp) <- c("Estimate", "Std. Error")
  rownames(temp) <- cnames
  temp[, 2]      <- sqrt(diag(cov))
 
  # Match the parameter names (there could be NAs in the vector of
  #   point estimates)
  if ((!is.null(cnames)) && (!is.null(names(parms)))) {
    parms <- parms[cnames]
  }
  temp[, 1] <- parms
 
  temp <- list(coefficients=temp, cov.scaled=cov)
  temp
 
} # END: setUpSummary

# Function to call for permutations
getPermutation <- function(fit0, nsub, perm.method=1) {

  # fit0           Base model fit
  # nsub
  # perm.method    1-3

  if (perm.method == 3) {
    # For gaussian family
    # Permute the residuals
    errors <- sample(fit0$residuals)

    # Add to linear predictor
    response <- fit0$linear.predictors + errors

    perm <- 1:nsub
  } else if (perm.method == 2) {
    # For binomial family
    perm     <- 1:nsub
    response <- rbinom(nsub, 1, fit0$fitted.values)    
  } else {
    perm     <- sample(1:nsub)
    response <- fit0$y
  }

  list(response=response, perm=perm)

} # END: getPermutation

# Function to call glm
callGLM <- function(y, X.main=NULL, X.int=NULL, int.vec=NULL,
                    family="binomial", prefix="SNP_", retX=TRUE,
                    retY=TRUE, inc.int.vec=1, int.vec.base=0) {

  # y           Response vector
  # X.main      Matrix of main effects (without intercept and int.vec)
  # X.int       Matrix for interactions
  # int.vec     Interaction vector or factor 
  # family
  # inc.int.vec 0 or 1 to include the interacting vector in the model

  temp <- getModelData(y, int.vec, X.main=X.main, X.int=X.int, 
                    prefix=prefix, inc.snp=inc.int.vec, 
                   snp.base=int.vec.base)

  y <- temp$y
  X <- temp$design

  if (is.null(prefix)) {
    fit <- glm(y ~ X-1, family=family,
            model=FALSE, x=retX, y=retY) 
  } else {
    fit <- glm(y ~ .-1, family=family, data=data.frame(X),
            model=FALSE, x=retX, y=retY)
  }

  fit

} # END: callGLM

# Function to return genotype counts
getGenoCounts <- function(snp, exclude=c(NA, NaN), check=1) {

  ret <- table(snp, exclude=exclude)
  if ((check) && (length(ret) < 3)) {
    temp        <- rep.int(0, times=3)
    names(temp) <- c("0", "1", "2")
    nm          <- names(ret)
    temp[nm]    <- ret
    ret         <- temp
  }
  ret

} # END: getGenoCounts

# Function to compute XBeta (linear.predictors) by matching the names
getXBeta <- function(X, beta) {

  cnames  <- intersect(colnames(X), names(beta))
  X       <- removeOrKeepCols(X, cnames, which=1)
  beta    <- beta[cnames]
  b2      <- beta
  dim(b2) <- c(length(beta), 1)
  ret     <- X %*% b2
  list(X=X, beta=beta, XBeta=ret)

} # END: getXBeta

# Function to return the model data for callGLM
getModelData <- function(y, snp, X.main=NULL, X.int=NULL, 
                    prefix=NULL, inc.snp=1, snp.base=0) {

  pflag <- !is.null(prefix)
  if (pflag) cnames <- colnames(X.main)

  # Append y and intercept to X.main
  X.main <- cbind(y, 1, X.main)
  if (pflag) colnames(X.main) <- c("y", "Intercept", cnames)

  mat <- getDesignMatrix(snp, X.main=X.main, X.int=X.int,
                   inc.snp=inc.snp, X.hasInt=1, prefix=prefix,
                   snp.base=snp.base)
 
  # Remove the response from mat
  y   <- mat[, 1]
  mat <- removeOrKeepCols(mat, 1, which=-1) 
  
  list(y=y, design=mat)

} # END: getModelData

# Function to return a design matrix.
# Missing values are automatically removed
getDesignMatrix <- function(snp, X.main=NULL, X.int=NULL,
                   inc.snp=1, X.hasInt=0, prefix=NULL,
                   snp.base=0) {
  
  # snp         Vector or factor. If a factor, see snp.base
  # X.main      Matrix for main effects
  # X.int       Matrix for interactions with snp
  # inc.snp     0 or 1 to include snp in the model
  # X.hasInt    0 or 1 if X.main has an intercept column
  # prefix      NULL or snp prefix for variable names
  #             If set to NULL, the default variable names
  #             from model.matrix will be kept.
  # snp.base    Baseline category for snp that will be left
  #             out of the returned matrix

  # Input matrices should have column names !!!
  # If column names for X.int contain a colon, then there will be a problem
  # Watch for X.main = NULL

  intFlag <- !is.null(X.int)
  if (intFlag) {
    intNames   <- colnames(X.int)
    intIds     <- grep(":", intNames, extended=FALSE)
    intIdsFlag <- length(intIds)
    if (intIdsFlag) {
      stop("ERROR: X.int column names cannot contain a colon (:)")
      intNames <- gsub(":", ".", intNames, extended=FALSE)
      colnames(X.int) <- intNames
    } 
  }
  pFlag   <- !is.null(prefix)
  snpFlag <- !is.null(snp)
  if (snpFlag) {
    facFlag <- is.factor(snp)
  } else {
    facFlag <- 0
  } 
  if (!snpFlag) {
    inc.snp <- 0
    intFlag <- 0
  }
 
  # An intercept will always be included
  if (!X.hasInt) X.main <- cbind(rep(1, times=length(snp)), X.main)

  if (intFlag) {
    #colnames(X.int) <- NULL
    mat <- model.matrix(~X.main + snp*X.int - 1 - X.int)
  } else {
    if (snpFlag) {
      mat <- model.matrix(~X.main + snp - 1)
    } else {
      mat <- model.matrix(~X.main - 1)
    }
  }

  # Set the column names
  if (pFlag) {
    assign <- attributes(mat)$assign
    max    <- max(assign)
    cnames <- colnames(mat)

    # Main effects
    ids  <- assign == 1
    temp <- cnames[ids]
    # Start from pos 7 (X.main is 6 chars)
    cnames[ids] <- substring(temp, 7)

    # Intercept column
    cnames[1] <- "Intercept"

    # SNP
    if (max > 1) {
      ids  <- assign == 2
      temp <- cnames[ids]
      # Start from pos 4 (snp is 3 chars)
      cnames[ids] <- paste(prefix, substring(temp, 4), sep="")

      # Interactions
        if (max > 2) {
          if (facFlag) {
            string <- "_"
          } else {
            string <- ""
          }

          ids  <- assign == 3
          temp <- cnames[ids]
          nint <- length(intNames)

          # Remove the string "snp"
          temp <- substring(temp, 4)
          temp <- unlist(strsplit(temp, ":", fixed=TRUE))
          n    <- length(temp)
  
          # temp is a character vector, every 2 consecutive elements
          # was from the same list  
          even <- seq(from=2, to=n, by=2)
          temp[even] <- substring(temp[even], 6)
          odd <- seq(from=1, to=n-1, by=2)
          temp[odd] <- paste(prefix, temp[odd], sep="")        
          if (nint == 1) {
            temp[even] <- intNames
          } 
          cnames[ids] <- paste(temp[odd], string, temp[even], sep="")

        } # END: if (max > 2)

    } # END: if (max > 1)

    colnames(mat) <- cnames

  } # END: if (pFlag)

  # Remove columns if needed
  if (inc.snp) {
    if (facFlag) {

      # Make sure baseline column is not in the model
      if (pFlag) {
        var <- paste(prefix, snp.base, sep="")
        if (intFlag) {
          var <- c(var, paste(prefix, snp.base, "_", intNames, sep=""))
        }
      } else {
        var <- paste("snp", snp.base, sep="")
        if (intFlag) {
          var <- c(var, paste(var, ":X.int", intNames, sep=""))
        }
      }
      temp <- var %in% colnames(mat)
      if (any(temp)) mat <- removeOrKeepCols(mat, var[temp], which=-1)
    }
  } else {
    # Get the snp column numbers
    ids      <- attributes(mat)$assign == 2
    if (any(ids)) {
      snp.cols <- (1:length(ids))[ids]
      mat <- removeOrKeepCols(mat, snp.cols, which=-1)
    } 
  }

  mat

} # END: getDesignMatrix

# Function to compute an effects table and standard errors
effects.init <- function(parms, cov, var1, var2, levels1, levels2, 
                     base1=0, base2=0, int.var=NULL, effects=1,
                     sep1="_", base2.name="baseline") {

  # parms       Vector of parameter estimates
  # var1        Vector of parameter names for one of the variables. If the
  #             length of this vector is > 1, then it is assumed that var1
  #             is a categorical variable.
  # var2        
  # levels1     Levels for var1 (snp)
  # levels2     Levels for var2 Can be NULL for a categorical variable
  # base1
  # base2
  # int.var     For no interaction, set to NULL. Otherwise a var1 x var2 
  #             matrix of interaction parameter names. int.var can also
  #             be a vector if length(var1) = 1 or length(var2) = 1.
  #             The order of this matrix must match the order of var1 and
  #             var2.
  #             The default is NULL.
  # effects     1 or 2  1 = joint, 2 = stratified
  # sep1        String to separate var1 with its levels
  #             The default is "_".

  int.flag <- !is.null(int.var)
  nv1      <- length(var1)
  nv2      <- length(var2)
  contv1   <- nv1 == 1  
  contv2   <- nv2 == 1
  joint    <- effects == 1

  # Check the number of interaction variables
  if (int.flag) {
    if (length(int.var) != nv1*nv2) {
      stop("ERROR with int.var")
    }
  }

  # Get the levels of the continuous variables, otherwise we assume
  #  the values are 0-1. 
  if (contv1) {
    nlev1   <- length(levels1)
    p1      <- rep(parms[var1], times=nlev1)
    names1  <- paste(var1, levels1, sep=sep1)
    v1      <- rep(var1, times=nlev1)
  } else {
    levels1 <- c(0, rep.int(1, times=nv1))
    nlev1   <- length(levels1)
    p1      <- c(0, parms[var1])
    base1   <- 0
    names1  <- c("baseline", var1)
    v1      <- c(var1[1], var1)
  }
  if (contv2) {
    nlev2   <- length(levels2)
    p2      <- rep(parms[var2], times=nlev2)
    names2  <- paste(var2, levels2, sep="_")
    v2      <- rep(var2, times=nlev2)
  } else {
    levels2 <- c(0, rep.int(1, times=nv2))
    nlev2   <- length(levels2)
    p2      <- c(0, parms[var2])
    base2   <- 0
    names2  <- c(base2.name, var2) 
    v2      <- c(var2[1], var2)
  }

  # Initialize the matrix of interaction parms
  pint <- matrix(data=0, nrow=nlev1, ncol=nlev2)

  # Get the interaction parm
  if (int.flag) {
    # Remove dimension of int.var if <= 1 categorical var
    if (sum(contv1 + contv2) != 0) dim(int.var) <- NULL 

    if (contv1 && contv2) {
      pint[,] <- parms[int.var]
      int.var <- matrix(data=int.var, nrow=nlev1, ncol=nlev2)
    } else {
      if (!contv1 && contv2) {
        temp <- c(0, parms[int.var])
        for (i in 1:nlev2) pint[, i] <- temp
        int.var <- rep(c(int.var[1], int.var), times=nlev2)
        dim(int.var) <- c(nlev1, nlev2)
      } else if (contv1 && !contv2) {
        temp <- c(0, parms[int.var])
        for (i in 1:nlev1) pint[i, ] <- temp
        int.var <- rep(c(int.var[1], int.var), each=nlev1)
        dim(int.var) <- c(nlev1, nlev2)
      } else {
        # Both are categorical
        pint[1, ] <- 0
        for (i in 2:nlev1) {
          pint[i, ] <- c(0, parms[int.var[i-1,]])
        }
       
        # Add baselines to int.var
        temp <- int.var
        int.var <- matrix(data="", nrow=nlev1, ncol=nlev2)
        int.var[2:nlev1, 2:nlev2] <- temp

        # We can assign anything to these 
        int.var[, 1] <- temp[1,1]
        int.var[1, ] <- temp[1,1]
      }
    }
  } 

  # Get the baseline value(s)
  if (joint) {
    base <- base1*p1[1] + base2*p2[1] + base1*base2*pint[1, 1]
    base <- exp(base)
    base <- rep(base, times=nlev2)
  } else {
    base <- rep(NA, times=nlev2)
    for (j in 1:nlev2) {
      temp <- base1*p1[1] + levels2[j]*p2[j] + base1*levels2[j]*pint[1, j]
      base[j] <- exp(temp) 
    }
  }

  # Initialize matrix of effects 
  eff <- matrix(data=NA, nrow=nlev1, ncol=nlev2)
  colnames(eff) <- names2
  rownames(eff) <- names1
  # Loop over the levels
  for (i in 1:nlev1) {
    for (j in 1:nlev2) {
      temp <- levels1[i]*p1[i] + levels2[j]*p2[j] + 
              levels1[i]*levels2[j]*pint[i, j]
      temp <- exp(temp)
      #if ((joint) || (levels1[i] != base1)) temp <- temp/base[j]    
      #eff[i, j] <- temp
      eff[i, j] <- temp/base[j]
    }
  }

  # Initialize matrix for standard errors
  se <- matrix(data=NA, nrow=nlev1, ncol=nlev2)
  colnames(se) <- names2
  rownames(se) <- names1

  # Get base levels
  sebase1 <- rep(base1, times=nlev1)
  if (joint) {
    sebase2 <- rep(base2, times=nlev2)
  } else {
    sebase2 <- levels2
  }

  # Loop over the levels
  for (i in 1:nlev1) {
    for (j in 1:nlev2) {
      b1 <- sebase1[i]
      b2 <- sebase2[j]

      # Get vector of variables and coefficients
      vars <- c(v1[i], v2[j])
      coef <- c(levels1[i] - b1, levels2[j] - b2)

      if (int.flag) {
        vars <- c(vars, int.var[i, j])
        coef <- c(coef, levels1[i]*levels2[j] - b1*b2)
      }

      # SE
      se[i, j] <- sqrt(getLinearComb.var(vars, cov, coef=coef))
    }
  }
 
  logEffects <- log(eff)
  lower95    <- exp(logEffects - 1.96*se)
  upper95    <- exp(logEffects + 1.96*se)  

  list(effects=eff, lower95=lower95, upper95=upper95, 
       logEffects=logEffects, logEffects.se=se)

} # END: effects.init

# Function to compute an effects table and standard errors from
#  the snp.logistic output
snp.effects <- function(fit, var, var.levels=c(0,1), method=NULL) {

  # fit      Output from snp.logistic or snp.matched
  # var      Name of variable to get effects with snp
  # var.levels   (For continuous var) First level is assumed to be
  #          the baseline level 
  #          The default is NULL.

  if (length(var) != 1) stop("Only 1 variable can be specified")

  # Determine the input object
  temp <- class(fit)
  if (temp == "snp.logistic") {
    which   <- 1
    snp     <- fit$model.info$snpName
    methods <- c("UML", "CML", "EB")
    cnames  <- colnames(fit$UML$cov)
    if (is.null(fit$UML)) return(NULL)
  } else if (temp == "snp.matched") {
    which   <- 2
    snp     <- fit$model.info$snp.vars
    methods <- c("CLR", "CCL", "HCL")
    for (m in methods) {
      temp <- fit[[m, exact=TRUE]]
      if (!is.null(temp)) {
        cnames <- colnames(temp$cov)
        break
      }
    }
  } else {
    stop("fit must be of class snp.logistic or snp.matched")
  }

  if (!is.null(method)) {
    temp <- methods %in% method
    methods <- methods[temp]
    if (!length(methods)) stop("Incorrect method")
  }

  levels <- var.levels
  sep1   <- "_"
  nsnp   <- length(snp)

  if (!(var %in% fit$model.info$main.vars)) {
    stop("var must be a main effect variable")
  }
  
  if (var %in% fit$model.info$factors) {
    
    facFlag    <- 1
    levels     <- levels(fit$model.info$data[, var])
    temp2      <- paste(var, "_", levels, sep="") 
  
    # Get the variables in the model fit
    temp       <- temp2 %in% cnames
    var2       <- temp2[temp]
    base2.name <- temp2[!temp] 
    temp       <- length(base2.name)
    if (!temp) base2.name <- "baseline"
    if (temp > 1) base2.name <- base2.name[1]
  } else {
    facFlag    <- 0
    var2       <- var
    base2.name <- NULL
  }

  levels1 <- 0:2
  if (is.null(levels)) {
    if (is.factor(fit$model.info$data[, var])) {
      base2 <- 0
    } else {
      stop("levels must be specified for a continuous var")
    }
  } else {
    base2 <- levels[1]
  }
  
  # Initialize the return list
  ret   <- list()

  intFlag <- 0
  int.var <- NULL
  if (var %in% fit$model.info$int.vars) intFlag <- 1

  for (var1 in snp) {
    if (intFlag) int.var <- paste(var1, ":", var2, sep="")
  
    for (method in methods) {
      temp2 <- fit[[method]]
      if (is.null(temp2)) next
      tlist <- list()
      eff1 <- effects.init(temp2$parms, temp2$cov, var1, var2, 
                 levels1, levels, base1=0, base2=base2,
                 int.var=int.var, effects=1, sep1=sep1,
                 base2.name=base2.name)
      eff2 <- effects.init(temp2$parms, temp2$cov, var1, var2, 
                 levels1, levels, base1=0, base2=base2,
                 int.var=int.var, effects=2, sep1=sep1,
                 base2.name=base2.name)
      eff3 <- effects.init(temp2$parms, temp2$cov, var2, var1, 
                 levels, levels1, base1=base2, base2=0,
                 int.var=int.var, effects=2, sep1=sep1,
                 base2.name="0")

      # Set attributes  
      attr(eff1, "var1")    <- var1 
      attr(eff1, "var2")    <- var2 
      attr(eff1, "levels1") <- levels1 
      attr(eff1, "levels2") <- levels 
      attr(eff2, "var1")    <- var1
      attr(eff2, "var2")    <- var2
      attr(eff2, "levels1") <- levels1
      attr(eff2, "levels2") <- levels
      attr(eff3, "var1")    <- var2
      attr(eff3, "var2")    <- var1
      attr(eff3, "levels1") <- levels
      attr(eff3, "levels2") <- levels1

      temp <- list(JointEffects=eff1, StratEffects=eff2, StratEffects.2=eff3)
      #temp <- list(JointEffects=eff1, StratEffects=eff2)

      class(temp) <- "snp.effects.method" 

      if (nsnp == 1) {
        ret[[method]] <- temp
      } else {
        tlist[[method]] <- temp
      }      
    }
    if (nsnp != 1) ret[[var1]] <- tlist
  }
  class(ret) <- "snp.effects"

  ret

} # END: snp.effects

# Function to compute the variance of a linear combination of parms
getLinearComb.var <- function(vars, cov, coef=NULL) {

  # vars       Vector variable names or indices in cov
  # cov        Covariance matrix 
  # coef       Coefficients for vars. The order must be the same
  #            The default is all coefficients are 1
 
  n <- length(vars)
  if (is.null(coef)) coef <- rep(1, times=n)
  if (n != length(coef)) stop("ERROR with parms and/or coef")

  sum <- 0
  # Sum up the variances
  for (i in 1:n) {
    sum <- sum + coef[i]*coef[i]*cov[vars[i], vars[i]]
  }
  if (n == 1) return(sum)

  # Sum up the covariances
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      sum <- sum + 2*coef[i]*coef[j]*cov[vars[i], vars[j]]
    }
  }

  sum

} # END: getLinearComb.var

# Function to compute odds ratios and standard errors from
# log-odds ratios and se. 
getORfromLOR <- function(data, lor.var, lor.se.var,
                         or="OR", or.se="OR.SE") {

  data[, or] <- exp(data[, lor.var])

  # Change the standard errors
  temp <- data[, lor.se.var]
  data[, or.se] <- temp*data[, or]

  data

} # END: getORfromLOR

# Function to compute confidence intervals
getCI <- function(data, var, var.se, se=1, 
                  lower="LOWER", upper="UPPER") {
 
  zcrit <- qnorm(0.05/2, lower.tail=FALSE)
  if (se) {
    temp  <- zcrit*data[, var.se]
  } else {
    temp  <- zcrit*sqrt(data[, var.se])
  }

  data[, lower] <- data[, var] - temp
  data[, upper] <- data[, var] + temp 

  data

} # END: getCI

# Function to compute normal pvalues
pvalue.normal <- function(test, sided=2) {
  
  sided*pnorm(test, lower.tail=FALSE)

} # END: pvalue.normal

# Function to perform an unadjusted analysis based on genotype
#  frequency counts for cases and controls
unadjustedGLM.counts <- function(file.list, op=NULL) {

  #################################################################
  # file.list      List of type file.list with additional fields:
  #  caseCounts    Variables for the genotype frequencies 0, 1, 2
  #                among the cases.
  #                No default
  #  controlCounts Variables for the genotype frequencies 0, 1, 2
  #                among the controls.
  #                No default.
  #  covar         Covariate in the model.
  #                The default is c(0, 1, 2)
  #  caseOrder     Order for caseCounts in terms of covar
  #                The default is c(1, 2, 3)
  #  controlOrder  Order for controlCounts in terms of covar
  #                The default is c(1, 2, 3)
  #  caseSep       The default is "/"
  #  controlSep    The default is "/"
  #################################################################
  # op            List with names:
  #  outfile      The default is NULL
  #  copyVars     Variables to copy to the output data set
  #################################################################

  file.list <- default.list(file.list,
       c("file", "file.type", "header", "delimiter", "caseCounts",
         "controlCounts", "covar", "caseOrder", "controlOrder",
         "caseSep", "controlSep"),
       list("ERROR", 3, 1, "\t", "ERROR", "ERROR", c(0,1,2),
            1:3, 1:3, "/", "/"),
       error=c(1,0,0,0,1,1,0,0,0,0,0)
       )
  
  covar  <- file.list$covar
  nn     <- length(covar)
  v1     <- file.list$caseCounts
  v0     <- file.list$controlCounts
  n1     <- length(v1)
  n0     <- length(v0)
  order1 <- file.list$caseOrder
  order0 <- file.list$controlOrder
  if (n1 == nn) {
    v1    <- v1[order1]
    flag1 <- 1
  } else {
    flag1 <- 0
    sep1  <- file.list$caseSep
  }
  if (n0 == nn) {
    v0    <- v0[order0]
    flag0 <- 1
  } else {
    flag0 <- 0
    sep0  <- file.list$controlSep
  }

  # Read in the data
  x <- loadData(file.list$file, file.list)

  vars <- getListName(op, "copyVars")
  if (!is.null(vars)) {
    x <- removeOrKeepCols(x, c(vars, v1, v0), which=1)
  } else {
    x <- removeOrKeepCols(x, c(v1, v0), which=1)
  }
  x  <- unfactor.all(x)
  nr <- nrow(x)

  # Add new variables
  newVars <- c("beta", "se", "test", "pvalue")
  for (var in newVars) x[, var] <- NA

  mat <- matrix(data=NA, nrow=nn, ncol=2)
  for (i in 1:nr) {
    mat[] <- NA

    # Case counts
    if (flag1) {
      mat[, 1] <- as.numeric(unlist(x[i, v1]))
    } else {
      mat[, 1] <- as.numeric(getVecFromStr(x[i, v1], delimiter=sep1))
    }
    
    # Control counts
    if (flag0) {
      mat[, 2] <- as.numeric(unlist(x[i, v0]))
    } else {
      mat[, 2] <- as.numeric(getVecFromStr(x[i, v0], delimiter=sep0))
    }
    temp <- try(glm(mat~covar, family=binomial), silent=TRUE)
    if (class(temp)[1] != "try-error") {
      temp <- summary(temp)$coefficients
      if (nrow(temp) == 2) x[i, newVars] <- temp[2,]
    } 
  }

  temp <- getListName(op, "outfile")
  if (!is.null(temp)) {
    write.table(x, file=temp, sep="\t", row.names=FALSE, quote=FALSE)
  }

  x

} # END: unadjustedGLM.counts

# Function to compute the inflation factor
inflationFactor <- function(tests, squared=0) {

  # tests    Vector of Z-test statistics or squared Z-test
  #          statistics
  # squared  0 or 1 if the tests are already squared

  i1  <- qchisq(0.5, df=1)
  if (!squared) tests <- tests*tests
  i2  <- median(tests, na.rm=TRUE)
  ret <- i2/i1
  ret

} # END: inflationFactor

# Function to compute genomic control adjusted p-values
GC.adj.pvalues <- function(tests, pvals=NULL, ifac=NULL) {

  # tests    Vector of Z-test statistics
  # ifac     NULL or the inflation factor to use
  #          The default is NULL

  test2 <- tests*tests
  if (is.null(pvals)) pvals <- 2*pnorm(abs(tests), lower.tail=FALSE)
  if (is.null(ifac)) ifac <- inflationFactor(test2, squared=1)
  ret  <- pchisq(test2/ifac, df=1, lower.tail=FALSE)
  ret

} # END: GC.adj.pvalues

# Function to swap 2 columns of a covariance matrix
swap2cols.cov <- function(mat, col1, col2, errorCheck=0) {

  # mat         Covariance matrix
  # col1        Column name or number
  # col2        Column name or number
  # errorCheck  0 or 1. If set to 1, then the matrix mat must
  #             have both row and column names for the error
  #             check to be done.
  #             The default is 0.
 
  if (col1 == col2) return(mat)

  nr <- nrow(mat)
  nc <- ncol(mat)
  if (nr != nc) stop("ERROR in swap2cols.cov: mat is not a square matrix")

  # Get row and column names
  cnames <- colnames(mat)
  rnames <- rownames(mat)
  cflag  <- !is.null(cnames)
  rflag  <- !is.null(rnames)

  # Get column numbers if variable names are passed in
  cflag1 <- is.character(col1)
  cflag2 <- is.character(col2)
  if (cflag1 || cflag2) {
    if (!cflag) stop("ERROR in swap2col.cov: mat must have column names")
    col1 <- (1:nr)[temp == col1]
    col2 <- (1:nr)[temp == col2]
  }

  # Let col1 be the smaller column
  temp <- col1
  if (col2 < col1) {
    col1 <- col2
    col2 <- temp
  }

  if (errorCheck && cflag && rflag) {
    mat0 <- mat
  } else {
    errorCheck <- 0
  }

  # Save col1 
  save <- mat[, col1]

  # Change columns col1 and col2
  if (col1 > 1) {
    vec <- 1:(col1-1) 
    mat[vec, col1] <- mat[vec, col2]
    mat[vec, col2] <- save[vec]
  }
  mat[col1, col1] <- mat[col2, col2]
  if (col1+1 < col2) {
    vec <- (col1+1):(col2-1)
    mat[vec, col1] <- mat[vec, col2]
    mat[vec, col2] <- save[vec]
  }
  mat[col2, col1] <- mat[col1, col2]
  mat[col2, col2] <- save[col1]
  if (col2 < nr) {
    vec <- (col2+1):nr
    mat[vec, col1] <- mat[vec, col2]
    mat[vec, col2] <- save[vec]
  }

  # Change rows col1 and col2
  vec  <- c(col1, col2)
  temp <- (1:nr)[-vec] 
  for (i in vec) {
    for (j in temp) mat[i, j] <- mat[j, i]
  }

  # Change row/col names
  if (cflag) {
    temp          <- cnames[col1]
    cnames[col1]  <- cnames[col2]
    cnames[col2]  <- temp
    colnames(mat) <- cnames
  }
  if (rflag) {
    temp          <- rnames[col1]
    rnames[col1]  <- rnames[col2]
    rnames[col2]  <- temp
    rownames(mat) <- rnames
  }

  # Error check
  if (errorCheck) {
    for (i in rnames) {
      for (j in cnames) {
        if (mat[i, j] != mat0[i, j]) {
          stop("ERROR in swap2cols: with the error check")
        }
      }
    }
  } 

  mat

} # END: swap2cols.cov

# Function to perform test for heterogeneity (logistic regression only)
heterTest <- function(data, X.vars, group.var, snp.var, op=NULL) {

  # data       Data frame
  # X.vars     Variables to be adjusted for
  # group.var  Variable that defines the groups 
  # snp.var    Name of the SNP variable
  # op         List with names:
  #   print    0 or 1 to print model summaries
  #            The default is 0
  #   levels   The levels of group.var that are to be used.
  #            The default is NULL so that all levels of group.var
  #            will be used.

  op <- default.list(op, c("print"), list(0))
  print <- op$print

  levels <- getListName(op, "levels")
  if (is.null(levels)) levels <- unique(data[, group.var])
  nlevels <- length(levels)
  
  # Variables in the model
  X.vars <- unique(c(X.vars, snp.var))

  nTest <- 0
  minP  <- 9999
  for (i in 1:(nlevels - 1)) {
    for (j in (i+1):nlevels) {
      # Get the subgroups to be included in the model
      vec  <- c(levels[i], levels[j])
      temp <- data[, group.var] %in% vec
      temp[is.na(temp)] <- FALSE   
      dat2 <- data[temp, ]

      # Define the design matrix
      X <- dat2[, X.vars]
      
      # Define the response. Set one group to 1
      y <- rep.int(0, times=nrow(X))
      temp <- dat2[, group.var] == vec[1]
      y[temp] <- 1

      nTest <- nTest + 1
      fit <- try(glm(y ~ ., data=X, family="binomial"), silent=TRUE)
      if ("try-error" %in% class(fit)) next
      if (!fit$converged) next

      s <- summary(fit)
      if (print) {
        print(vec)
        print(table(dat2[, group.var]))
        print(s)
      }
      temp <- s$coefficients
      pval <- temp[snp.var, 4]
      minP <- min(minP, pval)    
    }
  }  

  minP <- min(1, minP*nTest)

  minP

} # heterTest

# Function to compute frequency counts for a vector of intervals.
# The returned object will be a data frame of frequency counts for
# the partitioned intervals or if data is passed in, then the returned
# object will be data with a new column.
freqCounts.var <- function(vec, intervals, leftEndClosed=1, data=NULL, 
                           newVar="newVar", newVarCats=NULL) {
 
  # vec              Numeric vector
  #                  No default
  # intervals        Numeric vector of intervals. Ex: c(0, 10, 20, 50, 200)
  #                  No default
  # leftEndClosed    0 or 1 Set to 1 for intervals [a , b).
  #                  Note: if a=b, then the interval is [a, a] and the next
  #                  interval will be (a, b) 
  # data             Data frame to be returned with the new variable newVar
  #                  on it with the disjoint categories.
  #                  The default is NULL
  # newVar           The name of the new variable if data is passed in
  #                  The default is "newVar".  
  # newVarCats       Vector of categories for newVar
  #                  The default is NULL

  intervals <- sort(intervals)
  intervals <- c(-Inf, intervals, Inf)
  n         <- length(intervals) - 1
  ret <- data.frame(rep.int(0, times=n+1))
  colnames(ret) <- "FREQ"
  if (leftEndClosed) {
    left  <- "["
    right <- ")"
    lop   <- ">="
    rop   <- "<"
  } else {
    left  <- "("
    right <- "]"
    lop   <- ">"
    rop   <- "<="
  }

  dFlag <- !is.null(data)
  if (dFlag) {
    data[, newVar] <- "MISSING"
    # Check for integers
    temp <- (vec == as.integer(vec))
    temp[is.na(temp)] <- TRUE
    if (all(temp)) {
      intFlag <- 1
    } else {
      intFlag <- 0
    }
    catFlag <- !is.null(newVarCats)
  } 

  rnames <- NULL
  flag   <- 0
  for (i in 1:n) {
    a <- intervals[i]
    b <- intervals[i+1]
    
    if (a == b) {
      text  <- paste("(vec == ", a, ")", sep="") 
      rtemp <- paste("[", a, ", ", b, "]", sep="")
      dstr  <- as.character(a)
      flag <- 1
    } else {
      if (flag) {
        # Previous had a == b, so left interval should be open to obtain disjoint sets
        text  <- paste("(vec", ">", a, ") & (vec", rop, b, ")", sep="")
        rtemp <- paste("(", a, ", ", b, right, sep="")
      } else {
        text  <- paste("(vec", lop, a, ") & (vec", rop, b, ")", sep="")
        rtemp <- paste(left, a, ", ", b, right, sep="")
      } 
      flag <- 0
      if (dFlag) {
        if (a == -Inf) {
          if (leftEndClosed) {
            dstr <- paste("lt", b, sep="")
          } else {
            dstr <- paste("lteq", b, sep="")
          }
        } else if (b == Inf) {
          if ((intFlag) & (!leftEndClosed)) a <- a + 1
          dstr <- paste(a, "plus", sep="")
        } else {
          if (intFlag) {
            if (!leftEndClosed) {
              a <- a + 1
            } else {
              b <- b - 1
            }
          }
          dstr <- paste(a, "to", b, sep="")
        }
      }
    }
    # Get the logical vector
    temp <- eval(parse(text=text))
    temp[is.na(temp)] <- FALSE
    ret[i, "FREQ"] <- sum(temp, na.rm=TRUE)
    rnames <- c(rnames, rtemp)

    if ( (dFlag) & (any(temp)) ){ 
      if (catFlag) dstr <- newVarCats[i]
      data[temp, newVar] <- dstr
    }
  }
  
  # Count the number of missing
  rnames <- c(rnames, "NA")
  ret[n+1, "FREQ"] <- sum(is.na(vec))
  rownames(ret) <- rnames

  if (dFlag) {
    print(ret)
    return(data)
  }
  ret

} # freqCounts.var

# Function to standardize a continuous vector
standardize.z <- function(vec) {

  mu <- mean(vec, na.rm=TRUE)
  se <- sqrt(var(vec, na.rm=TRUE))
  ret <- (vec - mu)/se
  ret 

} # END: standardize.z

# Function to create a design matrix
dsgnMat <- function(data, vars, facVars, removeInt=1) {

  # data        Data frame
  # vars        Character vector of variable names or a formula
  # facVars     Character vector of factor names

  if (is.null(vars)) return(list(designMatrix=NULL, newVars=NULL)) 

  # See if vars is a character string containing a formula
  if ((length(vars) == 1) && (substr(vars, 1, 1) == "~")) {
    vars <- as.formula(vars)
  }

  # Determine if vars is a formula
  if ("formula" %in% class(vars)) {
    # Get the design matrix
    design <- model.matrix(vars, data=data)

    # Remove the intercept, if needed
    newVars <- colnames(design)
    if (removeInt) {
      if (newVars[1] == "(Intercept)") {
        design  <- removeOrKeepCols(design, 1, which=-1)
        newVars <- newVars[-1]
      }
    }

    return(list(designMatrix=design, newVars=newVars))    
  }

  design  <- removeOrKeepCols(data, vars, which=1)
  newVars <- NULL
  if (!is.null(facVars)) {
    temp <- vars %in% facVars
    if (any(temp)) {
      temp    <- vars[temp]
      temp    <- createDummy(design, vars=temp)
      design  <- temp$data
      newVars <- temp$newVars
    }
  } 
  design <- as.matrix(design)

  # Check for constant variables
  design <- checkForConstantVar(design, msg=1)$data

  if (!removeInt) {
    # Add intercept
    cnames <- colnames(design)
    design <- cbind(1, design)
    colnames(design) <- c("Intercept", cnames)
  }

  # Make sure matrix is numeric
  d <- dim(design)
  cnames <- colnames(design)
  design <- as.numeric(design)
  dim(design) <- d
  colnames(design) <- cnames

  list(designMatrix=design, newVars=newVars)

} # END: dsgnMat

# Function to return genotype stats
getGenoStats <- function(vec, MAF=1, freqCounts=1) {

  # Vec must be numeric coded as 0-1-2 or NA

  if (MAF) {
    MAF2 <- getMAF(vec)
  } else {
    MAF2 <- NULL
  }
  if (freqCounts) {
    counts <- getGenoCounts(vec)
  } else {
    counts <- NULL
  }
  n.miss   <- sum(is.na(vec))
  len      <- length(vec)
  missRate <- n.miss/len
  n        <- len - n.miss

  list(MAF=MAF2, n.miss=n.miss, freqCounts=counts, missRate=missRate, n=n)

} # END: getGenoStats


##### 5/17/2012: just added "genofiles" ########

myExtractGenotype.gwas2 = function(myList,phefile,gwasID,yvar,genofiles){

     ########### [1.1] Make pheno.list, snp.list to be used for mergePhonoGeno ###########

     pheno.list = snp.list=NULL
     pheno.list = list(file=phefile, id.var=gwasID, cc.var=yvar)
     pheno.list

      #> pheno.list
      #$file
      #[1] "/data/home/hans4/ap/bladder/BC.pheno_1221_2010.txt"
      #
      #$id.var
      #[1] "GWAS_ID"
      #
      #$cc.var
      #[1] "CASECONTROL_CASE"

      #> myList[1:10,]
      #   datNum  num
      #1      34 2784
      #2      38  131
      #3      60 1080
      #4     196 2052
      #5     138 2910
      #6      66 2466

     for(j in 1:nrow(myList)){

         datNum=myList[j,"datNum"]
         num=myList[j,"num"]

         genofile=genofiles[datNum]
         genofile   #[1] "/data/gwas/lung/part/lung_and_PLCO2_part_2.ldat.gz"

          Start = num+1
          End = num+1
          #snp.list = list(file=genofile,start.vec=Start,stop.vec=1+End)
          snp.list = list(file=genofile,start.vec=Start,stop.vec=End)
          snp.list
          #> snp.list
          #$file
          #[1] "/data/gwas/bladder/part/bladder_build8_subjects_part_34.ldat.gz"
          #
          #$start.vec
          #[1] 2784
          #
          #$stop.vec
          #[1] 2785

          #########[1.2] Combine geno + pheno #################

          myDat0 = mergePhenoGeno(snp.list, pheno.list, op=list(which=0))
          ncol(myDat0)
          dim(myDat0)
          colnames(myDat0)
          # [1] "GWAS_ID"              "STUDY"                "CASECONTROL_CASE"
          # [4] "SEX_FEMALE"           "stratum_ASTURIAS"     "stratum_BARCELONA"
          # [7] "stratum_ELCHE"        "stratum_TENERIFE"     "stratum_VALLES"
          #[10] "stratum_MAINE"        "stratum_VERMONT"      "stratum_ATBC"
          #[13] "stratum_PLCO"         "AGE_CAT_BASE_lt50"    "AGE_CAT_BASE_50_54"
          #[16] "AGE_CAT_BASE_55_59"   "AGE_CAT_BASE_65_69"   "AGE_CAT_BASE_70_74"
          #[19] "AGE_CAT_BASE_75p"     "DNA_SOURCE_BUCCAL"    "cig_cat_CURRENT"
          #[22] "cig_cat_FORMER"       "cig_cat_MISSING"      "cig_cat_NEVER"
          #[25] "cig_cat_OCCASIONAL"   "cig_cat_NEVER_OCC"    "cigever_REGULAR"
          #[28] "PLCO.cigever_REGULAR" "RS1014971"            "RS11892031"
          #[31] "RS2294008"            "RS401681"             "RS710521"
          #[34] "RS798766"             "RS8102137"            "RS9642880"
          #[37] "RS1495741"            "CIG_CAT"              "PLCO.FORMER"
          #[40] "PLCO.CURRENT"         "SEX_MALE"             "AGE_CAT_BASE_60_65"
          #[43] "AGE"                  "AGE2"                 "rs11924987"
          #[46] "rs6797523"

          snp=myDat0[,ncol(myDat0)]
          snpName=colnames(myDat0)[ncol(myDat0)]

          if(j==1) { outdat=myDat0 }
          if(j>1) { outdat=cbind(outdat,snp) ; colnames(outdat)[ncol(outdat)]=snpName }

     }##nd of j

      col
      outdat


}#end of










### if it only has one level of covs, then remove it
myReduce=function(xx){   ### this gives the names of columns with only one level

      #xx=COVS
      #> colnames(xx)
      # [1] "stratum_ASTURIAS"     "stratum_BARCELONA"    "stratum_ELCHE"
      # [4] "stratum_TENERIFE"     "stratum_VALLES"       "stratum_MAINE"
      # [7] "stratum_VERMONT"      "stratum_ATBC"         "stratum_PLCO"
      #[10] "age_cat_base_50_54"   "age_cat_base_55_59"   "age_cat_base_60_64"
      #[13] "age_cat_base_65_69"   "age_cat_base_70_74"   "age_cat_base_75plus"
      #[16] "gender_FEMALE"        "dna_source_BUCCAL"    "PLCO.cig_cat_FORMER"
      #[19] "PLCO.cig_cat_CURRENT"
      #> table(xx[,19])
      #
      #   0
      #2012
      
      x=xx[,1] ;x[1:5]=NA
      mySmall=function(x) length(levels(as.factor(as.character(x))))
      
      indic.one = apply(xx,2,mySmall)==1
      
      cnames=colnames(xx)[indic.one] # names of columns witn only one level
      xx2=xx[,!indic.one]
      
      list(dat=xx2,cnames=cnames,indic=indic.one)
      

}#end of


  #### 5/2/2012: put Z score and others...######
  #### 5/6/2010 take out arguments for  Estimate Std. Error    z value     Pr(>|z|)
 
  #bNames="Estimate"
  #sdName="Std. Error"
  #pName="Pr(>|z|)"
  
  myOR.CI4=function(xx,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",zName="z value",pval=F){

        #> xx=full$coef
        #> xx
        #                Estimate Std. Error    z value     Pr(>|z|)
        #(Intercept) -0.750026027 0.14114317 -5.3139379 1.072812e-07
        #geno         0.168774184 0.04662121  3.6201160 2.944709e-04
        #race1        0.297560013 0.10316615  2.8842796 3.923103e-03
        #sex1        -0.073478696 0.07416867 -0.9906972 3.218334e-01
        #age_cat1    -0.113698897 0.10863623 -1.0466020 2.952832e-01
        #age_cat2    -0.138700275 0.10581389 -1.3107945 1.899272e-01
        #age_cat3    -0.309324107 0.10913423 -2.8343453 4.591968e-03
        #age_cat4    -0.263627342 0.11893800 -2.2165106 2.665655e-02
        #age_cat5    -0.263389923 0.23039374 -1.1432165 2.529487e-01
        #smoking1     0.135488073 0.23874607  0.5674986 5.703754e-01
        #smoking2    -0.002006010 0.08506604 -0.0235818 9.811862e-01
        #smoking3    -0.373251080 0.08498864 -4.3917760 1.124285e-05
        #bmi1         0.169979498 0.08502958  1.9990631 4.560153e-02
        #bmi2         0.303177546 0.09434804  3.2133953 1.311756e-03
        #bmi3         0.579843176 0.11796522  4.9153741 8.861307e-07
        #everhbp1     0.583651325 0.06945305  8.4035372 4.332274e-17

        #tm1=as.vector(xx[,"Estimate"])
        #tm2=as.vector(xx[,"Std. Error"])
        #tm3=as.vector(xx[,"Pr(>|z|)"])
        
        tm1=as.vector(xx[,bName])
        tm2=as.vector(xx[,sName])
        tm3=as.vector(xx[,pName])
        tm4=as.vector(xx[,zName])
                
        OR=exp(tm1)
        CI0=cbind(tm1-1.96*tm2,tm1+1.96*tm2)
        CI=round(exp(CI0),2)

        #> OR
        # [1] 0.4723543 1.1838528 1.3465692 0.9291559 0.8925267 0.8704889 0.7339429
        # [8] 0.7682598 0.7684422 1.1450955 0.9979960 0.6884923 1.1852806 1.3541549
        #[15] 1.7857584 1.7925718
        #> CI
        #           [,1]      [,2]
        # [1,] 0.3581990 0.6228900
        # [2,] 1.0804705 1.2971269
        # [3,] 1.1000486 1.6483350
        # [4,] 0.8034428 1.0745392
        # [5,] 0.7213535 1.1043182
        # [6,] 0.7074449 1.0711094
        # [7,] 0.5926050 0.9089902
        # [8,] 0.6085076 0.9699518
        # [9,] 0.4892109 1.2070530
        #[10,] 0.7171615 1.8283801
        #[11,] 0.8447323 1.1790670
        #[12,] 0.5828480 0.8132853
        #[13,] 1.0033270 1.4002314
        #[14,] 1.1255315 1.6292173
        #[15,] 1.4171267 2.2502808
        #[16,] 1.5644328 2.0539799
        
        ans=cbind(OR=round(OR,3),CI1=CI[,1],CI2=CI[,2],beta=tm1,sd=tm2,Z=round(tm4,3))
      	#ans=cbind(OR=OR,CI1=CI[,1],CI2=CI[,2],beta=tm1,sd=tm2,Z=tm4)
          
        
        if(pval==T){
        ans=cbind(OR=round(OR,2),CI1=CI[,1],CI2=CI[,2],beta=round(tm1,3),sd=round(tm2,3),pval2=tm3,Z=round(tm4,2))
        
        }
        
        
        ans

   }# end of myOR.CI

            myMAF=function(xx,indic.control){  # minAll
                  
                  xx2=xx[indic.control]
                  tb=table(xx2)
                  tb
                  #xx
                  #   0    1    2     aa    aA     AA
                  #3898 1659  169
                  N=2*length(xx2)

                  #maj = (tb[1]*2 + tb[2])/N
                  mino = (tb[3]*2 + tb[2])/N
                  as.vector(mino)

           }#3nd


runAssoc=function(Y,Xs,COVS0,PRINT=T){
 
           
           ############# [1] Remove variables that only have one level ########################
           
           if(is.null(COVS0)==F) COVS=myReduce(COVS0)$dat
           
           
           ############# [2] Base Model #########################################################
             
            if(is.null(COVS0)==F) mdat=data.frame(Y=Y,COVS)
            if(is.null(COVS0)==T) mdat=data.frame(Y=Y)
            
            lm.base0=glm(Y~.,data=mdat,family=binomial(link="logit"))
            lm.base1=summary(lm.base0)
            lm.base2=myOR.CI4(lm.base1$coef,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",zName="z value",pval=T)  
            row.names(lm.base2)= row.names(lm.base1$coef)        
           
           ############# [2] Run the analysis #################################################
           
           emp = rep(NA,9) ; names(emp) = c("pval","OR","CI1","CI2","beta","sd","MAF","geno.case","geno.control")
           
           for(j in 1:ncol(Xs)){
           
                X1=Xs[,j]
                snp=colnames(Xs)[j]
                
                if(is.null(COVS0)==F)  try(tm1 <-runAssoc.small(Y,X1,COVS) )
                if(is.null(COVS0)==T)  try(tm1 <-runAssoc.small(Y,X1,COVS=NULL) )
                
                
                if(is.null(tm1)==T )    tm1=emp  
                
                #> tm1
                #               pval                  OR                 CI1                 CI2 
                #"0.150535820987201"  "1.15113030553701" "0.950146944507789"  "1.39462741840651" 
                #               beta                  sd                 MAF           geno.case 
                # "0.14074433405041" "0.097899470392951" "0.114173228346457"      "998//254//18" 
                #       geno.control 
                #     "880//234//30" 

                 out=c(snp=snp,tm1)
                
                
                if(j==1) myout=out
                if(j>1) myout=rbind(myout,out)
                
                if(PRINT==T) print(j)
                    
          }#end of j       
          
          myout[,c("OR","CI1","CI2")]        
                        
          
          row.names(myout)=1:nrow(myout)
          OUT=list(pval=myout, baseModel=lm.base1, baseModel2=lm.base2)
          OUT          
                     
   
}# end of logistic.single.models


  runAssoc.small=function(Y,X1,COVS){


      if(is.null(COVS)==F) mdat=data.frame(Y=Y,X1=X1,COVS)
      if(is.null(COVS)==T) mdat=data.frame(Y=Y,X1=X1)
      
      
      indic.control=(Y==0)

      lm0=NULL
      try(lm0<-glm(Y~X1+.,data=mdat,family=binomial(link="logit")))

      if(is.null(lm0)==F){

          lm1=summary(lm0)
          tt1=myOR.CI4(lm1$coef,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",zName="z value",pval=T)
          tt2=tt1[row.names(lm1$coef)=="X1",c("OR","CI1","CI2","beta","sd","pval2")]
          names(tt2)[length(tt2)]="pval"

      }#end of

      if(is.null(lm0)==T) {tt2=c(rep(NA,6));names(tt2)= c("OR","CI1","CI2","beta","sd","pval") }



      ######## maf ###########################################

      maf=myMAF(xx=X1,indic.control=indic.control)
      tb=table(X1,Y)
      #> tb
      #   Y
      #X1    0   1
      #  0 726 649
      #  1 467 422
      #  2  77  74
      case=paste(tb[,1],collapse="//") #[1] "726.467.77"
      contr=paste(tb[,2],collapse="//")

      ans = c(tt2[6],tt2[-6],MAF=round(maf,2),geno.case=case,geno.control=contr)
      ans


  }#end of


 ########## 5/23/2013: this allows no X1 in the model
 

runAssoc.baseOnly=function(Y,COVS0){
 
           
           ############# [1] Remove variables that only have one level ########################
           
           COVS=myReduce(COVS0)$dat
           
           
           ############# [2] Base Model #########################################################
             
            mdat=data.frame(Y=Y,COVS)
            lm.base0=glm(Y~.,data=mdat,family=binomial(link="logit"))
            lm.base1=summary(lm.base0)
            lm.base2=myOR.CI4(lm.base1$coef,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",zName="z value",pval=T)  
            row.names(lm.base2)= row.names(lm.base1$coef)        
           
               OUT = list(OR=lm.base2, beta=lm.base1,fit=lm.base0)
               OUT          
                     
   
}# end of logistic.single.models

     
########## 5/10/2012: I changed argument and internal function
############ 5/13/2010: pull out pvalfiles argument #####################################




myStratAnalysis3=function(stratVar,stname,Y,Xs,COVS,myhead="run",doPlot=F,YLIM,CEX.LAB,CEX.AXIS,PCH,COLs,XLAB,YLAB,CEX,snp.name,cex.snp,writeFile=F,all.strat=T,strat=NULL){

          doThis=F
          if(doThis==T){
          
                  stratVar = dat0[,"smoking"]; stname="smoking"
                  doPlot=T;YLIM=c(0,6);CEX.LAB=1.5;CEX.AXIS=1.5;CEX.MAIN=2;PCH=16 ; CEX=2; cex.snp = 0.7
                  XLAB="";YLAB="-log10(p-value)";COLs=c("dark green","blue","orange","red");snp.name=T;writeFile=T          
                  all.strat=F;strat="1"
                  myhead="run"
          }#end of  if(doThis==T){
          
          ######## [0] making table ####
          
          tb=table(stratVar)
          tb
          #  1   2   3   4
          #592 590 722 706

          if(all.strat==T){   nm=names(tb)   }
          if(all.strat==F) { nm=strat }

          OUT=rep(list(list()),length(nm))
          names(OUT)=nm

          fileNames=paste(myhead,".",stname,".strat.",nm,".csv",sep="")
         
          ############## [1] run for each strata ####################################

          for(j in 1:length(nm)){
                
                
                out0=NULL
                indic = (stratVar==nm[j])  &   is.na(stratVar)==F                
                LEV = nm[j]                
                yy = Y[indic]
                xx = Xs[indic,]
                covs = COVS[indic,]
          
                try(out0<-runAssoc(Y=yy,Xs=xx,COVS=covs,PRINT=T))   
                OUT[[j]]=out0[[1]]
                 #> colnames(out0[[1]])
                # [1] "snp"          "pval"         "OR"           "CI1"          "CI2"         
                # [6] "beta"         "sd"           "MAF"          "geno.case"    "geno.control"

               ####### write an output #####
               
               if(writeFile==T) write.csv(out0[[1]],fileNames[j])

               if(doPlot==T){ 
  
                      xx=out0[[1]][,1:2]
                      MAIN=paste("Association Straritifed by ", stname, "\n Stratum = ", j,sep="") 
                      
                      if(snp.name==F){ 
                                          
                          plot(1:nrow(xx),-log10(as.numeric(xx[,"pval"])),pch=PCH,cex=CEX,col=COLs[j],ylim=YLIM,cex.lab=CEX.LAB,cex.axis=CEX.AXIS,cex.main=CEX.MAIN,xlab=XLAB,ylab=YLAB,main=MAIN)
                          abline(h=-log10(0.05/nrow(xx)),lty="dotted",col="red",lwd=2)
                          abline(h=c(0:10),lty="dotted",col="blue")
                                    
                      }# end of 
                      
                      if(snp.name==T){ 
                          
                          par(mar=c(7,4,4,3))
                          plot(1:nrow(xx),-log10(as.numeric(xx[,"pval"])),pch=PCH,axes=F,cex=CEX,col=COLs[j],ylim=YLIM,cex.lab=CEX.LAB,cex.axis=CEX.AXIS,cex.main=CEX.MAIN,xlab=XLAB,ylab=YLAB,main=MAIN)
                          axis(1,at=1:nrow(xx),labels=xx[,"snp"],col.axis="black",las=3,font.lab=3,cex.axis=cex.snp)
                          axis(2)
                          abline(h=-log10(0.05/nrow(xx)),lty="dotted",col="red",lwd=2)                      
                          abline(h=c(0:10),lty="dotted",col="blue")
                  
                      
                      }#end of 
               
                  }# end of doPlot
 
          }# end of j all strata end


          OUT


}# end of myStrat


		myFlipGeno  = function(x){
	
				# > table(x)
				# x
				#    0    1    2 
				# 7112 3818  501 
			
				x2=x
				x2[x==2]=0
				x2[x==0]=2
			
				table(x2)
				# x2
				#    0    1    2 
				#  501 3818 7112 
				
				x2				
	
		}#


myPRS = function(betas, GENO){  # x is genotype: so that risk allele is 

		###### weighted sum across all risk alleles ####
		
		# > betas
		#  [1] 0.213 0.144 0.109 0.183 0.132 0.180 0.203 0.139 0.235 0.123 0.147 0.179 0.229 0.182 0.128 0.175 0.313 0.148 0.191 0.237 0.157 0.273 0.290 0.206
		
		# > GENO
		#      rs12441998 rs7736354 rs2153903 rs11580716 rs11626844 rs2809964 rs2561537 rs10516962 rs9308366 rs2837900 rs1032554 rs7259175 rs898456 rs16897883
		# [1,]          2         2         1          2          0         2         1          1         0         0         1         2       NA          0
		# [2,]          2         2         1          2          1         2         2          0         1         0         2         2       NA          0
		# [3,]          2         2         0          2          2         0         1          1         2         2         2         1       NA          0
		# [4,]          2         1         2          2          1         1         2          1         2         2         1         2       NA          0
		# [5,]          2         2         0          2          0         0         1          2         1         2         2         2       NA          1
		#      rs1010294 rs7852051 rs8040562 rs1352889 rs7244453 rs3886594 rs11631676 rs9509598 rs8066241 rs2332637
		# [1,]         0         1         0         0         2         0          1        NA        NA         1
		# [2,]         0         1         0         0         2         1          0        NA        NA         1
		# [3,]         1         1         0         0         2         0          1        NA        NA         0
		# [4,]         0         1         0         0         2         1          0        NA        NA         1
		# [5,]         0         0         0         0         0         0          0        NA        NA         0
		
		
		
		###### Change NA's to 0 in genotype since multiplication won't work for NA ####
		GENO2=GENO
		GENO2[is.na(GENO2)]=0
		# > GENO2
		#      rs12441998 rs7736354 rs2153903 rs11580716 rs11626844 rs2809964 rs2561537 rs10516962 rs9308366 rs2837900 rs1032554 rs7259175 rs898456 rs16897883
		# [1,]          2         2         1          2          0         2         1          1         0         0         1         2        0          0
		# [2,]          2         2         1          2          1         2         2          0         1         0         2         2        0          0
		# [3,]          2         2         0          2          2         0         1          1         2         2         2         1        0          0
		# [4,]          2         1         2          2          1         1         2          1         2         2         1         2        0          0
		# [5,]          2         2         0          2          0         0         1          2         1         2         2         2        0          1
		#      rs1010294 rs7852051 rs8040562 rs1352889 rs7244453 rs3886594 rs11631676 rs9509598 rs8066241 rs2332637
		# [1,]         0         1         0         0         2         0          1         0         0         1
		# [2,]         0         1         0         0         2         1          0         0         0         1
		# [3,]         1         1         0         0         2         0          1         0         0         0
		# [4,]         0         1         0         0         2         1          0         0         0         1
		# [5,]         0         0         0         0         0         0          0         0         0         0
		
		
		##### matrix multiplication ####
		
		GENO2 %*% betas
		as.vector(GENO2 %*% betas)

}#3nd of 		




			myPredict.lm = function(mydat.tr, mydat.val,covnames){  ## do linear regresion
			

					#covnames=c("CIG_CAT_CURRENT","CIG_CAT_FORMER",
				#				"STUDY_CPSII","STUDY_EAGLE","STUDY_PLCO",    # study                                       # gender
				#				"CIG_CAT_CURRENT.PLCO", "CIG_CAT_FORMER.PLCO",
				#				"EAGLE_EV2","PLCO_EV4",	"PLCO_EV5",	"ATBC_EV2")


						# 					covnames=c("CIG_CAT_CURRENT","CIG_CAT_FORMER",
						# 								"GENDER_FEMALE",
						# 								"AGE_CAT_51to55","AGE_CAT_56to60","AGE_CAT_61to65","AGE_CAT_66to70","AGE_CAT_71to75", "AGE_CAT_75p", # age
						# 								"STUDY_CPSII","STUDY_EAGLE","STUDY_PLCO",    # study                                       # gender
						# 								"CIG_CAT_CURRENT.PLCO", "CIG_CAT_FORMER.PLCO",
						# 								"EAGLE_EV2","PLCO_EV4",	"PLCO_EV5",	"ATBC_EV2")

					############# (1) Prediction #########################################
					
					Y=mydat.tr[,"CASECONTROL_CODE"]
					
					COVS = mydat.tr[,covnames]

					PRS = mydat.tr[,"PRS"]
					# > summary(PRS)
					#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   1.907   3.690   4.324   4.330   4.976   6.542
		
					lm1=lm(PRS~Y+., data=COVS) 	
					summary(lm1)		
		
					# Coefficients: (5 not defined because of singularities)
					#                      Estimate Std. Error t value Pr(>|t|)    
					# (Intercept)           4.87523    0.03239 150.519  < 2e-16 ***
					# Y                     0.20289    0.01138  17.823  < 2e-16 ***
					# CIG_CAT_CURRENT      -0.10160    0.02022  -5.025 5.14e-07 ***
					# CIG_CAT_FORMER       -0.04981    0.01841  -2.705  0.00684 ** 
					# GENDER_FEMALE        -0.02480    0.01644  -1.508  0.13153    
					# AGE_CAT_51to55        0.02387    0.02831   0.843  0.39909    
					# AGE_CAT_56to60       -0.01015    0.02709  -0.375  0.70801    
					# AGE_CAT_61to65       -0.02764    0.02690  -1.027  0.30434    
					# AGE_CAT_66to70       -0.03399    0.02784  -1.221  0.22217    
					# AGE_CAT_71to75       -0.02071    0.02965  -0.699  0.48485    
					# AGE_CAT_75p          -0.02216    0.03399  -0.652  0.51451    
					# STUDY_CPSII           0.11374    0.02113   5.383 7.55e-08 ***
					# STUDY_EAGLE          -1.23239    0.01643 -75.017  < 2e-16 ***
					# STUDY_PLCO                 NA         NA      NA       NA    
					# CIG_CAT_CURRENT.PLCO       NA         NA      NA       NA    
					# CIG_CAT_FORMER.PLCO        NA         NA      NA       NA    
					# EAGLE_EV2             0.17676    0.48899   0.361  0.71775    
					# PLCO_EV4                   NA         NA      NA       NA    
					# PLCO_EV5                   NA         NA      NA       NA    
					# ATBC_EV2              0.62654    0.50622   1.238  0.21587    
					# ---
					# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
					# 
					# Residual standard error: 0.4882 on 8234 degrees of freedom
					# Multiple R-squared:  0.6338,	Adjusted R-squared:  0.6332 
					# F-statistic:  1018 on 14 and 8234 DF,  p-value: < 2.2e-16
	

					# > summary(lm1)$coef
					#                    Estimate Std. Error    t value     Pr(>|t|)
					# (Intercept)      4.86208003 0.02091532 232.465061 0.000000e+00
					# Y                0.19924122 0.01129007  17.647478 1.910839e-68
					# CIG_CAT_CURRENT -0.09260910 0.01976425  -4.685688 2.834937e-06
					# CIG_CAT_FORMER  -0.04389931 0.01781925  -2.463590 1.377563e-02
					# STUDY_CPSII      0.09472955 0.02016568   4.697563 2.675453e-06
					# STUDY_EAGLE     -1.24815791 0.01507733 -82.783757 0.000000e+00
	
					betas.lm = summary(lm1)$coef[,1]		
					# > betas.lm
					#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER     STUDY_CPSII     STUDY_EAGLE 
					#      4.86208003      0.19924122     -0.09260910     -0.04389931      0.09472955     -1.24815791 



					###############(2) Validation ###########################################

					Y2=mydat.val[,"CASECONTROL_CODE"]

					covnames2=c("CASECONTROL_CODE",covnames)

					
					COVS2 = as.matrix(mydat.val[,covnames2])

					DAT=COVS2

	
					myPred.small = function(betas.lm, DAT){
							# > betas.lm
							#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER 
							#      4.86208003      0.19924122     -0.09260910     -0.04389931 
							# > DAT[1:5,]
							#   CIG_CAT_CURRENT CIG_CAT_FORMER STUDY_CPSII STUDY_EAGLE
							# 1               1              0           0           0
							# 2               1              0           0           0
							# 3               1              0           0           0
							# 4               1              0           0           0
							# 5               1              0           0           0

							DAT2=cbind(int=rep(1,nrow(DAT)), DAT)
							DAT2[1:5,]
							#   int CIG_CAT_CURRENT CIG_CAT_FORMER STUDY_CPSII STUDY_EAGLE
							# 1   1               1              0           0           0
							# 2   1               1              0           0           0
							# 3   1               1              0           0           0
							# 4   1               1              0           0           0
							# 5   1               1              0           0           0
	
							prs = as.vector(DAT2 %*% betas.lm)
							#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
							#   4.769   4.818   4.862   4.896   5.017   5.061 
							
							prs
				  }# end of myPred			
							

				prs.val = myPred.small(betas.lm, DAT)
				summary(prs.val)
				
				prs.real = mydat.val[,"PRS"]
				summary(prs.real)
				#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
				#   0.940   3.297   3.682   3.816   4.194   6.497 
	
	
	
				plot(prs.real, prs.val, xlim=c(0,7), ylim=c(0,7),pch=16)
				abline(0,1,col="red")
	
					par(mfrow=c(1,2))
					med1=median(prs.real)
					med2=median(prs.val)
					hist(prs.real, main=paste("Polygenic Risk Score (Observed Data) \n median=", round(med1,2), sep=""),
											nclass=10, prob=T,col="blue", xlab="Polygenic Risk Score",ylim=c(0,0.5),xlim=c(0,7)) 		
					hist(prs, main=paste("Polygenic Risk Score (Prediction) \n median=", round(med2,2), sep=""),
											nclass=10, prob=T,col="red", xlab="Polygenic Risk Score",ylim=c(0,0.5),xlim=c(0,7)) 
								
					
					
				print("Predicted PRS")
				print(summary(prs.val))				

				print("Observed PRS")
				print(summary(prs.real))				
			
			
		}#end of myPredic
		
		
		
		
		
		
	
# model=NULL
# model$fitted.values=prob.Q4
# model$y = prs.real

myROC=function (model, line.col = "red", auc.coords=c(.3,.1), grid = TRUE, grid.col = "blue", MAIN="", ...) {
   
   ######### make firsttable #############
   
    firsttable <- table(model$fitted.values, model$y)
	# > firsttable[1:20,]
	#                      
	#                       0 1
	#   0.00925087949022891 1 0
	#   0.00939018249797395 1 0
	#   0.00951867933483535 1 0
	#   0.00958150474149071 0 1
	#   0.00969549844314346 1 0
	#   0.00975581349021637 1 0
	#   0.00982924684431823 1 0
	#   0.0098350787655727  1 0
	#   0.00987811678466243 1 0
	#   0.00991407036852422 1 0
	#   0.00993360034564743 1 0
	#   0.00995426534575047 1 0
	#   0.00997226487183328 1 0
	#   0.00997298235353573 0 1
	#   0.0100126426963955  1 0
	#   0.0100171727529577  1 0
	#   0.0100219959202232  1 0
	#   0.0100310273465087  1 0
	#   0.0100882808549596  1 0
	#   0.0100900013622442  1 0
   
 
	colnames(firsttable) <- c("Non-diseased", "Diseased")
	rownames(firsttable) <- substr(rownames(firsttable),  1, 6)
	# firsttable[1:5,]
	#         
	#          Non-diseased Diseased
	#   0.0092            1        0
	#   0.0093            1        0
	#   0.0095            1        0
	#   0.0095            0        1
	#   0.0096            1        0
   
      firsttable1 <- cbind(as.numeric(rownames(firsttable)), firsttable)
		# > firsttable1[1:10,]
		#               Non-diseased Diseased
		# 0.0092 0.0092            1        0
		# 0.0093 0.0093            1        0
		# 0.0095 0.0095            1        0
		# 0.0095 0.0095            0        1
		# 0.0096 0.0096            1        0
		# 0.0097 0.0097            1        0
		# 0.0098 0.0098            1        0
		# 0.0098 0.0098            1        0
		# 0.0098 0.0098            1        0
		# 0.0099 0.0099            1        0
		# 
   
   
    rownames(firsttable1) <- rep("", nrow(firsttable1))
    colnames(firsttable1)[1] <- "predicted.prob"
    firsttable <- firsttable1[, 2:3]
    
		# > firsttable[1:10,]
		#  Non-diseased Diseased
		#             1        0
		#             1        0
		#             1        0
		#             0        1
		#             1        0
		#             1        0
		#             1        0
		#             1        0
		#             1        0
		#             1        0
    
    
    
    
    
    
    secondtable <- firsttable
    colnames(secondtable) <- c("1-Specificity", "Sensitivity")
    
    
    for (i in 1:length(secondtable[, 1])) {
        
        
        secondtable[i, 1] <- (sum(firsttable[, 1]) - sum(firsttable[(1:i), 1]))/sum(firsttable[, 1])
        secondtable[i, 2] <- (sum(firsttable[, 2]) - sum(firsttable[(1:i),  2]))/sum(firsttable[, 2])
        rownames(secondtable)[i] <- paste(">", rownames(secondtable)[i])
    }
    
    
	#  > secondtable[1:10,]
	#    1-Specificity Sensitivity
	# >      0.9990826   1.0000000
	# >      0.9981651   1.0000000
	# >      0.9972477   1.0000000
	# >      0.9972477   0.9984051
	# >      0.9963303   0.9984051
	# >      0.9954128   0.9984051
	# >      0.9944954   0.9984051
	# >      0.9935780   0.9984051
   
    
    
    secondtable <- rbind((c(1, 1)), secondtable)
    colnames(secondtable) <- c("1-Specificity", "Sensitivity")
    #model.des <- deparse(logistic.model$formula)

    rownames(secondtable)=c("0",firsttable1[,1])



    auc <- 0
    for (i in 1:(nrow(secondtable) - 1)) {
        auc <- auc + (secondtable[i, 1] - secondtable[(i + 1), 
            1]) * 0.5 * (secondtable[i, 2] + secondtable[(i + 
            1), 2])
    }
    
    

            plot(secondtable[, 1], secondtable[, 2], xlab = "1-Specificity", 
                ylab = "Sensitivity", xlim = (c(0, 1)), ylim = (c(0, 
                  1)), asp = 1, col = line.col, type = "l",  main=MAIN,  cex.axis=1.5, cex.lab=1.5, cex.main=1.5,...)
         
            lines(x = c(0, 1), y = c(0, 1), lty = 2, col = "blue")
            if (grid) {
                abline(v = 0, lty = 2, col = grid.col)
                abline(v = 0.2, lty = 2, col = grid.col)
                abline(v = 0.4, lty = 2, col = grid.col)
                abline(v = 0.6, lty = 2, col = grid.col)
                abline(v = 0.8, lty = 2, col = grid.col)
                abline(v = 1, lty = 2, col = grid.col)
                abline(h = 0, lty = 2, col = grid.col)
                abline(h = 0.2, lty = 2, col = grid.col)
                abline(h = 0.4, lty = 2, col = grid.col)
                abline(h = 0.6, lty = 2, col = grid.col)
                abline(h = 0.8, lty = 2, col = grid.col)
                abline(h = 1, lty = 2, col = grid.col)
            }
            auclabel <- paste("AUC =", round(auc, 
                3))
            if (!is.null(auc.coords)) {
                text(x = auc.coords[1], y = auc.coords[2], pos = 4, 
                  labels = auclabel, cex=2,...)
            }
        
        
        
    
    list(auc = auc, predicted.table = firsttable1, 
        diagnostic.table = secondtable)
        
        
        
  }      

   
   
			
			
lroc=function(logistic.model, graph = TRUE, add = FALSE, Main = "", line.col = "red", auc.coords = c(.3,.1), grid = TRUE, grid.col = "blue", ...) {
    if (add) {
        #title <- FALSE
    }
    if (length(grep("cbind", names(model.frame(logistic.model)))) > 
        0) {
        firsttable1 <- cbind(logistic.model$fitted.values, model.frame(logistic.model)[, 
            1][, 2:1])
        firsttable1 <- firsttable1[order(firsttable1[, 1]), ]
    }
    else {
        if (length(grep("(weights)", names(model.frame(logistic.model)))) > 
            0) {
            firsttable <- xtabs(as.vector(model.frame(logistic.model)[, 
                ncol(model.frame(logistic.model))]) ~ logistic.model$fitted.values + 
                logistic.model$y)
        }
        else {
            firsttable <- table(logistic.model$fitted.values, 
                logistic.model$y)
        }
        colnames(firsttable) <- c("Non-diseased", "Diseased")
        rownames(firsttable) <- substr(rownames(firsttable), 
            1, 6)
        firsttable1 <- cbind(as.numeric(rownames(firsttable)), 
        firsttable)
    }
   
   
    rownames(firsttable1) <- rep("", nrow(firsttable1))
    colnames(firsttable1)[1] <- "predicted.prob"
    firsttable <- firsttable1[, 2:3]
    
    secondtable <- firsttable
    
    for (i in 1:length(secondtable[, 1])) {
        secondtable[i, 1] <- (sum(firsttable[, 1]) - sum(firsttable[(1:i), 
            1]))/sum(firsttable[, 1])
        secondtable[i, 2] <- (sum(firsttable[, 2]) - sum(firsttable[(1:i), 
            2]))/sum(firsttable[, 2])
        rownames(secondtable)[i] <- paste(">", rownames(secondtable)[i])
    }
    secondtable <- rbind((c(1, 1)), secondtable)
    colnames(secondtable) <- c("1-Specificity", "Sensitivity")
    model.des <- deparse(logistic.model$formula)
    auc <- 0
    for (i in 1:(nrow(secondtable) - 1)) {
        auc <- auc + (secondtable[i, 1] - secondtable[(i + 1), 
            1]) * 0.5 * (secondtable[i, 2] + secondtable[(i + 
            1), 2])
    }
    if (graph) {
        if (!add) {
            plot(secondtable[, 1], secondtable[, 2], xlab = "1-Specificity", 
                ylab = "Sensitivity", xlim = (c(0, 1)), ylim = (c(0, 
                  1)), asp = 1, col = line.col, type = "l", main=Main,...)
            #if (title) {
                #title(main = title)
            #}
            lines(x = c(0, 1), y = c(0, 1), lty = 2, col = "blue")
            if (grid) {
                abline(v = 0, lty = 2, col = grid.col)
                abline(v = 0.2, lty = 2, col = grid.col)
                abline(v = 0.4, lty = 2, col = grid.col)
                abline(v = 0.6, lty = 2, col = grid.col)
                abline(v = 0.8, lty = 2, col = grid.col)
                abline(v = 1, lty = 2, col = grid.col)
                abline(h = 0, lty = 2, col = grid.col)
                abline(h = 0.2, lty = 2, col = grid.col)
                abline(h = 0.4, lty = 2, col = grid.col)
                abline(h = 0.6, lty = 2, col = grid.col)
                abline(h = 0.8, lty = 2, col = grid.col)
                abline(h = 1, lty = 2, col = grid.col)
            }
            auclabel <- paste("AUC =", round(auc, 
                3))
            if (!is.null(auc.coords)) {
                text(x = auc.coords[1], y = auc.coords[2], pos = 4, 
                  labels = auclabel, cex=2, ...)
            }
        }
        else {
            lines(secondtable[, 1], secondtable[, 2], col = line.col, 
                ...)
        }
    }
    list(model.description = model.des, auc = auc, predicted.table = firsttable1, 
        diagnostic.table = secondtable)
}


######## polygenic risk score prediction model using logistic regression model after converting PRS to binary using Q4 threshold 

myPredict.Q4 = function(DAT,covnames,Q4,PERCENT,MYSEED,excludeD=F){  ## do linear regresion
			
					# > covnames
					#  [1] "CIG_CAT_CURRENT"      "CIG_CAT_FORMER"       "STUDY_CPSII"          "STUDY_EAGLE"          "STUDY_PLCO"           "CIG_CAT_CURRENT.PLCO"
					#  [7] "CIG_CAT_FORMER.PLCO"  "EAGLE_EV2"            "PLCO_EV4"             "PLCO_EV5"             "ATBC_EV2"            


					#########(0) Dividing training/ validation data ##########################
					mydat.tr=mydat.val=NULL
					
					set.seed(MYSEED)
					#PERCENT = 0.85  # train data percent
	
					rows.train = sample(1:nrow(DAT), floor(nrow(DAT)*PERCENT), replace=F)
					#rows.train = sample(1:nrow(mydat0), 8012, replace=F)
		
					# > rows.train[1:5]
					# [1]  8252 10024  8709 10140  5224
					# > length(rows.train)
					# [1] 8012
					rows.val = (1:nrow(DAT))[-rows.train] 
					# > length(rows.val)
					# [1] 3434
					# > rows.val[1:10]
					#  [1]  1  3  8 10 14 17 18 21 23 25
		
					#indic.train = mydat0[,"STUDY"]!="PLCO"
					#sum(indic.train)
					#indic.valid = !indic.train

					### training data set
					# > sum(indic.train)
					# [1] 8249

					### 
					# > sum(indic.valid)
					# [1] 3197


					mydat.tr = DAT[rows.train,]
					mydat.val = DAT[rows.val,]
					# > dim(mydat.tr)
					# [1] 8249  271
					# > dim(mydat.val)
					#[1] 3197  271		
		

					#Q4 = 4.53
					#Q4=3	



					############# (1) Prediction #########################################
					
					Y=mydat.tr[,"CASECONTROL_CODE"]
					
					COVS = mydat.tr[,covnames]

					PRS = mydat.tr[,"PRS"] > Q4
					# > summary(PRS)
					#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   1.907   3.690   4.324   4.330   4.976   6.542
					
					table(PRS)
		
					if(excludeD==F) lm1=glm(PRS~Y+., data=COVS,family=binomial(link="logit")) 	
						
					
					if(excludeD==T)	lm1=glm(PRS~., data=COVS,family=binomial(link="logit"))
					print(summary(lm1))
		
					# Coefficients: (5 not defined because of singularities)
					#                      Estimate Std. Error t value Pr(>|t|)    
					# (Intercept)           4.87523    0.03239 150.519  < 2e-16 ***
					# Y                     0.20289    0.01138  17.823  < 2e-16 ***
					# CIG_CAT_CURRENT      -0.10160    0.02022  -5.025 5.14e-07 ***
					# CIG_CAT_FORMER       -0.04981    0.01841  -2.705  0.00684 ** 
					# GENDER_FEMALE        -0.02480    0.01644  -1.508  0.13153    
					# AGE_CAT_51to55        0.02387    0.02831   0.843  0.39909    
					# AGE_CAT_56to60       -0.01015    0.02709  -0.375  0.70801    
					# AGE_CAT_61to65       -0.02764    0.02690  -1.027  0.30434    
					# AGE_CAT_66to70       -0.03399    0.02784  -1.221  0.22217    
					# AGE_CAT_71to75       -0.02071    0.02965  -0.699  0.48485    
					# AGE_CAT_75p          -0.02216    0.03399  -0.652  0.51451    
					# STUDY_CPSII           0.11374    0.02113   5.383 7.55e-08 ***
					# STUDY_EAGLE          -1.23239    0.01643 -75.017  < 2e-16 ***
					# STUDY_PLCO                 NA         NA      NA       NA    
					# CIG_CAT_CURRENT.PLCO       NA         NA      NA       NA    
					# CIG_CAT_FORMER.PLCO        NA         NA      NA       NA    
					# EAGLE_EV2             0.17676    0.48899   0.361  0.71775    
					# PLCO_EV4                   NA         NA      NA       NA    
					# PLCO_EV5                   NA         NA      NA       NA    
					# ATBC_EV2              0.62654    0.50622   1.238  0.21587    
					# ---
					# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
					# 
					# Residual standard error: 0.4882 on 8234 degrees of freedom
					# Multiple R-squared:  0.6338,	Adjusted R-squared:  0.6332 
					# F-statistic:  1018 on 14 and 8234 DF,  p-value: < 2.2e-16
	
					#library(epicalc)
					
					
					
					
					ro1=lroc(lm1, auc.coords=c(.3,.1),Main="ROC curve for training data")
					# > ro1[[3]][8000:8020,]
					#  predicted.prob Non-diseased Diseased
					#          0.8747            1        0
					#          0.8747            0        1
					#          0.8747            0        1
					#          0.8747            0        1
					#          0.8747            1        0
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            1        0
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
				
					# > ro1$auc
					# [1] 0.910922

					# > summary(lm1)$coef
					#                    Estimate Std. Error    t value     Pr(>|t|)
					# (Intercept)      4.86208003 0.02091532 232.465061 0.000000e+00
					# Y                0.19924122 0.01129007  17.647478 1.910839e-68
					# CIG_CAT_CURRENT -0.09260910 0.01976425  -4.685688 2.834937e-06
					# CIG_CAT_FORMER  -0.04389931 0.01781925  -2.463590 1.377563e-02
					# STUDY_CPSII      0.09472955 0.02016568   4.697563 2.675453e-06
					# STUDY_EAGLE     -1.24815791 0.01507733 -82.783757 0.000000e+00
	
					betas.lm = summary(lm1)$coef[,1]		
					betas.lm
					#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER     STUDY_CPSII     STUDY_EAGLE 
					#      4.86208003      0.19924122     -0.09260910     -0.04389931      0.09472955     -1.24815791 



					###############(2) Validation ###########################################

					#Y2=mydat.val[,"CASECONTROL_CODE"]

					prs.real = ( mydat.val[,"PRS"] > Q4)*1
					summary(prs.real)
					#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   0.940   3.297   3.682   3.816   4.194   6.497 


					covnames2=c("CASECONTROL_CODE",covnames)
					if(excludeD==T) covnames2=covnames
					covnames2

					
					COVS2 = as.matrix(mydat.val[,covnames2])

					DAT=COVS2

	
					myPred.small.Q4 = function(betas.lm, DAT){
							# > betas.lm
							#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER 
							#      4.86208003      0.19924122     -0.09260910     -0.04389931 
							# > DAT[1:5,]
							#   CIG_CAT_CURRENT CIG_CAT_FORMER STUDY_CPSII STUDY_EAGLE
							# 1               1              0           0           0
							# 2               1              0           0           0
							# 3               1              0           0           0
							# 4               1              0           0           0
							# 5               1              0           0           0

							DAT2=cbind(int=rep(1,nrow(DAT)), DAT)
							DAT2[1:5,]
							#   int CIG_CAT_CURRENT CIG_CAT_FORMER STUDY_CPSII STUDY_EAGLE
							# 1   1               1              0           0           0
							# 2   1               1              0           0           0
							# 3   1               1              0           0           0
							# 4   1               1              0           0           0
							# 5   1               1              0           0           0
							# > cbind(names(betas.lm),colnames(DAT2))
							#       [,1]                   [,2]                  
							#  [1,] "(Intercept)"          "int"                 
							#  [2,] "Y"                    "CASECONTROL_CODE"    
							#  [3,] "CIG_CAT_CURRENT"      "CIG_CAT_CURRENT"     
							#  [4,] "CIG_CAT_FORMER"       "CIG_CAT_FORMER"      
							#  [5,] "STUDY_EAGLE"          "STUDY_EAGLE"         
							#  [6,] "STUDY_ATBC"           "STUDY_ATBC"          
							#  [7,] "STUDY_PLCO"           "STUDY_PLCO"          
							#  [8,] "CIG_CAT_CURRENT.PLCO" "CIG_CAT_CURRENT.PLCO"
							#  [9,] "CIG_CAT_FORMER.PLCO"  "CIG_CAT_FORMER.PLCO" 

	
							product = as.vector(DAT2 %*% betas.lm)
							
							
							prob = exp(product)/(1+exp(product))
							
							prob
							
				  }# end of myPred			
							

				#### predicted probablilty ##
				
				prob.Q4 = myPred.small.Q4(betas.lm, DAT)
				# > prob.Q4[1:10]
				#  [1] 0.10255710 0.10929332 0.11118991 0.16942529 0.07785262 0.09702257 0.26473337 0.09733927 0.07771870 0.09509142
				
				model=NULL
				model$fitted.values=prob.Q4
				model$y = prs.real

				out2=myROC(model, line.col = "red", MAIN="ROC curve for validation data")
				# > out2[[1]]
				# [1] 0.9109638
				# > out2[[2]][1:10,]
				#  predicted.prob Non-diseased Diseased
				#          0.0092            1        0
				#          0.0093            1        0
				#          0.0095            1        0
				#          0.0095            0        1
				#          0.0096            1        0
				#          0.0097            1        0
				#          0.0098            1        0
				#          0.0098            1        0
				#          0.0098            1        0
				#          0.0099            1        0

				#          0.8852            0        1
				#          0.8852            0        1
				#          0.8855            0        1
				#          0.8863            0        1
				#          0.8874            0        1
				#          0.8887            0        1
				#          0.9126            6       62
				#          0.9369            0        6


				# > out2[[3]][1:10,]
				# > out2[[3]][1:10,]
				#        1-Specificity Sensitivity
				# 0          1.0000000   1.0000000
				# 0.0092     0.9990826   1.0000000
				# 0.0093     0.9981651   1.0000000
				# 0.0095     0.9972477   1.0000000
				# 0.0095     0.9972477   0.9984051
				# 0.0096     0.9963303   0.9984051
				# 0.0097     0.9954128   0.9984051
				# 0.0098     0.9944954   0.9984051
				# 0.0098     0.9935780   0.9984051
				# 0.0098     0.9926606   0.9984051

				# > colSums(out2[[2]])
				# predicted.prob   Non-diseased       Diseased 
				#       437.4157      1090.0000       627.0000 
				# 
				# > 626/627
				# [1] 0.9984051  --> first ensitivity drop by identifying 1 disease with given threshold
				
				#(1087:1090)/1090
				#[1] 0.9972477 0.9981651 0.9990826 1.0000000
			
			
	
	
	
				#plot(prs.real, prob.Q4, xlim=c(0,7), ylim=c(0,7),pch=16)
				#abline(0,1,col="red")
	
				#	par(mfrow=c(1,2))
			#		med1=median(prs.real)
			#		med2=median(prs.val)
			#		hist(prs.real, main=paste("Polygenic Risk Score (Observed Data) \n median=", round(med1,2), sep=""),
			#								nclass=10, prob=T,col="blue", xlab="Polygenic Risk Score",ylim=c(0,0.5),xlim=c(0,7)) 		
			#		hist(prs, main=paste("Polygenic Risk Score (Prediction) \n median=", round(med2,2), sep=""),
			#								nclass=10, prob=T,col="red", xlab="Polygenic Risk Score",ylim=c(0,0.5),xlim=c(0,7)) 
								
					
					
				#print("Predicted PRS")
				#print(summary(prs.val))				

				#print("Observed PRS")
				#print(summary(prs.real))
				
				ans=list(roc.train=ro1, roc.val=out2,lm.train=summary(lm1))
				ans				
			
			
		}#end of myPredic
			



					mySimGeno.Freq = function(D, snpnames, freq.case, freq.control,SEED ){
					
							# > snpnames
							#  [1] "rs12441998" "rs7736354"  "rs2153903"  "rs11580716" "rs11626844" "rs2809964"  "rs2561537"  "rs10516962" "rs9308366"  "rs2837900"  "rs1032554"  "rs7259175"  "rs898456"  
							# [14] "rs16897883" "rs1010294"  "rs7852051"  "rs8040562"  "rs1352889"  "rs7244453"  "rs3886594"  "rs11631676" "rs9509598"  "rs8066241"  "rs2332637" 

					        #freq.case
					        #freq.control
							#  [1] 0.80621710 0.74805716 0.62772625 0.94246678 0.68914515 0.63512158 0.51316119 0.40386062 0.76660817 0.84256706 0.92140887 0.89558787 0.39709200 0.10253196 0.26121835 0.24617699
							# [17] 0.09100025 0.17197293 0.76284783 0.13336676 0.16370018 0.45537729 0.48295312 0.19115066
							# > 					        freq.control
							#  [1] 0.7716878 0.7247079 0.6082525 0.9260502 0.6674124 0.6031568 0.3948546 0.3742232 0.7319165 0.8196619 0.9085260 0.8803132 0.2427293 0.0866269 0.2465822 0.2244594 0.0764355 0.1580910
							# [19] 0.7445936 0.1183197 0.1480239 0.2834949 0.3017649 0.1670395
					        
					        set.seed(SEED)
					        myseeds = ceiling(runif(2*length(snpnames), 1234,999999));myseeds
					        
					        for(j in 1:length(snpnames)){
					        
					        	snp=snpnames[j]
					        	freq1=freq.case[j]
					        	freq2=freq.control[j]

					        	### case genotype frequency ##
					        	freq.geno1 = c((1-freq1)^2,  2*freq1*(1-freq1) ,freq1^2)
					        	freq.geno1  #[1] 0.0634752 0.3769353 0.5595895

					        	### control genotype frequency ##
					        	freq.geno2 = c((1-freq2)^2,  2*freq2*(1-freq2) ,freq2^2)
					        	freq.geno2  #[1] 0.05212646 0.35237148 0.59550205
					        	
					        	
					        	xx1 = rep(NA, length(D))
					        	
					        	set.seed(myseeds[2*j-1]) # 1, 3, 5,....
						        xx1[D==1] = sample(c(0,1,2), size=sum(D==1), replace=T, prob=freq.geno1)  
						        
						        set.seed(myseeds[2*j])
						        xx1[D==0]= sample(c(0,1,2), size=sum(D==0), replace=T, prob=freq.geno2)
								# > table(xx1[D==1])/sum(table(xx1[D==1]))
								#          0          1          2 
								# 0.03422274 0.30684455 0.65893271 
								# 
								# > freq.geno1
								# [1] 0.03755181 0.31246218 0.64998601	
								
								if(j==1) geno.val = xx1
								if(j>1) geno.val = cbind(geno.val,xx1)				        	
					        	
					        	
					        }#end of 
					        
					        colnames(geno.val)=snpnames
							# > dim(geno.val)
							# [1] 3434   24
							# > geno.val[1:10,]
							#       rs12441998 rs7736354 rs2153903 rs11580716 rs11626844 rs2809964 rs2561537 rs10516962 rs9308366 rs2837900 rs1032554 rs7259175 rs898456 rs16897883 rs1010294 rs7852051 rs8040562
							#  [1,]          2         1         1          2          2         2         0          0         1         2         2         2        0          0         1         0         1
							#  [2,]          2         1         1          2          2         2         1          1         2         2         2         2        0          0         2         1         0
							#  [3,]          1         1         1          2          2         2         1          1         2         2         2         1        0          0         1         0         0
							#  [4,]          2         1         0          2          2         1         0          1         2         1         2         2        0          1         0         0         0
							#  [5,]          2         2         1          2          2         2         0          1         2         1         2         2        0          0         0         1         0
							#  [6,]          2         2         0          2          2         1         1          1         2         2         1         2        2          0         0         1         0
							#  [7,]          1         1         2          2          0         0         0          0         1         2         2         2        0          1         0         0         0
							#  [8,]          1         1         0          2          0         2         1          2         2         1         2         1        0          0         1         0         1
							#  [9,]          1         2         2          2          0         0         1          1         2         2         1         2        0          0         0         1         0
							
							geno.val
					
					}#en of mySimGeno.Freq.small

  
  
  

#### 9/3/2013: do not use Disease status for predicting PRS #################  
#### 6/27/2013: this includes fiting PRS for training set again (not using the one fitted using WHOLE data ####

myPredict.Q4_3 = function(DAT,covnames1,covnames2,PERCENT,MYSEED,excludeD=F){  ## do linear regresion
			
					# > covnames1
					# [1] "CIG_CAT_CURRENT"      "CIG_CAT_FORMER"       "STUDY_EAGLE"          "STUDY_ATBC"           "STUDY_PLCO"           "CIG_CAT_CURRENT.PLCO"
					# [7] "CIG_CAT_FORMER.PLCO" 
					# > covnames2
					# [1] "CIG_CAT_CURRENT" "CIG_CAT_FORMER"

					#########(0) Dividing training/ validation data ##########################
					
					mydat.tr=mydat.val=NULL
					
					set.seed(MYSEED)
					#PERCENT = 0.85  # train data percent
	
					rows.train = sample(1:nrow(DAT), floor(nrow(DAT)*PERCENT), replace=F)
					#rows.train = sample(1:nrow(mydat0), 8012, replace=F)
		
					# > rows.train[1:5]
					# [1]  8252 10024  8709 10140  5224
					# > length(rows.train)
					# [1] 8012
					rows.val = (1:nrow(DAT))[-rows.train] 
					# > length(rows.val)
					# [1] 3434
					# > rows.val[1:10]
					#  [1]  1  3  8 10 14 17 18 21 23 25
		
					#indic.train = mydat0[,"STUDY"]!="PLCO"
					#sum(indic.train)
					#indic.valid = !indic.train

					### training data set
					# > sum(indic.train)
					# [1] 8249

					### 
					# > sum(indic.valid)
					# [1] 3197


					mydat.tr = DAT[rows.train,]
					mydat.val = DAT[rows.val,]
					# > dim(mydat.tr)
					# [1] 8249  271
					# > dim(mydat.val)
					#[1] 3197  271		
		

					#Q4 = 4.53
					#Q4=3




					####### [1] For training data: Fit polygenic logistic regression for: D ~ SNP1 + SNP2 + ... + covariates..  #####################
									# (i) estimating beta1, beta2, beta3,...
									# (ii) for getting PRS (and Q4?) 
			
					#geno.tr = mydat.tr[,snpNames]


					Y=mydat.tr[,"CASECONTROL_CODE"]
					COVS2=mydat.tr[,c(snpNames,covnames1)]

					################# [1.1] run polygenic logistic regression:D ~ SNP1 + SNP2 + ... + covariates..  ###

				out.tr=runAssoc.baseOnly(Y,COVS0=COVS2)
					out.tr[[1]]
					#                                OR   CI1   CI2   beta        sd        pval2      Z
					# (Intercept)                  0.00  0.00  0.00 -8.146     0.625 8.126920e-39 -13.03
					# rs12441998                   1.22  1.07  1.39  0.197     0.068 3.828691e-03   2.89
					# rs7736354                    1.16  1.02  1.32  0.152     0.066 2.077233e-02   2.31
					# rs2153903                    1.10  0.98  1.23  0.093     0.058 1.115643e-01   1.59
					# rs11580716                   1.22  0.96  1.54  0.197     0.120 1.010578e-01   1.64
					# rs11626844                   1.16  1.03  1.30  0.145     0.058 1.315047e-02   2.48
					# rs2809964                    1.17  1.04  1.31  0.157     0.058 7.086593e-03   2.69
					# rs2561537                    1.21  1.08  1.35  0.192     0.057 6.688150e-04   3.40
					# rs10516962                   1.16  1.04  1.30  0.148     0.057 8.915793e-03   2.62
					# rs9308366                    1.29  1.14  1.47  0.256     0.065 7.898964e-05   3.95
					# rs2837900                    1.14  0.98  1.33  0.135     0.078 8.495908e-02   1.72
					# rs1032554                    1.07  0.87  1.32  0.066     0.106 5.366455e-01   0.62
					# rs7259175                    1.28  1.08  1.51  0.245     0.086 4.345771e-03   2.85
					# rs898456                     1.24  1.10  1.41  0.218     0.065 7.364902e-04   3.38
					# rs16897883                   1.17  0.98  1.39  0.154     0.090 8.819536e-02   1.70
					# rs1010294                    1.18  1.03  1.35  0.164     0.070 1.929365e-02   2.34

					###### [1.2] get betas ###

					betas= out.tr[[1]][snpNames,"beta"]
					betas
					# rs12441998  rs7736354  rs2153903 rs11580716 rs11626844  rs2809964  rs2561537 rs10516962  rs9308366  rs2837900  rs1032554  rs7259175   rs898456 rs16897883  rs1010294  rs7852051 
					#      0.197      0.152      0.093      0.197      0.145      0.157      0.192      0.148      0.256      0.135      0.066      0.245      0.218      0.154      0.164      0.146 
					#  rs8040562  rs1352889  rs7244453  rs3886594 rs11631676  rs9509598  rs8066241  rs2332637 
					#      0.281      0.081      0.165      0.302      0.152      0.187      0.266      0.214 

					###### [1.3] get polygenic risk score ###

					GENO = as.matrix(mydat.tr[,snpNames])
					# > all(colnames(GENO)==names(betas))
					# [1] TRUE					
					prs.tr=myPRS(betas, GENO)
					summary(prs.tr)
					#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   1.142   3.988   4.528   4.590   5.206   6.992 
					summary(prs.tr[Y==1])
					#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   2.392   4.220   4.820   4.819   5.424   6.992 

					summary(prs.tr[Y==0])
					#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   1.142   3.803   4.297   4.359   4.888   6.706 
					
					#out.tr$prs=prs.tr
					# > names(out.tr)
					# [1] "OR"   "beta" "prs"

					Q4 = as.vector(summary(prs.tr[Y==0])["3rd Qu."])
					Q4
					#[1] 4.888



					PRS.high.tr = (prs.tr >=Q4)*1
					table(PRS.high.tr)
					# PRS.high.tr
					#    0    1 
					# 5100 2912  					
 					
 					
 					################# [1.2] Prediction model for PRS ~ D + covs  ###############
					
					COVS=data.frame(mydat.tr[,covnames1])
					
					lm1=glm(PRS.high.tr ~., data=COVS,family=binomial(link="logit")) 
		
					#if(excludeD==F) lm1=glm(PRS.high.tr ~., data=COVS,family=binomial(link="logit")) 	
					#if(excludeD==T)	lm1=glm(PRS.high.tr~., data=COVS,family=binomial(link="logit"))
					print(summary(lm1))
					# > summary(lm1)
					# 
					# Call:
					# glm(formula = PRS.high.tr ~ Y + ., family = binomial(link = "logit"), 
					#     data = COVS)
					# 
					# Deviance Residuals: 
					#     Min       1Q   Median       3Q      Max  
					# -2.2305  -0.4858  -0.2897   0.6079   2.6654  
					# 
					# Coefficients:
					#                      Estimate Std. Error z value Pr(>|z|)    
					# (Intercept)           1.04481    0.12437   8.401  < 2e-16 ***
					# Y                     1.35602    0.06681  20.297  < 2e-16 ***
					# CIG_CAT_CURRENT      -0.35735    0.16600  -2.153 0.031343 *  
					# CIG_CAT_FORMER       -0.37292    0.14185  -2.629 0.008564 ** 
					# STUDY_EAGLE          -4.19494    0.12133 -34.573  < 2e-16 ***
					# STUDY_ATBC           -0.44871    0.15097  -2.972 0.002958 ** 
					# STUDY_PLCO           -3.68633    0.21703 -16.985  < 2e-16 ***
					# CIG_CAT_CURRENT.PLCO  0.91385    0.25835   3.537 0.000404 ***
					# CIG_CAT_FORMER.PLCO   0.93675    0.23861   3.926 8.64e-05 ***
					# ---
					# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
					# 
					# (Dispersion parameter for binomial family taken to be 1)
					# 
					#     Null deviance: 10501.8  on 8011  degrees of freedom
					# Residual deviance:  6563.7  on 8003  degrees of freedom
					# AIC: 6581.7
					# 
					# Number of Fisher Scoring iterations: 5	
					#library(epicalc)
					
					lms=summary(lm1)
					lm1.or=myOR.CI4(lms$coef,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",zName="z value",pval=T)  
            		row.names(lm1.or)= row.names(lms$coef)        
					lm1.or
					#                        OR  CI1  CI2   beta    sd         pval2      Z
					# (Intercept)          2.84 2.23 3.63  1.045 0.124  4.443340e-17   8.40
					# Y                    3.88 3.40 4.42  1.356 0.067  1.362615e-91  20.30
					# CIG_CAT_CURRENT      0.70 0.51 0.97 -0.357 0.166  3.134279e-02  -2.15
					# CIG_CAT_FORMER       0.69 0.52 0.91 -0.373 0.142  8.563590e-03  -2.63
					# STUDY_EAGLE          0.02 0.01 0.02 -4.195 0.121 6.342418e-262 -34.57
					# STUDY_ATBC           0.64 0.47 0.86 -0.449 0.151  2.957536e-03  -2.97
					# STUDY_PLCO           0.03 0.02 0.04 -3.686 0.217  1.053966e-64 -16.99
					# CIG_CAT_CURRENT.PLCO 2.49 1.50 4.14  0.914 0.258  4.042453e-04   3.54
					# CIG_CAT_FORMER.PLCO  2.55 1.60 4.07  0.937 0.239  8.644692e-05   3.93
          
					
					
			out.pred = list(or=lm1.or, lm=summary(lm1), Q4=Q4)#, prs=prs.tr,D=Y)
					
					
					ro1=lroc(lm1, auc.coords=c(.3,.1),Main="ROC curve for training data")
					# > ro1[[3]][8000:8020,]
					#  predicted.prob Non-diseased Diseased
					#          0.8747            1        0
					#          0.8747            0        1
					#          0.8747            0        1
					#          0.8747            0        1
					#          0.8747            1        0
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            1        0
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
				
					# > ro1$auc
					# [1] 0.910922

					# > summary(lm1)$coef
					#                    Estimate Std. Error    t value     Pr(>|t|)
					# (Intercept)      4.86208003 0.02091532 232.465061 0.000000e+00
					# Y                0.19924122 0.01129007  17.647478 1.910839e-68
					# CIG_CAT_CURRENT -0.09260910 0.01976425  -4.685688 2.834937e-06
					# CIG_CAT_FORMER  -0.04389931 0.01781925  -2.463590 1.377563e-02
					# STUDY_CPSII      0.09472955 0.02016568   4.697563 2.675453e-06
					# STUDY_EAGLE     -1.24815791 0.01507733 -82.783757 0.000000e+00
	
					betas.lm0 = summary(lm1)$coef[,1]		
					betas.lm0
					#          (Intercept)                    Y      CIG_CAT_CURRENT       CIG_CAT_FORMER          STUDY_EAGLE           STUDY_ATBC           STUDY_PLCO 
					#            1.0448082            1.3560196           -0.3573513           -0.3729172           -4.1949418           -0.4487150           -3.6863339 
					# CIG_CAT_CURRENT.PLCO  CIG_CAT_FORMER.PLCO 
					#            0.9138454            0.9367494 


					if(length(covnames1)==length(covnames2)) betas.lm=betas.lm0
					if(length(covnames1) > length(covnames2)) {
					
						# > covnames1
						# [1] "CIG_CAT_CURRENT"      "CIG_CAT_FORMER"       "STUDY_EAGLE"          "STUDY_ATBC"           "STUDY_PLCO"           "CIG_CAT_CURRENT.PLCO"
						# [7] "CIG_CAT_FORMER.PLCO" 
						# > covnames2
						# [1] "CIG_CAT_CURRENT" "CIG_CAT_FORMER" 
						# > names(betas.lm0)
						# [1] "(Intercept)"          "Y"                    "CIG_CAT_CURRENT"      "CIG_CAT_FORMER"       "STUDY_EAGLE"          "STUDY_ATBC"          
						# [7] "STUDY_PLCO"           "CIG_CAT_CURRENT.PLCO" "CIG_CAT_FORMER.PLCO" 
						
						betas.lm = betas.lm0[c("(Intercept)", "Y", covnames2)]
						betas.lm
						 #    (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER 
						#       1.0448082       1.3560196      -0.3573513      -0.3729172 
					
					
					}


					###############(2) Validation ###########################################

					#Y2=mydat.val[,"CASECONTROL_CODE"]
					GENO.val = as.matrix(mydat.val[,snpNames])
					prs.val.obs = myPRS(betas, GENO.val)
					
					#prs.real = ( mydat.val[,"PRS"] > Q4)*1
					#summary(prs.real)
					#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   0.940   3.297   3.682   3.816   4.194   6.497 

					PRS.high.val.obs = (prs.val.obs >=Q4)*1
					# > table(PRS.high.val.obs)
					# PRS.high.val.obs
					#    0    1 
					# 2221 1213 					
					
					#covnames2b=c("CASECONTROL_CODE",covnames2)
					#if(excludeD==T) covnames2b=covnames
					
					covnames2b=covnames2
					covnames2b

					
					COVS2 = as.matrix(mydat.val[,covnames2b])

					DAT2=COVS2

	
					myPred.small.Q4 = function(betas.lm, DAT){
							# > betas.lm
							#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER 
							#      4.86208003      0.19924122     -0.09260910     -0.04389931 
							# > DAT[1:5,]
							#   CIG_CAT_CURRENT CIG_CAT_FORMER STUDY_CPSII STUDY_EAGLE
							# 1               1              0           0           0
							# 2               1              0           0           0
							# 3               1              0           0           0
							# 4               1              0           0           0
							# 5               1              0           0           0

							DAT2=cbind(int=rep(1,nrow(DAT)), DAT)
							DAT2[1:5,]
							#   int CIG_CAT_CURRENT CIG_CAT_FORMER STUDY_CPSII STUDY_EAGLE
							# 1   1               1              0           0           0
							# 2   1               1              0           0           0
							# 3   1               1              0           0           0
							# 4   1               1              0           0           0
							# 5   1               1              0           0           0
							# > cbind(names(betas.lm),colnames(DAT2))
							#       [,1]                   [,2]                  
							#  [1,] "(Intercept)"          "int"                 
							#  [2,] "Y"                    "CASECONTROL_CODE"    
							#  [3,] "CIG_CAT_CURRENT"      "CIG_CAT_CURRENT"     
							#  [4,] "CIG_CAT_FORMER"       "CIG_CAT_FORMER"      
							#  [5,] "STUDY_EAGLE"          "STUDY_EAGLE"         
							#  [6,] "STUDY_ATBC"           "STUDY_ATBC"          
							#  [7,] "STUDY_PLCO"           "STUDY_PLCO"          
							#  [8,] "CIG_CAT_CURRENT.PLCO" "CIG_CAT_CURRENT.PLCO"
							#  [9,] "CIG_CAT_FORMER.PLCO"  "CIG_CAT_FORMER.PLCO" 

	
							product = as.vector(DAT2 %*% betas.lm)
							
							
							prob = exp(product)/(1+exp(product))
							
							prob
							
				  }# end of myPred			
							

				#### predicted probablilty ##
				
				prob.Q4 = myPred.small.Q4(betas.lm, DAT2)
				# > prob.Q4[1:10]
				#  [1] 0.10255710 0.10929332 0.11118991 0.16942529 0.07785262 0.09702257 0.26473337 0.09733927 0.07771870 0.09509142
				
				model=NULL
				model$fitted.values=prob.Q4
				model$y = PRS.high.val.obs

				out.val=myROC(model, line.col = "red", MAIN="ROC curve for validation data")
				# > out2[[1]]
				# [1] 0.9109638
				# > out2[[2]][1:10,]
				#  predicted.prob Non-diseased Diseased
				#          0.0092            1        0
				#          0.0093            1        0
				#          0.0095            1        0
				#          0.0095            0        1
				#          0.0096            1        0
				#          0.0097            1        0
				#          0.0098            1        0
				#          0.0098            1        0
				#          0.0098            1        0
				#          0.0099            1        0

				#          0.8852            0        1
				#          0.8852            0        1
				#          0.8855            0        1
				#          0.8863            0        1
				#          0.8874            0        1
				#          0.8887            0        1
				#          0.9126            6       62
				#          0.9369            0        6


				# > out2[[3]][1:10,]
				# > out2[[3]][1:10,]
				#        1-Specificity Sensitivity
				# 0          1.0000000   1.0000000
				# 0.0092     0.9990826   1.0000000
				# 0.0093     0.9981651   1.0000000
				# 0.0095     0.9972477   1.0000000
				# 0.0095     0.9972477   0.9984051
				# 0.0096     0.9963303   0.9984051
				# 0.0097     0.9954128   0.9984051
				# 0.0098     0.9944954   0.9984051
				# 0.0098     0.9935780   0.9984051
				# 0.0098     0.9926606   0.9984051

				# > colSums(out2[[2]])
				# predicted.prob   Non-diseased       Diseased 
				#       437.4157      1090.0000       627.0000 
				# 
				# > 626/627
				# [1] 0.9984051  --> first ensitivity drop by identifying 1 disease with given threshold
				
				#(1087:1090)/1090
				#[1] 0.9972477 0.9981651 0.9990826 1.0000000
			
			
	
	
	
				#plot(prs.real, prob.Q4, xlim=c(0,7), ylim=c(0,7),pch=16)
				#abline(0,1,col="red")
	
				#	par(mfrow=c(1,2))
			#		med1=median(prs.real)
			#		med2=median(prs.val)
			#		hist(prs.real, main=paste("Polygenic Risk Score (Observed Data) \n median=", round(med1,2), sep=""),
			#								nclass=10, prob=T,col="blue", xlab="Polygenic Risk Score",ylim=c(0,0.5),xlim=c(0,7)) 		
			#		hist(prs, main=paste("Polygenic Risk Score (Prediction) \n median=", round(med2,2), sep=""),
			#								nclass=10, prob=T,col="red", xlab="Polygenic Risk Score",ylim=c(0,0.5),xlim=c(0,7)) 
								
					
					
				#print("Predicted PRS")
				#print(summary(prs.val))				

				#print("Observed PRS")
				#print(summary(prs.real))
				
				ans=list(out.tr=out.tr, out.pred=out.pred,roc.train=ro1, roc.val=out.val)
				ans				
			
			
}#end of myPredic myPredict.Q4_3


#### 6/27/2013: this includes fiting PRS for training set again (not using the one fitted using WHOLE data ####

myPredict.Q4_2 = function(DAT,covnames1,covnames2,PERCENT,MYSEED,excludeD=F){  ## do linear regresion
			
					# > covnames1
					# [1] "CIG_CAT_CURRENT"      "CIG_CAT_FORMER"       "STUDY_EAGLE"          "STUDY_ATBC"           "STUDY_PLCO"           "CIG_CAT_CURRENT.PLCO"
					# [7] "CIG_CAT_FORMER.PLCO" 
					# > covnames2
					# [1] "CIG_CAT_CURRENT" "CIG_CAT_FORMER"

					#########(0) Dividing training/ validation data ##########################
					
					mydat.tr=mydat.val=NULL
					
					set.seed(MYSEED)
					#PERCENT = 0.85  # train data percent
	
					rows.train = sample(1:nrow(DAT), floor(nrow(DAT)*PERCENT), replace=F)
					#rows.train = sample(1:nrow(mydat0), 8012, replace=F)
		
					# > rows.train[1:5]
					# [1]  8252 10024  8709 10140  5224
					# > length(rows.train)
					# [1] 8012
					rows.val = (1:nrow(DAT))[-rows.train] 
					# > length(rows.val)
					# [1] 3434
					# > rows.val[1:10]
					#  [1]  1  3  8 10 14 17 18 21 23 25
		
					#indic.train = mydat0[,"STUDY"]!="PLCO"
					#sum(indic.train)
					#indic.valid = !indic.train

					### training data set
					# > sum(indic.train)
					# [1] 8249

					### 
					# > sum(indic.valid)
					# [1] 3197


					mydat.tr = DAT[rows.train,]
					mydat.val = DAT[rows.val,]
					# > dim(mydat.tr)
					# [1] 8249  271
					# > dim(mydat.val)
					#[1] 3197  271		
		

					#Q4 = 4.53
					#Q4=3




					####### [1] For training data: Fit polygenic logistic regression for: D ~ SNP1 + SNP2 + ... + covariates..  #####################
									# (i) estimating beta1, beta2, beta3,...
									# (ii) for getting PRS (and Q4?) 
			
					#geno.tr = mydat.tr[,snpNames]


					Y=mydat.tr[,"CASECONTROL_CODE"]
					COVS2=mydat.tr[,c(snpNames,covnames1)]

					################# [1.1] run polygenic logistic regression:D ~ SNP1 + SNP2 + ... + covariates..  ###

				out.tr=runAssoc.baseOnly(Y,COVS0=COVS2)
					out.tr[[1]]
					#                                OR   CI1   CI2   beta        sd        pval2      Z
					# (Intercept)                  0.00  0.00  0.00 -8.146     0.625 8.126920e-39 -13.03
					# rs12441998                   1.22  1.07  1.39  0.197     0.068 3.828691e-03   2.89
					# rs7736354                    1.16  1.02  1.32  0.152     0.066 2.077233e-02   2.31
					# rs2153903                    1.10  0.98  1.23  0.093     0.058 1.115643e-01   1.59
					# rs11580716                   1.22  0.96  1.54  0.197     0.120 1.010578e-01   1.64
					# rs11626844                   1.16  1.03  1.30  0.145     0.058 1.315047e-02   2.48
					# rs2809964                    1.17  1.04  1.31  0.157     0.058 7.086593e-03   2.69
					# rs2561537                    1.21  1.08  1.35  0.192     0.057 6.688150e-04   3.40
					# rs10516962                   1.16  1.04  1.30  0.148     0.057 8.915793e-03   2.62
					# rs9308366                    1.29  1.14  1.47  0.256     0.065 7.898964e-05   3.95
					# rs2837900                    1.14  0.98  1.33  0.135     0.078 8.495908e-02   1.72
					# rs1032554                    1.07  0.87  1.32  0.066     0.106 5.366455e-01   0.62
					# rs7259175                    1.28  1.08  1.51  0.245     0.086 4.345771e-03   2.85
					# rs898456                     1.24  1.10  1.41  0.218     0.065 7.364902e-04   3.38
					# rs16897883                   1.17  0.98  1.39  0.154     0.090 8.819536e-02   1.70
					# rs1010294                    1.18  1.03  1.35  0.164     0.070 1.929365e-02   2.34

					###### [1.2] get betas ###

					betas= out.tr[[1]][snpNames,"beta"]
					betas
					# rs12441998  rs7736354  rs2153903 rs11580716 rs11626844  rs2809964  rs2561537 rs10516962  rs9308366  rs2837900  rs1032554  rs7259175   rs898456 rs16897883  rs1010294  rs7852051 
					#      0.197      0.152      0.093      0.197      0.145      0.157      0.192      0.148      0.256      0.135      0.066      0.245      0.218      0.154      0.164      0.146 
					#  rs8040562  rs1352889  rs7244453  rs3886594 rs11631676  rs9509598  rs8066241  rs2332637 
					#      0.281      0.081      0.165      0.302      0.152      0.187      0.266      0.214 

					###### [1.3] get polygenic risk score ###

					GENO = as.matrix(mydat.tr[,snpNames])
					# > all(colnames(GENO)==names(betas))
					# [1] TRUE					
					prs.tr=myPRS(betas, GENO)
					summary(prs.tr)
					#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   1.142   3.988   4.528   4.590   5.206   6.992 
					summary(prs.tr[Y==1])
					#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   2.392   4.220   4.820   4.819   5.424   6.992 

					summary(prs.tr[Y==0])
					#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   1.142   3.803   4.297   4.359   4.888   6.706 
					
					#out.tr$prs=prs.tr
					# > names(out.tr)
					# [1] "OR"   "beta" "prs"

					Q4 = as.vector(summary(prs.tr[Y==0])["3rd Qu."])
					Q4
					#[1] 4.888



					PRS.high.tr = (prs.tr >=Q4)*1
					table(PRS.high.tr)
					# PRS.high.tr
					#    0    1 
					# 5100 2912  					
 					
 					
 					################# [1.2] Prediction model for PRS ~ D + covs  ###############
					
					COVS=data.frame(mydat.tr[,covnames1])
		
					if(excludeD==F) lm1=glm(PRS.high.tr ~ Y+., data=COVS,family=binomial(link="logit")) 	
					if(excludeD==T)	lm1=glm(PRS.high.tr~., data=COVS,family=binomial(link="logit"))
					print(summary(lm1))
					# > summary(lm1)
					# 
					# Call:
					# glm(formula = PRS.high.tr ~ Y + ., family = binomial(link = "logit"), 
					#     data = COVS)
					# 
					# Deviance Residuals: 
					#     Min       1Q   Median       3Q      Max  
					# -2.2305  -0.4858  -0.2897   0.6079   2.6654  
					# 
					# Coefficients:
					#                      Estimate Std. Error z value Pr(>|z|)    
					# (Intercept)           1.04481    0.12437   8.401  < 2e-16 ***
					# Y                     1.35602    0.06681  20.297  < 2e-16 ***
					# CIG_CAT_CURRENT      -0.35735    0.16600  -2.153 0.031343 *  
					# CIG_CAT_FORMER       -0.37292    0.14185  -2.629 0.008564 ** 
					# STUDY_EAGLE          -4.19494    0.12133 -34.573  < 2e-16 ***
					# STUDY_ATBC           -0.44871    0.15097  -2.972 0.002958 ** 
					# STUDY_PLCO           -3.68633    0.21703 -16.985  < 2e-16 ***
					# CIG_CAT_CURRENT.PLCO  0.91385    0.25835   3.537 0.000404 ***
					# CIG_CAT_FORMER.PLCO   0.93675    0.23861   3.926 8.64e-05 ***
					# ---
					# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
					# 
					# (Dispersion parameter for binomial family taken to be 1)
					# 
					#     Null deviance: 10501.8  on 8011  degrees of freedom
					# Residual deviance:  6563.7  on 8003  degrees of freedom
					# AIC: 6581.7
					# 
					# Number of Fisher Scoring iterations: 5	
					#library(epicalc)
					
					lms=summary(lm1)
					lm1.or=myOR.CI4(lms$coef,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",zName="z value",pval=T)  
            		row.names(lm1.or)= row.names(lms$coef)        
					lm1.or
					#                        OR  CI1  CI2   beta    sd         pval2      Z
					# (Intercept)          2.84 2.23 3.63  1.045 0.124  4.443340e-17   8.40
					# Y                    3.88 3.40 4.42  1.356 0.067  1.362615e-91  20.30
					# CIG_CAT_CURRENT      0.70 0.51 0.97 -0.357 0.166  3.134279e-02  -2.15
					# CIG_CAT_FORMER       0.69 0.52 0.91 -0.373 0.142  8.563590e-03  -2.63
					# STUDY_EAGLE          0.02 0.01 0.02 -4.195 0.121 6.342418e-262 -34.57
					# STUDY_ATBC           0.64 0.47 0.86 -0.449 0.151  2.957536e-03  -2.97
					# STUDY_PLCO           0.03 0.02 0.04 -3.686 0.217  1.053966e-64 -16.99
					# CIG_CAT_CURRENT.PLCO 2.49 1.50 4.14  0.914 0.258  4.042453e-04   3.54
					# CIG_CAT_FORMER.PLCO  2.55 1.60 4.07  0.937 0.239  8.644692e-05   3.93
          
					
					
			out.pred = list(or=lm1.or, lm=summary(lm1), Q4=Q4)#, prs=prs.tr,D=Y)
					
					
					ro1=lroc(lm1, auc.coords=c(.3,.1),Main="ROC curve for training data")
					# > ro1[[3]][8000:8020,]
					#  predicted.prob Non-diseased Diseased
					#          0.8747            1        0
					#          0.8747            0        1
					#          0.8747            0        1
					#          0.8747            0        1
					#          0.8747            1        0
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            1        0
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8748            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
					#          0.8749            0        1
				
					# > ro1$auc
					# [1] 0.910922

					# > summary(lm1)$coef
					#                    Estimate Std. Error    t value     Pr(>|t|)
					# (Intercept)      4.86208003 0.02091532 232.465061 0.000000e+00
					# Y                0.19924122 0.01129007  17.647478 1.910839e-68
					# CIG_CAT_CURRENT -0.09260910 0.01976425  -4.685688 2.834937e-06
					# CIG_CAT_FORMER  -0.04389931 0.01781925  -2.463590 1.377563e-02
					# STUDY_CPSII      0.09472955 0.02016568   4.697563 2.675453e-06
					# STUDY_EAGLE     -1.24815791 0.01507733 -82.783757 0.000000e+00
	
					betas.lm0 = summary(lm1)$coef[,1]		
					betas.lm0
					#          (Intercept)                    Y      CIG_CAT_CURRENT       CIG_CAT_FORMER          STUDY_EAGLE           STUDY_ATBC           STUDY_PLCO 
					#            1.0448082            1.3560196           -0.3573513           -0.3729172           -4.1949418           -0.4487150           -3.6863339 
					# CIG_CAT_CURRENT.PLCO  CIG_CAT_FORMER.PLCO 
					#            0.9138454            0.9367494 


					if(length(covnames1)==length(covnames2)) betas.lm=betas.lm0
					if(length(covnames1) > length(covnames2)) {
					
						# > covnames1
						# [1] "CIG_CAT_CURRENT"      "CIG_CAT_FORMER"       "STUDY_EAGLE"          "STUDY_ATBC"           "STUDY_PLCO"           "CIG_CAT_CURRENT.PLCO"
						# [7] "CIG_CAT_FORMER.PLCO" 
						# > covnames2
						# [1] "CIG_CAT_CURRENT" "CIG_CAT_FORMER" 
						# > names(betas.lm0)
						# [1] "(Intercept)"          "Y"                    "CIG_CAT_CURRENT"      "CIG_CAT_FORMER"       "STUDY_EAGLE"          "STUDY_ATBC"          
						# [7] "STUDY_PLCO"           "CIG_CAT_CURRENT.PLCO" "CIG_CAT_FORMER.PLCO" 
						
						betas.lm = betas.lm0[c("(Intercept)", "Y", covnames2)]
						betas.lm
						 #    (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER 
						#       1.0448082       1.3560196      -0.3573513      -0.3729172 
					
					
					}


					###############(2) Validation ###########################################

					#Y2=mydat.val[,"CASECONTROL_CODE"]
					GENO.val = as.matrix(mydat.val[,snpNames])
					prs.val.obs = myPRS(betas, GENO.val)
					
					#prs.real = ( mydat.val[,"PRS"] > Q4)*1
					#summary(prs.real)
					#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
					#   0.940   3.297   3.682   3.816   4.194   6.497 

					PRS.high.val.obs = (prs.val.obs >=Q4)*1
					# > table(PRS.high.val.obs)
					# PRS.high.val.obs
					#    0    1 
					# 2221 1213 					
					
					covnames2b=c("CASECONTROL_CODE",covnames2)
					if(excludeD==T) covnames2b=covnames
					covnames2b

					
					COVS2 = as.matrix(mydat.val[,covnames2b])

					DAT2=COVS2

	
					myPred.small.Q4 = function(betas.lm, DAT){
							# > betas.lm
							#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER 
							#      4.86208003      0.19924122     -0.09260910     -0.04389931 
							# > DAT[1:5,]
							#   CIG_CAT_CURRENT CIG_CAT_FORMER STUDY_CPSII STUDY_EAGLE
							# 1               1              0           0           0
							# 2               1              0           0           0
							# 3               1              0           0           0
							# 4               1              0           0           0
							# 5               1              0           0           0

							DAT2=cbind(int=rep(1,nrow(DAT)), DAT)
							DAT2[1:5,]
							#   int CIG_CAT_CURRENT CIG_CAT_FORMER STUDY_CPSII STUDY_EAGLE
							# 1   1               1              0           0           0
							# 2   1               1              0           0           0
							# 3   1               1              0           0           0
							# 4   1               1              0           0           0
							# 5   1               1              0           0           0
							# > cbind(names(betas.lm),colnames(DAT2))
							#       [,1]                   [,2]                  
							#  [1,] "(Intercept)"          "int"                 
							#  [2,] "Y"                    "CASECONTROL_CODE"    
							#  [3,] "CIG_CAT_CURRENT"      "CIG_CAT_CURRENT"     
							#  [4,] "CIG_CAT_FORMER"       "CIG_CAT_FORMER"      
							#  [5,] "STUDY_EAGLE"          "STUDY_EAGLE"         
							#  [6,] "STUDY_ATBC"           "STUDY_ATBC"          
							#  [7,] "STUDY_PLCO"           "STUDY_PLCO"          
							#  [8,] "CIG_CAT_CURRENT.PLCO" "CIG_CAT_CURRENT.PLCO"
							#  [9,] "CIG_CAT_FORMER.PLCO"  "CIG_CAT_FORMER.PLCO" 

	
							product = as.vector(DAT2 %*% betas.lm)
							
							
							prob = exp(product)/(1+exp(product))
							
							prob
							
				  }# end of myPred			
							

				#### predicted probablilty ##
				
				prob.Q4 = myPred.small.Q4(betas.lm, DAT2)
				# > prob.Q4[1:10]
				#  [1] 0.10255710 0.10929332 0.11118991 0.16942529 0.07785262 0.09702257 0.26473337 0.09733927 0.07771870 0.09509142
				
				model=NULL
				model$fitted.values=prob.Q4
				model$y = PRS.high.val.obs

				out.val=myROC(model, line.col = "red", MAIN="ROC curve for validation data")
				# > out2[[1]]
				# [1] 0.9109638
				# > out2[[2]][1:10,]
				#  predicted.prob Non-diseased Diseased
				#          0.0092            1        0
				#          0.0093            1        0
				#          0.0095            1        0
				#          0.0095            0        1
				#          0.0096            1        0
				#          0.0097            1        0
				#          0.0098            1        0
				#          0.0098            1        0
				#          0.0098            1        0
				#          0.0099            1        0

				#          0.8852            0        1
				#          0.8852            0        1
				#          0.8855            0        1
				#          0.8863            0        1
				#          0.8874            0        1
				#          0.8887            0        1
				#          0.9126            6       62
				#          0.9369            0        6


				# > out2[[3]][1:10,]
				# > out2[[3]][1:10,]
				#        1-Specificity Sensitivity
				# 0          1.0000000   1.0000000
				# 0.0092     0.9990826   1.0000000
				# 0.0093     0.9981651   1.0000000
				# 0.0095     0.9972477   1.0000000
				# 0.0095     0.9972477   0.9984051
				# 0.0096     0.9963303   0.9984051
				# 0.0097     0.9954128   0.9984051
				# 0.0098     0.9944954   0.9984051
				# 0.0098     0.9935780   0.9984051
				# 0.0098     0.9926606   0.9984051

				# > colSums(out2[[2]])
				# predicted.prob   Non-diseased       Diseased 
				#       437.4157      1090.0000       627.0000 
				# 
				# > 626/627
				# [1] 0.9984051  --> first ensitivity drop by identifying 1 disease with given threshold
				
				#(1087:1090)/1090
				#[1] 0.9972477 0.9981651 0.9990826 1.0000000
			
			
	
	
	
				#plot(prs.real, prob.Q4, xlim=c(0,7), ylim=c(0,7),pch=16)
				#abline(0,1,col="red")
	
				#	par(mfrow=c(1,2))
			#		med1=median(prs.real)
			#		med2=median(prs.val)
			#		hist(prs.real, main=paste("Polygenic Risk Score (Observed Data) \n median=", round(med1,2), sep=""),
			#								nclass=10, prob=T,col="blue", xlab="Polygenic Risk Score",ylim=c(0,0.5),xlim=c(0,7)) 		
			#		hist(prs, main=paste("Polygenic Risk Score (Prediction) \n median=", round(med2,2), sep=""),
			#								nclass=10, prob=T,col="red", xlab="Polygenic Risk Score",ylim=c(0,0.5),xlim=c(0,7)) 
								
					
					
				#print("Predicted PRS")
				#print(summary(prs.val))				

				#print("Observed PRS")
				#print(summary(prs.real))
				
				ans=list(out.tr=out.tr, out.pred=out.pred,roc.train=ro1, roc.val=out.val)
				ans				
			
			
}#end of myPredic myPredict.Q4_2
			




	
        myCateg.quantile = function(x,indic,qts){   # this makes only do

               X=as.numeric(as.character(x))


               ####### get quantiles from the selected group ####
               x2=X[indic]

               qts2=quantile(x2,qts,na.rm=T)
               qts2
               #> qts2
               #     25%      50%      75%
               #1473.148 1945.450 2531.065

               #x2.tmp = x2[is.na(x2)==F]

               if (all(c(0.25,0.5,0.75) %in% qts)==T){

                     indic1 = X <= qts2[1]
                     indic2 = qts2[1] < X & X <= qts2[2]
                     indic3 = qts2[2] < X & X <= qts2[3]
                     indic4 = qts2[3] < X

                     sum(indic1,na.rm=T)+sum(indic2,na.rm=T)+sum(indic3,na.rm=T)+sum(indic4,na.rm=T)
                     #>sum(indic1,na.rm=T)+sum(indic2,na.rm=T)+sum(indic3,na.rm=T)+sum(indic4,na.rm=T)
                     #[1] 1212
                     #> sum(!is.na(x2))
                     #[1] 1212

                    x.cat = X

                    x.cat[indic1] = 1
                    x.cat[indic2] = 2
                    x.cat[indic3] = 3
                    x.cat[indic4] = 4

                    #### overlap the values with NA <-- they used to be NA

                    x.cat[is.na(X)==T] = NA
                    table(x.cat)
                    #x.cat
                    #  1   2   3   4
                    #303 303 303 303

                    checkThis=F

                    if(checkThis==T){

                        cbind(X,x.cat)[1:20,]
                        #           x2 x.cat
                        # [1,] 1289.40     1
                        # [2,] 1089.99     1
                        # [3,] 1442.11     1
                        # [4,]      NA    NA
                        # [5,] 1487.22     2
                        # [6,] 2703.93     4
                        # [7,] 2256.53     3
                        # [8,] 1036.38     1
                        # [9,] 1484.15     2
                        #[10,] 1244.08     1
                        #[11,]  627.99     1
                        #[12,] 2670.87     4
                        #[13,] 2042.88     3
                        #[14,] 1059.90     1
                        #[15,] 1282.77     1
                        #[16,] 2103.57     3
                        #[17,] 1988.73     3
                        #[18,] 2262.47     3
                        #[19,] 1052.24     1
                        #[20,] 3473.69     4

                        plot(x.cat,X)
                        abline(h=qts2)
                        table(x.cat)
                        #x.cat
                        #  1   2   3   4
                        #303 303 303 303

                    }#end of checkThis

               }# end of if (all(qts ==c(0.25,0.5,0.75)==T){


               if (all(c(0.33,0.66) %in% qts)==T){

                     indic1 = X <= qts2[1]
                     indic2 = qts2[1] < X & X <= qts2[2]
                     indic3 = qts2[2] < X

                     sum(indic1,na.rm=T)+sum(indic2,na.rm=T)+sum(indic3,na.rm=T)
                     #>sum(indic1,na.rm=T)+sum(indic2,na.rm=T)+sum(indic3,na.rm=T)+sum(indic4,na.rm=T)
                     #[1] 1212
                     #> sum(!is.na(x2))
                     #[1] 1212

                    x.cat = X

                    x.cat[indic1] = 1
                    x.cat[indic2] = 2
                    x.cat[indic3] = 3
                    

                    #### overlap the values with NA <-- they used to be NA

                    x.cat[is.na(X)==T] = NA
                    table(x.cat)
                    #x.cat
                    #  1   2   3   4
                    #303 303 303 303

                    checkThis=F

                    if(checkThis==T){

                        cbind(X,x.cat)[1:20,]
                        #           x2 x.cat
                        # [1,] 1289.40     1
                        # [2,] 1089.99     1
                        # [3,] 1442.11     1
                        # [4,]      NA    NA
                        # [5,] 1487.22     2
                        # [6,] 2703.93     4
                        # [7,] 2256.53     3
                        # [8,] 1036.38     1
                        # [9,] 1484.15     2
                        #[10,] 1244.08     1
                        #[11,]  627.99     1
                        #[12,] 2670.87     4
                        #[13,] 2042.88     3
                        #[14,] 1059.90     1
                        #[15,] 1282.77     1
                        #[16,] 2103.57     3
                        #[17,] 1988.73     3
                        #[18,] 2262.47     3
                        #[19,] 1052.24     1
                        #[20,] 3473.69     4

                        plot(x.cat,X)
                        abline(h=qts2)
                        table(x.cat)
                        #x.cat
                        #  1   2   3   4
                        #303 303 303 303

                    }#end of checkThis

               }# end of if (all(qts ==c(0.25,0.5,0.75)==T){

                 x.cat


        }# end of myQuantile
        
#tt=myCateg.quantile(x,indic,qts)
#cbind(tt,x)[1:10]
#> cbind(tt,x)[1:10,]
#      tt  x
# [1,] "1" NA
# [2,] "1" "1289.4"
# [3,] "1" "1089.99"
# [4,] NA  "1442.11"
# [5,] "2" NA
# [6,] "4" "1487.22"
# [7,] "3" "2781.51"
# [8,] "1" NA
# [9,] "2" "2703.93"
#[10,] "1" "1489.75"
#> table(tt)
#tt
#  1   2   3   4
#652 644 639 648

  
  
#####  5/18/2010 it also create first reference dummy too.. it's your choice to have it or not
##### 4/26/2010: this deals with no level (it autoatically deals with it #######


myDummyVar3=function(mat,refer=F,SORT=T){  ### this will create dummy variables with given level
                       ######## the most common one is reference and will be omitted

        ans2=NULL
        
        ####### in case it's not matrix make it ######
        
        if(is.null(dim(mat))==T) {
        
           dim(mat)=c(length(mat),1)
           colnames(mat)="x"
           
        }# end of    

        
        
        for(u in 1:ncol(mat)){

              x=mat[,u]
              cname=colnames(mat)[u]
              cname
              
              
              #if(is.null(level)==T) {   # if no level is specified do it
                if(SORT==T){
                
                    tb=sort(table(x),decreasing=T)
                    tb
                    # 61 62 31 11 64 41 63 
                    #80 77 51 49 48 46 13
                    
                }# sort
                
                
                if(SORT==F){
                
                    tb=table(x)
                    tb
                    # 61 62 31 11 64 41 63 
                    #80 77 51 49 48 46 13
                                    
                
                }#sort
                
                level= names(tb)
                level
                #[1]  "61" "62" "31" "11" "64" "41" "63"
              
              #}# end of is.null
              
                #> level
              
              #[1] 0 1 2
              #mat
              #      [,1] [,2]
              # [1,]    0    0
              # [2,]    0    1
              # [3,]    0    2
              # [4,]    1    0
              # [5,]    1    1
              # [6,]    1    2
              # [7,]    2    0
              # [8,]    2    1
              # [9,]    2    2
              level2=level[-1] #[1] 1 2  --> skip dummy variable for first level
              if(refer==T) level2=level
              
              ans=matrix(NA,nrow=nrow(mat),ncol=(length(level2)))
              colnames(ans)=paste(cname,level2,sep=".")
              
              #> ans
              #      [,1] [,2] [,3] [,4]
              # [1,]   NA   NA   NA   NA
              # [2,]   NA   NA   NA   NA
              # [3,]   NA   NA   NA   NA

              #Start=seq(1,ncol(ans),by=length(level2)) #[1] 1 3


              #where=Start[u]:(Start[u]+length(level2)-1) #[1] 1 2

              for(m in 1:length(level2)){  # skip first level, wich is zero level "0"

                   #ans[,where[m]] = (x==level2[m])*1

                   ans[,m] = (x==level2[m])*1


              }# end of m loop
              
              #> cbind(x,ans)[40:52,]
              #       x x.62 x.31 x.11 x.64 x.41 x.63
              # [1,] 11    0    0    1    0    0    0
              # [2,] 11    0    0    1    0    0    0
              # [3,] 11    0    0    1    0    0    0
              # [4,] 11    0    0    1    0    0    0
              # [5,] 11    0    0    1    0    0    0
              # [6,] 11    0    0    1    0    0    0
              # [7,] 11    0    0    1    0    0    0
              # [8,] 11    0    0    1    0    0    0
              # [9,] 11    0    0    1    0    0    0
              #[10,] 11    0    0    1    0    0    0
              #[11,] 31    0    1    0    0    0    0
              #[12,] 31    0    1    0    0    0    0
              #[13,] 31    0    1    0    0    0    0


      if (u==1) ans2=ans
      if (u!=1) ans2=cbind(ans2,ans)
      

     }# end of u loop

      ans2

}# end of myDummyVar

#tt=myDummyVar2(mat,level=NULL)

#cbind(mat,tt)[40:60,]
#       x x.62 x.31 x.11 x.64 x.41 x.63
# [1,] 11    0    0    1    0    0    0
# [2,] 11    0    0    1    0    0    0
# [3,] 11    0    0    1    0    0    0
# [4,] 11    0    0    1    0    0    0
# [5,] 11    0    0    1    0    0    0
# [6,] 11    0    0    1    0    0    0
# [7,] 11    0    0    1    0    0    0
# [8,] 11    0    0    1    0    0    0
# [9,] 11    0    0    1    0    0    0
#[10,] 11    0    0    1    0    0    0
#[11,] 31    0    1    0    0    0    0
#[12,] 31    0    1    0    0    0    0
#[13,] 31    0    1    0    0    0    0
#[14,] 31    0    1    0    0    0    0
#[15,] 31    0    1    0    0    0    0
#[16,] 31    0    1    0    0    0    0
#[17,] 31    0    1    0    0    0    0
#[18,] 31    0    1    0    0    0    0
#[19,] 31    0    1    0    0    0    0
#[20,] 31    0    1    0    0    0    0
#[21,] 31    0    1    0    0    0    0


#> ttt=myDummyVar3(mat,refer=T,SORT=F)
#> ttt[1:5,]
#     x.1 x.2 x.3 x.4
#[1,]   0   0   0   1
#[2,]   0   0   0   1
#[3,]   0   0   0   1
#[4,]   0   0   0   1
#[5,]   0   0   0   1
#> ttt=myDummyVar3(mat,refer=T,SORT=T)
#> ttt[1:5,]
#     x.3 x.4 x.1 x.2
#[1,]   0   1   0   0
#[2,]   0   1   0   0
#[3,]   0   1   0   0
#[4,]   0   1   0   0
#[5,]   0   1   0   0
  


  myAssign.PRS.NLST = function(DATA){		
		
		
		#candx_cal_year
		#rand_year
		
		#CANDX_DAYS	Days from randomization to first diagnosis of lung cancer 
	
		#CANCYR  Study year associated with first confirmed lung cancer.
		#CAN_SCR	Result of screen associated with the first confirmed lung cancer diagnosis
		#CANC_RPT_LINK	Is the reported lung cancer linked to a positive screen?
		#CANC_FREE_DAYS	Days from randomization to date when participant was last known to be free from lung cancer
        # CIGSMOK:  Smoking status at baseline,   0: Former, 1:Current
		
		################ (1.0) screen detected LC or not #############################
		
		
		#CAN_SCR	Result of screen associated with the first confirmed lung cancer diagnosis
			# 0	No Cancer
			# 1	Positive Screen
			# 2	Negative Screen
			# 3	Missed Screen
			# 4	Post Screening
		

		#CANC_RPT_LINK	Is the reported lung cancer linked to a positive screen?
		#CANC_FREE_DAYS	Days from randomization to date when participant was last known to be free from lung cancer
		#DEATH_DAYS:      Days since randomization at death from most definitive source
		#DCFDEATHLC:	Death certificate lung cancer death
		#FINALDEATHLC:	Final lung cancer death (Combined best information: EVP supplemented with DCF)
        
        ##############[1.1] Get covariates to be used for prediction: disease status, smoking status, and gender #####
        
		# > colnames(DATA)
		#   [1] "center"             "cigssmok"           "gender"             "rndgroup"          
		#   [5] "age"                "pkyr"               "scr_lat0"           "scr_lat1"          
		#   [9] "scr_lat2"           "ELIG"               "loclhil"            "locllow"           
		#  [13] "loclup"             "locrhil"            "locrlow"            "locrmid"           
		#  [17] "locunk"             "locrup"             "locoth"             "locmed"            
		#  [21] "loclmsb"            "locrmsb"            "loccar"             "loclin"            
		#  [25] "de_type"            "de_grade"           "de_stag"            "lesionsize"        
		#  [29] "scr_res0"           "scr_res1"           "scr_res2"           "scr_iso0"          
		#  [33] "scr_iso1"           "scr_iso2"           "RACE"               "ETHNIC"            
		#  [37] "EDUCAT"             "proc0"              "proc1"              "proc2"             
		#  [41] "biop0"              "biop1"              "biop2"              "invas0"            
		#  [45] "invas1"             "invas2"             "proclc"             "bioplc"            
		#  [49] "invaslc"            "medcomp0"           "medcomp1"           "medcomp2"          
		#  [53] "medcomplc"          "cancyr"             "can_scr"            "DCFICD"            
		#  [57] "evp_revr"           "evp_direct"         "evpincomplete"      "deathstat"         
		#  [61] "evpsel"             "evpcert"            "treatlc"            "mra_stat0"         
		#  [65] "mra_stat1"          "mra_stat2"          "no_proc_reas0"      "no_proc_reas1"     
		#  [69] "no_proc_reas2"      "study"              "conflc"             "evpdeath"          
		#  [73] "canc_rpt_link"      "canc_rpt_source"    "WDLOST"             "CONTACTSTATUS"     
		#  [77] "pid"                "smokeage"           "age_quit"           "smokeday"          
		#  [81] "smokeyr"            "dcfdeathlc"         "finaldeathlc"       "ndicd"             
		#  [85] "evpsent"            "deathcutoff"        "hist_group"         "rand_year"         
		#  [89] "scr_cal_year0"      "scr_days0"          "scr_age0"           "scr_cal_year1"     
		#  [93] "scr_days1"          "scr_age1"           "scr_cal_year2"      "scr_days2"         
		#  [97] "scr_age2"           "candx_cal_year"     "candx_days"         "candx_age"         
		# [101] "fup_cal_year"       "fup_days"           "fup_age"            "death_cal_year"    
		# [105] "death_days"         "death_age"          "canc_free_cal_year" "canc_free_days"    
		# [109] "canc_free_age"     
		
		
		### (1) Lung cancer status ###

		
			#tt1=as.numeric(DATA[,"cancyr"])
			# > table(tt1)
			# tt1
			#   0   1   2   3   4   5   6   7 
			# 176 111 146  47  67  65  22   1 
			# > sum(tt1[is.na(tt1)==F])
			# [1] 1276		
			#CANDX_DAYS
			#tt2=as.numeric(DATA[,"candx_days"])
			#  length(tt2[!is.na(tt2)])
			# [1] 635
			# > length(tt1[!is.na(tt1)])
			# [1] 635
		
		indic.LC =  (is.na(as.numeric(DATA[,"cancyr"]))==F)*1
			sum(indic.LC)
			#[1] 635
			
	   ### (2) smoking status: current, former #######		
		
			# > table(DATA[,"cigssmok"])
			# 
			#    0    1 
			# 8427 7343 
			
	   indic.former = 	(DATA[,"cigssmok"]==0)*1
	   indic.current = 	(DATA[,"cigssmok"]==1)*1
			# > sum(indic.former)
			# [1] 8427
			# > sum(indic.current)
			# [1] 7343	   
	   
	   	
	   ### (3) gender ##########
			# > table(DATA[,"gender"])
			# 
			#     1 
			# 15770 	   
	   indic.female = (DATA[,"gender"]==2)*1	


 


      ### (4) make a matrix ###########
		# > betas
		#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER   GENDER_FEMALE 
		#      -1.2165470       0.9561218       0.4949310      -0.2146295      -0.2304710 

      
      covmat = cbind(inter=rep(1,nrow(DATA)),  Y=indic.LC, Current=indic.current, Former=indic.former, Gender.female=indic.female)
			# > covmat[1:10,]
			#       inter Y Current Former Gender.female
			#  [1,]     1 0       1      0             0
			#  [2,]     1 0       0      1             0
			#  [3,]     1 0       0      1             0
			#  [4,]     1 0       1      0             0
			#  [5,]     1 0       0      1             0
			#  [6,]     1 0       1      0             0
			#  [7,]     1 0       0      1             0
			#  [8,]     1 0       0      1             0
			#  [9,]     1 0       1      0             0
			# [10,]     1 0       1      0             0
     
     ####(5) Get probability for highPRS ###
	
	#set.seed(SEED[1])
	
	product = as.vector(covmat %*% betas.lm)
	prob = matrix(exp(product)/(1+exp(product)),ncol=1)
		# > prob[1:5,]
		# [1] 0.3270372 0.1929154 0.1929154 0.3270372 0.1929154
	
	#### (6) simulate highPRS given probability ###########
	x=prob[1,] #[1] 0.327037
	pred.small = function(x){   sample(c(0,1), 1, prob=c(1-x, x))  }
	
	
	#set.seed(SEED[2])
	      
	highPRS = apply(prob, 1, pred.small)
		# > table(highPRS)
		# highPRS
		#     0     1 
		# 11560  4210 
		
		
	########### add LC death info ##########
	
	death = as.numeric(DATA[,"finaldeathlc"]) 
		#> death[1:10]
		# [1] NA NA NA NA NA NA NA NA  0 NA

		table(death)
		# lcdeath2
		#    0    1 
		# 1060  334		
    indic.LC.death = (death==1)*1
    indic.OCM.death = (death==0)*1
		# > sum(indic.LC.death,na.rm=T);sum(indic.OCM.death,na.rm=T)
		# [1] 302
		# [1] 1012
		
		
		
	DATA.big = cbind(DATA,covmat,highPRS=highPRS, indic.LC.death=indic.LC.death, indic.OCM.death = indic.OCM.death)
	
	DATA.big	


}#end of 



		
  myAssign.PRS.NLST = function(DATA){		
		
		
		#candx_cal_year
		#rand_year
		
		#CANDX_DAYS	Days from randomization to first diagnosis of lung cancer 
	
		#CANCYR  Study year associated with first confirmed lung cancer.
		#CAN_SCR	Result of screen associated with the first confirmed lung cancer diagnosis
		#CANC_RPT_LINK	Is the reported lung cancer linked to a positive screen?
		#CANC_FREE_DAYS	Days from randomization to date when participant was last known to be free from lung cancer
        # CIGSMOK:  Smoking status at baseline,   0: Former, 1:Current
		
		################ (1.0) screen detected LC or not #############################
		
		
		#CAN_SCR	Result of screen associated with the first confirmed lung cancer diagnosis
			# 0	No Cancer
			# 1	Positive Screen
			# 2	Negative Screen
			# 3	Missed Screen
			# 4	Post Screening
		

		#CANC_RPT_LINK	Is the reported lung cancer linked to a positive screen?
		#CANC_FREE_DAYS	Days from randomization to date when participant was last known to be free from lung cancer
		#DEATH_DAYS:      Days since randomization at death from most definitive source
		#DCFDEATHLC:	Death certificate lung cancer death
		#FINALDEATHLC:	Final lung cancer death (Combined best information: EVP supplemented with DCF)
        
        ##############[1.1] Get covariates to be used for prediction: disease status, smoking status, and gender #####
        
		# > colnames(DATA)
		#   [1] "center"             "cigssmok"           "gender"             "rndgroup"          
		#   [5] "age"                "pkyr"               "scr_lat0"           "scr_lat1"          
		#   [9] "scr_lat2"           "ELIG"               "loclhil"            "locllow"           
		#  [13] "loclup"             "locrhil"            "locrlow"            "locrmid"           
		#  [17] "locunk"             "locrup"             "locoth"             "locmed"            
		#  [21] "loclmsb"            "locrmsb"            "loccar"             "loclin"            
		#  [25] "de_type"            "de_grade"           "de_stag"            "lesionsize"        
		#  [29] "scr_res0"           "scr_res1"           "scr_res2"           "scr_iso0"          
		#  [33] "scr_iso1"           "scr_iso2"           "RACE"               "ETHNIC"            
		#  [37] "EDUCAT"             "proc0"              "proc1"              "proc2"             
		#  [41] "biop0"              "biop1"              "biop2"              "invas0"            
		#  [45] "invas1"             "invas2"             "proclc"             "bioplc"            
		#  [49] "invaslc"            "medcomp0"           "medcomp1"           "medcomp2"          
		#  [53] "medcomplc"          "cancyr"             "can_scr"            "DCFICD"            
		#  [57] "evp_revr"           "evp_direct"         "evpincomplete"      "deathstat"         
		#  [61] "evpsel"             "evpcert"            "treatlc"            "mra_stat0"         
		#  [65] "mra_stat1"          "mra_stat2"          "no_proc_reas0"      "no_proc_reas1"     
		#  [69] "no_proc_reas2"      "study"              "conflc"             "evpdeath"          
		#  [73] "canc_rpt_link"      "canc_rpt_source"    "WDLOST"             "CONTACTSTATUS"     
		#  [77] "pid"                "smokeage"           "age_quit"           "smokeday"          
		#  [81] "smokeyr"            "dcfdeathlc"         "finaldeathlc"       "ndicd"             
		#  [85] "evpsent"            "deathcutoff"        "hist_group"         "rand_year"         
		#  [89] "scr_cal_year0"      "scr_days0"          "scr_age0"           "scr_cal_year1"     
		#  [93] "scr_days1"          "scr_age1"           "scr_cal_year2"      "scr_days2"         
		#  [97] "scr_age2"           "candx_cal_year"     "candx_days"         "candx_age"         
		# [101] "fup_cal_year"       "fup_days"           "fup_age"            "death_cal_year"    
		# [105] "death_days"         "death_age"          "canc_free_cal_year" "canc_free_days"    
		# [109] "canc_free_age"     
		
		
		### (1) Lung cancer status ###

		
			#tt1=as.numeric(DATA[,"cancyr"])
			# > table(tt1)
			# tt1
			#   0   1   2   3   4   5   6   7 
			# 176 111 146  47  67  65  22   1 
			# > sum(tt1[is.na(tt1)==F])
			# [1] 1276		
			#CANDX_DAYS
			#tt2=as.numeric(DATA[,"candx_days"])
			#  length(tt2[!is.na(tt2)])
			# [1] 635
			# > length(tt1[!is.na(tt1)])
			# [1] 635
		
		indic.LC =  (is.na(as.numeric(DATA[,"cancyr"]))==F)*1
			sum(indic.LC)
			#[1] 635
			
	   ### (2) smoking status: current, former #######		
		
			# > table(DATA[,"cigssmok"])
			# 
			#    0    1 
			# 8427 7343 
			
	   indic.former = 	(DATA[,"cigssmok"]==0)*1
	   indic.current = 	(DATA[,"cigssmok"]==1)*1
			# > sum(indic.former)
			# [1] 8427
			# > sum(indic.current)
			# [1] 7343	   
	   
	   	
	   ### (3) gender ##########
			# > table(DATA[,"gender"])
			# 
			#     1 
			# 15770 	   
	   indic.female = (DATA[,"gender"]==2)*1	


 


      ### (4) make a matrix ###########
		# > betas
		#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER   GENDER_FEMALE 
		#      -1.2165470       0.9561218       0.4949310      -0.2146295      -0.2304710 

      
      covmat = cbind(inter=rep(1,nrow(DATA)),  Y=indic.LC, Current=indic.current, Former=indic.former, Gender.female=indic.female)
			# > covmat[1:10,]
			#       inter Y Current Former Gender.female
			#  [1,]     1 0       1      0             0
			#  [2,]     1 0       0      1             0
			#  [3,]     1 0       0      1             0
			#  [4,]     1 0       1      0             0
			#  [5,]     1 0       0      1             0
			#  [6,]     1 0       1      0             0
			#  [7,]     1 0       0      1             0
			#  [8,]     1 0       0      1             0
			#  [9,]     1 0       1      0             0
			# [10,]     1 0       1      0             0
     
     ####(5) Get probability for highPRS ###
	
	#set.seed(SEED[1])
	
	product = as.vector(covmat %*% betas.lm)
	prob = matrix(exp(product)/(1+exp(product)),ncol=1)
		# > prob[1:5,]
		# [1] 0.3270372 0.1929154 0.1929154 0.3270372 0.1929154
	
	#### (6) simulate highPRS given probability ###########
	x=prob[1,] #[1] 0.327037
	pred.small = function(x){   sample(c(0,1), 1, prob=c(1-x, x))  }
	
	
	#set.seed(SEED[2])
	      
	highPRS = apply(prob, 1, pred.small)
		# > table(highPRS)
		# highPRS
		#     0     1 
		# 11560  4210 
		
		
	########### add LC death info ##########
	
	death = as.numeric(DATA[,"finaldeathlc"]) 
		#> death[1:10]
		# [1] NA NA NA NA NA NA NA NA  0 NA

		table(death)
		# lcdeath2
		#    0    1 
		# 1060  334		
    indic.LC.death = (death==1)*1
    indic.OCM.death = (death==0)*1
		# > sum(indic.LC.death,na.rm=T);sum(indic.OCM.death,na.rm=T)
		# [1] 302
		# [1] 1012
		
		
		
	DATA.big = cbind(DATA,covmat,highPRS=highPRS, indic.LC.death=indic.LC.death, indic.OCM.death = indic.OCM.death)
	
	DATA.big	


}#end of 




  ######## duplicate simulations for controling randomness #####
  
  myAssign.PRS.NLST2 = function(DATA,duplicate=T,times=50){		
		
		
		#candx_cal_year
		#rand_year
		
		#CANDX_DAYS	Days from randomization to first diagnosis of lung cancer 
	
		#CANCYR  Study year associated with first confirmed lung cancer.
		#CAN_SCR	Result of screen associated with the first confirmed lung cancer diagnosis
		#CANC_RPT_LINK	Is the reported lung cancer linked to a positive screen?
		#CANC_FREE_DAYS	Days from randomization to date when participant was last known to be free from lung cancer
        # CIGSMOK:  Smoking status at baseline,   0: Former, 1:Current
		
		################ (1.0) screen detected LC or not #############################
		
		
		#CAN_SCR	Result of screen associated with the first confirmed lung cancer diagnosis
			# 0	No Cancer
			# 1	Positive Screen
			# 2	Negative Screen
			# 3	Missed Screen
			# 4	Post Screening
		

		#CANC_RPT_LINK	Is the reported lung cancer linked to a positive screen?
		#CANC_FREE_DAYS	Days from randomization to date when participant was last known to be free from lung cancer
		#DEATH_DAYS:      Days since randomization at death from most definitive source
		#DCFDEATHLC:	Death certificate lung cancer death
		#FINALDEATHLC:	Final lung cancer death (Combined best information: EVP supplemented with DCF)
        
        ##############[1.1] Get covariates to be used for prediction: disease status, smoking status, and gender #####
        
		# > colnames(DATA)
		#   [1] "center"             "cigssmok"           "gender"             "rndgroup"          
		#   [5] "age"                "pkyr"               "scr_lat0"           "scr_lat1"          
		#   [9] "scr_lat2"           "ELIG"               "loclhil"            "locllow"           
		#  [13] "loclup"             "locrhil"            "locrlow"            "locrmid"           
		#  [17] "locunk"             "locrup"             "locoth"             "locmed"            
		#  [21] "loclmsb"            "locrmsb"            "loccar"             "loclin"            
		#  [25] "de_type"            "de_grade"           "de_stag"            "lesionsize"        
		#  [29] "scr_res0"           "scr_res1"           "scr_res2"           "scr_iso0"          
		#  [33] "scr_iso1"           "scr_iso2"           "RACE"               "ETHNIC"            
		#  [37] "EDUCAT"             "proc0"              "proc1"              "proc2"             
		#  [41] "biop0"              "biop1"              "biop2"              "invas0"            
		#  [45] "invas1"             "invas2"             "proclc"             "bioplc"            
		#  [49] "invaslc"            "medcomp0"           "medcomp1"           "medcomp2"          
		#  [53] "medcomplc"          "cancyr"             "can_scr"            "DCFICD"            
		#  [57] "evp_revr"           "evp_direct"         "evpincomplete"      "deathstat"         
		#  [61] "evpsel"             "evpcert"            "treatlc"            "mra_stat0"         
		#  [65] "mra_stat1"          "mra_stat2"          "no_proc_reas0"      "no_proc_reas1"     
		#  [69] "no_proc_reas2"      "study"              "conflc"             "evpdeath"          
		#  [73] "canc_rpt_link"      "canc_rpt_source"    "WDLOST"             "CONTACTSTATUS"     
		#  [77] "pid"                "smokeage"           "age_quit"           "smokeday"          
		#  [81] "smokeyr"            "dcfdeathlc"         "finaldeathlc"       "ndicd"             
		#  [85] "evpsent"            "deathcutoff"        "hist_group"         "rand_year"         
		#  [89] "scr_cal_year0"      "scr_days0"          "scr_age0"           "scr_cal_year1"     
		#  [93] "scr_days1"          "scr_age1"           "scr_cal_year2"      "scr_days2"         
		#  [97] "scr_age2"           "candx_cal_year"     "candx_days"         "candx_age"         
		# [101] "fup_cal_year"       "fup_days"           "fup_age"            "death_cal_year"    
		# [105] "death_days"         "death_age"          "canc_free_cal_year" "canc_free_days"    
		# [109] "canc_free_age"     
		
		
		### (1) Lung cancer status ###

		
			#tt1=as.numeric(DATA[,"cancyr"])
			# > table(tt1)
			# tt1
			#   0   1   2   3   4   5   6   7 
			# 176 111 146  47  67  65  22   1 
			# > sum(tt1[is.na(tt1)==F])
			# [1] 1276		
			#CANDX_DAYS
			#tt2=as.numeric(DATA[,"candx_days"])
			#  length(tt2[!is.na(tt2)])
			# [1] 635
			# > length(tt1[!is.na(tt1)])
			# [1] 635
		
		indic.LC =  (is.na(as.numeric(DATA[,"cancyr"]))==F)*1
			sum(indic.LC)
			#[1] 635
			
	   ### (2) smoking status: current, former #######		
		
			# > table(DATA[,"cigssmok"])
			# 
			#    0    1 
			# 8427 7343 
			
	   indic.former = 	(DATA[,"cigssmok"]==0)*1
	   indic.current = 	(DATA[,"cigssmok"]==1)*1
			# > sum(indic.former)
			# [1] 8427
			# > sum(indic.current)
			# [1] 7343	   
	   
	   	
	   ### (3) gender ##########
			# > table(DATA[,"gender"])
			# 
			#     1 
			# 15770 	   
	   indic.female = (DATA[,"gender"]==2)*1	


		########### add LC death info ##########
	
		death = as.numeric(DATA[,"finaldeathlc"]) 
			#> death[1:10]
			# [1] NA NA NA NA NA NA NA NA  0 NA

			table(death)
			# lcdeath2
			#    0    1 
			# 1060  334		
		indic.LC.death = (death==1)*1
		indic.OCM.death = (death==0)*1
			# > sum(indic.LC.death,na.rm=T);sum(indic.OCM.death,na.rm=T)
			# [1] 302
			# [1] 1012
 


      ### (4) make a matrix ###########
		# > betas
		#     (Intercept)               Y CIG_CAT_CURRENT  CIG_CAT_FORMER   GENDER_FEMALE 
		#      -1.2165470       0.9561218       0.4949310      -0.2146295      -0.2304710 

      
      covmat = cbind(inter=rep(1,nrow(DATA)),  Y=indic.LC, Current=indic.current, Former=indic.former, Gender.female=indic.female)
			# > covmat[1:10,]
			#       inter Y Current Former Gender.female
			#  [1,]     1 0       1      0             0
			#  [2,]     1 0       0      1             0
			#  [3,]     1 0       0      1             0
			#  [4,]     1 0       1      0             0
			#  [5,]     1 0       0      1             0
			#  [6,]     1 0       1      0             0
			#  [7,]     1 0       0      1             0
			#  [8,]     1 0       0      1             0
			#  [9,]     1 0       1      0             0
			# [10,]     1 0       1      0             0
			
			
			
	DATA.NLST = cbind(DATA,covmat, indic.LC.death=indic.LC.death, indic.OCM.death = indic.OCM.death)
			
     
     ####(5) Get probability for highPRS ###
	
	#set.seed(SEED[1])
	rows = 1:nrow(DATA)
	N=times*nrow(DATA)
	
	if(duplicate==F){   index = 1:nrow(DATA) }
	
	
	if(duplicate==T){   index=sample(rows,N, replace=T)  }
	
			# > index[1:10]
			#  [1]  9036 10956  8344   135  5740  6811  5787 12903 12802  9372
			# > table(index)[1:10]
			# index
			#  1  2  3  4  5  6  7  8  9 10 
			#  6  9  6  3  5  6  7  7  4  5 
	
	product = as.vector(covmat[index,] %*% betas.lm)
	prob = matrix(exp(product)/(1+exp(product)),ncol=1)
		# > prob[1:5,]
		# [1] 0.3270372 0.1929154 0.1929154 0.3270372 0.1929154
	
	#### (6) simulate highPRS given probability ###########
	x=prob[1,] #[1] 0.327037
	pred.small = function(x){   sample(c(0,1), 1, prob=c(1-x, x))  }
	
	
	#set.seed(SEED[2])
	      
	highPRS = apply(prob, 1, pred.small)
		# > table(highPRS)
		# highPRS
		#     0     1 
		# 11560  4210 
		
	DATA.NLST.PRS = cbind(DATA.NLST[index,], highPRS=highPRS)	
		
		
		
	ans=list(DATA.nlst =DATA.NLST, DATA.nlst.prs = DATA.NLST.PRS)
	ans

}#end of 




  ######### 9/3/2013: PLCO version from NLST ####
  ######## duplicate simulations for controling randomness #####
  
  myAssign.PRS.PLCO = function(DATA,betas,duplicate=T,times=50){		
		

				############## PLCO real data #########################################

				#cstatusl_cat_09t13
						#Does the participant have confirmed or reported lung cancer? Censored for 09/T13 cut-off.							
						# -1="Cancer before randomization" 0="No cancer"
						# 1="Confirmed cancer"
						# 2="Reported cancer" 3="Unconfirmable likely cancer" 4="Erroneous report of cancer" 
						#5="Reported Cancer (Pending CDCC)" 14="Confirmed non-target cancer"
				 
				
				#exitdays_incl_cnf09t13_adj  
							# Days from
							# Randomization Until
							# Exit Date for Lung Incidence
				#exitdays_incl_cnf09t13_unadj  
							# Days from
							# Randomization Until
							# Exit Date for Lung Incidence
				#intstatl_cat
							#Interval Status			
							# A categorized version of intstatl. Censored for 09/T13 cut-off.
							# 0="No cancer" 1="Control with cancer" 2="Never Screened" 3="Post-Screening" 4="Interval"
							# 5="Screen Dx"							
				# del_type											
						# Lung Cancer Histologic type											
						# .F="No Form"
						# 2="Squamous Cell Carcinoma" 3="Spindle Cell Carcinoma"
						# 4="Small Cell Carcinoma"
						# 5="Intermediate Cell Carcinoma" 7="Adenocarcinoma"
						# 8="Acinar Adenocarcinoma"
						# 9="Papillary Adenocarcinoma" 10="Bronchioalveolar Adenocarcinoma" 11="Adenocarcinoma w/Mucus Formation" 12="Large Cell Carcinoma"
						# 13="Giant Cell Carcinoma"
						# 14="Clear Cell Carcinoma" 15="Adenosquamous Carcinoma" 18="Adenoid Cystic Carcinoma" 30="Non-small cell (recoded)" 31="Carcinoma NOS (recoded)" 32="Mixed small and non-small cell (recoded)"
						# 33="Neuroendocrine NOS (recoded)" 98="not available"											
											
							
				#lungcancer_09t13	
						# Does Participant have Confirmed Lung Cancer
						# Confirmed lung cancer indicator. Censored for 09/T13 cut-off.
						# (format: ynf)
						# 0="No" 1="Yes"

				#cstatusl_cat_09t13
						# 09T13 Cancer Status
						# Does the participant have confirmed or reported lung cancer? Censored for 09/T13 cut-off.
						# -1="Cancer before randomization" 0="No cancer"
						# 1="Confirmed cancer"
						# 2="Reported cancer" 3="Unconfirmable likely cancer" 4="Erroneous report of cancer" 5="Reported Cancer (Pending CDCC)" 14="Confirmed non-target cancer"

				# cig_stat
						# Cigarette Smoking Status
						# Derived from questions 10 and 12. Set to .A if answers contradict or if one answer is missing but the other isn't
						# .A="Ambiguous"
						# .F="No Form"
						# .M="Missing"
						# 0="Never Smoked Cigarettes" 1="Current Cigarette Smoker" 2="Former Cigarette Smoker"

				# f_dthl_09t13
							# Is Lung Cancer the Final Cause of Death?
							# Cause of death flag that considers best available information (DCF, CDQ). Censored for 09/T13 cut-off.
							# (format: ynf)
							# 0="No" 1="Yes"

				# lung_fh
						# Family History of Lung Cancer?
						# Derived from question 21
						# .F="No Form"
						# .M="Missing"
						# 0="No"
						# 1="Yes, immediate family member" 9="Possibly - relative or cancer type not clear"


        
        ##############[1.1] Get covariates to be used for prediction: disease status, smoking status, and gender #####
        
		# > colnames(DATA)
		#   [1] "X"                            "EDUCAT"                       "MARITAL"                      "OCCUPAT"                      "PIPE"                        
		#   [6] "CIGAR"                        "XRAY"                         "ASPPD"                        "IBUPPD"                       "CURHORM"                     
		#  [11] "THORM"                        "emphys_f"                     "bronchit_f"                   "bq_returned"                  "race7"                       
		#  [16] "hispanic_f"                   "fh_cancer"                    "asp_f"                        "ibup_f"                       "horm_f"                      
		#  [21] "lung_fh"                      "smoked_f"                     "smokea_f"                     "rsmoker_f"                    "ssmokea_f"                   
		#  [26] "cigpd_f"                      "filtered_f"                   "cig_stat"                     "cig_stop_rnd"                 "cig_years_rnd"               
		#  [31] "pack_years_rnd"               "smkl15yr_rnd"                 "sct_flag_rnd"                 "bmi_curr"                     "bmi_curc"                    
		#  [36] "weight_f"                     "height_f"                     "copd_f"                       "cig_per_day"                  "iid"                         
		#  [41] "bq_compdays"                  "sqxo_cig_stat"                "sqxbq_cig_stat"               "sqxbq_cig_change"             "sqx_smk100"                  
		#  [46] "sqx_smokea"                   "sqx_nvr_smk_reg"              "sqx_smk_lgt"                  "sqx_smk_menthol"              "sqx_smk30days"               
		#  [51] "sqx_amt_smk"                  "sqx_smk_wake"                 "sqx_smk_crave1"               "sqx_smk_crave2"               "sqx_smk_crave3"              
		#  [56] "sqx_smk_crave4"               "sqx_smk_plan_quit"            "sqx_smk_try_quit"             "sqx_smk_quit_dur"             "sqx_smk_age_quit"            
		#  [61] "sqx_nicotine_gum"             "sqx_nicotine_patch"           "sqx_nicotine_oth"             "sqx_zyban"                    "sqx_smk_aid_none"            
		#  [66] "sqx_advise_quit"              "sqx_smk_exp_child"            "sqx_smk_exp_adult"            "sqx_smk_exp_work"             "sqx_worry_lungca"            
		#  [71] "sqx_risk_lungca"              "sqx_smk_HIS"                  "sqx_smk_BEHV"                 "sqx_smk_FTND"                 "sqx_years_quit"              
		#  [76] "sqx_cig_years"                "sqxbq_pack_years_bqplussqx"   "sqxo_pack_years"              "GENDER"                       "RNDGROUP"                    
		#  [81] "center"                       "age"                          "agelevel"                     "icdtop_l"                     "delsumm"                     
		#  [86] "delval"                       "icdbeh_l"                     "icdgrd_l"                     "icdmor_l"                     "del_type"                    
		#  [91] "lung_t"                       "lung_n"                       "lung_m"                       "lclinstg"                     "del_stag"                    
		#  [96] "lclin_t"                      "lclin_n"                      "lclin_m"                      "lung_type"                    "intstatl_cat"                
		# [101] "fin_smcl"                     "is_small_cell"                "xry_result0"                  "xry_result1"                  "xry_result2"                 
		# [106] "xry_result3"                  "curative_pneuml"              "curative_wsll"                "curative_chemol"              "curative_radl"               
		# [111] "curative_cyberl"              "primary_trtl_NSC"             "primary_trtl_small"           "xry_prot"                     "ph_lung_bq_status"           
		# [116] "ph_lung_sqx_status"           "ph_lung_bq"                   "ph_lung_sqx"                  "ph_lung_rnd_status"           "ph_lung_rnd"                 
		# [121] "conf_orderl"                  "primary_trtl_NSC_days"        "primary_trtl_small_days"      "xry_days0"                    "xry_days1"                   
		# [126] "xry_days2"                    "xry_days3"                    "has_xry0"                     "as_mass0"                     "as_mass_cnt0"                
		# [131] "as_mass_loc0"                 "as_mass_loc_pos0"             "as_mass_loc_side0"            "as_nodule0"                   "as_nodule_cnt0"              
		# [136] "as_nodule_loc0"               "as_nodule_loc_pos0"           "as_nodule_loc_side0"          "has_xry1"                     "as_mass1"                    
		# [141] "as_mass_cnt1"                 "as_mass_loc1"                 "as_mass_loc_pos1"             "as_mass_loc_side1"            "as_nodule1"                  
		# [146] "as_nodule_cnt1"               "as_nodule_loc1"               "as_nodule_loc_pos1"           "as_nodule_loc_side1"          "has_xry2"                    
		# [151] "as_mass2"                     "as_mass_cnt2"                 "as_mass_loc2"                 "as_mass_loc_pos2"             "as_mass_loc_side2"           
		# [156] "as_nodule2"                   "as_nodule_cnt2"               "as_nodule_loc2"               "as_nodule_loc_pos2"           "as_nodule_loc_side2"         
		# [161] "has_xry3"                     "as_mass3"                     "as_mass_cnt3"                 "as_mass_loc3"                 "as_mass_loc_pos3"            
		# [166] "as_mass_loc_side3"            "as_nodule3"                   "as_nodule_cnt3"               "as_nodule_loc3"               "as_nodule_loc_pos3"          
		# [171] "as_nodule_loc_side3"          "exitdays_incl_cnf09t13_unadj" "adjust_ind_incl"              "adjust_ind_mort"              "exitdays_mort_09t13_unadj"   
		# [176] "exitdays_incl_cnf09t13_adj"   "exitstat_incl_cnf09t13"       "exitdays_mort_09t13_adj"      "exitstat_mort_09t13"          "dobyear"                     
		# [181] "cstatusl_cat_09t13"           "lungcancer_09t13"             "candxl_pre0xry"               "is_dead_09t13"                "f_dthl_09t13"                
		# [186] "has_nrf_09t13"                "elig_bq"                      "elig_sqx"                     "old_xryresult0"               "old_xryresult1"              
		# [191] "old_xryresult2"               "old_xryresult3"              
		
		### (1) Lung cancer status ###

		
		# > tt1=DATA[,"cstatusl_cat_09t13"]
		# > table(tt1)
		# tt1
		#     0     1     2     3     4     5    14 
		# 36974   684    15    20    72    16    22 
		
		#tt2=DATA[,"lungcancer_09t13"]
		# > tt2[1:10]
		#  [1] "0" "0" "0" "0" "0" "0" "0" "0" "0" "0"
		# > table(tt2)
		# tt2
		#     0     1 
		# 37130   673 
		indic.LC = as.numeric(DATA[,"cstatusl_cat_09t13"])==1
		sum(indic.LC)
		#[1] 684
			
			
	   ### (2) smoking status: current, former #######	
	   ## 0="Never Smoked Cigarettes" 1="Current Cigarette Smoker" 2="Former Cigarette Smoker"	
		
			# > table(DATA[,"cig_stat"])
		# 
		#     0     1     2 
		# 21235  3671 12897 
			
	   indic.former = 	(DATA[,"cig_stat"]==2)*1
	   indic.current = 	(DATA[,"cig_stat"]==1)*1
			sum(indic.former)
			# [[1] 12897
			sum(indic.current)
			# [1] 3671	   
	   
	   	
	   ### (3) gender ##########
			# > table(DATA[,"gender"])
			# 
			#     1 
			# 15770 	   
	   indic.female = (DATA[,"GENDER"]=="F")*1	
	   
	   
	   #### (4) Family history ########
		# > table(DATA[,"lung_fh"])
		# 
		#     0     1     9     M 
		# 32530  4287   756   230 
	   
	   indic.fhl = DATA[,"lung_fh"]==1
	   sum(indic.fhl)


		########### add LC death info ##########
	
		death = as.numeric(DATA[,"f_dthl_09t13"]) 
			#> death[1:10]
			# [1] NA NA NA NA NA NA NA NA  0 NA

			table(death)
			#     0     1 
			# 37372   431 
		indic.LC.death = (death==1)*1
		#indic.OCM.death = (death==0)*1
			# > sum(indic.LC.death,na.rm=T);sum(indic.OCM.death,na.rm=T)
			# [1] 302
			# [1] 1012
 


      ### (4) make a matrix ###########
		# > betas
		#         (Intercept)     CIG_CAT_CURRENT      CIG_CAT_FORMER FAMILY_HX_LUNG_CODE       GENDER_FEMALE 
		#          -1.1047176           1.0072395           0.1976378          -0.2478984          -0.1343584 

      
      covmat = cbind(inter=rep(1,nrow(DATA)),  Current=indic.current, Former=indic.former, FH=indic.fhl,Gender.female=indic.female)
		#      inter Current Former FH Gender.female
		# [1,]     1       0      0  0             1
		# [2,]     1       0      1  0             1
		# [3,]     1       0      0  0             1
		# [4,]     1       0      1  0             1
		# [5,]     1       0      0  0             1
			
			
	DATA.PLCO = cbind(DATA,covmat, indic.LC.death=indic.LC.death, indic.LC=indic.LC)
			
     
     ####(5) Get probability for highPRS ###
	
	#set.seed(SEED[1])
	rows = 1:nrow(DATA)
	N=times*nrow(DATA)
	
	if(duplicate==F){   index = 1:nrow(DATA) }
	
	
	if(duplicate==T){   index=sample(rows,N, replace=T)  }
	
			# > index[1:10]
			#  [1]  9036 10956  8344   135  5740  6811  5787 12903 12802  9372
			# > table(index)[1:10]
			# index
			#  1  2  3  4  5  6  7  8  9 10 
			#  6  9  6  3  5  6  7  7  4  5 
	
	product = as.vector(covmat[index,] %*% betas.lm)
	prob = matrix(exp(product)/(1+exp(product)),ncol=1)
		# > prob[1:5,]
		# [1] 0.3270372 0.1929154 0.1929154 0.3270372 0.1929154
	
	#### (6) simulate highPRS given probability ###########
	x=prob[1,] #[1] 0.327037
	pred.small = function(x){   sample(c(0,1), 1, prob=c(1-x, x))  }
	
	
	#set.seed(SEED[2])
	      
	highPRS = apply(prob, 1, pred.small)
		# > table(highPRS)
		# highPRS
		#     0     1 
		# 11560  4210 
		
	DATA.PLCO.PRS = cbind(DATA.PLCO[index,], highPRS=highPRS)	
		
		
		
	ans=list(DATA.plco =DATA.PLCO, DATA.plco.prs = DATA.PLCO.PRS)
	ans

}#end of 


CombineStratified=function(file1,file2,file3){


			dat.cu=as.matrix(read.table(file1,sep="\t",header=T))
			dim(dat.cu)
			#[1] 25 14
			# > dat.cu[1:5,]
			#           Exact.Trait Ancestry  Region                Reported.Genes  Mapped.Genes        SNP N.case N.control
			# 1 Squamous cell carc.       AS 12q23.1 SLC17A8, NR1H4, SCYL2, GAS2L3             - rs12296850   3415      2335
			# 2         Lung cancer       AS 5p15.33                   TERT, hTERT          TERT  rs2736100   3413      2329
			# 3         Lung cancer       AS 6p21.32                  HLA class II             -  rs2395185   3407      2323
			# 4         Lung cancer       AS 6p21.32                  HLA class II             - rs28366298   2162      1294
			# 5         Lung cancer       AS  6q22.1                  ROS1, DCBLD1 ROS1 - DCBLD1  rs9387478   3414      2324
			#      Pvalue         Beta         SE        OR OR.95CI.L OR.95CI.U
			# 1 0.9096650 -0.009497884 0.08371044 0.9905471 0.8406580  1.167161
			# 2 0.2979251 -0.041389802 0.03976373 0.9594551 0.8875178  1.037223
			# 3 0.8505236  0.008733529 0.04634392 1.0087718 0.9211793  1.104693
			# 4 0.5830829 -0.034943941 0.06366325 0.9656595 0.8523791  1.093995
			# 5 0.1610621 -0.056428602 0.04026276 0.9451340 0.8734157  1.022741
			dat.fo=as.matrix(read.table(file2,sep="\t",header=T))
			dim(dat.fo)
			dat.fo[1:5,]
			#           Exact.Trait Ancestry  Region                Reported.Genes  Mapped.Genes        SNP N.case N.control
			# 1 Squamous cell carc.       AS 12q23.1 SLC17A8, NR1H4, SCYL2, GAS2L3             - rs12296850   1927      2001
			# 2         Lung cancer       AS 5p15.33                   TERT, hTERT          TERT  rs2736100   1927      1985
			# 3         Lung cancer       AS 6p21.32                  HLA class II             -  rs2395185   1925      1977
			# 4         Lung cancer       AS 6p21.32                  HLA class II             - rs28366298   1927      2001
			# 5         Lung cancer       AS  6q22.1                  ROS1, DCBLD1 ROS1 - DCBLD1  rs9387478   1922      1968
			#       Pvalue         Beta         SE        OR OR.95CI.L OR.95CI.U
			# 1 0.16404391 -0.131676475 0.09462244 0.8766246 0.7282314 1.0552561
			# 2 0.01847407 -0.109892726 0.04664423 0.8959302 0.8176545 0.9816995
			# 3 0.14291171 -0.075655592 0.05164085 0.9271355 0.8378873 1.0258900
			# 4 0.05004650 -0.106782043 0.05449269 0.8987215 0.8076813 1.0000236
			# 5 0.96364874 -0.002147651 0.04712313 0.9978547 0.9098196 1.0944081

			dat.ne=as.matrix(read.table(file3,sep="\t",header=T))
			dim(dat.ne)

	### just take the subset 15 snps


	row.names(dat.cu)=dat.cu[,"SNP"]
	row.names(dat.fo)=dat.fo[,"SNP"]
	row.names(dat.ne)=dat.ne[,"SNP"]
 
	dat.cu = dat.cu[snp.orders,]  
	dat.fo = dat.fo[snp.orders,]  
	dat.ne = dat.ne[snp.orders,]  
 
			#### stack them together ###
			
			nSNPs=nrow(dat.cu)

			for(j in 1:nSNPs){

				tm1=c("Current",dat.cu[j,])
				tm2=c("Former",dat.fo[j,])
				tm3=c("Never",dat.ne[j,])
	
				tm=rbind(tm3,tm2,tm1)
				if(j==1) out=NULL
				out=rbind(out,tm)	


			}#3nd of 

			colnames(out)[1]="Smoking"
			row.names(out)=1:nrow(out)
			# > out[1:5,]
			#   Smoking   Exact.Trait           Ancestry Region    Reported.Genes                  Mapped.Genes SNP          N.case N.control Pvalue         Beta           
			# 1 "Never"   "Squamous cell carc." "AS"     "12q23.1" "SLC17A8, NR1H4, SCYL2, GAS2L3" "-"          "rs12296850" "350"  "1379"    "0.086994635"  " 0.311999801" 
			# 2 "Former"  "Squamous cell carc." "AS"     "12q23.1" "SLC17A8, NR1H4, SCYL2, GAS2L3" "-"          "rs12296850" "1927" "2001"    "1.640439e-01" "-0.131676475" 
			# 3 "Current" "Squamous cell carc." "AS"     "12q23.1" "SLC17A8, NR1H4, SCYL2, GAS2L3" "-"          "rs12296850" "3415" "2335"    "9.096650e-01" "-0.0094978836"
			# 4 "Never"   "Lung cancer"         "AS"     "5p15.33" "TERT, hTERT"                   "TERT"       "rs2736100"  "350"  "1365"    "0.001841117"  "-0.301010699" 
			# 5 "Former"  "Lung cancer"         "AS"     "5p15.33" "TERT, hTERT"                   "TERT"       "rs2736100"  "1927" "1985"    "1.847407e-02" "-0.109892726" 
			#   SE           OR          OR.95CI.L   OR.95CI.U  
			# 1 "0.18229945" "1.3661544" "0.9557039" "1.9528830"
			# 2 "0.09462244" "0.8766246" "0.7282314" "1.0552561"
			# 3 "0.08371044" "0.9905471" "0.8406580" "1.1671613"
			# 4 "0.09664094" "0.7400699" "0.6123649" "0.8944069"
			# 5 "0.04664423" "0.8959302" "0.8176545" "0.9816995"
			# > colnames(out)
			#  [1] "Smoking"        "Exact.Trait"    "Ancestry"       "Region"         "Reported.Genes" "Mapped.Genes"   "SNP"            "N.case"         "N.control"     
			# [10] "Pvalue"         "Beta"           "SE"             "OR"             "OR.95CI.L"      "OR.95CI.U"     
			# > 

			mycolnames=c("SNP","Reported.Genes","Region","Exact.Trait","Ancestry","Smoking","N.case","N.control","OR","Pvalue","OR.95CI.L","OR.95CI.U")
			out2=out[,mycolnames]
			# > out2[1:3,]
			#   SNP          Mapped.Genes Reported.Genes                  Region    Exact.Trait           Ancestry Smoking   N.case N.control OR          Pvalue         OR.95CI.L  
			# 1 "rs12296850" "-"          "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "Squamous cell carc." "AS"     "Never"   "350"  "1379"    "1.3661544" "0.086994635"  "0.9557039"
			# 2 "rs12296850" "-"          "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "Squamous cell carc." "AS"     "Former"  "1927" "2001"    "0.8766246" "1.640439e-01" "0.7282314"
			# 3 "rs12296850" "-"          "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "Squamous cell carc." "AS"     "Current" "3415" "2335"    "0.9905471" "9.096650e-01" "0.8406580"
			#   OR.95CI.U  
			# 1 "1.9528830"
			# 2 "1.0552561"
			# 3 "1.1671613"


			CI=paste(round(as.numeric(out2[,"OR.95CI.U"]),2),"-",round(as.numeric(out2[,"OR.95CI.L"]),2),sep="")
			out3=cbind(out2[,1:10],CI=CI)
			# > out3[1:4,]
			#   SNP          Reported.Genes                  Region    Exact.Trait           Ancestry Smoking   N.case N.control OR          Pvalue         CI         
			# 1 "rs12296850" "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "Squamous cell carc." "AS"     "Never"   "350"  "1379"    "1.3661544" "0.086994635"  "1.95-0.96"
			# 2 "rs12296850" "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "Squamous cell carc." "AS"     "Former"  "1927" "2001"    "0.8766246" "1.640439e-01" "1.06-0.73"
			# 3 "rs12296850" "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "Squamous cell carc." "AS"     "Current" "3415" "2335"    "0.9905471" "9.096650e-01" "1.17-0.84"
			# 4 "rs2736100"  "TERT, hTERT"                   "5p15.33" "Lung cancer"         "AS"     "Never"   "350"  "1365"    "0.7400699" "0.001841117"  "0.89-0.61"
   
			#### rename columns ####
			
			# > table(out3[,"Exact.Trait"])
			# 
			# Lung adenocarcinoma         Lung cancer               NSCLC    NSCLC - Survival     SCLC - Survival Squamous cell carc. 
			#                  24                  39                   3                   3                   3                   3 
			tm1=gsub("Lung adenocarcinoma","AD",out3[,"Exact.Trait"])
			tm2=gsub("Lung cancer","All",tm1)
			tm3=gsub("Squamous cell carc.","SQ",tm2)
			tm4=gsub("NSCLC - Survival","NSCLC.Surv",tm3)
			tm5=gsub("SCLC - Survival","SCLC.Surv",tm4)
			# > table(tm5)
			# tm5
			#         AD        All      NSCLC NSCLC.Surv  SCLC.Surv         SQ 
			#         24         39          3          3          3          3 

			out3[,"Exact.Trait"]=tm5
			out3[,"OR"]=round(as.numeric(out3[,"OR"]),2)
			# out3[1:5,]
			# > out3[1:5,]
			#   SNP          Reported.Genes                  Region    Exact.Trait Ancestry Smoking   N.case N.control OR     Pvalue         CI         
			# 1 "rs12296850" "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "SQ"        "AS"     "Never"   "350"  "1379"    "1.37" "0.086994635"  "1.95-0.96"
			# 2 "rs12296850" "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "SQ"        "AS"     "Former"  "1927" "2001"    "0.88" "1.640439e-01" "1.06-0.73"
			# 3 "rs12296850" "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "SQ"        "AS"     "Current" "3415" "2335"    "0.99" "9.096650e-01" "1.17-0.84"
			# 4 "rs2736100"  "TERT, hTERT"                   "5p15.33" "All"       "AS"     "Never"   "350"  "1365"    "0.74" "0.001841117"  "0.89-0.61"
			# 5 "rs2736100"  "TERT, hTERT"                   "5p15.33" "All"       "AS"     "Former"  "1927" "1985"    "0.9"  "1.847407e-02" "0.98-0.82"
 
 			### insert empty cells ###
 			
			 dels=(1:nrow(out3))[-(3*(1:nSNPs)-2)]
			 # > dels[1:10]
			#  [1]  2  3  5  6  8  9 11 12 14 15

			out3[dels, 1:5]=""
			out3[1:5,]
			# > out3[1:5,]
			#   SNP          Reported.Genes                  Region    Exact.Trait Ancestry Smoking   N.case N.control OR     Pvalue         CI         
			# 1 "rs12296850" "SLC17A8, NR1H4, SCYL2, GAS2L3" "12q23.1" "SQ"        "AS"     "Never"   "350"  "1379"    "1.37" "0.086994635"  "1.95-0.96"
			# 2 ""           ""                              ""        ""          ""       "Former"  "1927" "2001"    "0.88" "1.640439e-01" "1.06-0.73"
			# 3 ""           ""                              ""        ""          ""       "Current" "3415" "2335"    "0.99" "9.096650e-01" "1.17-0.84"
			# 4 "rs2736100"  "TERT, hTERT"                   "5p15.33" "All"       "AS"     "Never"   "350"  "1365"    "0.74" "0.001841117"  "0.89-0.61"
			# 5 ""           ""                              ""        ""          ""       "Former"  "1927" "1985"    "0.9"  "1.847407e-02" "0.98-0.82"
			
			out3
 

}#end of 


            match.order=function(x,y){   ### wait! y doesn't need to be exact same length!  it can be length(x) <= length(y)

                ####### this match y to x ----> gives the row numbers of y so that y[ans] will be exactly the same as x
                ####### they should be unique
                
                ans=NULL
                
                #> cbind(x,y[1:length(x)])
                #      x   y
                #       x   y
                # [1,] "a" "c"
                # [2,] "b" "h"
                # [3,] "c" "d"
                # [4,] "d" "g"
                # [5,] "e" "f"
                # [6,] "f" "a"
                # [7,] "g" "j"
                # [8,] "h" "i"
                # [9,] "i" "b"
                #[10,] "j" "e"


                x2=x
                dim(x2)=c(length(x),1)

                #xx=x[1]
                mymatch=function(xx,y){

                   tm=(1:length(y))[y %in% xx]
                   if(length(tm)==0) tm=NA  # doesn't exist then NA
                   tm
                   
                }# end of mymatch
                mymatch(x[1],y)

                #ans0=sapply(x2,mymatch,y=y)
                
                ans0=apply(x2,1,mymatch,y=y)
                
                ans=as.vector(ans0)

                ans2 = list(ans.na= ans,ans.no.na=ans[!is.na(ans)])
                ans2
                # [1]  6  9  1  3 10  5  4  2  8  7

                #> cbind(x,y[ans])
                #      x
                # [1,] "a" "a"
                # [2,] "b" "b"
                # [3,] "c" "c"
                # [4,] "d" "d"
                # [5,] "e" "e"
                # [6,] "f" "f"
                # [7,] "g" "g"
                # [8,] "h" "h"
                # [9,] "i" "i"
                #[10,] "j" "j"
                
                ######## example ####################
                #
                #x=c("a","b","c")
                #y=c("d","a","e","b")
                #
                #ans=match.order(x,y)  ### wait! y doesn't need to be exact same length!  it can be length(x) <= length(y)
                
                #> ans[[1]]
                #[1] 2 4 NA  --> "c" doens't exist in y
                #ans[[2]]
                #[1] 2 4
                #> x
                #[1] "a" "b" "c"
                #> y[ans[[2]]]
                #[1] "a" "b"
                #

            }# end of mmatch order
            
            
#x=c("a","b","c")
#y=c("d","a","e","b")
#
#ans=match.order(x,y)  ### wait! y doesn't need to be exact same length!  it can be length(x) <= length(y)

#> ans[[1]]
#[1] 2 4 NA  --> "c" doens't exist in y
#ans[[2]]
#[1] 2 4
#> x
#[1] "a" "b" "c"
#> y[ans[[2]]]
#[1] "a" "b"
# source("source.gwas.functions.R")





