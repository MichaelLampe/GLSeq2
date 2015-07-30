source ("GLSeq.Util.R")
# Allows us to give a warning if no counting protocol ended up being run.
occured <- FALSE
# Clears count.comm
count.comm <- ""
################################################
#RSEM Counting Protocol
################################################
#
# RSEM only works with NON-GAPPED aligners
# 
# Decided to hardcode in the Gapped-Alignment only.
# Might change this back in the future, but I think it is worthwhile as it inhibits some pretty weird errors.
################################################
if ("RSEM" %in% cAlgor){
  if (alignment == "alignment"){
    if (aAlgor == "Bowtie" || aAlgor == "Bowtie2"){
      occured <- TRUE
      setwd(base.dir)
      source ("GLSeq.RsemCount.R")
    } else{
      add.to.logs("RSEM was indicated as a counting option, but a supported alignment protocol was not used.  Not counting using RSEM.",log.file)
    }
  } else{
    # In case the user is just counting
    occured <- TRUE
    setwd(base.dir)
    source ("GLSeq.RsemCount.R")
  }
}
################################################
#HTSeq Counting Protocol
#
# Note: HTSeq tends to take pretty long, so I like that it is last.
# Then you have some files to check out while waiting for it...
################################################
if ("HTSeq" %in% cAlgor){
  occured <- TRUE
  setwd(base.dir)
  source ("GLSeq.HTSeq.R")
}
################################################
#Cufflinks Counting Protocol
################################################
if ("Cufflinks" %in% cAlgor){
  occured <- TRUE
  setwd(base.dir)
  source ("GLSeq.Cufflinks.R")
}
################################################
#FeatureCounts Counting Protocol
################################################
if ("FeatureCounts" %in% cAlgor){
  occured <- TRUE
  setwd(base.dir)
  # Quick fix in case there is a problem with the resName for some reason.
  # The main reason this condition occurs is if the user is only counting the
  # file and not aligning immediately prior to counting
  if (is.null(this.resName)) this.resName <- text.add
  ##########
  # Make sure GTF file is copied (It is important for Feature Counts)
  ##########
  ref.dir <- paste(base.dir, rGenome, sep="")
  refCopy <- paste("cd ", ref.dir, " && cp ",refGFFname," ",dest.dir, sep="")
  system(refCopy)
  #
  countfile <- paste(this.resName,"FeatureCounts","counts", sep=".") 
  # Feature counts is actually an RScript.  So instead of sourcing it, we simply prepare the program to run another R command.
  script <- paste( "Rscript GLSeq.FeatureCounts.R",countable.sam,rGenome,refGFFname,dest.dir,this.resName,paired.end)
  if (count.comm != "") count.comm <- paste(count.comm,"&& ","cd ", base.dir, "&& ",script,"&&","cd ",dest.dir,"&&","mv",countfile,destDirFeatureCountsCount) 
  if (count.comm == "") count.comm <- paste("cd", base.dir, "&&",script,"&&","cd ",dest.dir,"&&","mv",countfile,destDirFeatureCountsCount)
}

if (!occured){
  warning("No counting protocol was initiated even though you indicated counting should occur. Please make sure you have supplied a supported counting protocol.")
}