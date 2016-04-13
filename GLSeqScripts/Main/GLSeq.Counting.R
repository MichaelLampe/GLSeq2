source ("GLSeq.Util.R")
# Allows us to give a warning if no counting protocol ended up being run.
occured <- FALSE
# Clears count.comm
# We'll enclose the counting in paren and then background each counting process to parallelize them
count.comm <- "("

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
  if (!occured){
    ref.dir <- paste(base.dir, rGenome, sep="")
    ref.dir <- trailDirCheck(ref.dir)
    refCopy <- paste("cp ", refGFFname," ",dest.dir, sep="")
    printOrExecute(refCopy,Condor)
  }
    occured <- TRUE
    setwd(base.dir)
    source ("GLSeq.RsemCount.R")
    count.comm <- paste(count.comm,"&")
}

################################################
#HTSeq Counting Protocol
#
# Note: HTSeq tends to take pretty long, so I like that it is last.
# Then you have some files to check out while waiting for it...
################################################
if ("HTSeq" %in% cAlgor){
  if (!occured){
    ref.dir <- paste(base.dir, rGenome, sep="")
    ref.dir <- trailDirCheck(ref.dir)
    refCopy <- paste("cp ",refGFFname," ",dest.dir, sep="")
    printOrExecute(refCopy,Condor)
  }
  occured <- TRUE
  setwd(base.dir)
  source ("GLSeq.HTSeq.R")
  count.comm <- paste(count.comm,"&")
}

################################################
#Cufflinks Counting Protocol
################################################
if ("Cufflinks" %in% cAlgor){
  if (!occured){
    ref.dir <- paste(base.dir, rGenome, sep="")
    ref.dir <- trailDirCheck(ref.dir)
    refCopy <- paste("cp ",refGFFname," ",dest.dir, sep="")
    printOrExecute(refCopy,Condor)
  }
  occured <- TRUE
  setwd(base.dir)
  source ("GLSeq.Cufflinks.R")
  count.comm <- paste(count.comm,"&")
}

################################################
#FeatureCounts Counting Protocol
################################################
if ("FeatureCounts" %in% cAlgor){
  if (!occured){
    ref.dir <- paste(base.dir, rGenome, sep="")
    ref.dir <- trailDirCheck(ref.dir)
    refCopy <- paste("cp ",refGFFname," ",dest.dir, sep="")
    printOrExecute(refCopy,Condor)
  }
  occured <- TRUE
  setwd(base.dir)
  # Quick fix in case there is a problem with the resName for some reason.
  # The main reason this condition occurs is if the user is only counting the
  # file and not aligning immediately prior to counting
  if (is.null(this.resName)) this.resName <- paste(dest.dir,text.add)
  ##########
  # Make sure GTF file is copied (It is important for Feature Counts)
  ##########
  #
  countfile <- paste(this.resName,"FeatureCounts","counts", sep=".")
  # Feature counts is actually an RScript.  So instead of sourcing it, we simply prepare the program to run another R command.
  script <- paste("Rscript",
                  paste(base.dir,"GLSeq.FeatureCounts.R",sep=""),
                  countable.sam,
                  destDirFeatureCountsCount,
                  refGFFname,
                  dest.dir,
                  this.resName,
                  paired.end,
                  idAttr,
                  FeatureCountsSpecialOptions)
  counts.summary <- paste(this.resName,".FeatureCounts.summary.txt",sep="")
  counts.data <- paste(this.resName,".FeatureCounts.counts.csv",sep="")
  counts.stats <- paste(this.resName,".FeatureCounts.stats.csv",sep="")
  counts.annotation <- paste(this.resName,".FeatureCounts.annotations.csv",sep="")
  #
  move.command <- paste("mv",counts.summary,destDirFeatureCountsCount)
  move.command <- paste(move.command,"&& mv",counts.data,destDirFeatureCountsCount)
  move.command <- paste(move.command,"&& mv",counts.stats,destDirFeatureCountsCount)
  move.command <- paste(move.command,"&& mv",counts.annotation,destDirFeatureCountsCount)
  if (count.comm != "") count.comm <- paste(count.comm,script,"&&",move.command)
  if (count.comm == "(") count.comm <- paste(script,"&&",move.command)
  count.comm <- paste(count.comm,"&")
}

################################################
#Rockhopper Counting Protocol
################################################
if (aAlgor != "Rockhopper"){
  if ("Rockhopper" %in% cAlgor){
    if (!occured) {
      ref.dir <- paste(base.dir, rGenome, sep="")
      ref.dir <- trailDirCheck(ref.dir)
      refCopy <- paste("cp ",refGFFname," ",dest.dir, sep="")
      printOrExecute(refCopy,Condor)
    }
    occured <- TRUE
    setwd(base.dir)
    source ("GLSeq.RockhopperCount.R")
    count.comm <- paste(count.comm,"&")
  }
}
# Close paren
count.comm <- paste(count.comm,")")
if (!occured){
  warning("No counting protocol was initiated even though you indicated counting should occur. Please make sure you have supplied a supported counting protocol.")
}