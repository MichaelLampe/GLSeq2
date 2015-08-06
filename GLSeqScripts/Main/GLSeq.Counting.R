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
  if (!occured){
    ref.dir <- paste(base.dir, rGenome, sep="")
    ref.dir <- trailDirCheck(ref.dir)
    refCopy <- paste("cp ",paste(ref.dir,refGFFname,sep="")," ",dest.dir, sep="")
    printOrExecute(refCopy,Condor)
  }
  if (alignment == "alignment"){
    if (aAlgor == "Bowtie" || aAlgor == "Bowtie2"){
      occured <- TRUE
      setwd(base.dir)
      source ("GLSeq.RsemCount.R")
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
  if (!occured){
    ref.dir <- paste(base.dir, rGenome, sep="")
    ref.dir <- trailDirCheck(ref.dir)
    refCopy <- paste("cp ",paste(ref.dir,refGFFname,sep="")," ",dest.dir, sep="")
    printOrExecute(refCopy,Condor)
  }
  occured <- TRUE
  setwd(base.dir)
  source ("GLSeq.HTSeq.R")
}
################################################
#Cufflinks Counting Protocol
################################################
if ("Cufflinks" %in% cAlgor){
  if (!occured){
    ref.dir <- paste(base.dir, rGenome, sep="")
    ref.dir <- trailDirCheck(ref.dir)
    refCopy <- paste("cp ",paste(ref.dir,refGFFname,sep="")," ",dest.dir, sep="")
    printOrExecute(refCopy,Condor)
  }
  occured <- TRUE
  setwd(base.dir)
  source ("GLSeq.Cufflinks.R")
}
################################################
#FeatureCounts Counting Protocol
################################################
if ("FeatureCounts" %in% cAlgor){
  if (!occured){
    ref.dir <- paste(base.dir, rGenome, sep="")
    ref.dir <- trailDirCheck(ref.dir)
    refCopy <- paste("cp ",paste(ref.dir,refGFFname,sep="")," ",dest.dir, sep="")
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
  script <- paste("Rscript",paste(base.dir,"GLSeq.FeatureCounts.R",sep=""),countable.sam,rGenome,refGFFname,dest.dir,this.resName,paired.end,idAttr)
  if (count.comm != "") count.comm <- paste(count.comm,"&& ",script,"&&","mv",countfile,destDirFeatureCountsCount)
  if (count.comm == "") count.comm <- paste(script,"&&","mv",countfile,destDirFeatureCountsCount)
}

if (!occured){
  warning("No counting protocol was initiated even though you indicated counting should occur. Please make sure you have supplied a supported counting protocol.")
}