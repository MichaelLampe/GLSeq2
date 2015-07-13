source ("GLSeq.Util.R")
occured <- FALSE
count.comm <- ""
################################################
#HTSeq Counting Protocol
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
  if (is.null(this.resName)) this.resName <- text.add
  ##########
  # Make sure GTF file is copied
  ##########
  ref.dir <- paste(base.dir, rGenome, sep="")
  refCopy <- paste("cd ", ref.dir, " && cp ",refGFFname," ",dest.dir, sep="")
  system(refCopy)
  #
  countfile <- paste(this.resName,"FeatureCounts","counts", sep=".") 
  script <- paste( "Rscript GLSeq.FeatureCounts.R",countable.sam,rGenome,refGFFname,dest.dir,this.resName,paired.end,sep=" ")
  if (count.comm != "") count.comm <- paste(count.comm," && ","cd ", base.dir, " && ",script," && ","cd ",dest.dir,"&&","mv",countfile,destDirFeatureCountsCount) 
  if (count.comm == "") count.comm <- paste("cd ", base.dir, " && ",script," && ","cd ",dest.dir,"&&","mv",countfile,destDirFeatureCountsCount)
}
################################################
#RSEM Counting Protocol
################################################
#
# RSEM only works with NON-GAPPED aligners
# 
# Decided to hardcode in the Gapped-Alignment only.
#
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
if (!occured){
  warning("No counting protocol was initiated.  Please make sure you have supplied a supported counting setting")
}

