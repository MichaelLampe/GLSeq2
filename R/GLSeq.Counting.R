source ("GLSeq.Util.R")
occured <- FALSE
################################################
#HTSeq Counting Protocol
################################################
if ("HTSeq" %in% cAlgor){
  occured <- TRUE
  setwd(base.dir)
  source ("GLSeq.HTSeq.R")
}
################################################
#FeatureCounts Counting Protocol
################################################
if ("FeatureCounts" %in% cAlgor){
  occured <- TRUE
  ##########
  # Make sure GTF file is copied
  ##########
  setwd(base.dir)
  ref.dir <- paste(base.dir, rGenome, sep="")
  refCopy <- paste("cd ", ref.dir, " && cp ",refGFFname," ",dest.dir, sep="")
  system(refCopy)
  #
  script <- paste( "Rscript GLSeq.FeatureCounts.R",countable.sam,rGenome,refGFFname,dest.dir,this.resName,paired.end,sep=" ")
  if (count.comm != "") count.comm <- paste(count.comm," && ","cd ", base.dir, " && ",script," && ","cd ",dest.dir) 
  if (count.comm == "") count.comm <- paste("cd ", base.dir, " && ", "Rscript GLSeq.FeatureCounts.R"," && ","cd ",dest.dir)
}
################################################
#RSEM Counting Protocol
################################################
#
# RSEM only works with NON-GAPPED aligners
# 
################################################
if ("RSEM" %in% cAlgor){
  occured <- TRUE
  setwd(base.dir)
  source ("GLSeq.RsemCount.R")
}
if (!occured){
  warning("No counting protocol was initiated.  Please make sure you have supplied a supported counting setting")
}

