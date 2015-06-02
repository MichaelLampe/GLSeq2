source ("GLSeq.Util.R")
setwd(base.dir)
occured <- FALSE
################################################
#BWA Alignment Protocol
################################################
if (qAlgor == "BWA"){
  occured <- TRUE
  source ("GLSeq.BWA.R")
}
################################################
#Bowtie and Bowtie2 Alignment Protocol
################################################
if (qAlgor == "Bowtie" || qAlgor == "Bowtie2"){
  occured <- TRUE
  source ("GLSeq.Bowtie.R")
}

################################################
#Cushaw w/ and w/o GPU Accel Alignment Protocol
################################################
if (qAlgor == "Cushaw"){
  occured <- TRUE
  source ("GLSeq.CUSHAW.R")
}
if (!occured){
  warning("No alignment protocol was initiated.  Please make sure you have supplied a supported alignment setting")
}
