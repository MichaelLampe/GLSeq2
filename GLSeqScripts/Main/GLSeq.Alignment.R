source ("GLSeq.Util.R")
setwd(base.dir)
occured <- FALSE
comm.stack.pool <- ""

################################################
#BWA Alignment Protocol
################################################
if (aAlgor == "BWA"){
  occured <- TRUE
  source ("GLSeq.BWA.R")
}

################################################
#Bowtie and Bowtie2 Alignment Protocol
################################################
if (aAlgor == "Bowtie" || aAlgor == "Bowtie2"){
  occured <- TRUE
  source ("GLSeq.Bowtie.R")
}

################################################
#Cushaw w/ and w/o GPU Accel Alignment Protocol
################################################
if (aAlgor == "Cushaw"){
  occured <- TRUE
  source ("GLSeq.CUSHAW.R")
}

################################################
#Cushaw w/ and w/o GPU Accel Alignment Protocol
################################################
if (aAlgor == "Cushaw_GPU"){
  occured <- TRUE
  source ("GLSeq.CUSHAW.R")
}

################################################
# TopHat Alignment Protocol
################################################
if (aAlgor == "TopHat"){
  occured <- TRUE
  source ("GLSeq.TopHat.R")
}

################################################
# Rockhopper Alignment Protocol
################################################
if (aAlgor == "Rockhopper"){
  occured <- TRUE
  source ("GLSeq.Rockhopper.R")
}
################################################
# STAR Alignment Protocol
################################################
if (aAlgor == "STAR"){
  occured <- TRUE
  source("GLSeq.STAR.R")
}
################################################
# HISAT Alignment Protocol
################################################
if (aAlgor == "HISAT"){
  occured <- TRUE
  source("GLSeq.HISAT.R")
}

if (!occured){
  warning("No alignment protocol was initiated.  Please make sure you have supplied a supported alignment setting")
}
