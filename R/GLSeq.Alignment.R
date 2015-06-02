source ("GLSeq.Util.R")
setwd(dest.dir)
################################################
#BWA Alignment Protocol
################################################
if (qAlgor == "BWA"){
  source ("GLSeq.BWA.R")
}
################################################
#Bowtie and Bowtie2 Alignment Protocol
################################################
if (qAlgor == "Bowtie" || qAlgor == "Bowtie2"){
  source ("GLSeq.Bowtie.R")
}

################################################
#Cushaw w/ and w/o GPU Accel Alignment Protocol
################################################
if (qAlgor == "Cushaw"){
  source ("GLSeq.Cushaw.R")
}
