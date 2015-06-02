source ("GLSeq.Util.R")
setwd(dest.dir)
################################################
#HTSeq Counting Protocol
################################################
if (cAlgor == "HTSeq"){
  source ("GLSeq.HTSeq.R")
}
################################################
#FeatureCounts Counting Protocol
################################################
if (cAlgor == "FeatureCounts"){
  source ("GLSeq.FeatureCounts.R")
}
################################################
#RSEM Counting Protocol
################################################
if (cAlgor == "RSEM"){
  source ("GLSeq.RSEM.R")
}
