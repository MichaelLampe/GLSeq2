#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Oleg Moskvin; info@scienceforever.com 
# Michael Lampe
# June 2015 
#########################################################
# Collects versioning of the various tools used here for reference.  
#########################################################


log.file <----- a
#
# Summarize alignment and versioning
# 
if (alignment == "alignment"){
  if (aAlgor == "BWA"){
    
  }
  if (aAlgor == "Bowtie"){
    
  }
  if (aAlgor == "Bowtie2"){
    
  }
  if (aAlgor == "CUSHAW"){
    if(GPU.accel){
      
    } 
  }
}
if (counting = "counting"){
  if ("FeatureCounts" %in% cAlgor){
    
  }
  if ("HTSeq" %in% cAlgor){
    
  }
  if ("RSEM" %in% cAlgor){
    
  }
}