source("GLSeq.Util.R")

####################################
# Copy genome indices to the destination dir:
####################################
copy.genome <- function(base.dir,rGenome,refFASTAname,dest.dir){
  if (is.null(base.dir) || is.null(rGenome) || is.null(refFASTAname) || is.null(dest.dir)) stop("Arguments should not be NULL")
  dest.dir <- trailDirCheck(dest.dir)
  base.dir <- trailDirCheck(base.dir)
  rGenome <- trailDirCheck(rGenome)
  ref.dir <- paste(base.dir, rGenome, sep="")
  indCopy <- paste("cp",paste(ref.dir,refFASTAname,sep=""),dest.dir)
  indCopy
}

countable.sam.name <- function(this.resName){
  if(is.null(this.resName)) stop("Arguments should not be NULL")
  countable.sam <- paste(this.resName, "countable.sam", sep=".")
  countable.sam
}

assign.resName <- function(name,text.add){
  if (is.null(name) || is.null(text.add)) stop("Arguments should not be NULL")
  this.resName <- paste(name, text.add, sep=".")
  this.resName
}

assign.name <- function(file,paired.end){
  if(is.null(file) || is.null(paired.end)) stop("Arguments should not be NULL")
  if(!is.logical(paired.end)) stop("Paired end must be a LOGICAL")
  name <- file
  if (paired.end){
    name <- substr(name,1,(nchar(name) - 5))
  } else{
    name <- substr(name,1,(nchar(name) - 3))
  }
  name
}