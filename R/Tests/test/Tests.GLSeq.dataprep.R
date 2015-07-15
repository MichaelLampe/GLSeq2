tests.load.dataFile <- function(){
  checkEquals("GLSeq.vars.test.rda",load.dataFile("test"))
  checkException(load.dataFile(NULL),"Arguments should not be NULL")
}

tests.check.presplit <- function(){
  checkEquals(TRUE,check.presplit(TRUE,TRUE))
  checkEquals(FALSE,check.presplit(FALSE,TRUE))
  checkEquals(FALSE,check.presplit(TRUE,FALSE))
  checkEquals(FALSE,check.presplit(FALSE,FALSE))
  
  checkException(check.presplit(NULL,NULL),"Arguments should not be NULL")
  checkException(check.presplit(TRUE,NULL),"Arguments should not be NULL")
  checkException(check.presplit(NULL,TRUE),"Arguments should not be NULL")
  checkException(check.presplit("a","a"),"Arguments must be LOGICALS")
  checkException(check.presplit(TRUE,"a"),"Arguments must be LOGICALS")
  checkException(check.presplit("a",TRUE),"Arguments must be LOGICALS")
}

tests.check.nStreamsDataPrep <- function(){
  checkEqualsNumeric(2,check.nStreamsDataPrep(rbind(1,2,3),2))
  checkEqualsNumeric(3,check.nStreamsDataPrep(rbind(1,2,3),4))
  checkEqualsNumeric(1,check.nStreamsDataPrep(rbind(1,2,3),1))
  
  checkException(check.nStreamsDataPrep(rbind("test1","test2"),NULL),"Arguments should not be NULL")
  checkException(check.nStreamsDataPrep(NULL, 3),"Arguments should not be NULL")
  checkException(check.nStreamsDataPrep(NULL,NULL),"Arguments should not be NULL")  
}

tests.check.ifFiles <- function(){
  checkEquals(NULL,check.ifFiles("a","a"))
  checkEquals(NULL,check.ifFiles(NULL,"a"))
  checkEquals(NULL,check.ifFiles("a",NULL))
  checkException(check.ifFiles(NULL,NULL),"Please check the list of raw FASTQ files and the contents of the raw directory \n Please remember that files must be in .fq or .fq.gz format \n")
}

tests.get.files.unzipped <- function(){
  checkException(get.files.unzipped(NULL),"Arguments should not be NULL")
}

tests.get.files.zipped <- function(){
  checkException(get.files.zipped(NULL),"Arguments should not be NULL")
}

tests.prepare.unzipped.file.names <- function(){
  checkEquals("test.fastq",prepare.unzipped.file.names("test.fastq"))
  checkEquals("test.fq",prepare.unzipped.file.names("test.fq"))
  checkException(get.files.unzipped(NULL),"Arguments should not be NULL")
}

tests.prepare.zipped.file.names <- function(){
  checkEquals("test.fastq",prepare.zipped.file.names("test.fastq.gz"))
  checkEquals("test.fq",prepare.zipped.file.names("test.fq.gz"))
  checkException(get.files.zipped(NULL),"Arguments should not be NULL")
}

tests.copy.files.to.dest.unzipped <- function(){
  checkEquals("cp /home/data/*.fq /home/runs/ ; cp /home/data/*.fastq /home/runs/",copy.files.to.dest.unzipped("/home/data/","/home/runs/"))
  checkException(copy.files.to.dest.unzipped(NULL,NULL),"Arguments should not be NULL")
  checkException(copy.files.to.dest.unzipped("a",NULL),"Arguments should not be NULL")
  checkException(copy.files.to.dest.unzipped(NULL,"a"),"Arguments should not be NULL")
}

tests.copy.files.to.dest.zipped <- function(){
  checkEquals("cp /home/data/*.gz /home/runs/",copy.files.to.dest.zipped("/home/data/","/home/runs/"))
  checkException(copy.files.to.dest.zipped(NULL,NULL),"Arguments should not be NULL")
  checkException(copy.files.to.dest.zipped("a",NULL),"Arguments should not be NULL")
  checkException(copy.files.to.dest.zipped(NULL,"a"),"Arguments should not be NULL")
}

tests.chunk.data.files.presplit <- function(){
  test1.1compare <- c(2)
  test1.2compare <- c(4)
  test1.3compare <- c(6)
  test1.input <- rbind("1","2","3","4","5","6","7")
  checkTrue(all(test1.1compare==unlist(chunk.data.files.presplit(test1.input,3)[1])))
  checkTrue(all(test1.2compare==unlist(chunk.data.files.presplit(test1.input,3)[2])))
  checkTrue(all(test1.3compare==unlist(chunk.data.files.presplit(test1.input,3)[3])))
  checkTrue(!is.null(unlist(chunk.data.files.presplit(test1.input,3)[2])))
  checkTrue(is.null(unlist(chunk.data.files.presplit(test1.input,3)[4])))
  
  test2.1compare <- c(2,4)
  test2.2compare <- c(6,8)
  test2.input <- rbind("1","2","3","4","5","6","7","8")
  checkTrue(all(test2.1compare==unlist(chunk.data.files.presplit(test2.input,2)[1])))
  checkTrue(all(test2.2compare==unlist(chunk.data.files.presplit(test2.input,2)[2])))
  checkTrue(!is.null(unlist(chunk.data.files.presplit(test2.input,2)[2])))
  checkTrue(is.null(unlist(chunk.data.files.presplit(test2.input,2)[3])))
  
  test3.1compare <- c(2)
  test3.2compare <- c(4)
  test3.input <- rbind("1","2","3","4")
  checkTrue(all(test3.1compare==unlist(chunk.data.files.presplit(test3.input,5)[1])))
  checkTrue(all(test3.2compare==unlist(chunk.data.files.presplit(test3.input,5)[2])))
  checkTrue(!is.null(unlist(chunk.data.files.presplit(test3.input,5)[2])))
  checkTrue(is.null(unlist(chunk.data.files.presplit(test3.input,5)[3])))
  
  
  checkException(chunk.data.files.unsplit(cbind("test1","test2"),NULL),"Arguments should not be NULL")
  checkException(chunk.data.files.unsplit(NULL, 3),"Arguments should not be NULL")
  checkException(chunk.data.files.unsplit(NULL,NULL),"Arguments should not be NULL")  
}

tests.chunk.data.files.unsplit <- function(){
  test1.1compare <- c(1,2)
  test1.2compare <- c(3,4,5)
  test1.3compare <- c(6,7)
  test1.input <- rbind("1","2","3","4","5","6","7")
  checkTrue(all(test1.1compare==unlist(chunk.data.files.unsplit(test1.input,3)[1])))
  checkTrue(all(test1.2compare==unlist(chunk.data.files.unsplit(test1.input,3)[2])))
  checkTrue(all(test1.3compare==unlist(chunk.data.files.unsplit(test1.input,3)[3])))
  
  test2.1compare <- c(1,2,3,4)
  test2.2compare <- c(5,6,7,8)
  test2.input <- rbind("1","2","3","4","5","6","7","8")
  checkTrue(all(test2.1compare==unlist(chunk.data.files.unsplit(test2.input,2)[1])))
  checkTrue(all(test2.2compare==unlist(chunk.data.files.unsplit(test2.input,2)[2])))
  
  test3.1compare <- c(1)
  test3.2compare <- c(2)
  test3.3compare <- c(3)
  test3.4compare <- c(4)
  test3.input <- rbind("1","2","3","4")
  checkTrue(all(test3.1compare==unlist(chunk.data.files.unsplit(test3.input,5)[1])))
  checkTrue(all(test3.2compare==unlist(chunk.data.files.unsplit(test3.input,5)[2])))
  checkTrue(all(test3.3compare==unlist(chunk.data.files.unsplit(test3.input,5)[3])))
  checkTrue(all(test3.4compare==unlist(chunk.data.files.unsplit(test3.input,5)[4])))
  
  checkException(chunk.data.files.unsplit(cbind("test1","test2"),NULL),"Arguments should not be NULL")
  checkException(chunk.data.files.unsplit(NULL, 3),"Arguments should not be NULL")
  checkException(chunk.data.files.unsplit(NULL,NULL),"Arguments should not be NULL")  
}

tests.copy.artificial.fq <- function(){
  checkEquals("cp /home/data/art.fq /home/runs/",copy.artificial.fq("/home/data/","art.fq","/home/runs/"))
  checkEquals("cp /home/data/art.fq /home/runs/",copy.artificial.fq("/home/data","art.fq","/home/runs/"))
  checkEquals("cp /home/data/art.fq /home/runs/",copy.artificial.fq("/home/data/","art.fq","/home/runs"))
  checkEquals("cp /home/data/art.fq /home/runs/",copy.artificial.fq("/home/data","art.fq","/home/runs"))
  checkException(copy.artificial.fq(NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(copy.artificial.fq("a","a",NULL),"Arguments should not be NULL")  
  checkException(copy.artificial.fq(NULL,"a","a"),"Arguments should not be NULL")  
  checkException(copy.artificial.fq("a",NULL,"a"),"Arguments should not be NULL")  
}

tests.get.file.base <- function(){
  checkEquals("test",get.file.base("test.fq"))
  checkEquals("test",get.file.base("test.fq.gz"))
  checkEquals("test",get.file.base("test.fastq"))
  checkEquals("test",get.file.base("test.fastq.gz"))
  checkEquals("fastRun",get.file.base("fastRun.fastq"))
  checkEquals("fqastRun",get.file.base("fqastRun.fq"))
  checkEquals("fqgzastRun",get.file.base("fqgzastRun.fq.gz"))
  checkEquals("fqgzastRun",get.file.base("fqgzastRun"))
  checkEquals("fqgzastRun",get.file.base("fqgzastRun.1.fq.gz"))
  checkEquals("fqgzastRun",get.file.base("fqgzastRun.1.fq"))
  checkEquals("test.fastq",get.file.base("test.fastq.2.fq"))
  checkException(get.file.base(NULL),"Arguments should not be NULL")
}

tests.copy.file <- function(){
  checkEquals("cp /home/data/test.fq /home/runs/",copy.file("/home/data/","test.fq","/home/runs/")) 
  checkEquals("cp /home/data/test.fq /home/runs/",copy.file("/home/data","test.fq","/home/runs/")) 
  checkEquals("cp /home/data/test.fq /home/runs/",copy.file("/home/data/","test.fq","/home/runs")) 
  checkEquals("cp /home/data/test.fq /home/runs/",copy.file("/home/data","test.fq","/home/runs"))
  checkEquals("cp /home/data/test.fq.gz /home/runs/",copy.file("/home/data/","test.fq.gz","/home/runs/")) 
  checkEquals("cp /home/data/test.fq.gz /home/runs/",copy.file("/home/data","test.fq.gz","/home/runs/")) 
  checkEquals("cp /home/data/test.fq.gz /home/runs/",copy.file("/home/data/","test.fq.gz","/home/runs")) 
  checkEquals("cp /home/data/test.fq.gz /home/runs/",copy.file("/home/data","test.fq.gz","/home/runs"))
  checkException(tests.copy.file(NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(tests.copy.file("a","a",NULL),"Arguments should not be NULL")
  checkException(tests.copy.file("a",NULL,"a"),"Arguments should not be NULL")
  checkException(tests.copy.file(NULL,"a","a"),"Arguments should not be NULL")
}

tests.unzip.gz.files <- function(){
  checkEquals("gunzip test.fq.gz",unzip.gz.files("test.fq.gz"))
  checkException(unzip.gz.files(NULL),"Arguments should not be NULL") 
}

tests.first.read.name <- function(){
  checkEquals("test.1.fq",first.read.name("test"))
  checkException(first.read.name(NULL),"Arguments should not be NULL")
}

tests.second.read.name <- function(){
  checkEquals("test.2.fq",second.read.name("test"))
  checkException(second.read.name(NULL),"Arguments should not be NULL")
}

tests.left.dirty.name <- function(){
  checkEquals("dirty.test.1.fq",left.dirty.name("test"))
  checkException(left.dirty.name(NULL),"Arguments should not be NULL")
}

tests.right.dirty.name <- function(){
  checkEquals("dirty.test.2.fq",right.dirty.name("test"))
  checkException(right.dirty.name(NULL),"Arguments should not be NULL")
}

tests.trimAssemble.PE <- function(){
  trimPath <- "/home/trim.jar"
  qScore <- "phred33"
  fqFile <- "result.fq"
  headcrop <- 12
  artificial.fq <- "test.art"
  trimMin <- 2
  
  test1.trimCommand <- paste("ILLUMINACLIP:test.art:2:30:10 HEADCROP:12 SLIDINGWINDOW:3:30 MINLEN:2")
  test1.command <- paste("java -jar /home/trim.jar PE -threads 20 -phred33 -trimlog result.1.fq.pairedtrim.log result.1.fq result.2.fq p.result.1.fq u.result.1.fq p.result.2.fq u.result.2.fq",test1.trimCommand)
  checkEquals(test1.command,trimAssemble.PE(fqFile,trimPath,qScore,headcrop,artificial.fq,trimMin))
  
  trimPath <- "/home/trim.jar/"
  checkException(trimAssemble.PE(fqFile,trimPath,qScore,headcrop,artificial.fq,trimMin),"Invalid path to Trimmomatic JAR file")
  
  checkException(trimAssemble.PE(NULL,NULL,NULL,NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(trimAssemble.PE("NULL","NULL","NULL",NULL,"NULL",NULL),"Arguments should not be NULL")
  checkException(trimAssemble.PE(NULL,"NULL",NULL,"NULL",NULL,NULL),"Arguments should not be NULL")
}

tests.split.unsplit.files.PE <- function(){
  checkEquals("cat /home/runs/test.fq | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[1234]/; @i = 0 if @i == 8' > test.1.fq && cat /home/runs/test.fq | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[5678]/; @i = 0 if @i == 8' > test.2.fq",split.unsplit.files.PE("/home/runs/","test.fq"))
  checkEquals("cat /home/runs/test.fq | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[1234]/; @i = 0 if @i == 8' > test.1.fq && cat /home/runs/test.fq | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[5678]/; @i = 0 if @i == 8' > test.2.fq",split.unsplit.files.PE("/home/runs","test.fq"))
  checkException(split.unsplit.files.PE(NULL,NULL),"Arguments should not be NULL")
  checkException(split.unsplit.files.PE("a",NULL),"Arguments should not be NULL")
  checkException(split.unsplit.files.PE(NULL,"a"),"Arguments should not be NULL")
}

tests.file.shuffle.PE <- function(){
  checkEquals("mv test.1.fq dirty.test.1.fq && mv test.2.fq dirty.test.2.fq && mv p.test.1.fq test.1.fq && mv p.test.2.fq test.2.fq",file.shuffle.PE("test.fq"))
  checkException(file.shuffle.PE(NULL),"Arguments should not be NULL")
}

tests.preQualityCheck.PE <- function(){
  checkEquals("/home/fastqcpath dirty.test.1.fq dirty.test.2.fq",preQualityCheck.PE("/home/fastqcpath","test.fq"))
  checkException(preQualityCheck.PE(NULL,NULL),"Arguments should not be NULL")
  checkException(preQualityCheck.PE("a",NULL),"Arguments should not be NULL")
  checkException(preQualityCheck.PE(NULL,"a"),"Arguments should not be NULL")
}

tests.postQualityCheck.PE <- function(){
  checkEquals("/home/fastqcpath test.1.fq test.2.fq",postQualityCheck.PE("/home/fastqcpath","test.fq"))
  checkException(postQualityCheck.PE(NULL,NULL),"Arguments should not be NULL")
  checkException(postQualityCheck.PE("a",NULL),"Arguments should not be NULL")
  checkException(postQualityCheck.PE(NULL,"a"),"Arguments should not be NULL")
}

tests.remove.unneeded.files <- function(){
  checkEquals("rm dirty.test.1.fq ; rm dirty.test.2.fq ; rm u.test.1.fq ; rm u.test.2.fq",remove.unneeded.files("test.fq"))
  checkException(remove.unneeded.files(NULL),"Arguments should not be NULL")
}

tests.dirty.name.SE <- function(){
  checkEquals("dirty.test.fq",dirty.name.SE("test"))
  checkException(dirty.name.SE(NULL),"Arguments should not be NULL")
}

tests.final.name.SE <- function(){
  checkEquals("test.fq",final.name.SE("test"))
  checkException(final.name.SE(NULL),"Arguments should not be NULL")
}

tests.trimAssemble.SE <- function(){
  trimPath <- "/home/trim.jar"
  qScore <- "phred33"
  leftDirtyName <- "dirty.l"
  headcrop <- 12
  artificial.fq <- "test.art"
  trimMin <- 2
  
  test1.trimCommand <- paste("ILLUMINACLIP:test.art:2:30:10 HEADCROP:12 SLIDINGWINDOW:3:30 MINLEN:2")
  test1.command <- paste("java -jar /home/trim.jar SE -threads 20 -phred33 -trimlog dirty.l.pairedtrim.log dirty.l p.dirty.l.fq", test1.trimCommand)
  checkEquals(test1.command,trimAssemble.SE(leftDirtyName,trimPath,qScore,headcrop,artificial.fq,trimMin))
  
  trimPath <- "/home/trim.jar/"
  checkException(trimAssemble.SE(leftDirtyName,trimPath,qScore,headcrop,artificial.fq,trimMin),"Invalid path to Trimmomatic JAR file")
  
  checkException(trimAssemble.PE(NULL,NULL,,NULL,NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(trimAssemble.PE("NULL","NULL","NULL",NULL,"NULL",NULL),"Arguments should not be NULL")
  checkException(trimAssemble.PE(NULL,"NULL",NULL,"NULL",NULL,NULL),"Arguments should not be NULL")
}

tests.file.shuffle.SE <- function(){
  checkEquals("mv test.fq dirty.test.fq  &&  mv p.test.fq test.fq",file.shuffle.SE("test.fq"))
  checkException(file.shuffle.SE(NULL),"Arguments should not be NULL")
}

tests.preQualityCheck.SE <- function(){
  checkEquals("/home/fastqc dirty.test.fq",preQualityCheck.SE("/home/fastqc","test.fq"))
  checkException(preQualityCheck.SE(NULL,NULL),"Arguments should not be NULL")
}

tests.postQualityCheck.SE <- function(){
  checkEquals("/home/fastqc test.fq",preQualityCheck.SE("/home/fastqc","test.fq"))
  checkException(postQualityCheck.SE(NULL,NULL),"Arguments should not be NULL")
}

remove.unneeded.files.SE <- function(){
  checkEquals("rm dirty.test.fq ; rm u.test.fq",remove.unneeded.files("test.fq"))
  checkException(remove.unneeded.files(NULL),"Arguments should not be NULL")
}