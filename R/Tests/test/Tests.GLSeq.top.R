test.create.text.add <- function(){
  # Functions takes in two args
  # ExpID = Experiment Identification String
  # runAttempt = Integer number of runs done with that ExpId previously
  checkEquals("MyTextRun.1",create.text.add("MyTextRun",1))
  checkEquals("My_Text_Run.1",create.text.add("My Text Run",1))
  checkEquals("My_Text_Run.1",create.text.add("My%Text?Run",1))
  checkEquals("My____Text___Run.1",create.text.add("My,,,,Text,,^Run",1))
  checkEquals("My__Text_Run.1",create.text.add("My\n Text\nRun",1))
  checkEquals("6.5",create.text.add(6,5))
  checkException(create.text.add(NULL,NULL),"Arguments should not be NULL")
  checkException(create.text.add("a",NULL),"Arguments should not be NULL")
  checkException(create.text.add(NULL,"a"),"Arguments should not be NULL")
}

test.create.dest.dir <- function(){
  # Function takes in two args
  # Dest.Dir.Base = Directory that the overall folder will be created in
  # Text.add = The unique label of the directory to be created.
  # Util gives the trailDirCheck to us.
  checkEquals("/home/runs/MyTextRun.1/",create.dest.dir("/home/runs","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1/",create.dest.dir("/home/runs/","MyTextRun.1"))  
  checkException(create.dest.dir(NULL,NULL),"Arguments should not be NULL")
  checkException(create.dest.dir("a",NULL),"Arguments should not be NULL")
  checkException(create.dest.dir(NULL,"a"),"Arguments should not be NULL")
}

test.create.dest.dir.log <- function(dest.dir,text.add){
  # Function takes in two args
  # Dest.Dir = The folder containing all the run data.
  # Text.add = The unique label of the directory to be created.
  checkEquals("/home/runs/MyTextRun.1/MyTextRun.1.stat/",create.dest.dir.log("/home/runs/MyTextRun.1/","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1/MyTextRun.1.stat/",create.dest.dir.log("/home/runs/MyTextRun.1","MyTextRun.1"))
  checkException(create.dest.dir.log(NULL,NULL),"Arguments should not be NULL")
  checkException(create.dest.dir.log("a",NULL),"Arguments should not be NULL")
  checkException(create.dest.dir.log(NULL,"a"),"Arguments should not be NULL")
}

test.instantiate.logs <- function(){
  # Function takes in three args
  # DestDirLog = The directory where the log file should be written
  # Dest.dir = The destination directory as a backup if a unique file is not provided.
  # Text.add = The unique label of the file to be created
  # 
  checkEquals("/home/runs/MyTextRun.1/MyTextRun.1.RunLog.txt",instantiate.logs("/home/","MyTextRun.1","/home/runs/MyTextRun.1/"))
  checkEquals("/home/runs/MyTextRun.1/MyTextRun.1.RunLog.txt",instantiate.logs("/home/","MyTextRun.1","/home/runs/MyTextRun.1"))
  checkEquals("/home/MyTextRun.1.RunLog.txt",instantiate.logs("/home/","MyTextRun.1",NULL))
  checkEquals("/home/MyTextRun.1.RunLog.txt",instantiate.logs("/home","MyTextRun.1",NULL))
  checkEquals("/home/MyTextRun.1.RunLog.txt",instantiate.logs("/home/","MyTextRun.1"))
  checkEquals("/home/MyTextRun.1.RunLog.txt",instantiate.logs("/home","MyTextRun.1"))
  checkException(create.dest.dir.log(NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(create.dest.dir.log("a",NULL,NULL),"Arguments should not be NULL")
  checkException(create.dest.dir.log(NULL,"a",NULL),"Arguments should not be NULL")
  checkException(create.dest.dir.log("a",NULL,"a"),"Arguments should not be NULL")
}

test.create.run.directory <- function(){
  # Function takes in two args
  # Dest.dir = The destination directory of the folder
  # (Optional) = Log.file location
  checkEquals("mkdir /home/",create.run.directory("/home/"))
  checkEquals("mkdir /home/",create.run.directory("/home/"))
  checkEquals("mkdir /home/",create.run.directory("/home/","ok"))
  checkEquals("mkdir /home/",create.run.directory("/home/","ok"))
  checkEquals("mkdir /home/",create.run.directory("/home/",NULL))
  checkEquals("mkdir /home/",create.run.directory("/home/",NULL))
  checkException(create.run.directory(NULL),"Arguments should not be NULL")
}

test.create.log.directory <- function(){
  # Function takes in two args
  # DestDirLog = Log directory
  # (Optional) = Log file location
  checkEquals("mkdir /home/",create.run.directory("/home/"))
  checkEquals("mkdir /home/",create.run.directory("/home/"))
  checkEquals("mkdir /home/",create.run.directory("/home/","ok"))
  checkEquals("mkdir /home/",create.run.directory("/home/","ok"))
  checkEquals("mkdir /home/",create.run.directory("/home/",NULL))
  checkEquals("mkdir /home/",create.run.directory("/home/",NULL))
  checkException(create.run.directory(NULL),"Arguments should not be NULL")
}

test.create.HtSeq.folder <- function() {
  # Function takes in three args
  # Dest.Dir = Where the folder should be created
  # Text.Add = Unique name
  #(Optional) = Log File Location
  checkEquals("/home/runs/MyTextRun.1.HTSeq.Counting/",create.HtSeq.folder("/home/runs","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.HTSeq.Counting/",create.HtSeq.folder("/home/runs/","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.HTSeq.Counting/",create.HtSeq.folder("/home/runs","MyTextRun.1",NULL))
  checkEquals("/home/runs/MyTextRun.1.HTSeq.Counting/",create.HtSeq.folder("/home/runs/","MyTextRun.1","a"))
  checkException(create.HtSeq.folder(NULL,NULL),"Arguments should not be NULL")
  checkException(create.HtSeq.folder(NULL,"a"),"Arguments should not be NULL")
  checkException(create.HtSeq.folder("a",NULL),"Arguments should not be NULL")
}

test.create.FeatureCounts.folder <- function() {
  # Function takes in three args
  # Dest.Dir = Where the folder should be created
  # Text.Add = Unique name
  #(Optional) = Log File Location
  checkEquals("/home/runs/MyTextRun.1.FeatureCounts.Counting/",create.FeatureCounts.folder("/home/runs","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.FeatureCounts.Counting/",create.FeatureCounts.folder("/home/runs/","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.FeatureCounts.Counting/",create.FeatureCounts.folder("/home/runs","MyTextRun.1",NULL))
  checkEquals("/home/runs/MyTextRun.1.FeatureCounts.Counting/",create.FeatureCounts.folder("/home/runs/","MyTextRun.1","a"))
  checkException(create.FeatureCounts.folder(NULL,NULL),"Arguments should not be NULL")
  checkException(create.FeatureCounts.folder(NULL,"a"),"Arguments should not be NULL")
  checkException(create.FeatureCounts.folder("a",NULL),"Arguments should not be NULL")
}

test.create.RSEM.folder <- function() {
  # Function takes in three args
  # Dest.Dir = Where the folder should be created
  # Text.Add = Unique name
  #(Optional) = Log File Location
  checkEquals("/home/runs/MyTextRun.1.RSEM.Counting/",create.RSEM.folder("/home/runs","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.RSEM.Counting/",create.RSEM.folder("/home/runs/","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.RSEM.Counting/",create.RSEM.folder("/home/runs","MyTextRun.1",NULL))
  checkEquals("/home/runs/MyTextRun.1.RSEM.Counting/",create.RSEM.folder("/home/runs/","MyTextRun.1","a"))
  checkException(create.RSEM.folder(NULL,NULL),"Arguments should not be NULL")
  checkException(create.RSEM.folder(NULL,"a"),"Arguments should not be NULL")
  checkException(create.RSEM.folder("a",NULL),"Arguments should not be NULL")
}

test.create.Cufflinks.folder <- function() {
  # Function takes in three args
  # Dest.Dir = Where the folder should be created
  # Text.Add = Unique name
  #(Optional) = Log File Location
  checkEquals("/home/runs/MyTextRun.1.Cufflinks.Counting/",create.Cufflinks.folder("/home/runs","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.Cufflinks.Counting/",create.Cufflinks.folder("/home/runs/","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.Cufflinks.Counting/",create.Cufflinks.folder("/home/runs","MyTextRun.1",NULL))
  checkEquals("/home/runs/MyTextRun.1.Cufflinks.Counting/",create.Cufflinks.folder("/home/runs/","MyTextRun.1","a"))
  checkException(create.Cufflinks.folder(NULL,NULL),"Arguments should not be NULL")
  checkException(create.Cufflinks.folder(NULL,"a"),"Arguments should not be NULL")
  checkException(create.Cufflinks.folder("a",NULL),"Arguments should not be NULL")
}

test.create.Collect.folder <- function() {
  # Function takes in three args
  # Dest.Dir = Where the folder should be created
  # Text.Add = Unique name
  #(Optional) = Log File Location
  checkEquals("/home/runs/MyTextRun.1.Collect/",create.Collect.folder("/home/runs","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.Collect/",create.Collect.folder("/home/runs/","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.Collect/",create.Collect.folder("/home/runs","MyTextRun.1",NULL))
  checkEquals("/home/runs/MyTextRun.1.Collect/",create.Collect.folder("/home/runs/","MyTextRun.1","a"))
  checkException(create.Collect.folder(NULL,NULL),"Arguments should not be NULL")
  checkException(create.Collect.folder(NULL,"a"),"Arguments should not be NULL")
  checkException(create.Collect.folder("a",NULL),"Arguments should not be NULL")
}

test.create.oldRun.Collect.folder <- function() {
  # Function takes in three args
  # Dest.Dir = Where the folder should be created
  # Text.Add = Unique name
  #(Optional) = Log File Location
  checkEquals("/home/runs/MyTextRun.1.Collect/",create.oldRun.Collect.folder("/home/runs","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.Collect/",create.oldRun.Collect.folder("/home/runs/","MyTextRun.1"))
  checkEquals("/home/runs/MyTextRun.1.Collect/",create.oldRun.Collect.folder("/home/runs","MyTextRun.1",NULL))
  checkEquals("/home/runs/MyTextRun.1.Collect/",create.oldRun.Collect.folder("/home/runs/","MyTextRun.1","a"))
  checkException(create.oldRun.Collect.folder(NULL,NULL),"Arguments should not be NULL")
  checkException(create.oldRun.Collect.folder(NULL,"a"),"Arguments should not be NULL")
  checkException(create.oldRun.Collect.folder("a",NULL),"Arguments should not be NULL")
}

test.copy.attribute.file.to.dest <- function(){
  # Function takes in three args
  # Attribute file path
  # Destination for that attribute file
  # Log file location (Optional)
  checkEquals("cp /home/storage/attrFile.R /home/runs/",copy.attribute.file.to.dest("/home/storage/attrFile.R","/home/runs/"))
  checkEquals("cp /home/storage/attrFile.R /home/runs/",copy.attribute.file.to.dest("/home/storage/attrFile.R","/home/runs"))
  checkEquals("cp /home/storage/attrFile.R /home/runs/",copy.attribute.file.to.dest("/home/storage/attrFile.R","/home/runs/","a"))
  checkEquals("cp /home/storage/attrFile.R /home/runs/",copy.attribute.file.to.dest("/home/storage/attrFile.R","/home/runs","a"))
  checkException(copy.attribute.file.to.dest("a",NULL),"Arguments should not be NULL")
  checkException(copy.attribute.file.to.dest(NULL,"a"),"Arguments should not be NULL")
  checkException(copy.attribute.file.to.dest(NULL,NULL),"Arguments should not be NULL")
}

test.fqfiles.table.pe.assemble <- function(){
  # Function takes in one arg
  # Table of all the .fq files.
  test1 <- NULL
  test1 <- rbind(test1,c("test.1.fq","test.2.fq"))
  checkEqualsNumeric(dim(test1),dim(fqfiles.table.pe.assemble(c("test"))))
  
  test2 <- NULL
  test2 <- rbind(test2,c("test.1.fq","test.2.fq"))
  test2 <- rbind(test2,c("test2.1.fq","test2.2.fq"))
  checkEqualsNumeric(dim(test2),dim(fqfiles.table.pe.assemble(c("test","test2"))))
  
  checkException(fqfiles.table.pe.assemble(NULL),"Arguments should not be NULL")  
}

test.rename.preprocessed.files <- function(){
  checkException(rename.preprocessed.files(NULL),"Arguments should not be NULL")
}

test.copy.preprocessed.files <- function() {
  checkException(copy.preprocessed.files(NULL,"a"),"Arguments should not be NULL")
  checkException(copy.preprocessed.files("a",NULL),"Arguments should not be NULL")
  checkException(copy.preprocessed.files(NULL,NULL),"Arguments should not be NULL")
}

test.convert.file.list.to.table <- function(){
  checkException(convert.file.list.to.table(NULL,"a"),"Arguments should not be NULL")
  checkException(convert.file.list.to.table(TRUE,NULL),"Arguments should not be NULL")
  checkException(convert.file.list.to.table(NULL,NULL),"Arguments should not be NULL")
}

test.find.files.for.dataprep <- function(){
  checkException(find.files.for.dataprep("a",NULL),"Arguments should not be NULL")
  checkException(find.files.for.dataprep(NULL,TRUE),"Arguments should not be NULL")
  checkException(find.files.for.dataprep(NULL,NULL),"Arguments should not be NULL")
}

test.prepare.chunk.function <- function(){
  
  test1.1compare <- c(1,2)
  test1.2compare <- c(3,4,5)
  test1.3compare <- c(6,7)
  test1.input <- rbind("1","2","3","4","5","6","7")
  checkTrue(all(test1.1compare==unlist(prepare.chunk.function(test1.input,3)[1])))
  checkTrue(all(test1.2compare==unlist(prepare.chunk.function(test1.input,3)[2])))
  checkTrue(all(test1.3compare==unlist(prepare.chunk.function(test1.input,3)[3])))
  
  test2.1compare <- c(1,2,3,4)
  test2.2compare <- c(5,6,7,8)
  test2.input <- rbind("1","2","3","4","5","6","7","8")
  checkTrue(all(test2.1compare==unlist(prepare.chunk.function(test2.input,2)[1])))
  checkTrue(all(test2.2compare==unlist(prepare.chunk.function(test2.input,2)[2])))
  
  test3.1compare <- c(1)
  test3.2compare <- c(2)
  test3.3compare <- c(3)
  test3.4compare <- c(4)
  test3.input <- rbind("1","2","3","4")
  checkTrue(all(test3.1compare==unlist(prepare.chunk.function(test3.input,5)[1])))
  checkTrue(all(test3.2compare==unlist(prepare.chunk.function(test3.input,5)[2])))
  checkTrue(all(test3.3compare==unlist(prepare.chunk.function(test3.input,5)[3])))
  checkTrue(all(test3.4compare==unlist(prepare.chunk.function(test3.input,5)[4])))
  
  checkException(prepare.chunk.function(cbind("test1","test2"),NULL),"Arguments should not be NULL")
  checkException(prepare.chunk.function(NULL, 3),"Arguments should not be NULL")
  checkException(prepare.chunk.function(NULL,NULL),"Arguments should not be NULL")  
}

test.check.nStreams <- function(){
  checkEqualsNumeric(2,check.nStreams(rbind(1,2,3),2))
  checkEqualsNumeric(3,check.nStreams(rbind(1,2,3),4))
  checkEqualsNumeric(1,check.nStreams(rbind(1,2,3),1))
  
  checkException(check.nStreams(cbind("test1","test2"),NULL),"Arguments should not be NULL")
  checkException(check.nStreams(NULL, 3),"Arguments should not be NULL")
  checkException(check.nStreams(NULL,NULL),"Arguments should not be NULL")  
}

test.RSEM.finish <- function() {
  
  checkEquals(" wait && cd /home/runs/ && mv *.RSEM.counts.* /home/runs/test.RSEM.counting/ && mv *.index.* /home/runs/test.RSEM.counting/",RSEM.finish("","/home/runs/test.RSEM.counting/","/home/runs/"))
  checkEquals(" wait && cd /home/runs/ && mv *.RSEM.counts.* /home/runs/test.RSEM.counting/ && mv *.index.* /home/runs/test.RSEM.counting/",RSEM.finish("","/home/runs/test.RSEM.counting","/home/runs/"))
  checkEquals(" wait && cd /home/runs/ && mv *.RSEM.counts.* /home/runs/test.RSEM.counting/ && mv *.index.* /home/runs/test.RSEM.counting/",RSEM.finish("","/home/runs/test.RSEM.counting/","/home/runs"))
  checkEquals(" wait && cd /home/runs/ && mv *.RSEM.counts.* /home/runs/test.RSEM.counting/ && mv *.index.* /home/runs/test.RSEM.counting/",RSEM.finish("","/home/runs/test.RSEM.counting","/home/runs"))
  
  checkException(check.RSEM.finish(NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(check.RSEM.finish("a",NULL,NULL),"Arguments should not be NULL")
  checkException(check.RSEM.finish("a","a",NULL),"Arguments should not be NULL")
  checkException(check.RSEM.finish(NULL,"a",NULL),"Arguments should not be NULL")
  checkException(check.RSEM.finish(NULL,"a","a"),"Arguments should not be NULL")
  checkException(check.RSEM.finish(NULL,NULL,"a"),"Arguments should not be NULL")
}

test.run.data.prep <- function() {
  checkException(run.data.prep(NULL,NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(run.data.prep("a",NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(run.data.prep("a","a",NULL,NULL),"Arguments should not be NULL")
  checkException(run.data.prep("a","a","a",NULL),"Arguments should not be NULL")
  checkException(run.data.prep(NULL,"a","a","a"),"Arguments should not be NULL")
  checkException(run.data.prep(NULL,"a","a",NULL),"Arguments should not be NULL")
  checkException(run.data.prep("a",NULL,"a",NULL),"Arguments should not be NULL")
  checkException(run.data.prep("a",NULL,"a","a"),"Arguments should not be NULL")
  checkException(run.data.prep(NULL,NULL,"a",NULL),"Arguments should not be NULL")
  checkException(run.data.prep(NULL,NULL,"a","a"),"Arguments should not be NULL")
  checkException(run.data.prep(NULL,NULL,NULL,"a"),"Arguments should not be NULL")
}

test.cushaw.gpu.run <- function(){
  checkException(cushaw.gpu.run(NULL),"Arguments should not be NULL")
}

test.execute.comm.stack <- function(){
  checkException(execute.comm.stack(NULL),"Arguments should not be NULL")
}

test.previous.run.directories <- function(){
  checkEquals("/home/run/",previous.run.directories("/home/run/","test"))
  checkEquals("/home/run/",previous.run.directories("/home/run","test"))
  checkException(previous.run.directories(NULL,NULL),"Arguments should not be NULL")
  checkException(previous.run.directories(NULL,"a"),"Arguments should not be NULL")
  checkException(previous.run.directories("a",NULL),"Arguments should not be NULL")
}

test.previous.run.FeatureCounts <- function(){
  checkEquals("/home/run/test.FeatureCounts.counting/",previous.run.FeatureCounts("/home/run/","test"))
  checkEquals("/home/run/test.FeatureCounts.counting/",previous.run.FeatureCounts("/home/run","test"))
  checkException(previous.run.FeatureCounts(NULL,NULL),"Arguments should not be NULL")
  checkException(previous.run.FeatureCounts(NULL,"a"),"Arguments should not be NULL")
  checkException(previous.run.FeatureCounts("a",NULL),"Arguments should not be NULL")
}
test.previous.run.HTSeq <- function(){
  checkEquals("/home/run/test.HTSeq.counting/",previous.run.HTSeq("/home/run/","test"))
  checkEquals("/home/run/test.HTSeq.counting/",previous.run.HTSeq("/home/run","test"))
  checkException(previous.run.HTSeq(NULL,NULL),"Arguments should not be NULL")
  checkException(previous.run.HTSeq(NULL,"a"),"Arguments should not be NULL")
  checkException(previous.run.HTSeq("a",NULL),"Arguments should not be NULL")
}
test.previous.run.RSEM <- function(){
  checkEquals("/home/run/test.RSEM.counting/",previous.run.RSEM("/home/run/","test"))
  checkEquals("/home/run/test.RSEM.counting/",previous.run.RSEM("/home/run","test"))
  checkException(previous.run.RSEM(NULL,NULL),"Arguments should not be NULL")
  checkException(previous.run.RSEM(NULL,"a"),"Arguments should not be NULL")
  checkException(previous.run.RSEM("a",NULL),"Arguments should not be NULL")
}