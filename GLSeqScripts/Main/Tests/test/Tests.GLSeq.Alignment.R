tests.copy.genome <- function(){
  checkEquals("cd /home/base/rGenome && cp refFasta /home/dest/dir/",copy.genome("/home/base/","rGenome","refFasta","/home/dest/dir/"))
  checkEquals("cd /home/base/rGenome && cp refFasta /home/dest/dir/",copy.genome("/home/base","rGenome","refFasta","/home/dest/dir/"))
  checkEquals("cd /home/base/rGenome && cp refFasta /home/dest/dir/",copy.genome("/home/base/","rGenome","refFasta","/home/dest/dir"))
  checkEquals("cd /home/base/rGenome && cp refFasta /home/dest/dir/",copy.genome("/home/base","rGenome","refFasta","/home/dest/dir"))
  checkException(copy.genome(NULL,NULL,NULL,NULL),"Arguments should not be NULL")
  checkException(copy.genome("NULL","NULL","NULL",NULL),"Arguments should not be NULL")
  checkException(copy.genome("NULL","NULL",NULL,"NULL"),"Arguments should not be NULL")
  checkException(copy.genome("NULL",NULL,"NULL","NULL"),"Arguments should not be NULL")
  checkException(copy.genome(NULL,"NULL","NULL","NULL"),"Arguments should not be NULL")
}

tests.countable.sam.name <- function(){
  checkEquals("test.name.countable.sam",countable.sam.name("test.name"))
  checkException(countable.sam.name(NULL),"Arguments should not be NULL")
}

tests.assign.resName <- function(){
  checkEquals("test.add",assign.resName("test","add"))
  checkException(assign.resName(NULL,NULL),"Arguments should not be NULL")
  checkException(assign.resName("NULL",NULL),"Arguments should not be NULL")
  checkException(assign.resName(NULL,"NULL"),"Arguments should not be NULL")
}

tests.assign.name <- function(){
  checkEquals("test",assign.name("test.fq",FALSE))
  checkEquals("test",assign.name("test.1.fq",TRUE))
  checkEquals("test",assign.name("test.2.fq",TRUE))
  checkEquals("test",assign.name("test_1.fq",TRUE))
  checkEquals("test",assign.name("test_2.fq",TRUE))
  checkEquals("test.1",assign.name("test.1.fq",FALSE))
  checkEquals("test.2",assign.name("test.2.fq",FALSE))
  checkEquals("test_1",assign.name("test_1.fq",FALSE))
  checkEquals("test_2",assign.name("test_2.fq",FALSE))
  checkException(assign.name(NULL,NULL),"Arguments should not be NULL")
  checkException(assign.name("NULL",NULL),"Arguments should not be NULL")
  checkException(assign.name(NULL,"NULL"),"Arguments should not be NULL")
  checkException(assign.name("OK","OK"),"Paired end must be a LOGICAL")
}