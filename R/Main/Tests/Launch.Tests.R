library('RUnit')

source("GLSeq.Top.Functions.R")
source("GLSeq.Dataprep.Functions.R")

test.suite <- defineTestSuite("Tests.GLSeq.top",dirs = file.path("Tests/test"),testFileRegexp = '.R$')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)
