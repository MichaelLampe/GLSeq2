library('RUnit')

source("../Main/GLSeq.Top.Functions.R")

test.suite <- defineTestSuite("Tests.GLSeq.top",dirs = file.path("../Tests/test"),testFileRegexp = '.R$')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)