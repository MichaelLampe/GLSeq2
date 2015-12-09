__author__ = 'mlampe'

from sys import argv
from subprocess import Popen
from subprocess import PIPE
from sys import exit
import os
# This is a wrapper for the rockhopper program so that multichromosomal organisms are just as easy to run.

# Takes in a bunch of args from the Rscript that can be used to construct a rockhopper call
rock_path = argv[1].strip("\n")
destination_directory = argv[2].strip("\n")
rockhopper_files = argv[3].strip("\n")
unique_output_directory = argv[4].strip("\n")
args = argv[5].strip("\n")
# Gets the result of rhFnaConverter 
referenceFilesDir = destination_directory+"ReferenceGenome" + "/"
directories = ""
folders = os.listdir(referenceFilesDir)

for i, folder in enumerate(folders):
    folder = referenceFilesDir + folder
    if i == 0:
        directories = folder
    else:
        directories = directories + "," + folder

#programArgs = ""
#if (paired == "TRUE"):
#  if (libstrand == "NULL"):
#    programArgs = "-s false"
#  elif (libstrand == "F"):
#    programArgs = "-fr -s true"
#  elif (libstrand == "R"):
#    programArgs = "-rf -s true"
#else:
#  if (libstrand == "NULL"):
#    programArgs = "-s false"
#  elif (libstrand == "F"):
#    programArgs = "-c false"
#  elif (libstrand == "R"):
#    programArgs = "-c true"

programArgs = args.replace("?"," ")
# Uses that output in addition to the Rscripts message to generate a command
rock_options = " -o " + unique_output_directory + " -v true -SAM -TIME " + programArgs
rock_align = "java -Xmx4096m -cp " + rock_path + " Rockhopper "+ "-g " + directories + rock_options + " \"" + rockhopper_files + "\""
print(rock_align)
process = Popen(rock_align,shell=True)
out,err = process.communicate()

# Try to kill, otherwise assume dead
try:process.kill()
except:pass
# Exits with good error code
exit(0)