__author__ = 'mlampe'

from sys import argv
from subprocess import Popen
from subprocess import PIPE
from sys import exit
# This is a wrapper for the rockhopper program so that multichromosomal organisms are just as easy to run.

# Takes in a bunch of args from the Rscript that can be used to construct a rockhopper call
rock_path = argv[1].strip("\n")
destination_directory = argv[2].strip("\n")
rockhopper_files = argv[3].strip("\n")
unique_output_directory = argv[4].strip("\n")
# Gets the result of rhFnaConverter
referenceFiles = destination_directory+"referenceGenomes"
reference_genome_command = "cat " + referenceFiles
output = Popen(reference_genome_command,shell=True,stdout=PIPE)
directories, err = output.communicate()
print(directories)
print(err)
# Remove any new lines in the file that could cause problems
directories = directories.strip("\n")
# Try to kill, otherwise assume dead
try:output.kill()
except:pass
# Uses that output in addition to the Rscripts message to generate a command
rock_options = " -o " + unique_output_directory + " -e false -SAM -TIME"
rock_align = "java -Xmx1200m -cp " + rock_path + " Rockhopper "+ "-g " + directories + rock_options + " \"" + rockhopper_files + "\""
print(rock_align)
process = Popen(rock_align,shell=True)
out,err = process.communicate()

# Try to kill, otherwise assume dead
try:process.kill()
except:pass
# Exits with good error code
exit(0)