__author__ = 'mlampe'

import sys
import re
import os
from subprocess import Popen

def checkTrail(directory):
    pattern = "/$"
    check = re.search(pattern,directory)
    if (check == None):
        return (directory + "/")
    else:
        return (directory)
#if (len(sys.argv) <= 1):
   # sys.exit("Please add the appropriate command line arguments. The argument the directory containing the files and the name of the you wish to convert")

def getFileName(line,directory):
    name = line.replace(">","",1)
    pattern = "^[a-zA-Z0-9_.]+"
    name = re.search(pattern,name)
    return (directory + name.group(0) + ".fna")

def getDirName(line,directory):
    name = line.replace(">","",1)
    pattern = "^[a-zA-Z0-9_.]+"
    name = re.search(pattern,name)
    return checkTrail(directory + name.group(0))

def writeConvertedFile(file_name,header,dna_sequence):
    with open(new_file_name,"w") as create_file:
        print("Processed subsequence:" + header)
        create_file.write(header)
        create_file.write(dna_sequence)

directory = sys.argv[1]
raw_file_name = sys.argv[2]

directory = checkTrail(directory)

command_arg = ""
count = 0
file = directory + raw_file_name
with open(file,"r+b") as raw_file:
    new_file_dir = ""
    new_file_name = ""
    header = ""
    dna_sequence = ""
    for line in raw_file:
        if ">" in line:
            # Don't run this write the first time
            if header != "":
                writeConvertedFile(new_file_name,header,dna_sequence)
            line = line.replace("\n","")
            line = line.replace(">","")
            header = ">" + line + "\n"
            new_file_dir = getDirName(line,directory)
            try:
              os.mkdir(new_file_dir)
            # Assume this means that the directory is already there.
            except:
              pass
            if (command_arg == ""):
                command_arg = new_file_dir
            else:
                command_arg = command_arg + "," + new_file_dir
            new_file_name = getFileName(line,new_file_dir)
            dna_sequence = ""
            count = count + 1
        else:
            dna_sequence = dna_sequence + line
# Write the final file in
if (dna_sequence != ""):
    writeConvertedFile(new_file_name,header,dna_sequence)
    count = count + 1
    if (command_arg == ""):
        command_arg = new_file_dir
    else:
        command_arg = command_arg + "," + new_file_dir
os.remove(file)
create_formatted_file = "echo " + command_arg + " > referenceGenomes.txt"
Popen(create_formatted_file,shell=True)
sys.exit(0)
