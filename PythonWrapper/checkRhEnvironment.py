__author__ = 'Michael Lampe'

import sys
import shutil
import os


directoryToCheck = sys.argv[1]
if not directoryToCheck.endswith("/"):
    directoryToCheck += "/"
fna_ext = ".fna"
ptt_ext = ".ptt"

# Looks into the directory that is supposed to be checked
for directory in os.listdir(directoryToCheck):
    directory = directoryToCheck + directory
    fna, ptt = False, False

    # Goes through the files in a given directory looking for certain file types.
    for files in os.listdir(directory):
        files = directory + files
        if files.endswith(fna_ext):
            fna = True
        elif files.endswith(ptt_ext):
            ptt = True

    # Remove the directory file if it doen't contain both an FNA and a PTT file
    if not fna or not ptt:
        shutil.rmtree(directory)

exit(0)
