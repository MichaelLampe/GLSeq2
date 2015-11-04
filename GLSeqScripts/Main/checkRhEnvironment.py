__author__ = 'mlampe'

import sys
import shutil
import os

directoryToCheck = sys.argv[1]
if not directoryToCheck.endswith("/"):
    directoryToCheck += "/"
fna_ext = ".fna"
ptt_ext = ".ptt"

for directory in os.listdir(directoryToCheck):
    directory = directoryToCheck + directory
    fna = False
    ptt = False
    for files in os.listdir(directory):
        files = directory + files
        if files.endswith(fna_ext):
            fna = True
        elif files.endswith(ptt_ext):
            ptt = True

    # Remove if doesn't contain both ptt and fna file.
    if not fna or not ptt:
        shutil.rmtree(directory)

exit(0)
