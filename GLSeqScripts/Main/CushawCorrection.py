__author__ = 'mlampe'
import argparse
import subprocess
parser = argparse.ArgumentParser(description="Runs a GlSeq2 bioinformatic pipeline.\n  "
                                                 "To utilize this pipeline with a graphical interface, "
                                                 "try the command {0}\n".format("python -m glseq.glseq-gui"))

parser.add_argument('unclean_file',
                    type=str,
                    help="Unclean sam file to clean")

parser.add_argument('clean_file',
                    type=str,
                    help="Output clean file")

args = vars(parser.parse_args())
unclean_file = args["unclean_file"].strip()
clean_file = args["clean_file"].strip()

awk_command = ["!", "/", "\\", "t", "\\", "*", "\\", "t", "/"]
awk_command = "".join(awk_command)
full_command = "awk " + awk_command + " " + unclean_file + " > " + clean_file
subprocess.call(full_command, shell=True)