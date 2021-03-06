__original_author__ = 'Michael Lampe'
__author__ = 'Michael Lampe'
__version__ = '1.0.1'
__last_update__ = "Febuary 2016"

import os
import re
import sys

"""
Precompile regexes
"""
GENE_ID_REGEX = re.compile("(?<=gene_id \").*(?=\"; transcript_id)", re.IGNORECASE)

"""
File extension
"""
PTT_FILE_EXTENSION = ".ptt"

# User interaction stuff
def _print_help():
    print "HELP FOR GTF/GFF TO PTT FILE CONVERTER:\n"
    print "Originally written by " + __original_author__
    print "Last Updated by " + __author__ + ", " + __last_update__ + "\n"
    print "Version " + __version__
    print "\nObjective:"
    print "This program converts GTF and GFF files into PTT files with the limited amount of convertable information " \
          "available between the files for programs that require PTT files."
    print "\nUse:"
    print "GtfToPtt.py [Input GTF or GFF file] [Output Directory] [Advanced Options]"
    print "The output directory is the place where individual sequence names from the input GTF or GFF files will " \
          "be placed into individual folders that are either preexisting or created at run time."
    print "\nAdvanced Options:"
    print "Advanced options should be passed in the format: OPTION_COMMAND=OPTION_VALUE"
    print "\nCurrent advanced options:"
    print "OPTION_COMMAND = NUMBER_OF_INPUT_FILE_COLUMNS"
    print "DEFAULT_VALUE = 9\n"


def convert_gtf_to_ptt(gtf_file, reference_genome_destination, gtf_column_count=9):
    """
    Converts a gtf file into a ptt file based on the expected formatting

    :param gtf_file: The input GTF file to be converted
    :param reference_genome_destination: The destination of our converted file.
    :param gtf_column_count: The number of columns in the GTF file that is to be converted.
    :return:
    """

    """
    Append trailing slash
    """
    if not reference_genome_destination.endswith("/"):
        directory_name = reference_genome_destination + "/"

    """
    Format species name based on GTF/GFF file name.
    Remove gtf file ending and absolute file path stuff
    """
    species = re.search("[0-9a-zA-Z_.]*$", gtf_file, re.IGNORECASE).group(0)
    species = species.split(".")[0]

    """
    PTT File Format
    Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product

    Dictionary containing the name of each sequence found in the single gtf/gff file and the number of proteins that
    Have been added to these files.
    """
    files_data = dict()
    with open (gtf_file) as file:
        for line in file:
            div_in = line.split("\t")
            if not (len(div_in) < gtf_column_count):

                """
                The sequence name is the first item in a gtf/gff file
                If we don't have a matching key, we'll want to do some things to establish the given
                Sequence within our structure
                """
                if div_in[0] not in files_data:

                    """
                    Create folder
                    """
                    folder_name = directory_name + div_in[0] + "/"
                    new_ptt_file = folder_name + div_in[0]

                    """
                    Makes directory if not available
                    """
                    if not os.path.exists(os.path.dirname(folder_name)):
                        os.makedirs(os.path.dirname(folder_name))

                    """
                    Create file with header
                    Species, which chromosome
                    """
                    h_line_1 = species + ", " + div_in[0]

                    """
                    Number of proteins.  We don't know this when writing this so we just give it a dummy value
                    """
                    h_line_2 = str("<?> proteins")

                    """
                    Description of column meanings
                    """
                    h_line_3 = "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct"

                    """
                    Add to dictionary so this doesn't occur again
                    At index 0 we have the number of proteins (Which the first one will be 1)
                    At index 1 we have a list of all the lines.  This saves us from rewriting the entire file when we want
                    to finally write to a file.
                    """
                    files_data[div_in[0]] = [1, list()]

                    """
                    If we add these now we won't need to do an insert when we want to input the number of proteins
                    """
                    files_data[div_in[0]][1].append(h_line_1)
                    files_data[div_in[0]][1].append(h_line_2)
                    files_data[div_in[0]][1].append(h_line_3)

                ptt_file = directory_name + div_in[0] + "/" + div_in[0]         # File Name
                location = div_in[3] + ".." + div_in[4]                         # Location
                strand = div_in[6]                                              # Strand
                length = str(int(div_in[4]) - int(div_in[3]) - 2)               # Length of protein
                pid = "-"                                                       # PID Gene
                gene = re.search(GENE_ID_REGEX, line).group(0)                  # Gene
                syn = str(species) + "_" + str(files_data[div_in[0]][0])        # Synonym
                code = "-"                                                      # Code
                cog = "-"                                                       # COG
                product = "-"                                                   # Product
                files_data[div_in[0]][0] = files_data[div_in[0]][0] + 1         # Increase protein number counts

                ptt_file_line = location + "\t" + \
                                strand + "\t" + \
                                length + "\t" + \
                                pid + "\t" + \
                                gene + "\t" + \
                                syn + "\t" + \
                                code + "\t" + \
                                cog + "\t" + \
                                product

                files_data[div_in[0]][1].append(ptt_file_line)

    for key in files_data:
        """
        Fix the number of the proteins in the header of each file
        Then create the files.
        """
        files_data[key][1][1] = str(files_data[key][0] - 1) + " proteins"
        file_name = directory_name + key + "/" + key + PTT_FILE_EXTENSION

        with open(file_name, "w") as new_ptt_file:
            for line in files_data[key][1]:
                new_ptt_file.write(line + "\n")


if __name__ == "__main__":

    # Check args length
    if len(sys.argv) < 3:
        _print_help()

    gtf_file = ""
    directory_name = ""
    gtf_column_count = 9

    for i , command in enumerate(sys.argv):
        command = str(command)                  # Make sure it is a string
        if i == 0:
            pass                                # This is the script name, just ignore it
        elif i == 1:
            gtf_file = command                 # Get input GTF file
        elif i == 2:
            directory_name = command            # Directory to write files and folders to
        elif i > 2:                             # Advanced commands
            try:
                split_command = command.split("=")
                option_name, opt_value = split_command[0], split_command[1]

                command_executed = False
                if option_name is "NUMBER_OF_INPUT_FILE_COLUMNS":
                    gtf_column_count = int(opt_value)
                    command_executed = True

                """
                If user supplied a dud advanced command
                """
                if not command_executed:
                    print "ERROR: Advanced option command: " + option_name + " could not be matched with any " \
                        "advanced options.  Please check input or type GtfToPtt.py --help for additional help."
                    exit(4)
            # If values improperly linked
            except IndexError:
                print "ERROR: Advanced command was incorrectly used. Please input advanced commands in the format:" \
                        "\n COMMAND_NAME=COMMAND_VALUE\n and only pass them after the file name and directory name."
                exit(3)

    # Check if proper file type
    if not gtf_file.endswith("gtf") and not gtf_file.endswith("gff"):
        print "ERROR: Input file to be converted should have either gtf or gff file ending."
        exit(2)

    convert_gtf_to_ptt(gtf_file, directory_name, gtf_column_count)
    exit(0)