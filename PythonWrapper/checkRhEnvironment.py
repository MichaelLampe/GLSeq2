__author__ = 'Michael Lampe'

import shutil
import os


def rh_environment_check(check_directory):
    """
    Rockhopper aligner requires files where both PTT and FNA files are expected.
    This script is utilized to check individual directories and ask if they contain both of these file types.
    If both are not present, they are removed from the file system.

    :param check_directory:
    :return:
    """
    FNA_EXT = ".fna"
    PTT_EXT = ".ptt"

    if not check_directory.endswith("/"):  # Check trailing directory
        check_directory += "/"

    """
    Looks into the directory that is supposed to be checked
    """
    for directory in os.listdir(check_directory):
        directory = check_directory + directory
        has_fna, has_ptt = False, False

        """
        Goes through the files in a given directory looking for certain file types.
        """
        for files in os.listdir(directory):
            files = directory + files

            if files.endswith(FNA_EXT):
                has_fna = True
            elif files.endswith(PTT_EXT):
                has_ptt = True

        """
        Remove the directory file if it doesn't contain both an FNA and a PTT file
        """
        if not has_fna or not has_ptt:
            shutil.rmtree(directory)


if __name__ == "__main__":
    import sys

    """
    Get the single command line arg
    """
    directory_to_check = sys.argv[1]
    try:
        rh_environment_check(directory_to_check)
        exit(0)
    except:
        exit(1)
