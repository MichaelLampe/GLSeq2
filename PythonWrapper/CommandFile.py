__author__ = 'Michael Lampe'

from subprocess import Popen
import re

"""
Mkdir can sometimes cause errors if the user already had a directory made.
There is a case where try to copy both fastq and fq files, so let's ignore that here.
This will alleviate that.
We'll compile this outside so we don't have to compile a lot
"""
no_useful_error_messages = re.compile("mkdir (.*?)|cp (.*?)", re.IGNORECASE)

class CommandFile:
    file_count = 0

    """
    Gets all the important path variables from the server so that our shell scripts
    Have access to a normal user environment (pretty much)
    Putting this up here means we only hit the system with 4 echo calls instead of a bunch
    """
    def __init__(self, run_name, condor_path, node):
        # Keep track of how many command files we have
        CommandFile.file_count = CommandFile.file_count + 1
        # Just make the files based on the number of their command and the run name because it is easy.
        self.condor_path = condor_path
        self.run_name = run_name
        self.file_name = run_name + "_Node" + str(node.number) + ".sh"

        """
        Should be a list of commands that the class command has prepared.
        They will be strings that are just a single command
        """
        self.command_bundle = node.command.split(" && ")

    """
    Creates a template bash file with the given command node.
    """
    def generate_bash_file(self):
        file = self.condor_path + self.run_name + "/" + self.file_name

        # Need to pass a string if using shell=True
        with open(file,"w") as command:
            """
            Write the bash header
            """
            command.write("#!/bin/bash\n\n")

            for x in xrange(0,len(self.command_bundle)):
                self.command_bundle[x] = self.command_bundle[x].replace("\"","")
                """
                Write the command in with a new line
                """
                command.write(self.command_bundle[x].lstrip())
                command.write("\n")

                """
                Not all of the commands have useful error messages
                """
                result = re.search(no_useful_error_messages,self.command_bundle[x])
                if (result==None):
                    command.write(self.error_checking(x,self.command_bundle[x]))

            """
            Normal exit is 0, which means and worked.
            """
            command.write("\nexit 0")

        """
        Makes the newly created file an executable
        """
        Popen(self.create_executable(file))
        return self.file_name

    """
    Adjust the files such that it makes the file executable.
    """
    def create_executable(self,sh_file):
        command_list = list()
        command_list.append("chmod")

        # Allows for it to be executable on the file system
        command_list.append("+x")
        command_list.append(sh_file)
        return command_list

    """
    Adds an error checking script to each command in the shell script
    Based on checking exit codes and if they are 0 or not.
    """
    def error_checking(self,command_number,command):
        error_check = "\n"
        error_check += "if [[ $? -ne 0 ]]; then\n"
        error_check += "  echo \"Step " + str(command_number) + " failed. The command was " + command + " Exiting now.\" \n"
        error_check += "  exit " + str(command_number + 1) + "\n"
        error_check += "fi\n"
        error_check += "\n"
        return error_check
