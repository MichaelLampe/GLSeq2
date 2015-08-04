__author__ = 'mlampe'

from subprocess import Popen
class CommandFile:
    file_count = 0
    def __init__(self,run_name,command_bundle,parent):
        # Keep track of how many command files we have
        CommandFile.file_count = CommandFile.file_count + 1
        # Just make the files based on the number of their command and the run name because it is easy.
        self.run_name = run_name
        self.file_name = run_name + "_" + str(CommandFile.file_count) + ".sh"
        # Should be a list of commands that the class command has prepared.
        # They will be strings that are just a single command
        self.command_bundle = command_bundle.split(" && ")
        #
        self.parent = parent

    def generate_bash_file(self):
        file = self.run_name + "/" + self.file_name
        with open(file,"w") as command:
            # Bash header
            command.write("#!/bin/bash\n\n")
            # Index through the commands writing them as we go w/ error checking
            for x in range(0,len(self.command_bundle)):
                # lstrip removes the leading white space which could be annoying
                self.command_bundle[x] = self.command_bundle[x].replace("\"","")
                command.write(self.command_bundle[x].lstrip())
                command.write(self.error_checking(x,self.command_bundle[x]))
        # Makes the newly created file an executable
        Popen(self.create_executable(file))
        return self.file_name

    def create_executable(self,sh_file):
        command_list = list()
        command_list.append("chmod")
        # Allows for it to be executable on the file system
        command_list.append("+x")
        command_list.append(sh_file)
        return command_list

    # Adds an error checking script to each command in the shell script
    def error_checking(self,command_number,command):
        error_check = "\n"
        error_check = error_check + "if [[ $? -ne 0 ]]; then\n"
        error_check = error_check + "  echo \"Step " + str(command_number) + " failed. The command was " + command + " Exiting now.\" \n"
        error_check = error_check + "  exit 1\n"
        error_check = error_check + "fi\n"
        error_check = error_check + "\n"
        return error_check
