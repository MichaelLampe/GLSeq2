__author__ = 'mlampe'

class CommandFile:
    file_count = 0
    def __init__(self,run_name,command_bundle):
        # Keep track of how many command files we have
        CommandFile.file_count = CommandFile.file_count + 1
        # Just make the files based on the number of their command and the run name because it is easy.
        self.file_name = run_name + "_" + str(CommandFile.file_count)
        # Should be a list of commands that the class command has prepared.
        # They will be strings that are just a single command
        self.command_bundle = command_bundle

    def generate_bash_file(self):
        with open(self.file_name,"w") as command:
            # Bash header
            command.write("#!/bin/bash\n\n")
            # Index through the commands writing them as we go w/ error checking
            for x in range(0,len(self.command_bundle)):
                command.write(self.command_bundle[x])
                command.write(self.error_checking(x,self.command_bundle[x]))

    def error_checking(self,command_number,command):
        error_check = "\n"
        error_check = error_check + "if [[ $? -ne 0]]; then\n"
        error_check = error_check + "  echo \"Step " + str(command_number) + " failed. The command was \"" + command + "\" Exiting now.\" \n"
        error_check = error_check + "  exit 1\n"
        error_check = error_check + "fi\n"
        error_check = error_check + "\n"
        return error_check
