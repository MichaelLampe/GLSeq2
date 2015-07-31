__author__ = 'mlampe'

import re
command = "command 1 && command2 & parallel1 && parallel2 && parallel3 &"

class CommandsProcessor:
    def __init__(self,commands):
        # The list of all commands broken into the order they were sent
        self.commands = commands
        # These are keywords that can be matched so we can better groups jobs
        self.group_keywords = {
            'mkdir',
            'cd',
            'mv',
            'cp',
            'date',
            'samtools'
            'rm'
        }
        self.parallelize_keywords = {
            'HTSeq.scripts.count'
            'GLSeq.FeatureCounts.R'
            'RSEM'
            'cufflinks'
        }

    # Just kinda sounds cool when you call it
    # This takes all the commands and divides it into the order they were run and tries
    # to paralleize commands run together and group smaller commands
    def handle(self):
        # All the commands in the order they should appear
        ordered_commands = list()
        for command in self.commands:
            # Parallel commands are divided into another list, so that order is retained, but
            # the ability to be run in parallel is easily identified.
            parallel_list = self.parallelize_commands(command)
            ordered_commands.append(parallel_list)
        grouped_ordered_commands = list()
        for command in ordered_commands:
            grouped_ordered_commands.append(self.group_commands(command))


    def parallelize_commands(self,command):
        command = command.split("&&")
        parallel_commands = list()
        single_command = list()
        for comm in command:
            if ("&" in comm):
                # Split on parallelizer
                parts = comm.split("&")
                # Add the last command of the first to that command
                if ("wait" not in parts[0]):
                    single_command.append(parts[0])
                # This command is finished, add it to the whole deal & remove empty strings
                single_command = filter(bool,single_command)
                parallel_commands.append(single_command)
                # Clear it
                single_command = list()
                # Add the first part to the new list
                if ("wait" not in parts[1]):
                    single_command.append(parts[1])
            else:
                if ("wait" not in comm):
                    single_command.append(comm)
        single_command = filter(bool,single_command)
        parallel_commands.append(single_command)
        # Remove empty strings
        parallel_commands = filter(bool,parallel_commands)
        return parallel_commands

    def group_commands(self,command):
        print command

# Command 1 needs to be run first
command1 = "command 1 && command2 & parallel1 && parallel2 && parallel3 &"
# Command 2 is run after
command2 = "comm2 && comm3 & wait"
together = list()
together.append(command1)
together.append(command2)
proc = CommandsProcessor(together)
proc.handle()