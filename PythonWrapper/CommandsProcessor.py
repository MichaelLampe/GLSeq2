__author__ = 'mlampe'

import re

class CommandsProcessor:
    def __init__(self,commands):
        # The list of all commands broken into the order they were sent
        self.commands = commands

    # This takes all the commands and divides it into the order they were run and tries
    # to paralleize commands run together and group smaller commands
    def handle(self):
        # All the commands in the order they should appear
        command_list = list()
        # Splits the commands by parallel processes (&)
        for command in self.commands:
            command = self.split_background(command)
            command = filter(bool,command)
            command_list.append(command)
        # Reassign and clear
        self.commands = command_list
        command_list = list()
        for command in self.commands:
            command = self.split_linked(command)
            command = filter(bool,command)
            command_list.append(command)
        return command_list

    # This breaks parts that have only one "&" into an array that will then be able to be run in parallel
    def split_background(self,command):
        re1='(?<!&)&(?!&)'	# Matches only a single & that cannot be preceded or succeeded by another &
        rg = re.compile(re1,re.IGNORECASE)
        # Would be nice to add counting split processing here.
        command_split = rg.split(command)
        return command_split

    # This breaks linked commands (&& or ;) into parts that are in sequence (In this same array)
    def split_linked(self,command):
        """
        I know this somewhat violates duck typing but it is an easy (And really readable) way to check if it is a
        single command.  Without this, iteration would cause individual string characters to be printed and the lack
        of iteration would cause a regex error because you can't regex a list type.
        """
        if type(command) is list:
            split_list = list()
            for c in command:
                re1 = '&&|;'
                rg = re.compile(re1,re.IGNORECASE)
                split_list.append(rg.split(c))
            return split_list
        else:
            re1 = '&&|;'
            rg = re.compile(re1,re.IGNORECASE)
            command_split = rg.split(command)
            return command_split