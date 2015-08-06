__author__ = 'mlampe'

import sys
import os
from Wrapper import GlSeqRun
from CommandsProcessor import CommandsProcessor
from CommandStack import Stack

# Command line args incoming~
glseq_path = sys.argv[1]
update_database = sys.argv[2]
prepare = sys.argv[3]
align = sys.argv[4]
count = sys.argv[5]
collect = sys.argv[6]
run_name = sys.argv[7]
protocol_id = sys.argv[8]
attribute_file_path = sys.argv[9]

current_run = GlSeqRun(glseq_path,update_database, prepare, align, count, collect, run_name, protocol_id, attribute_file_path,"TRUE")

# Initiates a run through the wrapper, and takes in all the output which is the various commands that would be run
commands = current_run.run()
# We get back a list of commands that are only cleaned up (Aka they could run if we just threw them at the shell)
# Let's parallize and break apart the commands in the CommandsProcessor
# Command 1 needs to be run first

# # Test Stuff
# command0 = "first"
# command3 = "mkdir ok"
# command4 = "mkdir ok && mkdir eh"
# command5 = "mkdir ok"
# command6 = "mkdir ok"
# command7 = "mkdir ok"
# command8 = "mkdir ok"
# command1 = "\"mkdir \"I love cookies\" && mkdir 2 && mkdir 3\""
# # Command 2 is run after
# command1point5 = "job1 & job2 & job3 &"
# command2 = "above && mkdir above2 ; mkdir above3"
# commands = list()
# commands.append(command0)
# commands.append(command3)
# commands.append(command4)
# commands.append(command5)
# commands.append(command6)
# commands.append(command7)
# commands.append(command8)
# commands.append(command1)
# commands.append(command1point5)
# commands.append(command2)
# run_name = "temp"
if not os.path.exists(os.path.dirname(run_name + "/")):
    os.makedirs(os.path.dirname(run_name + "/"))

processor = CommandsProcessor(commands)
# This will break apart the commands in groups that can be parallelized while still being ordered.
ordered_commands = processor.handle()
command_stack = Stack(ordered_commands)
# This is temp
#run_name = "temp"
command_stack.create_stack(run_name)
command_stack.create_dag_workflow(run_name)
command_stack.submit()

