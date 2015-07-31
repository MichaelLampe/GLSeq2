__author__ = 'mlampe'

import sys
from Wrapper import GlSeqRun
from CommandsProcessor import CommandsProcessor

# Command line args incoming~
update_database = sys.argv[1]
prepare = sys.argv[2]
align = sys.argv[3]
count = sys.argv[4]
collect = sys.argv[5]
run_name = sys.argv[6]
protocol_id = sys.argv[7]
attribute_file_path = sys.argv[8]

current_run = GlSeqRun(update_database, prepare, align, count, collect, run_name, protocol_id, attribute_file_path,"TRUE")

# Initiates a run through the wrapper, and takes in all the output which is the various commands that would be run
commands = current_run.run()
# We get back a list of commands that are only cleaned up (Aka they could run if we just threw them at the shell)
# Let's parallize and break apart the commands in the CommandsProcessor
processor = CommandsProcessor(commands)
processor.handle()

