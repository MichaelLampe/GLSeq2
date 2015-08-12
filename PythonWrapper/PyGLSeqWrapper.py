__author__ = 'mlampe'

import sys
import os
from Wrapper import GlSeqRun
from CommandsProcessor import CommandsProcessor
from CommandStack import Stack
from Command import Command
from CommandFile import CommandFile
# Command line args incoming~
def print_message():
    print("Welcome to the GLSeq Wrapper for Condor!")
    print("The wrapper takes in the same arguments as the traditional Rscript version of GLSeq")
    print("\n")
    print("python PyGLSeqWrapper.py [absolute path to GLSeq.top.R] [update/noupdate] [dataprep/nodataprep]")
    print("[alignmet/noalignment] [counting/nocounting] [collect/nocollect] [Name of Run] [Protocol ID (0 is default)]")
    print("[Attribute file path or directory containing multiple attribute files]")
    print("If running from an attribute file directory, the field [Name of Run] can be excluded (But doesn't need to be)")
    print("\n")
    print("These commands will instantiate a GLSeq run, develop a Dagman protocol for your run, generate unique shell scripts and logging, and finally submit your job to HTCondor.")

def trail_check(directory):
    if directory.endswith("/"):
        return directory
    else:
        directory = directory + "/"
        return directory

def clear_classes():
    Command.parallel_track = list()
    Command.job_count = 1
    CommandsProcessor.file_count = 0
    CommandFile.file_count = 0

if len(sys.argv) < 1:
    print("Insufficient arguments")
    print_message()
    exit(1)

# Makes all the scope stuff fine
run_name = sys.argv[7]
if (len(sys.argv) < 10):
    # We can exclude run name in directory runs.
    glseq_path = sys.argv[1]
    update_database = sys.argv[2]
    prepare = sys.argv[3]
    align = sys.argv[4]
    count = sys.argv[5]
    collect = sys.argv[6]
    protocol_id = sys.argv[7]
    attribute_file_path = sys.argv[8]
else:
    glseq_path = sys.argv[1]
    update_database = sys.argv[2]
    prepare = sys.argv[3]
    align = sys.argv[4]
    count = sys.argv[5]
    collect = sys.argv[6]
    run_name = sys.argv[7]
    protocol_id = sys.argv[8]
    attribute_file_path = sys.argv[9]


if os.path.isdir(attribute_file_path):
    attribute_file_path = trail_check(attribute_file_path)
    files = os.listdir(attribute_file_path)
    for file_number in range(0,len(files)):
        if files[file_number].endswith(".R"):
            # We'd rather have one fail and the rest run then one fail and none run.
            current_att_file = attribute_file_path + files[file_number]
            try:
                # We'll make the run name the attribute file name for easy of use.
                run_name = os.path.splitext(files[file_number])[0]
                # Replaces the run name with attribute file name for easy of use.
                # Creates a GLSeq wrapper run which will print out all the commands to be taken in by the python wrapper
                if not os.path.exists(os.path.dirname(run_name + "/")):
                    os.makedirs(os.path.dirname(run_name + "/"))
                #
                current_run = GlSeqRun(glseq_path,update_database, prepare, align, count, collect, run_name, protocol_id, current_att_file,"TRUE")
                # # Initiates a run through the wrapper, and takes in all the output which is the various commands that would be run
                commands = current_run.run()
                # Creates the directory for the Dagman shell scripts and log files
                # Trims and groups commands
                processor = CommandsProcessor(commands)
                # This will break apart the commands in groups that can be parallelized while still being ordered.
                ordered_commands = processor.handle()
                # Starts a command stack which will do the processing moving forward
                command_stack = Stack(ordered_commands)
                # Creates a new stack by a run name
                command_stack.create_stack(run_name)
                # Creates the dag workflow requested
                command_stack.create_dag_workflow(run_name)
                # Submits the workflow to condor.
                #command_stack.submit()
                command_stack.plot_graph(run_name)
                clear_classes()
            except:
                import traceback
                print(traceback.print_exc())
                print("Your attribute file " + current_att_file + " failed to run.")
                pass
else:
    # Creates the directory for the Dagman shell scripts and log files
    if not os.path.exists(os.path.dirname(run_name + "/")):
        os.makedirs(os.path.dirname(run_name + "/"))
    # Creates a GLSeq wrapper run which will print out all the commands to be taken in by the python wrapper
    current_run = GlSeqRun(glseq_path,update_database, prepare, align, count, collect, run_name, protocol_id, attribute_file_path,"TRUE")
    # # Initiates a run through the wrapper, and takes in all the output which is the various commands that would be run
    commands = current_run.run()
    # Trims and groups commands
    processor = CommandsProcessor(commands)
    # This will break apart the commands in groups that can be parallelized while still being ordered.
    ordered_commands = processor.handle()
    # Starts a command stack which will do the processing moving forward
    command_stack = Stack(ordered_commands)
    # Creates a new stack by a run name
    command_stack.create_stack(run_name)
    # Creates the dag workflow requested
    command_stack.create_dag_workflow(run_name)
    # Submits the workflow to condor.
    command_stack.submit()
# It seems to have worked
exit(0)