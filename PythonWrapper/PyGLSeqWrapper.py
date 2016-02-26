__author__ = 'Michael Lampe'

import sys, os
from IndividualGlSeqRun import GlSeqRun
from CommandAssembler import CommandProcessor
from CommandStack import Stack
from CommandFile import CommandFile

def print_message():
    print("Welcome to the GLSeq Wrapper for Condor!")
    print("The wrapper takes in the same arguments as the traditional Rscript version of GLSeq")
    print("\n")
    print("python PyGLSeqWrapper.py [absolute path to GLSeq.top.R] [update/noupdate] [dataprep/nodataprep]")
    print("[alignmet/noalignment] [counting/nocounting] [collect/nocollect] [Name of Run] [Protocol ID (0 is default)]")
    print("[Attribute file path or directory containing multiple attribute files]")
    print("If running from an attribute file directory, "
          "the field [Name of Run] can be excluded (But doesn't need to be)")
    print("\n")
    print("These commands will instantiate a GLSeq run, develop a "
          "Dagman protocol for your run, generate unique shell scripts"
          " and logging, and finally submit your job to HTCondor.")

def trail_check(directory):
    if not directory.endswith("/"):
        directory += "/"
    return directory

def clear_classes():
    CommandFile.file_count = 0

def start_run(glseq_instance):
    commands = glseq_instance.run()

    # Creates the directory for the Dagman shell scripts and log files
    # Trims and groups commands
    script_path = glseq_instance.glseq_path.split("/")
    script_path = script_path[0:len(script_path) - 1]
    script_path = "/".join(script_path) + "/"

    processor = CommandProcessor(commands, glseq_instance.run_name, script_path, glseq_instance.condor_path)
    # Puts the commands into an organized graph structure
    graph = processor.create_graph()
    # Starts a command stack which will do the processing moving forward
    command_stack = Stack(graph, glseq_instance.condor_path)

    # If MatPlotLib is available, create a nice graph.
    try:
        command_stack.plot_graph(glseq_instance.run_name)
    except:
        pass

    # Creates a new stack by a run name
    command_stack.create_stack(glseq_instance.run_name)
    # Submits the workflow to condor.
    command_stack.submit()
    clear_classes()

# Just catch if there is not enough arguments right away.
if len(sys.argv) < 1:
    print("Insufficient arguments")
    print_message()
    exit(1)

# What a command passes in
command = {
    "glseq_path" : "",
    "update_database" : "",
    "prepare" : "",
    "align" : "",
    "count": "",
    "collect": "",
    "attribute_file_path" : "",
    "run_name" : "",
    "protocol_id" : ""
}

# ['/home/GLBRCORG/mrlampe/_GlowCurrent/Scripts/PyGLSeqWrapper.py',
# 'Placeholder',
# 'nodataprep',
#  'alignment',
# 'counting',
# 'collect',
#  'RunningOneWindows',
#  '0',
#  '/home/GLBRCORG/mrlampe/ssh/RunningOneWindows.R']

# Based on the command line args, fill the command fields.
try:
    command["glseq_path"] = sys.argv[1]
    command["update_database"] = sys.argv[2]
    command["prepare"] = sys.argv[3]
    command["align"] = sys.argv[4]
    command["count"] = sys.argv[5]
    command["collect"] = sys.argv[6]

    if len(sys.argv) < 10:
        # We can exclude run name in directory runs.
        command["attribute_file_path"] = sys.argv[8]
        command["protocol_id"] = sys.argv[7]

    else:
        command["attribute_file_path"] = sys.argv[9]
        command["run_name"] = sys.argv[7]
        command["protocol_id"] = sys.argv[8]
except IndexError:
    print_message()
    print("There was an error in your commands.  Please check your command string")
    exit(1)

# If the file path provided is a directory, we'll use all the R files in it.
if os.path.isdir(command["attribute_file_path"]):
    attribute_file_path = trail_check(command["attribute_file_path"])
    files = os.listdir(attribute_file_path)

    for file_number, current_file in enumerate(files):
        if files[file_number].endswith(".R"):
            # We'd rather have one fail and the rest run then one fail and none run.
            current_att_file = attribute_file_path + current_file
            try:
                # We'll make the run name the attribute file name for easy of use.
                run_name = os.path.splitext(current_file)[0]
                # Replaces the run name with attribute file name for easy of use.
                # Creates a GLSeq wrapper run which will print out all the commands to be taken in by the python wrapper
                if not os.path.exists(os.path.dirname(run_name + "/")):
                    os.makedirs(os.path.dirname(run_name + "/"))

                current_run = GlSeqRun(command["glseq_path"],
                                       command["update_database"],
                                       command["prepare"],
                                       command["align"],
                                       command["count"],
                                       command["collect"],
                                       run_name,
                                       command["protocol_id"],
                                       command["attribute_file_path"],
                                       "TRUE")
                start_run(current_run)
            except Exception as e:
                print e
                print("Your attribute file " + current_att_file + " failed to run.")
                pass
else:
    # Creates a GLSeq wrapper run which will print out all the commands to be taken in by the python wrapper
    try:
        current_run = GlSeqRun(command["glseq_path"],
                               command["update_database"],
                               command["prepare"],
                               command["align"],
                               command["count"],
                               command["collect"],
                               command["run_name"],
                               command["protocol_id"],
                               command["attribute_file_path"],
                               "TRUE")
        start_run(current_run)
    except Exception as e:
        print "The run name is : " + command["run_name"]
        print e
        print("Your attribute file " + command["attribute_file_path"] + " failed to run.")
        exit(1)
exit(0)
