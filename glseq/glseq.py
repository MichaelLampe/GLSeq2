__author__ = 'mlampe'

import copy
import os
import argparse
from condor_glseq.IndividualGlSeqRun import GlSeqRun
from condor_glseq.CommandAssembler import CommandProcessor
from condor_glseq.CommandStack import Stack
from condor_glseq.CommandFile import CommandFile

class GlSeq():
    def trail_check(self, directory):
        if not directory.endswith("/"):
            directory += "/"
        return directory

    def clear_classes(self):
        CommandFile.file_count = 0

    def start_run(self, glseq_instance):
        commands = glseq_instance.run()

        if glseq_instance.condor == "":
            # Don't do all the condor stuff it not running as Condor.
            return

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
        self.clear_classes()

    def start_run_instance(self, command_line_args, error_message="0"):
        """
        Takes in command line args and makes run instances for either an individual attribute file or directory them.

        :param command_line_args: Dictionary of args.
        :return: 0 = Worked ; Anything else = error message
        """
        """
        Condor mode
        """
        if command_line_args.get("condor"):
            command_line_args["condor"] = "TRUE"
        else:
            command_line_args["condor"] = ""

        """
        Directory mode
        """
        if os.path.isdir(command_line_args["attribute_file_path"]):
            attribute_file_path = self.trail_check(command_line_args["attribute_file_path"])
            command_line_args["attribute_file_path"] = os.listdir(attribute_file_path)

            if len(command_line_args["attribute_file_path"]) <= 0:
                return "No files available in attribute file directory supplied."

            if len([file for file in command_line_args["attribute_file_path"] if file.endswith(".R")]) <= 0:
                return "No attributes files (Files ending with .R) available in attribute file directory supplied."

        else:
            command_line_args["attribute_file_path"] = [command_line_args["attribute_file_path"]]
            if not os.path.isfile(command_line_args["attribute_file_path"][0]):
                return "Unable to locate attribute file supplied.  Please check whether supplied file exists.  " \
                       "\n\nIf you meant to run a directory instead, we have determined that the directory you wanted does not exist as well." \
                       "\n\nRemember to supply the absolute path."

        """
        Check if glseq path supplied
        """
        if command_line_args["glseq_path"] is None:
            # Local version if pip installed
            path = __file__.split("/")
            path = path[:len(path)-1]
            string_path = ""
            for p in path:
                string_path += p + "/"

            absolute_path = os.path.abspath(string_path + "GLSeqScripts/Main/GLSeq.top.R")
            if os.path.isfile(absolute_path):
                command_line_args["glseq_path"] = absolute_path
            else:
                return "No copy of GLSeq.top.R supplied and no local copy installed. " \
                       "The command \'pip install glseq\' installs a local copy of this file.  " \
                       "You can then run this script as python -m glseq.glseq <options>"

        """
        Go through and start each run.
        """
        for file_number, current_file in enumerate(command_line_args["attribute_file_path"]):
            # Only try running .R files.
            current_attribute_file_path = command_line_args["attribute_file_path"][file_number]
            if current_file.endswith(".R"):
                if command_line_args.get("verbose"):
                    print "Attempting to run attribute file {0}".format(current_file)

                # Copy array so we don't overwrite anything important.
                current_commands = copy.deepcopy(command_line_args)
                current_commands["attribute_file_path"] =  current_attribute_file_path

                try:
                    # Switch run name if more than 1 file running or run name not provided.
                    if (len(command_line_args["attribute_file_path"]) > 1) or (current_commands["run_name"] == ""):
                        run_name = current_file.split("/")[-1].replace(".R","")
                        current_commands["run_name"] = run_name

                    # Make directory if not currently in existence.
                    if not os.path.exists(os.path.dirname(current_commands["run_name"] + "/")):
                        os.makedirs(os.path.dirname(current_commands["run_name"] + "/"))

                    # Instantiate run
                    current_run = GlSeqRun(current_commands["glseq_path"],
                                           current_commands["update_database"],
                                           current_commands["prepare"],
                                           current_commands["align"],
                                           current_commands["count"],
                                           current_commands["collect"],
                                           current_commands["run_name"],
                                           current_commands["protocol_id"],
                                           current_commands["attribute_file_path"],
                                           current_commands["condor"]
                                           )

                    # Start run
                    self.start_run(current_run)

                    if command_line_args.get("verbose"):
                        print "Successfully started run with attribute file {0}".format(current_commands["attribute_file_path"])

                except Exception as e:
                    # Overwrites good return code
                    if error_message == "0":
                        error_message = ""

                    # Adds to the error message for each error we see.
                    current_message = ""
                    current_message += "The run name is : " + command_line_args["run_name"]
                    current_message += "\n\n"
                    current_message += str(e)
                    current_message += "\n\n"
                    current_message += "Your attribute file " + \
                                     str(command_line_args["attribute_file_path"][file_number]) + \
                                     " failed to run.\n\n"

                    if command_line_args.get("verbose"):
                        print current_message

                    error_message += current_message

        return error_message

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs a GlSeq2 bioinformatic pipeline.\n  "
                                                 "To utilize this pipeline with a graphical interface, "
                                                 "try the command {0}\n".format("python -m glseq.glseq-gui"))

    parser.add_argument("--condor",
                       action="store_true",
                       help="Attempts to run the current job on a local Condor server.\n")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Outputs errors as they happen instead of at the end,"
                             " slightly more vocal about what is occuring.")

    parser.add_argument('update_database',
                        type=str,
                        choices=["noupdate", "update"],
                        help="Indicates whether, after the run, a database should be updated.  Currently unimplemented.\n")

    parser.add_argument('prepare',
                        type=str,
                        choices=["nodataprep", "dataprep"],
                        help="Indicates whether input FASTQ files should be prepared prior to being aligned.\n")

    parser.add_argument('align',
                        type=str,
                        choices=["noalignment","alignment"],
                        help="Indicates if alignment should happen on the input files.\n")

    parser.add_argument('count',
                        type=str,
                        choices=["nocounting","counting"],
                        help="Indicates if output SAM files, or previously prepared SAM files, should be counted.\n")

    parser.add_argument('collect',
                        type=str,
                        choices=["nocollect","collect"],
                        help="Indicates if, after the run, a summary step should be run to collect additional information.\n")

    parser.add_argument('attribute_file_path',
                        type=str,
                        help="The absolute path to the attribute file which will be used to instantiate the pipeline.  "
                             "If the --condor option is active, you can also provide a directory containing multiple "
                             "attribute files to instantiate runs for all files within the directory.\n")

    # Default for this is defined above to get absolute path.
    parser.add_argument('glseq_path',
                        type=str,
                        nargs="?",
                        help="Absolute path to the GlSeq scripts.  Uses the setup.py installed version if nothing provided.\n")

    parser.add_argument('run_name',
                        type=str,
                        nargs='?',
                        const="",
                        default="",
                        help="Name of run if running a single attribute file. "
                             " If running a directory of attribute files, this option will be ignored and runs "
                             "will be named by their attribute file names.\n")

    parser.add_argument('protocol_id',
                        type=int,
                        nargs='?',
                        const=0,
                        default=0,
                        help="Unused currently, but planned to be used to retrieve attributes "
                             "from a protocol table in a database.  If you are not using a database,"
                             " the default value of 0 is used.\n")
    r = GlSeq()
    return_message = r.start_run_instance(vars(parser.parse_args()))
    exit("\nReturn Message: " + return_message + "\n")