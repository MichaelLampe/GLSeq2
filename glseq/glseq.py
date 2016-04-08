__author__ = 'mlampe'

import sys
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
        if command_line_args.condor:
            command_line_args["condor"] = "TRUE"
        else:
            command_line_args["condor"] = ""

        """
        Directory mode
        """
        if os.path.isdir(command_line_args["attribute_file_path"]):
            attribute_file_path = self.trail_check(command_line_args["attribute_file_path"])
            command_line_args["attribute_file_path"] = os.listdir(attribute_file_path)
        else:
            command_line_args["attribute_file_path"] = [self.trail_check(command_line_args["attribute_file_path"])]

        """
        Go through and start each run.
        """
        for file_number, current_file in enumerate(command_line_args["attribute_file_path"]):
            # Only try running .R files.
            current_attribute_file_path = command_line_args["attribute_file_path"][file_number]
            if current_file.endswith(".R"):
                if command_line_args.verbose:
                    print "Attempting to run attribute file {0}".format(current_file)

                current_att_file = current_attribute_file_path + current_file

                # Copy array so we don't overwrite anything important.
                current_commands = copy.deepcopy(command_line_args)
                current_commands["attribute_file_path"] = current_att_file

                try:
                    # Switch run name if more than 1 file running or run name not provided.
                    if (len(command_line_args["attribute_file_path"]) > 1) or (current_commands["run_name"] == ""):
                        run_name = os.path.splitext(current_file)[0]
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

                    if command_line_args.verbose:
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
                                     command_line_args["attribute_file_path"] + \
                                     " failed to run.\n\n"

                    if command_line_args.verbose:
                        print current_message

                    error_message += current_message

        return error_message

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs a GlSeq2 bioinformatic pipeline.\n")

    parser.add_argument("--condor",
                       action="store_true",
                       help="Attempts to run the current job on a local Condor server.\n")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Outputs errors as they happen instead of at the end,"
                             " slightly more vocal about what is occuring.")

    parser.add_argument('glseq_path',
                        type=str,
                        help="Absolute path to the GlSeq scripts.\n")

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

    print parser.parse_args()
    r = GlSeq()
    return_message = r.start_run_instance(parser.parse_args())
    exit(return_message)