__author__ = 'Michael Lampe'

import subprocess
import os
# The goal of this class is to create a wrapper from within which we can run
# GLSeq and receive the output command


class GlSeqRun:
    def __init__(self, glseq_path, update_database, prepare_data, align,
                 count, collect, run_name, protocol_id, attribute_file_path, condor):

        self.glseq_path = glseq_path
        self.update = update_database
        self.prepare = prepare_data
        self.align = align
        self.count = count
        self.collect = collect
        self.run_name = run_name
        self.protocol_id = protocol_id
        self.attribute_file_path = attribute_file_path
        self.condor = condor

        self.condor_path = attribute_file_path.split("/")
        self.condor_path = self.condor_path[0:len(self.condor_path) - 1]
        self.condor_path = "/".join(self.condor_path) + "/"

    def run(self):
        # Scarcity server uses Python 2.6, so we need to do it this way.
        output = subprocess.Popen(self.define_run_args(), stdout=subprocess.PIPE)
        out, err = output.communicate()

        # If condor or not
        if self.condor != "":
            output = subprocess.Popen(self.define_run_args(), stdout=subprocess.PIPE)
            out, err = output.communicate()
            directory = self.condor_path + self.run_name
        else:
            # Run non condor as daemon.
            args = self.define_run_args()
            args.append("&")
            subprocess.Popen(args)
            return

        print directory
        # Create Directory
        if not os.path.exists(directory):
            os.makedirs(directory)

        output_file = directory + "/" + self.run_name + "_GLSeqOutput.txt"
        with open(output_file, 'w') as log_file:
            if out is not None:
                log_file.write("The output was:\n")
                log_file.write(out)
            else:
                log_file.write("There was no output.\n")
            if err is not None:
                log_file.write("The errors were:\n")
                log_file.write(err)
            else:
                log_file.write("There were no error messages.\n")
        try:
            output.kill()
        except:
            pass
        # The Rscript will give us a command output here.
        # This is just a big long string that we need to parse.
        return self.parse_command(out)

    def define_run_args(self):
        run = list()
        run.append("Rscript")
        run.append(self.glseq_path)
        run.append(self.update)
        run.append(self.prepare)
        run.append(self.align)
        run.append(self.count)
        run.append(self.collect)
        run.append(self.run_name)
        run.append(self.protocol_id)
        run.append(self.attribute_file_path)
        run.append(self.condor)
        return run

    def parse_command(self,command):
        initial_pattern = "[1] \""
        broken_command = command.split(initial_pattern)
        # First element will always be an empty string
        broken_command.pop(0)
        for x in range(len(broken_command)-1, -1, -1, ):
            # Remove quotes from end of the command
            broken_command[x] = broken_command[x].replace("\"", "")
            broken_command[x] = broken_command[x].replace("\n", "")
            if len(broken_command[x]) <= 0:
                broken_command.pop(x)
        # This broken apart array is what we will give to Command to construct things in CentralInterface
        return broken_command
