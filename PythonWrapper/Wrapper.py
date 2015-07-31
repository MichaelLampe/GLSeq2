__author__ = 'mlampe'

import subprocess
# The goal of this class is to create a wrapper from within which we can run
# GLSeq and receive the output command

class GlSeqRun:
    def __init__(self,update_database,prepare_data,align,count,collect,run_name,protocol_id,attribute_file_path,condor):
        self.update = update_database
        self.prepare = prepare_data
        self.align = align
        self.count = count
        self.collect = collect
        self.run_name = run_name
        self.protocol_id = protocol_id
        self.attribute_file_path = attribute_file_path
        self.condor = condor;

    def run(self):
        # Scarcity server uses Python 2.6, so we need to do it this way.
        output = subprocess.Popen(self.defineRunArgs(),stdout=subprocess.PIPE)
        out, err = output.communicate()
        # The Rscript will give us a command output here.
        # This is just a big long string that we need to parse.
        return self.parseCommand(out)

    def defineRunArgs(self):
        run = list()
        run.append("Rscript")
        run.append("GLSeq.top.R")
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

    def parseCommand(self,command):
        initial_pattern = "[1] \""
        broken_command = command.split(initial_pattern)
        # First element will always be an empty string
        broken_command.pop(0)
        for x in range (0,len(broken_command)):
            # Remove quotes from end of the command
            broken_command[x].replace("\"","")
        # This broken apart array is what we will give to Command to construct things in CentralInterface
        return broken_command




    #def createBashScript(self):