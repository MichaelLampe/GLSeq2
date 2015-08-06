__author__ = 'mlampe'

from Command import Command
# Pydagman can be found at https://github.com/brandentimm/pydagman
from pydagman.dagfile import Dagfile
from subprocess import Popen

class Stack:
    # This will take all the commands together and organize them into a command stack composed of command objects.
    def __init__(self,command_list):
        self.command_list = command_list
        self.command_stack = list()

        self.group_keywords = [
            # The space is important for keeping them as just the linux commands
            'date ',
            'mkdir ',
            'cd ',
            'mv ',
            'cp ',
            'date ',
            'samtools ',
            'rm ',
            'echo '
        ]

    def create_stack(self,run_name):
        # Connects steps that have a key word
        # The last step is always either the end of this parallel command
        # Or a larger step that is not a keyword in hopes of creating "larger" chunked steps.
        #
        # The format of the command_list is as follows:
        # [All individual commands received][Commands executed in parallel][The linked commands within each of those]
        # Individual Commands
        for z in range(len(self.command_list) - 1,-1,-1):
            # Parallel command
            self.concatenate_commands(self.command_list,z)
            for y in range (len(self.command_list[z]) -1,-1,-1):
                # Commands linked together
                for x in range (len(self.command_list[z][y]) - 1 ,-1,-1):
                    self.concatenate_linked_commands(self.command_list[z][y],x)
        for parent in range (0, len(self.command_list)):
            for step in self.command_list[parent]:
                command = Command(step,parent)
                self.command_stack.append(command)
        for command in self.command_stack:
            command.create_bash_files(run_name)
            self.create_submit_file(run_name)
            command.create_dag_jobs(run_name)

    def concatenate_commands(self,command,current_index):
        for word in self.group_keywords:
            if (word in command[current_index][0][0]):
                if (current_index != 0):
                    # Adds the previous command to the one ahead of it if it is in the list
                    command[current_index -1][0] += command[current_index][0]
                    # Then removes the command just duplicated from its original location
                    command[current_index].pop(0)
                    # Checks if the location is now empty and removes it if it is
                    if(len(command[current_index]) < 1):
                        command.pop(current_index)

    def concatenate_linked_commands(self,command,current_index):
        for word in self.group_keywords:
            if (word in command[current_index - 1]):
                if (current_index != 0):
                    # This links members of the same step/not commands printed individually
                    command[current_index - 1] = command[current_index-1] + " && " + command[current_index]
                    command.pop(current_index)
                    # Exit after you find one because that's all good
                    break

    def create_dag_workflow(self,run_name):
        mydag = Dagfile()

        # Takes care of parallelizing parent child relationships on the graph
        for x in range(0,len(Command.parallel_track)):
            if (x != 0):
                # Add all the previous jobs that are parents
                for y in range(0,len(Command.parallel_track[x-1])):
                    # Add all the children
                    for z in range(0,len(Command.parallel_track[x])):
                        child_job = Command.parallel_track[x][z].dag_jobs[0]
                        parent_job = Command.parallel_track[x-1][y].dag_jobs[len(Command.parallel_track[x-1][y].dag_jobs)-1]
                        child_job.add_parent(parent_job)
        # Put everything together in the Dagfile workflow
        for command in self.command_stack:
            for job in command.dag_jobs:
                mydag.add_job(job)
        self.dag_file = run_name + "/my_workflow.dag"
        mydag.save(self.dag_file)

    def create_submit_file(self,run_name):
        submit_file_name = str(run_name + "/" + run_name + ".submit")
        with open(str(submit_file_name),'w') as submit_file:
            submit_file.write("universe = vanilla\n")
            submit_file.write("executable = $(execute)\n")
            # Only run on linux
            submit_file.write("requirements = OpSys == \"LINUX\"\n")
            #submit_file.write("arguments = $(args)\n")
            submit_file.write("log = $(log)\n")
            submit_file.write("out = $(output)\n")
            submit_file.write("err = $(error)\n")
            submit_file.write("request_memory = $(mem)\n")
            submit_file.write("request_cpus = $(cpus)\n")
            submit_file.write("queue\n")
        self.submit_file = submit_file_name

    def submit(self):
        submit_command = list()
        submit_command.append("condor_submit_dag")
        submit_command.append(self.dag_file)
        Popen(submit_command)
