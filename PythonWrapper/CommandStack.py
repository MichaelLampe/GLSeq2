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
        for steps in self.command_list:
            for step in steps:
                for x in range (len(step) - 1 ,0,-1):
                    for word in self.group_keywords:
                        if (word in step[x - 1]):
                            if (x != 0):
                                step[x - 1] = step[x-1] + " && " + step[x]
                                step.pop(x)
                                # Exit after you find one because that's all good
                                break
        for parent in range (0, len(self.command_list)):
            for step in self.command_list[parent]:
                command = Command(step,parent)
                self.command_stack.append(command)

        for command in self.command_stack:
            command.create_bash_files(run_name)
            self.create_submit_file(run_name)
            command.create_dag_jobs(run_name)

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
            #submit_file.write("arguments = $(args)\n")
            submit_file.write("log = job.log\n")
            submit_file.write("out = job.out\n")
            submit_file.write("err = job.err\n")
            submit_file.write("request_memory = $(mem)\n")
            submit_file.write("request_cpus = $(cpus)\n")
            submit_file.write("queue\n")
        self.submit_file = submit_file_name

    def submit(self):
        submit_command = list()
        submit_command.append("condor_submit_dag")
        submit_command.append(self.dag_file)
        Popen(submit_command)

