__author__ = 'mlampe'

from CommandFile import CommandFile
# Pydagman can be found at https://github.com/brandentimm/pydagman
from pydagman.job import Job

class Command():
    parallel_track = list()
    job_count = 1
    def __init__(self,commands,parent_node):
        self.parent_node = parent_node
        Command.parallel_track.append(list())
        Command.parallel_track[self.parent_node].append(self)
        # A list of commands that can be run in parallel and sequentially
        self.commands = commands
        self.associated_bash_files = list()
        self.dag_jobs = list()

    # To string method I used for debugging of individual command nodes
    def toString(self):
        return("The command is " + str(self.commands) + "." + " This command has parent node of " + str(self.parent_node) + ".")

    # Creates the bash files that will be addedd to the job flow
    def create_bash_files(self,run_name):
        for command in self.commands:
            new_command = CommandFile(run_name,command,self.parent_node)
            bash_file = self.command_file = new_command.generate_bash_file()
            self.associated_bash_files.append(bash_file)

    # Sets up dag jobs.
    # Will need to add a controller in the future to optimize mem and cpu usage, as well as GPU.
    # This will really help when we are mass scheduling batch jobs.
    def create_dag_jobs(self,run_name):
        for subcommand in range(0,len(self.commands)):
            current_job = Job((run_name + "/" + run_name +".submit"),"JOB" + str(Command.job_count))
            current_job.add_var("mem","4G")
            current_job.add_var("cpus","2")
            current_job.add_var("execute","./" + run_name + "/" + self.associated_bash_files[subcommand])
            current_job.retry(2)
            current_job.pre_skip("1")
            # Still need to add parent interactions, which is done in the comm stack
            Command.job_count = Command.job_count + 1
            self.dag_jobs.append(current_job)
        for job in range(0,len(self.dag_jobs)):
            if (job != 0):
                self.dag_jobs[job].add_parent(self.dag_jobs[job-1])
