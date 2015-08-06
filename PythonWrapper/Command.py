__author__ = 'mlampe'

from CommandFile import CommandFile
# Pydagman can be found at https://github.com/brandentimm/pydagman
from pydagman.job import Job
import os
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
        # Creates the directories on the file system if they are not already created
        log_file_dir,out_file_dir, err_file_dir = create_condor_directories(run_name)
        for subcommand in range(0,len(self.commands)):
            #
            output_file_name = str(run_name) + "_JOB" + str(Command.job_count)
            current_job = Job((run_name + "/" + run_name +".submit"),"JOB" + str(Command.job_count))
            current_job.add_var("log",log_file_dir + output_file_name + "_LogFile.txt")
            current_job.add_var("output",out_file_dir + output_file_name  + "_OutputFile.txt")
            current_job.add_var("error",err_file_dir + output_file_name + "_ErrorFile.txt")
            self.assign_resources(current_job)
            current_job.add_var("execute","./" + run_name + "/" + self.associated_bash_files[subcommand])
            current_job.retry(1)
            current_job.pre_skip("1")
            # Still need to add parent interactions, which is done in the comm stack
            Command.job_count = Command.job_count + 1
            self.dag_jobs.append(current_job)
        for job in range(0,len(self.dag_jobs)):
            if (job != 0):
                self.dag_jobs[job].add_parent(self.dag_jobs[job-1])
    # Assign memory,cpu, and GPU resources to processes
    def assign_resources(self,job):
        # Scarce resource is 256 mb and 1 cpu and 0 gpus
        scarce_resource = [
            "256", # Memory
            "1", # CPUS
            "0" # GPUS
        ]
        # Low resource is 1GB and 1 cpu and 0 gpus
        low_resource = [
            "1GB", # Memory
            "1", # CPUS
            "0" # GPUS
        ]
        # Medium resource is 4GB and 2 cpus and 0 gpus
        medium_resource = [
            "4G",
            "4",
            "0"
        ]
        # High resource is 8GB and 8 cpu and 0 gpus
        high_resource = [
            "8G",
            "8",
            "0"
        ]
        # High resource is 8GB and 6 cpu and 0 gpus
        gpu_needed = [
            "8G",
            "6", # As long as we only use CUSHAW-GPU keep this here (Tested as the fastest GPU-CPU combo)
            "1"
        ]
        resources = medium_resource
        job.add_var("mem",resources[0])
        job.add_var("cpus",resources[1])
        job.add_var("gpus",resources[2])

# Statics
def create_condor_directories(run_name):
    # The names of all the directories
    log_file_dir =str(run_name) + "/" + "LogFile/"
    out_file_dir=str(run_name) + "/" + "OutputFile/"
    err_file_dir = str(run_name) + "/" + "ErrorFile/"
    # Makes sure that they are already made
    if not os.path.exists(os.path.dirname(log_file_dir)):
        os.makedirs(os.path.dirname(log_file_dir))
    if not os.path.exists(os.path.dirname(out_file_dir)):
        os.makedirs(os.path.dirname(out_file_dir))
    if not os.path.exists(os.path.dirname(err_file_dir)):
        os.makedirs(os.path.dirname(err_file_dir))
    # Returns all of their paths
    return log_file_dir,out_file_dir,err_file_dir