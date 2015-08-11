__author__ = 'mlampe'

from CommandFile import CommandFile
# Pydagman can be found at https://github.com/brandentimm/pydagman
from pydagman.job import Job
import os
import re
from re import search, compile
from CondorGrapher import Node
import CommandStack
class Command():
    parallel_track = list()
    job_count = 1
    def __init__(self,commands,parent_node):
        self.parent_node = parent_node
        Command.parallel_track.append(list())
        Command.parallel_track[self.parent_node].append(self)
        # A list of commands that can be run in parallel and sequentially
        self.commands = commands
        self.dag_jobs = list()
        self.associated_bash_files = list()

    # To string method I used for debugging of individual command nodes
    def toString(self):
        return("The command is " + str(self.commands) + "." + " This command has parent node of " + str(self.parent_node) + ".")

    # Creates the bash files that will be addedd to the job flow
    def create_bash_files(self,run_name):
        # Generates an individual BASH file.  This can be a group of commands or just one.
        # Each loop here is an individual parallel job at this step
        for command in self.commands:
            new_command = CommandFile(run_name,command,self.parent_node)
            bash_file = self.command_file = new_command.generate_bash_file()
            self.associated_bash_files.append(bash_file)

    # Sets up dag jobs.
    # Will need to add a controller in the future to optimize mem and cpu usage, as well as GPU.
    # This will really help when we are mass scheduling batch jobs.
    def create_dag_jobs(self,run_name,graph,command_stack):
        # Creates the directories on the file system if they are not already created
        log_file_dir,out_file_dir, err_file_dir = create_condor_directories(run_name)
        for subcommand in range(0,len(self.commands)):
            mem, cpus, gpus = self.assign_resources(self.commands[subcommand])
            output_file_name = str(run_name) + "_JOB" + str(Command.job_count)
            if (int(gpus) >= 1):
                current_job = Job((run_name + "/" + run_name +".gpu.submit"),"JOB" + str(Command.job_count))
            else:
                current_job = Job((run_name + "/" + run_name +".submit"),"JOB" + str(Command.job_count))
            current_job.add_var("log",log_file_dir + output_file_name + "_LogFile.txt")
            current_job.add_var("output",out_file_dir + output_file_name  + "_OutputFile.txt")
            current_job.add_var("error",err_file_dir + output_file_name + "_ErrorFile.txt")
            # Creates a new node for visualization
            Node(current_job.name,mem,cpus,gpus,command_stack)
            current_job.add_var("mem",mem)
            current_job.add_var("cpus",cpus)
            current_job.add_var("gpus",gpus)
            current_job.add_var("execute","./" + run_name + "/" + self.associated_bash_files[subcommand])
            current_job.retry(1)
            current_job.pre_skip("1")
            # Still need to add parent interactions, which is done in the comm stack
            Command.job_count = Command.job_count + 1
            self.dag_jobs.append(current_job)
        for job in range(0,len(self.dag_jobs)):
            if (job != 0):
                graph.add_edge(self.dag_jobs[job].name,self.dag_jobs[job-1].name)
                self.dag_jobs[job].add_parent(self.dag_jobs[job-1])

    # Assign memory,cpu, and GPU resources to processes
    def assign_resources(self,command):
        def found_match(pattern,command):
            x = search(pattern,command)
            return x!=None

        # Scarce resource is 256 mb and 1 cpu and 0 gpus
        scarce_resource = [
            "256", # Memory
            "1", # CPUS
            "-1" # GPUS
            # return scarce_resource[0],scarce_resource[1],scarce_resource[2]
        ]
        # Low resource is 1GB and 1 cpu and 0 gpus
        low_resource = [
            "1G", # Memory
            "1", # CPUS
            "-1" # GPUS
            # return low_resource[0],low_resource[1],low_resource[2]
        ]
        # Medium resource is 4GB and 2 cpus and 0 gpus
        medium_resource = [
            "4G",
            "4",
            "-1"
            # return medium_resource[0],medium_resource[1],medium_resource[2]
        ]
        # High resource is 8GB and 8 cpu and 0 gpus
        high_resource = [
            "8G",
            "8",
            "-1"
            # return high_resource[0],high_resource[1],high_resource[2]
        ]
        # Likely will never be used, but it is just another case
        extreme_resource = [
            "16G",
            "16"
            "-1"
            # return extreme_resource[0],extreme_resource[1],extreme_resource[2]
        ]
        # High Ram uses 32 GB of Ram and 6 cores
        # Primarily planned for use with Samtools sort
        high_ram = [
            "32G",
            "6",
            "-1"
        ]
        gpu_needed = [
            "8G",
            "6", # As long as we only use CUSHAW-GPU keep this here (Tested as the fastest GPU-CPU combo)
            "1"
            # return gpu_needed[0],gpu_needed[1],gpu_needed[2]
        ]

        # Default
        default = medium_resource[0],medium_resource[1],medium_resource[2]

        #
        # The ^(.*?) ____ (?=[ ]) regex is there so that we only match the first word
        # As we sometimes use the aligners/counters as run names and that can
        # cause a lot of things to ask for a lot of resources
        cushaw_gpu = compile("^(.*?)cushaw2-gpu(?=[ ])",re.IGNORECASE)
        if found_match(cushaw_gpu,command):
            return gpu_needed[0],gpu_needed[1],gpu_needed[2]
        # This is often the slowest step, so lets give it a lot of RAM to work with
        samtools_sort = compile("samtools sort",re.IGNORECASE)
        if found_match(samtools_sort,command):
            return high_ram[0],high_ram[1],high_ram[2]
        rsem = compile("^(.*?)rsem(?=[ ])",re.IGNORECASE)
        # RSEM uses sort and also takes a long time.
        if found_match(rsem,command):
            return high_ram[0],high_ram[1],high_ram[2]
        # High Resource
        cushaw = compile("^(.*?)cushaw(?=[ ])",re.IGNORECASE)
        if found_match(cushaw,command):
            return high_resource[0],high_resource[1],high_resource[2]
        bwa = compile("^(.*?)bwa(?=[ ])",re.IGNORECASE)
        if found_match(bwa,command):
            return high_resource[0],high_resource[1],high_resource[2]
        bowtie = compile("^(.*?)bowtie(?=[ ])",re.IGNORECASE)
        if found_match(bowtie,command):
            return high_resource[0],high_resource[1],high_resource[2]
        tophat = compile("^(.*?)tophat(?=[ ])",re.IGNORECASE)
        if found_match(tophat,command):
            return high_resource[0],high_resource[1],high_resource[2]
        cufflinks = compile("^(.*?)cufflinks(?=[ ])",re.IGNORECASE)
        if found_match(cufflinks,command):
            return high_resource[0],high_resource[1],high_resource[2]
        # Medium Resources
        fastqc = compile("^(.*?)fastqc(?=[ ])",re.IGNORECASE)
        if found_match(fastqc,command):
            return medium_resource[0],medium_resource[1],medium_resource[2]
        bam2wig = compile("^(.*?)bam2wig(?=[ ])",re.IGNORECASE)
        if found_match(bam2wig,command):
            return medium_resource[0],medium_resource[1],medium_resource[2]
        picardTools = compile("^(.*?)picard-tools(?=[ ])",re.IGNORECASE)
        if found_match(picardTools,command):
            return medium_resource[0],medium_resource[1],medium_resource[2]
        trim = compile("^(.*?)trimmomatic(?=[ ])",re.IGNORECASE)
        if found_match(trim,command):
            return medium_resource[0],medium_resource[1],medium_resource[2]
        rockhopper = compile("^(.*?)rockhopper(?=[ ])",re.IGNORECASE)
        if found_match(rockhopper,command):
            return medium_resource[0],medium_resource[1],medium_resource[2]
        # bifxapps isn't an app, but is our app folder.
        # This should catch any absolute path someone introduces to try and assign it
        #  to a "medium resources" setting (Assuming the app is in bifxapps)
        bifxapps = compile("^(.*?)bifxapps(?=[ ])",re.IGNORECASE)
        if found_match(bifxapps,command):
            return medium_resource[0],medium_resource[1],medium_resource[2]
        # Low Resource
        samtools = compile("^(.*?)samtools(?=[ ])",re.IGNORECASE)
        if found_match(samtools,command):
            return low_resource[0],low_resource[1],low_resource[2]
        # Scarce Resource
        htseq = compile("^(.*?)HTSeq(?=[ ])",re.IGNORECASE)
        if found_match(htseq,command):
            return scarce_resource[0],scarce_resource[1],scarce_resource[2]
        featureCounts = compile("^(.*?)featurecounts(?=[ ])",re.IGNORECASE)
        if found_match(featureCounts,command):
            return scarce_resource[0],scarce_resource[1],scarce_resource[2]
        # If nothing else has returned by now, return whatever was set as the defautl value.
        return default
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