__author__ = 'Michael Lampe'

# Pydagman can be found at https://github.com/brandentimm/pydagman
from pydagman.dagfile import Dagfile
from pydagman.job import Job
from subprocess import Popen
from subprocess import PIPE
from CommandFile import CommandFile
import networkx as nx
import os

"""
Grabs a given file path based on a command
"""
def get_paths(command):
    # Grap the output of the echo $PATH command
    output = Popen(command, shell=True, stdout=PIPE)
    path_out, err = output.communicate()
    # Make sure connection is dead.
    try:
        output.kill()
    except:
        pass
    return path_out

"""
Setup the condor environment with any relevant language path.
"""
def get_environment():
    path = get_paths('echo $PATH').strip("\n")
    python_path = get_paths('echo $PYTHONPATH').strip("\n")
    r_path = get_paths('echo $R_LIBS').strip("\n")
    perl_path = get_paths('echo $PERL5LIB').strip("\n")
    environment = "PATH=" + path + ";"
    environment = environment + "PYTHONPATH=" + python_path + ";"
    environment = environment + "R_LIBS=" + r_path + ";"
    environment = environment + "PERL5LIB=" + perl_path
    return environment

"""
Creates the given directories for Condor
"""
def create_condor_directories(run_name, condor_path):
    # The names of all the directories
    run = condor_path + str(run_name) + "/"
    log_file_dir = run + "LogFile/"
    out_file_dir = run + "OutputFile/"
    err_file_dir = run + "ErrorFile/"

    # Makes sure that they are already made
    if not os.path.exists(os.path.dirname(log_file_dir)):
        os.makedirs(os.path.dirname(log_file_dir))
    if not os.path.exists(os.path.dirname(out_file_dir)):
        os.makedirs(os.path.dirname(out_file_dir))
    if not os.path.exists(os.path.dirname(err_file_dir)):
        os.makedirs(os.path.dirname(err_file_dir))

    # Returns all of their paths
    return log_file_dir, out_file_dir, err_file_dir


class Stack:
    def __init__(self, graph, condor_path):
        # Graph containing all the commands after they have been organized
        self.graph = graph

        self.condor_path = condor_path

        # Keeps track of the shell files created
        self.associated_bash_files = dict()
        # Keeps track of the dog jobs
        self.dag_jobs = dict()

    def plot_graph(self, run_name):
        self.graph.plot_graph(run_name, self.condor_path)

    def create_stack(self, run_name):
        self._create_bash_files(run_name)
        self.create_dag_jobs(run_name)
        for node in self.graph.nodes():
            parents = self.graph.get_parents(node)
            for parent in parents:
                self.dag_jobs[node].add_parent(self.dag_jobs[parent])
        self.create_dag_workflow(run_name)
        self.create_submit_file(run_name)

    def _create_bash_files(self, run_name):
        """
        Generates an individual BASH file.  This can be a group of commands or just one.
        Each loop here is an individual parallel job at this step
        """
        file_number = 1

        # Topological sort is great for ordering these files in the order they will execute.
        for node in nx.topological_sort(self.graph.G):
            new_command = CommandFile(run_name, self.condor_path, node)
            bash_file = new_command.generate_bash_file()
            self.associated_bash_files[node] = (file_number, bash_file)
            file_number += 1

    def create_dag_jobs(self, run_name):
        # Creates the directories on the file system if they are not already created
        self.log_file_dir, self.out_file_dir, self.err_file_dir = create_condor_directories(run_name, self.condor_path)

        # Topological sort is great for ordering these files in the order they will execute.
        for node in nx.topological_sort(self.graph.G):
            output_file_name = str(run_name) + "_JOB" + str(self.associated_bash_files[node][0])

            # If to run with a GPU or not (Moreso if to request one)
            if int(self.graph.G.node[node]['gpus']) >= 1:
                current_job = Job((self.condor_path + run_name + "/" + run_name + ".gpu.submit"), "JOB" +
                                  str(self.associated_bash_files[node][0]))
            else:
                current_job = Job((self.condor_path + run_name + "/" + run_name + ".submit"), "JOB" +
                                  str(self.associated_bash_files[node][0]))

            """
            Adding memory request and cpu request
            Because we create a bunch of shell files that are run, there are not args.
            """
            self.static_job_values(current_job, output_file_name)
            current_job.add_var("args", "")
            current_job.add_var("mem", self.graph.G.node[node]['mem'])
            current_job.add_var("cpus", self.graph.G.node[node]['cpus'])

            # If a GPU file, then add the number of GPUs requested
            if int(self.graph.G.node[node]['gpus']) >= 1:
                current_job.add_var("gpus", self.graph.G.node[node]['gpus'])

            # Run the created SH file

            current_job.add_var("execute", self.condor_path + run_name + "/" + self.associated_bash_files[node][1])
            # Still need to add parent interactions, which is done in the comm stack
            self.dag_jobs[node] = current_job

    def static_job_values(self, job, file_name):
        job.add_var("log", self.log_file_dir + file_name + "_LogFile.txt")
        job.add_var("output", self.out_file_dir + file_name  + "_OutputFile.txt")
        job.add_var("error", self.err_file_dir + file_name + "_ErrorFile.txt")
        job.retry(1)
        job.pre_skip("1")

    def create_dag_workflow(self, run_name):
        mydag = Dagfile()

        # Topological sort is great for ordering these files in the order they will execute.
        for node in nx.topological_sort(self.graph.G):
            mydag.add_job(self.dag_jobs[node])
        self.dag_file = self.condor_path + run_name + "/my_workflow.dag"
        mydag.save(self.dag_file)

    def create_submit_file(self, run_name):
        # Just create both submit files as this makes sure that the GPU condition is handled properly
        submit_file_name = str(self.condor_path + run_name + "/" + run_name + ".submit")
        environment = get_environment()
        with open(str(submit_file_name), 'w') as submit_file:
            submit_file.write("universe = vanilla\n")
            submit_file.write("executable = $(execute)\n")
            submit_file.write("arguments = $(args)\n")
            submit_file.write("log = $(log)\n")
            submit_file.write("out = $(output)\n")
            submit_file.write("err = $(error)\n")
            submit_file.write("request_memory = $(mem)\n")
            submit_file.write("request_cpus = $(cpus)\n")
            submit_file.write("environment = " + environment + "\n")
            submit_file.write("queue\n")

        self.submit_file = submit_file_name
        # This is a special submit file that allows GPUs
        submit_file_name = str(self.condor_path + run_name + "/" + run_name + ".gpu.submit")
        environment = get_environment()
        with open(str(submit_file_name), 'w') as submit_file:
            submit_file.write("universe = vanilla\n")
            submit_file.write("executable = $(execute)\n")
            submit_file.write("arguments = $(args)\n")
            submit_file.write("log = $(log)\n")
            submit_file.write("out = $(output)\n")
            submit_file.write("err = $(error)\n")
            submit_file.write("request_memory = $(mem)\n")
            submit_file.write("request_cpus = $(cpus)\n")
            submit_file.write("request_GPUs = $(gpus)\n")

            # This will only work on our (GLBRC) file systems.
            submit_file.write("Requirements = TARGET.FileSystemDomain == \"glbrc.org\"\n")
            submit_file.write("environment = " + environment + "\n")
            submit_file.write("queue\n")
        self.submit_file = submit_file_name

    def submit(self):
        submit_command = list()
        submit_command.append("condor_submit_dag")
        submit_command.append(self.dag_file)
        # Runs but doesn't wait for a connection so we can keep iterating.
        Popen(submit_command, stdin=None, stdout=None, stderr=None, close_fds=True)

