__author__ = 'Michael Lampe'
__version__ = "1.0.0"
__maintainer__ = "Michael Lampe"
__email__ = "MrLampe@Wisc.edu"
__status__ = "Development"

import re
import os
import sys
import numpy as np
try:
    import matplotlib
# Let's us run this without a display
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print("Unable to import matplotlib libraries - efficiency analyzer closing.")
    exit(1)

memory = re.compile("Memory ")
job_name = re.compile("JOB[0-9]*")

class ResourceSummarizer:
    def __init__(self, log_directory, run_name):
        self.log_directory = log_directory
        # Makes sure has an ending slash
        if not self.log_directory.endswith("/"):
            self.log_directory += "/"
        # Absolute file paths
        self.log_files = [(self.log_directory + file) for file in os.listdir(log_directory)]

    def summary(self):
        job_numbers = list()
        efficiency = dict()
        resources = dict()
        for log_file in self.log_files:
            # Get the string "JOB##"
            # It will be used as the key for all of this
            # And also the label for the bars.
            # The try is here in case we get a badly formed log file.
            # This may occur primarily if we run the efficiency analyzer as a member of the DAG, in which case its log
            #  file wouldn't be complete.
            try:
                job_number = self.parse_job_name(log_file)
                job_numbers.append(job_number)
                # Calculate the efficiency of the request vs allocation
                used, expected = self.analyze_log_file(log_file)
                resources[job_number] = (used, expected)
                efficiency[job_number] = self.resource_efficiency(used, expected)
            except:
                pass
        # After we have everything, plot it
        job_numbers = self.natural_sort(job_numbers)
        self.plot_visual(job_numbers, resources, efficiency)

    def resource_efficiency(self, used, expected):
        return float(used)/float(expected) * 100

    def parse_job_name(self, file_name):
        return re.search(job_name, file_name).group()

    def analyze_log_file(self, file_name):
        with open(file_name) as log_file:
            for line in log_file:
                if re.search(memory, line):
                    # Remove junk characters
                    line = re.sub("\t", "", line)
                    line = re.sub("\n", "", line)
                    # Remove all the spaces and convert the rest to an array
                    value = filter(bool, re.split(" ", line))
                    # Just pick out the numbers
                    value = [int(val) for val in value if val.isdigit()]
                    used, expected = value[0], value[1]
                    return used, expected

    def plot_visual(self, job_numbers, resources, efficiency):
        def determine_max_ur(list1, list2):
            a , b = max(list1), max(list2)
            if a>b:
                return a
            return b

        def determine_max_eff(list1):
            a = max(list1)
            if a < 100:
                return 100
            return a

        size = np.arange(len(job_numbers))
        figure = plt.figure(figsize=(20 ,20), dpi=600)

        # Subplot 1 - Usage and request side by side
        axis = figure.add_subplot(211)
        width = 0.25
        # Order the data by job number (Incoming job numbers are ordered)
        usage = [resources[job_number][0] for job_number in job_numbers]
        request = [resources[job_number][1] for job_number in job_numbers]
        # Assign data
        u = axis.bar(size, usage,width, color="blue")
        r = axis.bar(size+width, request, width,color="red")
        # Make the graph informative
        axis.legend((u[0], r[0]),("Usage", "Request"))
        axis.set_ylim(0, determine_max_ur(usage, request))
        axis.set_title("Usage of and Requested Memory")
        axis.set_ylabel("Memory (MB)")

        # Subplot 2 - Efficiency calculation by step
        axis = figure.add_subplot(212)
        # Order data by job number
        eff = [efficiency[job_number] for job_number in job_numbers]
        # Assign data
        axis.bar(size, eff, color="purple")
        # Make the graph informative
        axis.set_ylim(0, determine_max_eff(eff))
        axis.set_title("Percent of Memory Utilized based on Requested Resources")
        axis.set_ylabel("Usage of Memory Requested (%)")

        plt.savefig(self.run_name + "/ResourceUtilizationSummary")

    # These two functions i got from SO.
    # They naturally sort a string (Numbers give order)
    # http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    def natural_sort(self, list1):
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(list1, key=alphanum_key)

# Input args
if len(sys.argv) < 1:
    print "Error, insufficient arguments"
    exit(1)
else:
    ResourceSummarizer(sys.argv[1], sys.argv[2]).summary()
