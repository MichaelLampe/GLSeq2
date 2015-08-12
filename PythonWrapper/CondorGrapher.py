__author__ = 'mlampe'

import networkx as nx
import matplotlib
# Let's us run this without a display
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy

class Graph:
    def __init__(self):
        self.G = nx.Graph()
    def add_node(self,node,attributes):
        # Indexed by its number
        self.G.add_node(node.number,attributes)

    def add_edge(self,parent,child):
        self.G.add_edge(child,parent)

    def plot_graph(self,run_name):
        # Graph one plots the memory utilization of various steps
        plt.figure(1,figsize=(30,30))
        plt.axis("off")
        plt.title(run_name + "/" + "Condor Memory Utilization Scheme")
        pos=nx.graphviz_layout(self.G)
        nx.draw_networkx(self.G,pos,k=10,with_labels=False,node_size = 1500,node_color=self.memory_colors())
        nx.draw_networkx_labels(self.G,pos,labels=self.memory_labels())
        plt.savefig(run_name + "/" + "CondorMemoryUtilization.png")
        plt.clf()

        # Graph two plots the CPU utilization of various steps
        plt.figure(2,figsize=(30,30))
        plt.axis("off")
        plt.title("Condor CPU Utilization Scheme")
        pos=nx.graphviz_layout(self.G)
        nx.draw_networkx(self.G,pos,k=10,with_labels=False,node_size = 1500,node_color=self.cpu_colors())
        nx.draw_networkx_labels(self.G,pos,self.cpu_labels())
        plt.savefig(run_name + "/" + "CondorCpuUtilization.png")
        plt.clf()

        # Graph three plots the GPU utilization of various steps
        plt.figure(3,figsize=(30,30))
        plt.axis("off")
        plt.title("Condor GPU Utilization Scheme")
        pos=nx.graphviz_layout(self.G)
        nx.draw_networkx(self.G,pos,k=10,with_labels=False,node_size = 1500,node_color=self.gpu_colors())
        nx.draw_networkx_labels(self.G,pos,self.gpu_labels())
        plt.savefig(run_name + "/" + "CondorGpuUtilization.png")
        plt.clf()

    def memory_labels(self):
        labels={}
        for node in self.G:
            labels[node] = node + "\n" + self.G.node[node]['memory']
        return labels

    def memory_colors(self):
        memory_color = []
        for node in self.G:
            if (self.G.node[node]["memory"][-1] == "G"):
                memory = 1024*int(self.G.node[node]["memory"][:-1])
            else:
                # Assumes we will never use K and only have M (Mega) other than G
                memory = int(self.G.node[node]["memory"][:-1])
            # Might want to do something more clever here to get a better gradient
            if memory < 10:
                memory_color.append("#0000FF")
            elif memory < 100:
                memory_color.append("#1919FF")
            elif memory < 300:
                memory_color.append("#3333FF")
            elif memory < 800:
                memory_color.append("#4D4DFF")
            elif memory < 1500:
                memory_color.append("#6666FF")
            elif memory < 2000:
                memory_color.append("8080FF")
            elif memory < 3000:
                memory_color.append("#9999FF")
            elif memory < 4000:
                memory_color.append("#CCCCFF")
            elif memory < 5000:
                memory_color.append("#FFE6FF")
            elif memory < 6000:
                memory_color.append("#FFCCCC")
            elif memory < 7000:
                memory_color.append("#FF9999")
            elif memory < 8000:
                memory_color.append("#FF4D4D")
            elif memory < 9000:
                memory_color.append("#FF6666")
            elif memory < 10000:
                memory_color.append("#FF1919")
            else:
                memory_color.append("#FF0000")
        return memory_color

    def cpu_labels(self):
        labels={}
        for node in self.G:
            labels[node] = node + "\n" + self.G.node[node]['cpu']
        return labels

    def cpu_colors(self):
        cpu_color = []
        for node in self.G:
            cpu = int(self.G.node[node]["cpu"])
            if cpu < 1:
                cpu_color.append("#0000FF")
            elif cpu < 2:
                cpu_color.append("#3333FF")
            elif cpu < 3:
                cpu_color.append("8080FF")
            elif cpu < 4:
                cpu_color.append("#9999FF")
            elif cpu < 5:
                cpu_color.append("#FFE6FF")
            elif cpu < 6:
                cpu_color.append("#FFCCCC")
            elif cpu < 7:
                cpu_color.append("#FF9999")
            elif cpu < 8:
                cpu_color.append("#FF4D4D")
            elif cpu < 9:
                cpu_color.append("#FF6666")
            elif cpu < 10:
                cpu_color.append("#FF1919")
            else:
                cpu_color.append("#FF0000")
        return cpu_color

    def gpu_labels(self):
        labels={}
        for node in self.G:
            labels[node] = node + "\n" + self.G.node[node]['gpu']
        return labels

    def gpu_colors(self):
        gpu_color = []
        for node in self.G:
            gpu = self.G.node[node]["gpu"]
            # Tuple of RGB
            if gpu < 1:
                gpu_color.append("#FFE6FF")
            elif gpu == 1:
                gpu_color.append("#FF9999")
            elif gpu == 2:
                gpu_color.append("#FF4D4D")
            elif gpu == 3:
                gpu_color.append("#FF6666")
            else:
                gpu_color.append("#FF0000")
        return gpu_color

    def prepare_colors(self,attribute):
        colors = list()
        if (attribute == "memory"):
            for node in Node.nodes:
                # If ends with gigabyte
                if(node.memory.endswith("G")):
                    memory = node.memory[:-1]
                    memory = int(memory) * 1024
                    colors.append(memory)
                else:
                    memory = int(node.memory)
                    colors.append(memory)
            median = numpy.median(colors)
            for color in range (0,len(colors)):
                if (colors[color] > median + 1025):
                    colors[color] = "r"
                elif (colors[color] < median - 1025):
                    colors[color] = "b"
                else:
                    colors[color] = "y"
        elif (attribute == "cpu"):
            for node in Node.nodes:
                # If ends with gigabyte
                    cpus = int(node.cpu)
                    colors.append(cpus)
            median = numpy.median(colors)
            for color in range (0,len(colors)):
                if (colors[color] > median + 4):
                    colors[color] = "r"
                elif (colors[color] < median - 4):
                    colors[color] = "b"
                else:
                    colors[color] = "y"
        elif (attribute == "gpu"):
            for node in Node.nodes:
                # If ends with gigabyte
                    gpus = int(node.gpu)
                    colors.append(gpus)
            median = numpy.median(colors)
            for color in range (0,len(colors)):
                if (colors[color] > median + 1):
                    colors[color] = "r"
                elif (colors[color] < median - 1):
                    colors[color] = "b"
                else:
                    colors[color] = "y"
        return colors

class Node:
    nodes = list()
    def __init__(self,number,memory,cpu,gpu,command_stack):
        self.number = str(number)
        # self.memory = str(memory)
        # self.cpu = str(cpu)
        # self.gpu= str(gpu)
        attributes = {
            'memory' : memory,
            'cpu' : cpu,
            'gpu' : gpu
        }
        command_stack.graph.add_node(self,attributes)
        Node.nodes.append(self)
