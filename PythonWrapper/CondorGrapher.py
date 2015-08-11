__author__ = 'mlampe'

import networkx as nx
import matplotlib
# Let's us run this without a display
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import CommandStack
import numpy

class Graph:
    def __init__(self):
        self.G = nx.Graph()

    def add_node(self,node):
        # Indexed by its number
        self.G.add_node(node.number,Memory=node.memory,CPU=node.cpu,GPU=node.gpu)

    def add_edge(self,parent,child):
        self.G.add_edge(child,parent)

    def plot_graph(self,run_name):
        plt.figure(1,figsize=(15,30))
        limit = plt.axis("off")
        title= plt.title(run_name + "/" + "Condor Memory Utilization Scheme")
        pos=nx.spring_layout(self.G,scale=2)
        nx.draw_networkx(self.G,pos,with_labels=True,node_size = 1500,node_color=self.prepare_colors("memory"))
        plt.savefig(run_name + "/" + "CondorMemoryUtilization.png")
        plt.clf()
        plt.figure(2,figsize=(15,30))
        limit = plt.axis("off")
        title= plt.title("Condor CPU Utilization Scheme")
        pos=nx.spring_layout(self.G,scale=2)
        nx.draw_networkx(self.G,pos,with_labels=True,node_size = 1500,node_color=self.prepare_colors("cpu"))
        plt.savefig(run_name + "/" + "CondorCpuUtilization.png")
        plt.clf()
        plt.figure(3,figsize=(15,30))
        limit = plt.axis("off")
        title= plt.title("Condor GPU Utilization Scheme")
        pos=nx.spring_layout(self.G,scale=2)
        nx.draw_networkx(self.G,pos,with_labels=True,node_size = 1500,node_color=self.prepare_colors("gpu"))
        plt.savefig(run_name + "/" + "CondorGpuUtilization.png")
        plt.clf()

    def prepare_colors(self,attribute):
        colors = list()
        if (attribute == "memory"):
            for node in Node.Nodes:
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
            for node in Node.Nodes:
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
            for node in Node.Nodes:
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
    Nodes = list()
    def __init__(self,number,memory,cpu,gpu,command_stack):
        self.number = number
        self.memory = memory
        self.cpu = cpu
        self.gpu = gpu
        command_stack.graph.add_node(self)
        Node.Nodes.append(self)