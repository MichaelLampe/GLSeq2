__author__ = 'mlampe'

# Pygraph library
# Modified from https://github.com/pmatiello/python-graph
from pygraph.classes.digraph import digraph as graph


class GraphCondor:
    def __init__(self,dag_file_name):
        self.dag_file_name = dag_file_name
        self.condor_graph = graph()

    def read_connections(self):
        with open(self.dag_file_name,'r'):
            for line in file:
                if "PARENT" in line:
                    # Splits the line and the left space is the string "PARENT"
                    # with the right half being the value of that
                    split = line.split(" ")
                    parent = split[1]
                    child = split[3]
                    self.create_connection(parent,child)
    def create_connection(self,parent,child):
        # Adds the parent and child nodes if they don't exist already
        [self.condor_graph.add_node(n) for n in (parent,child) if not self.condor_graph.has_node(n)]
        # Adds the edge from parent to child
        self.condor_graph.add_edge(parent,child)

    # Still need to add visualization