__author__ = 'mlampe'

import networkx as nx
import matplotlib
# Let's us run this without a display
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re


# This is a wrapper for the networkx Digraph class
# This is the structure that will hold onto the commands data
class Graph:
    def __init__(self):
        # Digraph allows directionality which is key
        self.G = nx.DiGraph()

    def draw_graph(self):
        pos = nx.graphviz_layout(self.G)
        nx.draw_networkx(self.G, pos)

    def nodes(self):
        return self.G.nodes()

    def add_node(self, new_node, parents=None, attributes=None):
        # Indexed by its number
        # Only add an attribute if it is worth adding one, aka there is something there
        if parents is None:
            if attributes is None:
                self.G.add_node(new_node)
            else:
                self.G.add_node(new_node, attributes)
        else:
            if attributes is None:
                self.G.add_node(new_node)
            else:
                self.G.add_node(new_node, attributes)
            self.add_edge(parents, new_node)

    def get_root_nodes(self):
        # We'll only do this when we ask for it so we aren't doing this a ton
        root_nodes = [node for node, degree in self.G.nodes().in_degree().items() if degree == 0]
        return root_nodes

    def add_edge(self, parents, children):
        if type(parents) is not list:
            if type(children) is not list:
                self.G.add_edge(parents, children)
            else:
                for child in children:
                    self.G.add_edge(parents, child)
        else:
            if type(children) is not list:
                for parent in parents:
                    self.G.add_edge(parent, children)
            else:
                for parent in parents:
                    for child in children:
                        self.G.add_edge(parent, child)

    def remove_edge(self, parents, children):
        if type(parents) is not list:
            if type(children) is not list:
                self.G.remove_edge(parents, children)
            else:
                for child in children:
                    self.G.remove_edge(parents, child)
        else:
            if type(children) is not list:
                for parent in parents:
                    self.G.remove_edge(parent, children)
            else:
                for parent in parents:
                    for child in children:
                        self.G.remove_edge(parent, child)

    def insert_node(self, parents, new_node):
        if type(parents) is list:
            children = [child for parent in parents for child in self.get_children(parent)]
            children = sorted(set(children))
        else:
            children = self.get_children(parents)
        self.add_edge(parents, new_node)
        self.add_edge(new_node, children)
        self.remove_edge(parents, children)

    def remove_node(self, node):
        parents , children = self.get_parents(node) , self.get_children(node)
        # Links the gap that would be created by removing the node
        self.add_edge(parents,children)
        self.remove_edge(parents,node)
        self.remove_edge(node,children)
        self.G.remove_node(node)

    # Different from remove node as it makes no attempt to bridge the gap.
    def delete_node(self,node):
        self.G.remove_node(node)

    def get_parents(self, node):
        return self.G.predecessors(node)

    def get_children(self, node):
        return self.G.successors(node)

    # This will split the nodes such that they will all be at the same level after splitting
    def split_node_parallel(self, node, new_nodes):
        parents, children = self.get_parents(node), self.get_children(node)
        # Assign the parents and children of the initial node to all the newly created nodes
        for new_node in new_nodes:
            for parent in parents:
                self.G.add_edge(parent, new_node)
            for child in children:
                self.G.add_edge(new_node, child)
        # Remove the initial node so that it doesn't stay in the graph
        self.remove_node(node)

    # This will split the nodes such that the first new node is the parent of the next new node and so on
    def split_node_sequential(self, node, new_nodes):
        parents, children = self.get_parents(node), self.get_children(node)
        # Assign the parents and children of the initial node to all the newly created nodes
        for n in range(0, len(new_nodes)):
            # Head node
            if n == 0:
                for parent in parents:
                    self.G.add_edge(parent, new_nodes[n])
            # General case where we put the node prior as the parent of the current node
            elif n >= 1:
                self.G.add_edge(new_nodes[n - 1], new_nodes[n])
            # Tail node
            if n == len(new_nodes) - 1:
                for child in children:
                    self.G.add_edge(new_nodes[n], child)
        # Remove the initial node so that it doesn't stay in the graph
        self.remove_node(node)

    def merge_nodes(self,nodes, new_data):
        parents = [parent for node in nodes for parent in self.get_parents(node) if parent not in nodes]
        children = [child for node in nodes for child in self.get_children(node) if child not in nodes]
        # Setup the new node
        merged_node = Node(new_data)
        self.add_node(merged_node)
        self.add_edge(parents,merged_node)
        self.add_edge(merged_node,children)
        # Remove all the old nodes
        for node in nodes:
            try:
                self.delete_node(node)
            except:
                pass

    # Stuff for drawing usage stats
    def plot_graph(self, run_name):
        # Graph one plots the memory utilization of various steps
        plt.figure(1, figsize=(30, 30))
        plt.axis("off")
        plt.title(run_name + "/" + "Condor Memory Utilization Scheme")
        pos = nx.graphviz_layout(self.G)
        nx.draw_networkx(self.G, pos, k=10, with_labels=False, node_size=1500, node_color=self.memory_colors())
        nx.draw_networkx_labels(self.G, pos, labels=self.memory_labels())
        plt.savefig(run_name + "/" + "CondorMemoryUtilization.png")
        plt.clf()

        # Graph two plots the CPU utilization of various steps
        plt.figure(2, figsize=(30, 30))
        plt.axis("off")
        plt.title("Condor CPU Utilization Scheme")
        pos = nx.graphviz_layout(self.G)
        nx.draw_networkx(self.G, pos, k=10, with_labels=False, node_size=1500, node_color=self.cpu_colors())
        nx.draw_networkx_labels(self.G, pos, self.cpu_labels())
        plt.savefig(run_name + "/" + "CondorCpuUtilization.png")
        plt.clf()

        # Graph three plots the GPU utilization of various steps
        plt.figure(3, figsize=(30, 30))
        plt.axis("off")
        plt.title("Condor GPU Utilization Scheme")
        pos = nx.graphviz_layout(self.G)
        nx.draw_networkx(self.G, pos, k=10, with_labels=False, node_size=1500, node_color=self.gpu_colors())
        nx.draw_networkx_labels(self.G, pos, self.gpu_labels())
        plt.savefig(run_name + "/" + "CondorGpuUtilization.png")
        plt.clf()

    def memory_labels(self):
        labels = {}
        for node in self.G.nodes():
            labels[node] = str(node.number) + "\n" + self.G.node[node]['mem']
        return labels

    def memory_colors(self):
        memory_color = []
        for node in self.G.nodes():
            if self.G.node[node]['mem'][-1] == "G":
                memory = 1024 * int(self.G.node[node]['mem'][:-1])
            else:
                # Assumes we will never use K and only have M (Mega) other than G
                memory = int(self.G.node[node]['mem'][:-1])
            # Might want to do something more clever here to get a better gradient
            if memory < 10:
                memory_color.append("#0000FF")
            elif memory < 50:
                memory_color.append("#1919FF")
            elif memory < 100:
                memory_color.append("#3333FF")
            elif memory < 500:
                memory_color.append("#4D4DFF")
            elif memory < 1000:
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
        labels = {}
        for node in self.G.nodes():
            labels[node] = str(node.number) + "\n" + self.G.node[node]['cpus']
        return labels

    def cpu_colors(self):
        cpu_color = []
        for node in self.G.nodes():
            cpu = int(self.G.node[node]['cpus'])
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
        labels = {}
        for node in self.G.nodes():
            labels[node] = str(node.number) + "\n" + self.G.node[node]['gpus']
        return labels

    def gpu_colors(self):
        gpu_color = []
        for node in self.G.nodes():
            gpu = int(self.G.node[node]["gpus"])
            # Tuple of RGB
            if gpu < 1:
                gpu_color.append("#0000FF")
            elif gpu == 1:
                gpu_color.append("#FF9999")
            elif gpu == 2:
                gpu_color.append("#FF4D4D")
            elif gpu == 3:
                gpu_color.append("#FF6666")
            else:
                gpu_color.append("#FF0000")
        return gpu_color

# An individual unit of data for the graph.
class Node:
    node_count = 0
    # Keeps all the nodes in a list for easy access

    def __init__(self, command):
        # If this is the first node added and the only node class currently instantiated
        # Then it is the root node
        # Each node keeps track of a command
        self.command = command
        # The root node will be 0, and so on as we make more nodes.
        # This is not the actual count of how many nodes currently exist, but
        # the count of how many times we have created a node.  If a node is deleted this will just continue
        # to count up and not go back to fill that number
        self.number , Node.node_count = Node.node_count , Node.node_count + 1


# This does the parsing and processing of the incoming commands.
# All commands are initially loaded onto the graph, so this just uses the graph methods to correctly
# Split, merge, and add nodes and edges.
class CommandProcessor:
    group_keywords = [
        # The space is important for keeping them as just the linux commands
        'date ',
        'mkdir ',
        'cd ',
        'mv ',
        'cp ',
        'date',
        'rm ',
        'echo ',
        'wait',
    ]

    def __init__(self, commands):
        # The list of all commands broken into the order they were sent
        self.commands = commands
        self.graph = Graph()

    def create_graph(self):
        # A new graph structure
        nodes = list()
        for c in range(0, len(self.commands)):
            current_node = Node(self.commands[c])
            nodes.append(current_node)
            if c == 0:
                self.graph.add_node(current_node)
            else:
                self.graph.add_node(current_node, nodes[c - 1])
        # At each of these steps the number of nodes could change so we'll
        # Keep iterating over them individually for each of these commands
        for node in self.graph.nodes():
            self.split_background_proc(node)
        # This splits parallel processes within parallel process away into their own command set
        for node in self.graph.nodes():
            self.split_parallel_proc_away(node)
        for node in self.graph.nodes():
            self.split_parallel_proc(node)
        for node in self.graph.nodes():
            self.split_linked(node)
        for n in range(len(self.graph.nodes())-1,-1,-1):
            # Removes white space
            self.graph.nodes()[n].command = self.graph.nodes()[n].command.strip()
            # Get rid of any nodes that have been filled with whitespace
            current_command = self.graph.nodes()[n].command
            if current_command is "":
                self.graph.remove_node(self.graph.nodes()[n])
        for node in nx.topological_sort(self.graph.G):
            try:
                self.merge_trivial_commands(node)
            except nx.NetworkXError:
                # The node could have already been merged in which case this will just skip it
                pass
        for node in self.graph.nodes():
            mem, cpus, gpus = self.assign_resources(node.command)
            self.graph.G.node[node]['mem'] = mem
            self.graph.G.node[node]['cpus'] = cpus
            self.graph.G.node[node]['gpus'] = gpus
        # Renumber the nodes
        for i, node in enumerate(nx.topological_sort(self.graph.G),start=1):
            node.number = i
        return self.graph


    def split_parallel_proc_away(self, node):
        parallel = re.compile('(\(.*?(?<!&)&(?!&).*?\))')
        if re.search(parallel, node.command):
            # Remove parens
            commands = node.command
            commands = re.split("\(|\)", commands)
            for x in range(0, len(commands)):
                commands[x] = commands[x].strip()
            # Removes a junk split
            commands = [command for command in commands if command is not ""]
            new_nodes = [Node(command) for command in commands]
            self.graph.split_node_sequential(node, new_nodes)

    def split_parallel_proc(self, node):
        ampersand = re.compile('(?<!&)&(?!&)')
        if re.search(ampersand, node.command):
            command = node.command
            commands = re.split(ampersand, command)
            for x in range(0, len(commands)):
                commands[x] = commands[x].strip()
            commands = filter(bool, commands)
            new_nodes = [Node(command) for command in commands]
            self.graph.split_node_parallel(node, new_nodes)

    def split_background_proc(self, node):
        def filter_split(split_commands):
            background = re.compile('(?<!&)&(?!&)')
            # This undoes the split that was done by the regex split_background
            # Thus, we effectively have only split by single &'s that are outside of any parenthesis
            for c in range(len(split_commands) - 1, -1, -1):
                if split_commands[c] is not None:
                    if re.search(background, split_commands[c]):
                        if c < len(split_commands) - 1:
                            if split_commands[c + 1] is not None:
                                split_commands[c] = split_commands[c] + split_commands[c + 1]
                        if c > 0:
                            if split_commands[c - 1] is not None:
                                split_commands[c] = split_commands[c - 1] + split_commands[c]
            # Remove some junk characters to make figuring out if we have an end command easier
            split_commands = [command for command in split_commands if command is not "" or " "]
            # An end command is if there is a trailing, nonparallel job attached to some parallel jobs
            # All of the & symbols not within parenthesis leave the value "None" in their place, so we can use
            # Those as markers prior to removing them to figure out if there was an & at the end or not.
            end_command = ""
            for c in range(len(split_commands) - 1, -1, -1):
                if split_commands[c] is not None:
                    end_command = split_commands[c] + end_command
                    split_commands.pop(c)
                else:
                    break
            # If the last command is not None, that means we should do the last stuff sequentially AFTER the parallel stuff
            # Thus, we'll eat up the end until we find a None and place that as a node at the end
            # Here we remove the item that was after the split which is now a duplicate
            for c in range(len(split_commands) - 1, -1, -1):
                if split_commands[c] is not None:
                    if re.search(background, split_commands[c]):
                        if c < len(split_commands) - 1:
                            split_commands.pop(c + 1)
            # Here we remove the item that was before the split which is not also a duplicate
            # We need to move forward in this case otherwise we get stuck eating the entire sequence
            # But that also means we'll get an index error if we eat anything, so that's why we catch it here.
            for c in range(0, len(split_commands)):
                try:
                    if split_commands[c] is not None:
                        if re.search(background, split_commands[c]):
                            if c > 0:
                                split_commands.pop(c - 1)
                except IndexError:
                    pass
            # Cleans up a bunch of stuff that might be leftover
            cleaned_commands = filter(bool, split_commands)
            cleaned_commands = [command for command in cleaned_commands if command is not " "]
            return cleaned_commands, end_command

        # This picks up both the single &, but notices if they are within a parenthesis or not.
        # If the & is not within the parenthesis, it will be replaced with None when split, while the
        # in the parenthesis nothing will happen.
        split_background = re.compile('(\(.*?(?<!&)&(?!&).*?\))|(?<!&)&(?!&)')
        split_commands = re.split(split_background, node.command)
        if len(split_commands) > 1:
            split_commands, end_command = filter_split(split_commands)
            new_nodes = [Node(command) for command in split_commands]
            self.graph.split_node_parallel(node, new_nodes)
            if end_command is not "":
                self.graph.insert_node(new_nodes, Node(end_command))

    def merge_trivial_commands(self, node):
        # Breadth first search to find as many mergers as possible while best retaining order
        def bfs(node):
            queue , merge_nodes, new_command = [node] , [node] , node.command
            while queue:
                # Remove the first node and look at that
                current_node = queue.pop(0)
                # Now we look for the keywords
                # If we find it, add the current node to the nodes that should be merged together
                # Also add all of that node's children to the queue so that we also look at them
                # That way we'll merge the largest number of nodes at each merger
                # By isolating the graph to only allow one parent, we remove a ton of issues where we can
                # Introduce circular dependencies (Because then one parent becomes the parent of the other if we
                # had more than 1 parent)
                if len(self.graph.get_parents(current_node)) <= 1:
                    for word in CommandProcessor.group_keywords:
                        if word in current_node.command:
                            merge_nodes.append(current_node)
                            # We'll do this here to keep our graph class more general
                            new_command = new_command + " && " + current_node.command
                            queue = queue + self.graph.get_children(current_node)
                            break
            return merge_nodes , new_command
        for word in CommandProcessor.group_keywords:
            if word in node.command:
                merge_nodes , new_command = bfs(node)
                if new_command != node.command:
                    self.graph.merge_nodes(merge_nodes,new_command)
                break


    def split_linked(self, node):
        split_linked = re.compile("&&|;")
        split_commands = re.split(split_linked, node.command)
        # If something was actually split
        if len(split_commands) > 1:
            new_nodes = [Node(command) for command in split_commands]
            self.graph.split_node_sequential(node, new_nodes)

    def assign_resources(self, command):
        def found_match(pattern, command):
            x = re.search(pattern, command)
            return x != None
        # These defaults were initially used until I got a better handle on the memory requirements of each of the
        # individual protocols

        # Scarce resource is 256 mb and 1 cpu and 0 gpus
        scarce_resource = [
            "256",  # Memory
            "1",  # CPUS
            "0"  # GPUS
            # return scarce_resource[0],scarce_resource[1],scarce_resource[2]
        ]
        # Low resource is 1GB and 1 cpu and 0 gpus
        low_resource = [
            "1G",  # Memory
            "1",  # CPUS
            "0"  # GPUS
            # return low_resource[0],low_resource[1],low_resource[2]
        ]
        # Medium resource is 4GB and 2 cpus and 0 gpus
        medium_resource = [
            "4G",
            "4",
            "0"
            # return medium_resource[0],medium_resource[1],medium_resource[2]
        ]
        # High resource is 8GB and 8 cpu and 0 gpus
        high_resource = [
            "8G",
            "8",
            "0"
            # return high_resource[0],high_resource[1],high_resource[2]
        ]
        # Likely will never be used, but it is just another case
        extreme_resource = [
            "16G",
            "16"
            "0"
            # return extreme_resource[0],extreme_resource[1],extreme_resource[2]
        ]
        # High Ram uses 32 GB of Ram and 6 cores
        # Primarily planned for use with Samtools sort
        high_ram = [
            "32G",
            "6",
            "0"
        ]
        gpu_needed = [
            "8G",
            "6",  # As long as we only use CUSHAW-GPU keep this here (Tested as the fastest GPU-CPU combo)
            "1"
            # return gpu_needed[0],gpu_needed[1],gpu_needed[2]
        ]

        # Default
        default = medium_resource[0], medium_resource[1], medium_resource[2]
        trim = re.compile("trimmomatic(.*?).jar", re.IGNORECASE)
        if found_match(trim, command):
            return "12G", "20", "0"
        # The ^(.*?) ____ (?=[ ]) regex is there so that we only match the first word
        # As we sometimes use the aligners/counters as run names and that can
        # cause a lot of things to ask for a lot of resources
        cushaw_gpu = re.compile("^(.*?)cushaw2-gpu(?=[ ])", re.IGNORECASE)
        if found_match(cushaw_gpu, command):
            return gpu_needed[0], gpu_needed[1], gpu_needed[2]
        # This is often the slowest step, so lets give it a lot of RAM to work with
        samtools_sort = re.compile("samtools sort", re.IGNORECASE)
        if found_match(samtools_sort, command):
            return "18G", "6", "0"
        rsem_expression = re.compile("^(.*?)rsem-calculate-expression(?=[ ])", re.IGNORECASE)
        # RSEM uses sort and also takes a long time.
        # I've seen it range from 3GB up to 17GB on pretty much the same files.
        # We'll set it to 32 to be pretty safe when given even larger files.
        if found_match(rsem_expression, command):
            return "18G", "8", "0"
        # We can give low memory to other RSEM tools that we use.
        rsem = re.compile("^(.*?)rsem(?=[ ])", re.IGNORECASE)
        # RSEM uses sort and also takes a long time.
        if found_match(rsem, command):
            return "1G", "1", "0"
        # High Resource
        cushaw = re.compile("^(.*?)cushaw(?=[ ])", re.IGNORECASE)
        if found_match(cushaw, command):
            return "400M", "8", "0"
        bwa = re.compile("^(.*?)bwa(?=[ ])", re.IGNORECASE)
        if found_match(bwa, command):
            return "1G", "8", "0"
        bowtie = re.compile("^(.*?)bowtie(?=[ ])", re.IGNORECASE)
        if found_match(bowtie, command):
            return "1G", "4", "0"
        tophat = re.compile("^(.*?)tophat(?=[ ])", re.IGNORECASE)
        if found_match(tophat, command):
            return "4G", "8", "0"
        cufflinks = re.compile("^(.*?)cufflinks(?=[ ])", re.IGNORECASE)
        if found_match(cufflinks, command):
            return high_resource[0], high_resource[1], high_resource[2]
        # Medium Resources
        fastqc = re.compile("^(.*?)fastqc(?=[ ])", re.IGNORECASE)
        if found_match(fastqc, command):
            return "200M", "1", "0"
        bam2wig = re.compile("bam2wig(?=[ ])", re.IGNORECASE)
        if found_match(bam2wig, command):
            return medium_resource[0], medium_resource[1], medium_resource[2]
        picardTools = re.compile("picard-tools(?=[ ])", re.IGNORECASE)
        if found_match(picardTools, command):
            return medium_resource[0], medium_resource[1], medium_resource[2]
        rockhopper = re.compile("rockhopperWrapper.py", re.IGNORECASE)
        if found_match(rockhopper, command):
            return "4G", "1", "0"
        # bifxapps isn't an app, but is our app folder.
        # This should catch any absolute path someone introduces to try and assign it
        #  to a "medium resources" setting (Assuming the app is in bifxapps)
        bifxapps = re.compile("^(.*?)bifxapps(?=[ ])", re.IGNORECASE)
        if found_match(bifxapps, command):
            return medium_resource[0], medium_resource[1], medium_resource[2]
        # Low Resource
        samtools = re.compile("^(.*?)samtools(?=[ ])", re.IGNORECASE)
        if found_match(samtools, command):
            return "100M", "1", "0"

        rockhopper_converter = re.compile("rhFnaConverter.py",re.IGNORECASE)
        if found_match(rockhopper_converter,command):
            return "50M","1","0"
        # Scarce Resource
        htseq = re.compile("HTSeq.scripts.count(?=[ ])", re.IGNORECASE)
        if found_match(htseq, command):
            return "75M", "1", "0"
        featureCounts = re.compile("featurecounts.R(?=[ ])", re.IGNORECASE)
        if found_match(featureCounts, command):
            return scarce_resource[0], scarce_resource[1], scarce_resource[2]
        mkcopymv = re.compile("mkdir (.*?)|cp (.*?)|mv (.*?)", re.IGNORECASE)
        if found_match(mkcopymv, command):
            return "5M", "1", "0"
        gunzip = re.compile("gunzip (.*?)", re.IGNORECASE)
        if found_match(gunzip, command):
            return "100M", "1", "0"
        bash = re.compile("bash (.*?)", re.IGNORECASE)
        if found_match(bash, command):
            return "100M", "1", "0"
        # If nothing else has returned by now, return whatever was set as the defautl value.
        return default
