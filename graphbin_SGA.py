#!/usr/bin/python

"""graphbin_SGA.py: Improved binning of metagenomic contigs using SGA assembly graphs.

GraphBin is a metagenomic contig binning tool that makes use of the contig 
connectivity information from the assembly graph to bin contigs. It utilizes 
the binning result of an existing binning tool and a label propagation algorithm 
to correct mis-binned contigs and predict the labels of contigs which are 
discarded due to short length.

graphbin_SGA.py makes use of the assembly graphs produced by SGA (String Graph Assembler).
"""

import sys
import os
import subprocess
import csv
import operator
import time
import argparse

from igraph import *
from labelpropagation.labelprop import LabelProp

__author__ = "Vijini Mallawaarachchi, Anuradha Wickramarachchi, and Yu Lin"
__copyright__ = "Copyright 2019, GraphBin Project"
__credits__ = ["Benjamin Kaehler", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Prototype"


# Sample command
# -------------------------------------------------------------------
# python graphbin_SGA.py     --graph /path/to/graph_file.asqg
#                            --binned /path/to/binning_result.csv
#                            --output /path/to/output_folder
# -------------------------------------------------------------------

start_time = time.time()

# Setup argument parser
#-----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--graph", required=True, help="path to the assembly graph file")
ap.add_argument("--binned", required=True, help="path to the .csv file with the initial binning output from an existing tool")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--max_iteration", required=False, nargs='?', type=int, help="maximum number of iterations for label propagation algorithm. [default: 100]")
ap.add_argument("--diff_threshold", required=False, nargs='?', type=float, help="difference threshold for label propagation algorithm. [default: 0.1]")

args = vars(ap.parse_args())

assembly_graph_file = args["graph"]
n_contigs = 0
n_bins = 0
contig_bins_file = args["binned"]
output_path = args["output"]

max_iteration = 100
diff_threshold = 0.1


print("\nWelcome to GraphBin: Improved Binning of Metagenomic Contigs using Assembly Graphs")
print("This version of GraphBin makes use of the assembly graph produced by SGA which is based on the OLC (more recent string graph) approach.\n")

print("GraphBin started\n-----------------")


# Validate max_iteration and diff_threshold
#---------------------------------------------------
try:

    if args["max_iteration"] is not None:
        max_iteration = int(args["max_iteration"])

except:
    print("\nPlease enter a valid number for max_iterations")
    print("Exiting GraphBin...\n")
    sys.exit(2)

try:

    if args["diff_threshold"] is not None:
        diff_threshold = float(args["diff_threshold"])

except:
    print("\nPlease enter a valid number for diff_threshold")
    print("Exiting GraphBin...\n")
    sys.exit(2)


# Check if output folder exists
#---------------------------------------------------

# Handle for missing trailing forwardslash in output folder path
if output_path[-1:] != "/":
    output_path = output_path + "/"

# Create output folder if it does not exist
if not os.path.isdir(output_path):
    subprocess.run("mkdir -p "+output_path, shell=True)


print("Assembly graph file:", assembly_graph_file)
print("Existing binning output file:", contig_bins_file)
print("Final binning output file:", output_path)
print("Maximum number of iterations:", max_iteration)
print("Difference threshold:", diff_threshold)


# Get the number of bins from the initial binning result
#--------------------------------------------------------

try:
    all_bins_list = []

    with open(contig_bins_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            all_bins_list.append(row[1])
            
    bins_list = list(set(all_bins_list))
    bins_list.sort()

    n_bins = len(bins_list)
    print("Number of bins available in binning result:", n_bins)
except:
    print("\nPlease make sure that the correct path to the binning result file is provided and it is having the correct format")
    print("Exiting GraphBin...\n")
    sys.exit(2)


print("\nConstructing the assembly graph...")

# Get the links from the .asqg file
#-----------------------------------

links = []

try:
    # Get contig connections from .asqg file
    with open(assembly_graph_file) as file:
        line = file.readline()
        
        while line != "":

            # Count the number of contigs
            if "VT" in line:
                n_contigs += 1
            
            # Identify lines with link information
            elif "ED" in line:
                link = []
                strings = line.split("\t")[1].split()
                link.append(int(strings[0][7:]))
                link.append(int(strings[1][7:]))
                links.append(link)
            line = file.readline()

except:
    print("\nPlease make sure that the correct path to the assembly graph file is provided")
    print("Exiting GraphBin...\n")
    sys.exit(2)

node_count = n_contigs


## Construct the assembly graph
#-------------------------------

try:

    # Create the graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)

    for i in range(len(assembly_graph.vs)):
        assembly_graph.vs[i]["id"]= i
        assembly_graph.vs[i]["label"]= str(i)

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            assembly_graph.add_edge(link[0], link[1])

    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

except:
    print("\nPlease make sure that the correct path to the assembly graph file is provided")
    print("Exiting GraphBin...\n")
    sys.exit(2)


# Get initial binning result
#----------------------------

print("\nObtaining the initial binning result...")

bins = [[] for x in range(n_bins)]

try:
    with open(contig_bins_file) as contig_bins:
        readCSV = csv.reader(contig_bins, delimiter=',')
        for row in readCSV:
            bin_num = int(row[1])-1
            contig_num = int(row[0])
            # print(contig_num,bin_num)
            bins[bin_num].append(contig_num)

    print("\nInitial Binning result\n-----------------------")

    for i in range(n_bins):
        bins[i].sort()

except:
    print("\nPlease make sure that the correct path to the binning result file is provided and it is having the correct format")
    print("Exiting GraphBin...\n")
    sys.exit(2)


# Remove labels of ambiguous vertices
#-------------------------------------

def getClosestLabelledVertices(graph, node, binned_contigs):
    queu_l = [graph.neighbors(node, mode='ALL')]
    visited_l = [node]
    labelled = []

    while len(queu_l) > 0:
        active_level = queu_l.pop(0)
        is_finish = False
        visited_l += active_level

        for n in active_level:
            if n in binned_contigs:
                is_finish = True
                labelled.append(n)
        if is_finish:
            return labelled
        else:
            temp = []
            for n in active_level:
                temp += graph.neighbors(n, mode='ALL')
            temp = list(set(temp))
            temp2 = []

            for n in temp:
                if n not in visited_l:
                    temp2.append(n)
            if len(temp2) > 0:
                queu_l.append(temp2)
    return labelled


print("\nDetermining ambiguous vertices...")

remove_labels = []

neighbours_have_same_label_list = []

for b in range(n_bins):

    for i in bins[b]:

        my_bin = b

        dist = {}

        # Get set of closest labelled vertices with distance = 1
        closest_neighbours = assembly_graph.neighbors(i, mode=ALL)

        # Determine whether all the closest labelled vertices have the same label as its own
        neighbours_have_same_label = True
        
        neighbours_binned = False
        
        for neighbour in closest_neighbours:
            for k in range(n_bins):
                if neighbour in bins[k]:
                    neighbours_binned = True
                    if k != my_bin:
                        neighbours_have_same_label = False
                        break
                        
        if not neighbours_have_same_label:
            remove_labels.append(i)
        elif neighbours_binned:
            neighbours_have_same_label_list.append(i)

for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

# Further remove labels of ambiguous vertices
binned_contigs = []

for n in range(n_bins):
    binned_contigs = sorted(binned_contigs+bins[n])
            
for b in range(n_bins):

    for i in bins[b]:
        
        if i not in neighbours_have_same_label_list:

            my_bin = b

            # Get set of closest labelled vertices
            closest_neighbours = getClosestLabelledVertices(assembly_graph, i, binned_contigs)

            if len(closest_neighbours) > 0:

                # Determine whether all the closest labelled vertices have the same label as its own
                neighbours_have_same_label = True

                for neighbour in closest_neighbours:
                    for k in range(n_bins):
                        if neighbour in bins[k]:
                            if k != my_bin:
                                neighbours_have_same_label = False
                                break

                if not neighbours_have_same_label and i not in remove_labels:
                    remove_labels.append(i)

remove_labels.sort()

print("Removing labels of ambiguous vertices...")

# Remove labels of ambiguous vertices
for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

print("\nObtaining Refined Binning result...")
        

# Get vertices which are not isolated and not in components without any labels
#-----------------------------------------------------------------------------

print("\nDeteremining vertices which are not isolated and not in components without any labels...")

non_isolated = []

for i in range(node_count):
    
    if i not in non_isolated and i in binned_contigs:

        component = []
        component.append(i)
        length = len(component)
        neighbours = assembly_graph.neighbors(i, mode=ALL)

        for neighbor in neighbours:
            if neighbor not in component:
                component.append(neighbor)

        component = list(set(component))

        while length!= len(component):

            length = len(component)

            for j in component:

                neighbours = assembly_graph.neighbors(j, mode=ALL)

                for neighbor in neighbours:
                    if neighbor not in component:
                        component.append(neighbor)

        labelled = False
        for j in component:
            if j in binned_contigs:
                labelled = True
                break

        if labelled:
            for j in component:
                if j not in non_isolated:
                    non_isolated.append(j)

print("Number of non-isolated contigs:", len(non_isolated))


# Run label propagation
#-----------------------

data = []

for contig in range(node_count):
    
    # Consider vertices that are not isolated

    if contig in non_isolated:
        line = []
        line.append(contig)

        assigned = False

        for i in range(n_bins):
            if contig in bins[i]:
                line.append(i+1)
                assigned = True
        
        if not assigned:
            line.append(0)

        neighbours = assembly_graph.neighbors(contig, mode=ALL)

        neighs = []

        for neighbour in neighbours:
            n = []
            n.append(neighbour)
            n.append(1.0)
            neighs.append(n)

        line.append(neighs)

        data.append(line)

lp = LabelProp()

lp.load_data_from_mem(data)

print("\nStarting label propagation with eps="+str(diff_threshold)+" and max_iteration="+str(max_iteration))

ans = lp.run(diff_threshold, max_iteration, show_log=True, clean_result=False) 
ans.sort()

print("\nObtaining Label Propagation result...")

for l in ans:
    for i in range(n_bins):
        if l[1]==i+1 and l[0] not in bins[i]:
            bins[i].append(l[0])


elapsed_time = time.time() - start_time


# Remove labels of ambiguous vertices
#-------------------------------------

print("\nDetermining ambiguous vertices...")

remove_labels = []

for b in range(n_bins):

    for i in bins[b]:

        my_bin = b

        closest_neighbours = assembly_graph.neighbors(i, mode=ALL)

        # Determine whether all the closest labelled vertices have the same label as its own
        neighbours_have_same_label = True
        
        for neighbour in closest_neighbours:
            for k in range(n_bins):
                if neighbour in bins[k]:
                    if k != my_bin:
                        neighbours_have_same_label = False
                        break
                        
        if not neighbours_have_same_label:
            remove_labels.append(i)

remove_labels.sort()

print("Removing labels of ambiguous vertices...")

# Remove labels of ambiguous vertices
for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

print("\nObtaining the Final Refined Binning result...")

for i in range(n_bins):
    bins[i].sort()

# Print elapsed time for the process
print("\nElapsed time: ", elapsed_time, " seconds")


# Write result to output file
#-----------------------------

print("\nWriting the Final Binning result to file...")

output_bins = []

for i in range(node_count):
    for k in range(n_bins):
        if i in bins[k]:
            line = []
            line.append("contig-"+str(i+1))
            line.append(k+1)
            output_bins.append(line)

output_file = output_path + 'graphbin_output.csv'

with open(output_file, mode='w') as output_file:
    output_writer = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in output_bins:
        output_writer.writerow(row)

print("\nFinal binning results can be found at", output_file.name)


# Exit program
#--------------

print("\nThank you for using GraphBin!\n")