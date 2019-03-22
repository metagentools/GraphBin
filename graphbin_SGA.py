#!/usr/bin/python

"""graphbin_SGA.py: Improved binning of metagenomic contigs using SGA assembly graphs.

GraphBin is a metagenomic contig binning tool that makes use of the contig 
connectivity information from the assembly graph to bin contigs. It utilizes 
the binning result of an existing binning tool and a label propagation algorithm 
to correct mis-binned contigs and predict the labels of contigs which are 
discarded due to short length.

graphbin_SGA.py makes use of the assembly graphs produced by SGA (String Graph Assembler).
"""

import sys, getopt
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
#                            --n_contigs INT
#                            --binned /path/to/binning_result.csv
#                            --output /path/to/output_folder
# -------------------------------------------------------------------

start_time = time.time()

# Setup argument parser
#-----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--graph", required=True, help="path to the assembly graph file")
ap.add_argument("--n_contigs", required=True, help="number of contigs")
ap.add_argument("--binned", required=True, help="path to the .csv file with the initial binning output from an existing tool")
ap.add_argument("--output", required=True, help="path to the output folder")

args = vars(ap.parse_args())

assembly_graph_file = args["graph"]
n_contigs = int(args["n_contigs"])
n_bins = 0
contig_bins_file = args["binned"]
output_path = args["output"]

print("\nGraphBin started\n----------------")

print("Assembly graph file:", assembly_graph_file)
print("Total number of contigs available:", n_contigs)
print("Existing binning output file:", contig_bins_file)
print("Final binning output file:", output_path)


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


# Get the links from the .asqg file
#-----------------------------------

links = []

try:
    # Get contig paths from contigs.paths
    with open(assembly_graph_file) as file:
        line = file.readline()
        
        while line != "":
            
            # Identify lines with link information
            if "ED" in line:
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

bins = [[] for x in range(n_bins)]

try:
    with open(contig_bins_file) as contig_bins:
        readCSV = csv.reader(contig_bins, delimiter=',')
        for row in readCSV:
            bin_num = int(row[1])-1
            contig_num = int(row[0])
            # print(contig_num,bin_num)
            bins[bin_num].append(contig_num)

    print("\nInitial Binning result\n----------------")

    for i in range(n_bins):
        bins[i].sort()
        print("Bin", i+1, "-", len(bins[i]), ":\n", bins[i])
except:
    print("\nPlease make sure that the correct path to the binning result file is provided and it is having the correct format")
    print("Exiting GraphBin...\n")
    sys.exit(2)


# Remove labels of ambiguous vertices
#-------------------------------------

remove_labels = []

for b in range(n_bins):

    for i in bins[b]:

        my_bin = b

        dist = {}

        closest_neighbours = assembly_graph.neighbors(i, mode=ALL)

        # Determine whether all the closest neighboring vertices have the same label as its own
        neighbours_have_same_label = True
        
        for neighbour in closest_neighbours:
            for k in range(n_bins):
                if neighbour in bins[k]:
                    if k != my_bin:
                        neighbours_have_same_label = False
                        break
                        
        if not neighbours_have_same_label:
            remove_labels.append(i)

# Further remove labels of ambiguous vertices
for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

for b in range(n_bins):

    for i in bins[b]:

        my_bin = b

        dist = {}

        # Get the distant to all the vertices
        for j in range(node_count):
            dis = assembly_graph.shortest_paths_dijkstra(source=i, target=j, weights=None, mode=OUT)[0][0]
            if dis != 0:
                dist[j] = dis
        
        # Sort the vertices in the ascending order of their distance
        sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))

        closest_neighbours = []

        distances = [sys.maxsize for x in range(n_bins)]

        for element in sorted_dist:

            count_is_million = True

            for k in range(n_bins):
                if distances[k] == sys.maxsize:
                    count_is_million = False

            if not count_is_million:

                for h in range(n_bins):
                    if element[0] in bins[h] and distances[h] == sys.maxsize:
                        distances[h] = element[1]

        min_dist = sys.maxsize
        min_index = sys.maxsize

        for j in range(n_bins):

            if distances[j] < min_dist:
                min_dist = distances[j]
                min_index = j
        
        # Get the closest vertices
        for element in sorted_dist:
            if element[1] == min_dist:
                closest_neighbours.append(element[0])

        # Determine whether all the closest vertices have the same label as its own
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
print("\nRemove labels of contigs:", remove_labels)

# Remove labels of ambiguous vertices
for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

print("\nRefined Binning result\n----------------------")

for i in range(n_bins):
    print("Bin", i+1, "-", len(bins[i]), ":\n", bins[i])
        

# Run label propagation
#-----------------------

data = []

for contig in range(node_count):
    
    neighbours = assembly_graph.neighbors(contig, mode=ALL)
    
    # Consider vertices that are not isolated

    if len(neighbours) > 0:
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

print("\nStarting label propagation\n---------------------------")

ans = lp.run(0.00001, 100, show_log=True, clean_result=False) 
ans.sort()

for l in ans:
    for i in range(n_bins):
        if l[1]==i+1 and l[0] not in bins[i]:
            bins[i].append(l[0])

print("\nLabel Propagation result\n-------------------------")

for i in range(n_bins):
    print("Bin", i+1, "-", len(bins[i]), ":\n", bins[i])

elapsed_time = time.time() - start_time


# Remove labels of ambiguous vertices
#-------------------------------------

remove_labels = []

for b in range(n_bins):

    for i in bins[b]:

        my_bin = b

        closest_neighbours = assembly_graph.neighbors(i, mode=ALL)

        # Determine whether all the closest binned neighboring vertices have the same label as its own
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
print("\nRemove labels of contigs:", remove_labels)

# Remove labels of ambiguous vertices
for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)


# Print the final result
#------------------------

print("\nFinal Refined Binning result\n----------------------")

for i in range(n_bins):
    print("Bin", i+1, "-", len(bins[i]), ":\n", bins[i])

# Print elapsed time for the process
print("\nElapsed time: ", elapsed_time, " seconds")


# Write result to output file
#-----------------------------

output_bins = []

for i in range(node_count):
    for k in range(n_bins):
        if i in bins[k]:
            line = []
            line.append("NODE_"+str(i+1))
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