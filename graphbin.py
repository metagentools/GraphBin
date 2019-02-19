#!/usr/bin/python

"""graphbin.py: Improved binning of metagenomic contigs using assembly graphs."""

import sys, getopt
import csv
import operator
import time
import argparse

from igraph import *
from labelprop import LabelProp

__author__ = "Vijini Mallawaarachchi, Anuradha Wickramarachchi, and Yu Lin"
__copyright__ = "Copyright 2019, GraphBin Project"
__credits__ = ["Benjamin Kaehler", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Prototype"


# Sample command
# ------------------------------------------------------------------
# python main.py    --graph /path/to/graph_file.gfa
#                   --contigs /path/to/contigs_file.fasta
#                   --paths /path/to/paths_file.paths
#                   --binned /path/to/binning_result.csv
#                   --output /path/to/output_folder
# ------------------------------------------------------------------

start_time = time.time()

# Setup argument parser
#-----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--graph", required=True, help="path to the assembly graph file")
ap.add_argument("--contigs", required=True, help="path to the contigs.fasta file")
ap.add_argument("--paths", required=True, help="path to the contigs.paths file")
ap.add_argument("--binned", required=True, help="path to the .csv file with the initial binning output from an existing tool")
ap.add_argument("--output", required=True, help="path to the output folder")

args = vars(ap.parse_args())

assembly_graph_file = args["graph"]
contig_file = args["contigs"]
contig_paths = args["paths"]
n_bins = 0
contig_bins_file = args["binned"]
output_path = args["output"]

print("\nGraphBin started\n----------------")

print("Assembly graph file:", assembly_graph_file)
print("Contigs file:", contig_file)
print("Contig paths file:", contig_paths)
print("Existing binning output file:", contig_bins_file)
print("Final binning output file:", output_path)


# Get the number of bins from the initial binning result
#---------------------------------------------------

try:
    all_bins_list = []

    with open(contig_bins_file) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            all_bins_list.append(row[1])
            
    bins_list = list(set(all_bins_list))
    bins_list.sort()

    n_bins = len(bins_list)
except:
    print("\nPlease make sure that the correct path to the binning result file is provided and it is having the correct format")
    print("Exiting GraphBin...\n")
    sys.exit(2)


# Get contig paths from contigs.paths
#-------------------------------------

paths = []
links = []

try:
    with open(contig_paths) as file:
        name = file.readline()
        path = file.readline()
        
        while name != "" and path != "":
                
            while ";" in path:
                path = path[:-2]+","+file.readline()
                
            paths.append(path.split("\n")[0])
            
            name = file.readline()
            path = file.readline()
except:
    print("\nPlease make sure that the correct path to the contig paths file is provided")
    print("Exiting GraphBin...\n")
    sys.exit(2)

node_count = int(len(paths)/2)

print("\nTotal number of contigs available:", node_count)


## Construct the assembly graph
#-------------------------------

try:
    # Get links from assembly_graph_with_scaffolds.gfa
    with open(assembly_graph_file) as file:
        line = file.readline()
        
        while line != "":
            
            # Identify lines with link information
            if "L" in line:
                strings = line.split("\t")
                links.append(strings[1]+strings[2]+" "+strings[3]+strings[4])
            line = file.readline()

    # Create the graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)

    for i in range(len(assembly_graph.vs)):
        assembly_graph.vs[i]["id"]= i
        assembly_graph.vs[i]["label"]= str(i)

    # Iterate paths
    for i in range(len(paths)):
        segments = paths[i].split(",")
        start = segments[0]
        end = segments[len(segments)-1]
        
        new_links = []
        connections = []
        
        # Iterate links
        for link in links:
            link_list = link.split()
            
            if start in link_list[0]:
                new_links.append(link_list[1])
            elif start in link_list[1]:
                new_links.append(link_list[0])
            if end in link_list[0]:
                new_links.append(link_list[1])
            elif end in link_list[1]:
                new_links.append(link_list[0])
        
        # Determine connections
        for new_link in new_links:
            for j in range(len(paths)):
                if new_link in paths[j] and int(j/2) not in connections and int(j/2)!=int(i/2):
                    ind = int(j/2)
                    connections.append(ind)
        
        # Add connections in graph
        for connection in connections:
            assembly_graph.add_edge(int(i/2),connection)

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

        # Get distant to connected vertices
        for j in range(node_count):
            dis = assembly_graph.shortest_paths_dijkstra(source=i, target=j, weights=None, mode=OUT)[0][0]
            if dis != 0:
                dist[j] = dis
        
        # Sort the vertices in the ascending order of their distance
        sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))

        closest_neighbours = []
        
        # Get the closest neighboring vertices
        for element in sorted_dist:
            if element[1] == 1:
                closest_neighbours.append(element[0])

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

# Remove labels of ambiguous vertices
for i in remove_labels:
    for n in range(n_bins):
        if i in bins[n]:
            bins[n].remove(i)

for b in range(n_bins):

    for i in bins[b]:

        my_bin = b

        dist = {}

        # Get distant to connected vertices
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

        dist = {}

        # Get distant to connected vertices
        for j in range(node_count):
            dis = assembly_graph.shortest_paths_dijkstra(source=i, target=j, weights=None, mode=OUT)[0][0]
            if dis != 0:
                dist[j] = dis
        
        # Sort the vertices in the ascending order of their distance
        sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))

        closest_neighbours = []

        distances = [sys.maxsize for x in range(n_bins)]

        for element in sorted_dist:

            count_is_large = True

            for k in range(n_bins):
                if distances[k] == sys.maxsize:
                    count_is_large = False

            if not count_is_large:

                for h in range(n_bins):
                    if element[0] in bins[h] and distances[h] == sys.maxsize:
                        distances[h] = element[1]

        # Get the minimum distance to binned neighboring vertices
        min_dist = sys.maxsize
        min_index = sys.maxsize

        for j in range(n_bins):

            if distances[j] < min_dist:
                min_dist = distances[j]
                min_index = j
        
        # Get the closest binned neighboring vertices
        for element in sorted_dist:
            if element[1] == min_dist:
                closest_neighbours.append(element[0])

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