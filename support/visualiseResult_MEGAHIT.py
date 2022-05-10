#!/usr/bin/python3

"""visualiseResult_MEGAHIT.py: Visualise the binning result from on the MEGAHIT assembly graph.

Visualize the binning result by denoting coloured contigs in the assembly
graph according to their corresponding bins. You can visualise the initial
binning result obtained from an existing binning tool and the final binning
result obtained from GraphBin and compare.

"""

import sys
import os
import csv
import argparse
import re
import subprocess
import random

from igraph import *
from bidirectionalmap.bidirectionalmap import BidirectionalMap


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "GPL"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"


# Sample command
# -------------------------------------------------------------------
# python visualiseResult_MEGAHIT.py --initial /path/to/file_with_initial_binning_result
#                                   --final /path/to/file_with_final_GraphBin_binning_result
#                                   --graph /path/to/graph_file.gfa
#                                   --output /path/to/output_folder
#                                   --prefix prefix for output image files
#                                   --type type_of_the_image
#                                   --width width_of_image
#                                   --height height_of_image
#                                   --dpi dpi_value
# -------------------------------------------------------------------


# Setup argument parser
#-----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--initial", required=True, type=str, help="path to the file containing the initial binning result from an existing tool")
ap.add_argument("--final", required=True, type=str, help="path to the file containing the final GraphBin binning result")
ap.add_argument("--graph", required=True, type=str, help="path to the assembly graph file")
ap.add_argument("--output", required=True, type=str, help="path to the output folder")
ap.add_argument("--prefix", required=False, type=str, default='', help="prefix for the output image files")
ap.add_argument("--type", required=False, type=str, default='png', help="type of the image (jpg, png, eps, svg)")
ap.add_argument("--width", required=False, type=int, default=2000, help="width of the image in pixels")
ap.add_argument("--height", required=False, type=int, default=2000, help="height of the image in pixels")
ap.add_argument("--vsize", required=False, type=int, default=50, help="size of the vertices")
ap.add_argument("--lsize", required=False, type=int, default=8, help="size of the vertex labels")
ap.add_argument("--margin", required=False, type=int, default=50, help="margin of the figure")
ap.add_argument("--dpi", required=False, type=int, default=300, help="dpi value")

args = vars(ap.parse_args())

initial_binning_result = args["initial"]
final_binning_result = args["final"]
assembly_graph_file = args["graph"]
output_path = args["output"]
prefix = args["prefix"]
dpi = args["dpi"]
width = args["width"]
height = args["height"]
vsize = args["vsize"]
lsize = args["lsize"]
margin = args["margin"]
image_type = args["type"]

print("\nWelcome to binning result visualiser of GraphBin!")
print("This version of the visualiser makes use of the assembly graph produced by MEGAHIT which is based on the de Bruijn graph approach.\n")


# Validate prefix
#---------------------------------------------------
try:
    if prefix != '':
        if not prefix.endswith("_"):
            prefix = prefix+"_"
    else:
        prefix = ''

except:
    print("\nPlease enter a valid string for prefix")
    print("Exiting visualiseResult...\nBye...!\n")
    sys.exit(1)


# Format type if provided
#---------------------------------------------------
if image_type.startswith("."):
    image_type = image_type[1:]


# Check if output folder exists
#---------------------------------------------------

# Handle for missing trailing forwardslash in output folder path
if output_path[-1:] != "/":
    output_path = output_path + "/"

# Create output folder if it does not exist
if not os.path.isdir(output_path):
    subprocess.run("mkdir -p "+output_path, shell=True)

print("Assembly graph file:", assembly_graph_file)
print("Initial binning results file:", initial_binning_result)
print("Final binning results file:", final_binning_result)
print("Final output path:", output_path)
print("Image type:", image_type)
print("Width of image:", width, "pixels")
print("Height of image:", height, "pixels")
print("Size of the vertices:", vsize, "pt")
print("Size of the vertex labels:", lsize, "pt")
print("Size of the margin:", margin, "pt")


# Get the number of bins from the initial binning result
#---------------------------------------------------

try:
    all_bins_list = []
    n_bins = 0

    with open(initial_binning_result) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            all_bins_list.append(row[1])
            
    bins_list = list(set(all_bins_list))
    bins_list.sort()

    n_bins = len(bins_list)
    print("Number of bins available in initial binning result:", n_bins)
except:
    print("\nPlease make sure that the correct path to the initial binning result file is provided and it is having the correct format")
    print("Exiting visualiseResult...\nBye...!\n")
    sys.exit(1)


print("\nConstructing the assembly graph...")

## Construct the assembly graph
#-------------------------------

node_count = 0

links = []

my_map = BidirectionalMap()

try:

    # Get links from .gfa file
    with open(assembly_graph_file) as file:

        line = file.readline()

        while line != "":

            # Identify lines with link information
            if line.startswith("L"):
                link = []

                strings = line.split("\t")

                start_1 = 'NODE_'
                end_1 = '_length'

                link1 = int(re.search('%s(.*)%s' % (start_1, end_1), strings[1]).group(1))

                start_2 = 'NODE_'
                end_2 = '_length'

                link2 = int(re.search('%s(.*)%s' % (start_2, end_2), strings[3]).group(1))

                link.append(link1)
                link.append(link2)
                links.append(link)

            elif line.startswith("S"):
                strings = line.split()

                start = 'NODE_'
                end = '_length'

                contig_num = int(re.search('%s(.*)%s' % (start, end), strings[1]).group(1))

                my_map[node_count] = int(contig_num)

                node_count += 1
            
            line = file.readline()


    print("\nTotal number of contigs available:", node_count)

    contigs_map = my_map
    contigs_map_rev = my_map.inverse

    # Create graph
    assembly_graph = Graph()

    # Add vertices
    assembly_graph.add_vertices(node_count)

    # Create list of edges
    edge_list = []

    for i in range(node_count):
        assembly_graph.vs[i]["id"]= i
        assembly_graph.vs[i]["label"]= 'NODE_'+str(contigs_map[i])
        
    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contigs_map_rev[link[0]], contigs_map_rev[link[1]]))

    # Add edges to the graph 
    assembly_graph.add_edges(edge_list)
    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

except:
    print("\nPlease make sure that the correct path to the assembly graph file is provided.")
    print("Exiting visualiseResult...\nBye...!\n")
    sys.exit(1)


# Get initial binning result
#----------------------------

print("\nObtaining the initial binning result...")

bins = [[] for x in range(n_bins)]

try:
    with open(initial_binning_result) as contig_bins:
        readCSV = csv.reader(contig_bins, delimiter=',')
        for row in readCSV:
            start = 'NODE_'
            end = ''
            contig_num = contigs_map_rev[int(re.search('%s(.*)%s' % (start, end), row[0]).group(1))]
            
            bin_num = int(row[1])-1
            bins[bin_num].append(contig_num)

    for i in range(n_bins):
        bins[i].sort()

except:
    print("\nPlease make sure that the correct path to the binning result file is provided and it is having the correct format")
    print("Exiting visualiseResult...\nBye...!\n")
    sys.exit(1)


# Get list of colours according to number of bins
#-------------------------------------------------

print("\nPicking colours...")

def colours(n):
  ret = []
  r = int(random.random() * 256)
  g = int(random.random() * 256)
  b = int(random.random() * 256)
  step = 256 / n
  for i in range(n):
    r += step
    g += step
    b += step
    r = int(r) % 256
    g = int(g) % 256
    b = int(b) % 256

    ret.append('#%02x%02x%02x' % (r,g,b))

  return ret

my_colours = colours(n_bins)


# Visualise the initial assembly graph
#--------------------------------------

print("\nDrawing and saving the assembly graph with the initial binning result...")

initial_out_fig_name = output_path+prefix+"initial_binning_result."+image_type

node_colours = []

for i in range(node_count):
    no_bin = True
    for j in range(n_bins):
        if i in bins[j]:
            node_colours.append(my_colours[j])
            no_bin = False
    
    if no_bin:
        node_colours.append("grey")

assembly_graph.vs["color"] = node_colours

visual_style = {}

# Set bbox and margin
visual_style["bbox"] = (width,height)
visual_style["margin"] = margin

# Set vertex size
visual_style["vertex_size"] = vsize

# Set vertex lable size
visual_style["vertex_label_size"] = lsize

# Don't curve the edges
visual_style["edge_curved"] = False

# Set the layout
my_layout = assembly_graph.layout_fruchterman_reingold()
visual_style["layout"] = my_layout

# Plot the graph
plot(assembly_graph, initial_out_fig_name, **visual_style)


# Get the final GraphBin binning result
#---------------------------------------

print("\nObtaining the final GraphBin binning result...")

bins = [[] for x in range(n_bins)]

try:
    with open(final_binning_result) as contig_bins:
        readCSV = csv.reader(contig_bins, delimiter=',')
        for row in readCSV:
            start = 'NODE_'
            end = ''
            contig_num = contigs_map_rev[int(re.search('%s(.*)%s' % (start, end), row[0]).group(1))]
            
            bin_num = int(row[1])-1
            bins[bin_num].append(contig_num)

    for i in range(n_bins):
        bins[i].sort()

except:
    print("\nPlease make sure that the correct path to the final binning result file is provided and it is having the correct format")
    print("Exiting visualiseResult...\nBye...!\n")
    sys.exit(1)


# Visualise the final assembly graph
#------------------------------------

print("\nDrawing and saving the assembly graph with the final GraphBin binning result...")

final_out_fig_name = output_path+prefix+"final_GraphBin_binning_result."+image_type

node_colours = []

for i in range(node_count):
    no_bin = True
    for j in range(n_bins):
        if i in bins[j]:
            node_colours.append(my_colours[j])
            no_bin = False
    
    if no_bin:
        node_colours.append("grey")

assembly_graph.vs["color"] = node_colours

visual_style = {}

# Set bbox and margin
visual_style["bbox"] = (width,height)
visual_style["margin"] = margin

# Set vertex size
visual_style["vertex_size"] = vsize

# Set vertex lable size
visual_style["vertex_label_size"] = lsize

# Don't curve the edges
visual_style["edge_curved"] = False

# Set the layout
visual_style["layout"] = my_layout

# Plot the graph
plot(assembly_graph, final_out_fig_name, **visual_style)

print("\nVisualization of the initial binning results can be found at", initial_out_fig_name)
print("Visualization of the final GraphBin binning results can be found at", final_out_fig_name)


# Exit program
#--------------

print("\nThank you for using visualiseResult for GraphBin!\n")