#!/usr/bin/python

"""visualiseResult_SGA.py: Visualise the binning result from on the SGA assembly graph.

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
from base64 import b16encode

__author__ = "Vijini Mallawaarachchi, Anuradha Wickramarachchi, and Yu Lin"
__copyright__ = "Copyright 2019, GraphBin Project"
__credits__ = ["Benjamin Kaehler", "Gavin Huttley"]
__license__ = "GPL"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"


# Sample command
# -------------------------------------------------------------------
# python visualiseResult.py     --initial /path/to/file_with_initial_binning_result
#                               --final /path/to/file_with_final_GraphBin_binning_result
#                               --graph /path/to/graph_file.asqg
#                               --output /path/to/output_folder
#                               --prefix prefix for output image files
#                               --type type_of_the_image
#                               --width width_of_image
#                               --height height_of_image
#                               --dpi dpi_value
# -------------------------------------------------------------------


# Setup argument parser
#-----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--initial", required=True, help="path to the file containing the initial binning result from an existing tool")
ap.add_argument("--final", required=True, help="path to the file containing the final GraphBin binning result")
ap.add_argument("--graph", required=True, help="path to the assembly graph file")
ap.add_argument("--output", required=True, help="path to the output folder")
ap.add_argument("--prefix", required=False, nargs='?', help="prefix for the output image files")
ap.add_argument("--type", required=False, nargs='?', help="type of the image (jpg, png, eps, svg)")
ap.add_argument("--width", required=False, nargs='?', help="width of the image in pixels")
ap.add_argument("--height", required=False, nargs='?', help="height of the image in pixels")
ap.add_argument("--dpi", required=False, nargs='?', help="dpi value")

args = vars(ap.parse_args())

initial_binning_result = args["initial"]
final_binning_result = args["final"]
assembly_graph_file = args["graph"]
output_path = args["output"]
prefix = "my"
dpi = 300
width = 1000
height = 1000
image_type = "png"

print("\nWelcome to binning result visualiser of GraphBin!")
print("This version of the visualiser makes use of the assembly graph produced by SGA which is based on the OLC (more recent string graph) approach.\n")

# Obtain prefix and type if provided
#---------------------------------------------------
if args["prefix"] is not None:
    prefix = args["prefix"]

if args["type"] is not None:
    if args["type"].startswith("."):
        image_type = args["type"][1:]
    else:
        image_type = args["type"]


# Validate height, width and dpi
#---------------------------------------------------
try:

    if args["width"] is not None:
        width = int(args["width"])

except:
    print("\nPlease enter a valid number for width")
    print("Exiting visualiseResult.py...\nBye...!\n")
    sys.exit(2)

try:

    if args["height"] is not None:
        height = int(args["height"])

except:
    print("\nPlease enter a valid number for height")
    print("Exiting visualiseResult.py...\nBye...!\n")
    sys.exit(2)

try:

    if args["dpi"] is not None:
        dpi = int(args["dpi"])

except:
    print("\nPlease enter a valid dpi value")
    print("Exiting visualiseResult.py...\nBye...!\n")
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
print("Initial binning results file:", initial_binning_result)
print("Final binning results file:", final_binning_result)
print("Final output path:", output_path)
print("Image type:", image_type)
print("Width of image:", width, "pixels")
print("Height of image:", height, "pixels")


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
    print("Exiting visualiseResult.py...\nBye...!\n")
    sys.exit(2)


print("\nConstructing the assembly graph...")

# Get the links from the .asqg file
#-----------------------------------

links = []
n_contigs = 0

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
    print("Exiting visualiseResult.py...\nBye...!\n")
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
    print("Exiting visualiseResult.py...\nBye...!\n")
    sys.exit(2)


# Get initial binning result
#----------------------------

print("\nObtaining the initial binning result...")

bins = [[] for x in range(n_bins)]

try:
    with open(initial_binning_result) as contig_bins:
        readCSV = csv.reader(contig_bins, delimiter=',')
        for row in readCSV:
            bin_num = int(row[1])-1
            contig_num = int(row[0])
            bins[bin_num].append(contig_num)

    for i in range(n_bins):
        bins[i].sort()

except:
    print("\nPlease make sure that the correct path to the initial binning result file is provided and it is having the correct format")
    print("Exiting visualiseResult.py...\nBye...!\n")
    sys.exit(2)


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

initial_out_fig_name = output_path+prefix+"_initial_binning_result."+image_type

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
visual_style["margin"] = 17

# Set vertex colours
# visual_style["vertex_color"] = '#fecc5c'

# # Scale vertices based on degree
# outdegree = g.outdegree()
visual_style["vertex_size"] = 25

# Set vertex lable size
visual_style["vertex_label_size"] = 8

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
            bin_num = int(row[1])-1
            contig_num = int(row[0][7:])
            bins[bin_num].append(contig_num)

    for i in range(n_bins):
        bins[i].sort()

except:
    print("\nPlease make sure that the correct path to the final binning result file is provided and it is having the correct format")
    print("Exiting visualiseResult.py...\nBye...!\n")
    sys.exit(2)


# Visualise the final assembly graph
#------------------------------------

print("\nDrawing and saving the assembly graph with the final GraphBin binning result...")

final_out_fig_name = output_path+prefix+"_final_GraphBin_binning_result."+image_type

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
visual_style["margin"] = 17

# Set vertex colours
# visual_style["vertex_color"] = '#fecc5c'

# # Scale vertices based on degree
# outdegree = g.outdegree()
visual_style["vertex_size"] = 25

# Set vertex lable size
visual_style["vertex_label_size"] = 8

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