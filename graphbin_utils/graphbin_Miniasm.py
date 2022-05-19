#!/usr/bin/env python3

"""graphbin_Miniasm.py: Refined binning of metagenomic contigs using Miniasm assembly graphs.

GraphBin is a metagenomic contig binning tool that makes use of the contig 
connectivity information from the assembly graph to bin contigs. It utilizes 
the binning result of an existing binning tool and a label propagation algorithm 
to correct mis-binned contigs and predict the labels of contigs which are 
discarded due to short length.

graphbin_Miniasm.py makes use of the assembly graphs produced by Miniasm long read assembler.
"""

import sys
import os
import subprocess
import csv
import time
import logging

from igraph import *
from Bio import SeqIO

from graphbin_utils.labelpropagation.labelprop import LabelProp
from graphbin_utils.bidirectionalmap.bidirectionalmap import BidirectionalMap
from graphbin_utils.graphbin_Func import getClosestLabelledVertices
from graphbin_utils.graphbin_Options import PARSER


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.6"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"


# Sample command
# -------------------------------------------------------------------
# python graphbin_Miniasm.py    --graph /path/to/graph_file.asqg
#                               --binned /path/to/binning_result.csv
#                               --output /path/to/output_folder
# -------------------------------------------------------------------


def run(args):
    # Setup logger
    #-----------------------

    logger = logging.getLogger("GraphBin %s" % __version__)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    start_time = time.time()

    assembly_graph_file = args.graph
    contigs_file = args.contigs
    contig_bins_file = args.binned
    output_path = args.output
    prefix = args.prefix
    delimiter = args.delimiter
    max_iteration = args.max_iteration
    diff_threshold = args.diff_threshold
    MIN_BIN_COUNT = 10


    # Setup output path for log file
    #---------------------------------------------------

    fileHandler = logging.FileHandler(output_path+"/"+prefix+"graphbin.log")
    fileHandler.setLevel(logging.INFO)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info("Welcome to GraphBin: Refined Binning of Metagenomic Contigs using Assembly Graphs.")
    logger.info("This version of GraphBin makes use of the assembly graph produced by Miniasm.")

    logger.info("Assembly graph file: "+assembly_graph_file)
    logger.info("Existing binning output file: "+contig_bins_file)
    logger.info("Final binning output file: "+output_path)
    logger.info("Maximum number of iterations: "+str(max_iteration))
    logger.info("Difference threshold: "+str(diff_threshold))

    logger.info("GraphBin started")


    # Get the number of bins from the initial binning result
    #--------------------------------------------------------

    n_bins = 0

    try:
        all_bins_list = []

        with open(contig_bins_file) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=delimiter)
            for row in readCSV:
                all_bins_list.append(row[1])
            
        bins_list = list(set(all_bins_list))
        bins_list.sort()

        n_bins = len(bins_list)
        logger.info("Number of bins available in the binning result: "+str(n_bins))
    
    except BaseException as err:
        logger.error(f"Unexpected {err=}, {type(err)=}")
        logger.error("Please make sure that the correct path to the binning result file is provided and it is having the correct format.")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)


    logger.info("Constructing the assembly graph")

    # Get the links from the .gfa file
    #-----------------------------------

    my_map = BidirectionalMap()

    node_count = 0

    nodes = []

    links = []

    try:
        # Get contig connections from .gfa file
        with open(assembly_graph_file) as file:
            
            for line in file.readlines():
                line = line.strip()

                # Count the number of contigs
                if line.startswith("S"):
                    strings = line.split("\t")
                    my_node = strings[1]
                    my_map[node_count] = my_node
                    nodes.append(my_node)
                    node_count += 1

                # Identify lines with link information
                elif line.startswith("L"):
                
                    link = []
                    strings = line.split("\t")
                
                    if strings[1] != strings[3]:
                        start = strings[1]
                        end = strings[3]
                        link.append(start)
                        link.append(end)
                        links.append(link)

    except BaseException as err:
        logger.error(f"Unexpected {err=}, {type(err)=}")
        logger.error("Please make sure that the correct path to the assembly graph file is provided.")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    contigs_map = my_map
    contigs_map_rev = my_map.inverse

    logger.info("Total number of contigs available: "+str(node_count))


    ## Construct the assembly graph
    #-------------------------------

    try:

        # Create the graph
        assembly_graph = Graph()

        # Create list of edges
        edge_list = []

        # Add vertices
        assembly_graph.add_vertices(node_count)

        # Name vertices
        for i in range(len(assembly_graph.vs)):
            assembly_graph.vs[i]["id"]= i
            assembly_graph.vs[i]["label"]= str(contigs_map[i])

        # Iterate links
        for link in links:
            # Remove self loops
            if link[0] != link[1]:
                # Add edge to list of edges
                edge_list.append((contigs_map_rev[link[0]], contigs_map_rev[link[1]]))

        # Add edges to the graph
        assembly_graph.add_edges(edge_list)
        assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    except BaseException as err:
        logger.error(f"Unexpected {err=}, {type(err)=}")
        logger.error("Please make sure that the correct path to the assembly graph file is provided.")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    logger.info("Total number of edges in the assembly graph: "+str(len(edge_list)))


    # Get initial binning result
    #----------------------------

    logger.info("Obtaining the initial binning result")

    bins = [[] for x in range(n_bins)]

    try:
        with open(contig_bins_file) as contig_bins:
            readCSV = csv.reader(contig_bins, delimiter=delimiter)
            for row in readCSV:
                contig_num = contigs_map_rev[row[0]]
            
                bin_num = bins_list.index(row[1])
                bins[bin_num].append(contig_num)

    except BaseException as err:
        logger.error(f"Unexpected {err=}, {type(err)=}")
        logger.error("Please make sure that you have provided the correct assembler type and the correct path to the binning result file in the correct format.")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)


    # Remove labels of ambiguous vertices
    #-------------------------------------


    logger.info("Determining ambiguous vertices")

    remove_by_bin = {}

    remove_labels = []

    neighbours_have_same_label_list = []

    for b in range(n_bins):

        for i in bins[b]:

            my_bin = b

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
                if my_bin in remove_by_bin:
                    if len(bins[my_bin]) - len(remove_by_bin[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.append(i)
                        remove_by_bin[my_bin].append(i)
                else:
                    if len(bins[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.append(i)
                        remove_by_bin[my_bin] = [i]

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
                        if my_bin in remove_by_bin:
                            if len(bins[my_bin]) - len(remove_by_bin[my_bin]) >= MIN_BIN_COUNT:
                                remove_labels.append(i)
                                remove_by_bin[my_bin].append(i)
                        else:
                            if len(bins[my_bin]) >= MIN_BIN_COUNT:
                                remove_labels.append(i)
                                remove_by_bin[my_bin] = [i]

    logger.info("Removing labels of ambiguous vertices")

    # Remove labels of ambiguous vertices
    for i in remove_labels:
        for n in range(n_bins):
            if i in bins[n]:
                bins[n].remove(i)

    logger.info("Obtaining refined binning result")


    # Get vertices which are not isolated and not in components without any labels
    #-----------------------------------------------------------------------------

    logger.info("Deteremining vertices which are not isolated and not in components without any labels")

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

    logger.info("Number of non-isolated contigs: "+str(len(non_isolated)))


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

    # Check if initial binning result consists of contigs belonging to multiple bins

    multiple_bins = False
    
    for item in data:
        if type(item[1]) is int and type(item[2]) is int:
            multiple_bins = True
            break

    if multiple_bins:
        logger.error("Initial binning result consists of contigs belonging to multiple bins. Please make sure that each contig in the initial binning result belongs to only one bin.")
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)
    
    # Label propagation

    lp = LabelProp()

    lp.load_data_from_mem(data)

    logger.info("Starting label propagation with eps="+str(diff_threshold)+" and max_iteration="+str(max_iteration))

    ans = lp.run(diff_threshold, max_iteration, show_log=True, clean_result=False)

    logger.info("Obtaining Label Propagation result")

    for l in ans:
        for i in range(n_bins):
            if l[1]==i+1 and l[0] not in bins[i]:
                bins[i].append(l[0])


    # Remove labels of ambiguous vertices
    #-------------------------------------

    logger.info("Determining ambiguous vertices")

    remove_by_bin = {}
    
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
                if my_bin in remove_by_bin:
                    if len(bins[my_bin]) - len(remove_by_bin[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.append(i)
                        remove_by_bin[my_bin].append(i)
                else:
                    if len(bins[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.append(i)
                        remove_by_bin[my_bin] = [i]

    logger.info("Removing labels of ambiguous vertices")

    # Remove labels of ambiguous vertices
    for i in remove_labels:
        for n in range(n_bins):
            if i in bins[n]:
                bins[n].remove(i)

    elapsed_time = time.time() - start_time

    logger.info("Obtaining the Final Refined Binning result")

    final_bins = {}

    for i in range(n_bins):
        for contig in bins[i]:
            final_bins[contig] = bins_list[i]

    # Print elapsed time for the process
    logger.info("Elapsed time: "+str(elapsed_time)+" seconds")


    # Write result to output file
    #-----------------------------

    logger.info("Writing the Final Binning result to file")

    output_bins = []

    output_bins_path = output_path + prefix + "bins/"
    output_file = output_path + prefix + 'graphbin_output.csv'

    if not os.path.isdir(output_bins_path):
        subprocess.run("mkdir -p "+output_bins_path, shell=True)

    bin_files = {}

    for bin_name in set(final_bins.values()):
        bin_files[bin_name] = open(
            output_bins_path + prefix + "bin_" + bin_name + ".fasta", 'w+')

    for n, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
    
        contig_num = contigs_map_rev[record.id]

        if contig_num in final_bins:
            bin_files[final_bins[contig_num]].write(
                f'>{str(record.id)}\n{str(record.seq)}\n')

    # Close output files
    for c in set(final_bins.values()):
        bin_files[c].close()

    for b in range(len(bins)):

        for contig in bins[b]:
            line = []
            line.append(str(contigs_map[contig]))
            line.append(bins_list[b])
            output_bins.append(line)
    
    with open(output_file, mode='w') as out_file:
        output_writer = csv.writer(out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in output_bins:
            output_writer.writerow(row)
    
    logger.info("Final binning results can be found in "+str(output_bins_path))


    unbinned_contigs = []

    for i in range(node_count):
        if i in remove_labels or i not in non_isolated:
            line = []
            line.append(str(contigs_map[i]))
            unbinned_contigs.append(line)

    if len(unbinned_contigs)!=0:
        unbinned_file = output_path + prefix + 'graphbin_unbinned.csv'

        with open(unbinned_file, mode='w') as out_file:
            output_writer = csv.writer(out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
            for row in unbinned_contigs:
                output_writer.writerow(row)

        logger.info("Unbinned contigs can be found at "+unbinned_file)


    # Exit program
    #--------------

    logger.info("Thank you for using GraphBin! Bye...!")

    logger.removeHandler(fileHandler)
    logger.removeHandler(consoleHeader)


def main():
    # Setup argument parser
    #-----------------------
    ap = PARSER
    args = ap.parse_args()
    run(args)


if __name__ == "__main__":
    main()
