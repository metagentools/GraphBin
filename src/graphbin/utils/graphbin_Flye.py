#!/usr/bin/env python3

"""graphbin_Flye.py: Refined binning of metagenomic contigs using Flye assembly graphs.

GraphBin is a metagenomic contig binning tool that makes use of the contig 
connectivity information from the assembly graph to bin contigs. It utilizes 
the binning result of an existing binning tool and a label propagation algorithm 
to correct mis-binned contigs and predict the labels of contigs which are 
discarded due to short length.

graphbin_Flye.py makes use of the assembly graphs produced by Flye long read assembler.
"""

import csv
import logging
import os
import subprocess
import sys
import time

from collections import defaultdict

from cogent3.parse.fasta import MinimalFastaParser
from igraph import *

from graphbin.utils.bidirectionalmap.bidirectionalmap import BidirectionalMap
from graphbin.utils.graphbin_Func import getClosestLabelledVertices
from graphbin.utils.graphbin_Options import PARSER
from graphbin.utils.labelpropagation.labelprop import LabelProp


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.6.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"


def run(args):
    # Setup logger
    # -----------------------
    logger = logging.getLogger("GraphBin %s" % __version__)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    consoleHeader = logging.StreamHandler()
    consoleHeader.setFormatter(formatter)
    consoleHeader.setLevel(logging.INFO)
    logger.addHandler(consoleHeader)

    start_time = time.time()

    assembly_graph_file = args.graph
    contigs_file = args.contigs
    contig_paths = args.paths
    contig_bins_file = args.binned
    output_path = args.output
    prefix = args.prefix
    delimiter = args.delimiter
    max_iteration = args.max_iteration
    diff_threshold = args.diff_threshold
    MIN_BIN_COUNT = 10

    # Setup output path for log file
    # ---------------------------------------------------

    fileHandler = logging.FileHandler(output_path + "/" + prefix + "graphbin.log")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    logger.info(
        "Welcome to GraphBin: Refined Binning of Metagenomic Contigs using Assembly Graphs."
    )
    logger.info(
        "This version of GraphBin makes use of the assembly graph produced by Flye which is a long reads assembler based on the de Bruijn graph approach."
    )

    logger.info("Assembly graph file: " + assembly_graph_file)
    logger.info("Existing binning output file: " + contig_bins_file)
    logger.info("Contig paths file: " + contig_paths)
    logger.info("Final binning output file: " + output_path)
    logger.info("Maximum number of iterations: " + str(max_iteration))
    logger.info("Difference threshold: " + str(diff_threshold))

    logger.info("GraphBin started")

    # Get the number of bins from the initial binning result
    # --------------------------------------------------------

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
        logger.info("Number of bins available in the binning result: " + str(n_bins))

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the binning result file is provided and it is having the correct format."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    logger.info("Constructing the assembly graph")

    # Get contig names
    # -----------------------------------

    contig_names = BidirectionalMap()

    contig_num = 0

    with open(contig_paths, "r") as file:

        for line in file.readlines():

            if not line.startswith("#"):
                name = line.strip().split()[0]
                contig_names[contig_num] = name
                contig_num += 1

    contig_names_rev = contig_names.inverse

    # Get the paths and edges
    # -----------------------------------

    paths = {}
    segment_contigs = {}

    try:

        with open(contig_paths) as file:

            for line in file.readlines():

                if not line.startswith("#"):

                    strings = line.strip().split()

                    contig_name = strings[0]

                    path = strings[-1]
                    path = path.replace("*", "")

                    if path.startswith(","):
                        path = path[1:]

                    if path.endswith(","):
                        path = path[:-1]

                    segments = path.rstrip().split(",")

                    contig_num = contig_names_rev[contig_name]

                    if contig_num not in paths:
                        paths[contig_num] = segments

                    for segment in segments:

                        if segment != "":

                            if segment not in segment_contigs:
                                segment_contigs[segment] = set([contig_num])
                            else:
                                segment_contigs[segment].add(contig_num)

        links_map = defaultdict(set)

        # Get links from assembly_graph.gfa
        with open(assembly_graph_file) as file:

            for line in file.readlines():
                line = line.strip()

                # Identify lines with link information
                if line.startswith("L"):
                    strings = line.split("\t")

                    f1, f2 = "", ""

                    if strings[2] == "+":
                        f1 = strings[1][5:]
                    if strings[2] == "-":
                        f1 = "-" + strings[1][5:]
                    if strings[4] == "+":
                        f2 = strings[3][5:]
                    if strings[4] == "-":
                        f2 = "-" + strings[3][5:]

                    links_map[f1].add(f2)
                    links_map[f2].add(f1)

        # Create list of edges
        edge_list = []

        for i in paths:
            segments = paths[i]

            new_links = []

            for segment in segments:

                my_segment = segment
                my_segment_num = ""

                my_segment_rev = ""

                if my_segment.startswith("-"):
                    my_segment_rev = my_segment[1:]
                    my_segment_num = my_segment[1:]
                else:
                    my_segment_rev = "-" + my_segment
                    my_segment_num = my_segment

                if my_segment in links_map:
                    new_links.extend(list(links_map[my_segment]))

                if my_segment_rev in links_map:
                    new_links.extend(list(links_map[my_segment_rev]))

                if my_segment in segment_contigs:
                    for contig in segment_contigs[my_segment]:
                        if i != contig:
                            # Add edge to list of edges
                            edge_list.append((i, contig))

                if my_segment_rev in segment_contigs:
                    for contig in segment_contigs[my_segment_rev]:
                        if i != contig:
                            # Add edge to list of edges
                            edge_list.append((i, contig))

                if my_segment_num in segment_contigs:
                    for contig in segment_contigs[my_segment_num]:
                        if i != contig:
                            # Add edge to list of edges
                            edge_list.append((i, contig))

            for new_link in new_links:
                if new_link in segment_contigs:
                    for contig in segment_contigs[new_link]:
                        if i != contig:
                            # Add edge to list of edges
                            edge_list.append((i, contig))

                if new_link.startswith("-"):
                    if new_link[1:] in segment_contigs:
                        for contig in segment_contigs[new_link[1:]]:
                            if i != contig:
                                # Add edge to list of edges
                                edge_list.append((i, contig))

        node_count = len(contig_names_rev)

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the assembly graph file is provided."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    logger.info("Total number of contigs available: " + str(node_count))

    ## Construct the assembly graph
    # -------------------------------

    try:

        # Create the graph
        assembly_graph = Graph()

        # Add vertices
        assembly_graph.add_vertices(node_count)

        # Name vertices
        for i in range(len(assembly_graph.vs)):
            assembly_graph.vs[i]["id"] = i
            assembly_graph.vs[i]["label"] = str(contig_names[i])

        # Add edges to the graph
        assembly_graph.add_edges(edge_list)
        assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the assembly graph file is provided."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    logger.info("Total number of edges in the assembly graph: " + str(len(edge_list)))

    # Get initial binning result
    # ----------------------------

    logger.info("Obtaining the initial binning result")

    bins = [[] for x in range(n_bins)]

    try:
        with open(contig_bins_file) as contig_bins:
            readCSV = csv.reader(contig_bins, delimiter=delimiter)
            for row in readCSV:
                contig_num = contig_names_rev[row[0]]

                bin_num = bins_list.index(row[1])
                bins[bin_num].append(contig_num)

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that you have provided the correct assembler type and the correct path to the binning result file in the correct format."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    # Remove labels of ambiguous vertices
    # -------------------------------------

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
        binned_contigs = sorted(binned_contigs + bins[n])

    for b in range(n_bins):

        for i in bins[b]:

            if i not in neighbours_have_same_label_list:

                my_bin = b

                # Get set of closest labelled vertices
                closest_neighbours = getClosestLabelledVertices(
                    assembly_graph, i, binned_contigs
                )

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
                            if (
                                len(bins[my_bin]) - len(remove_by_bin[my_bin])
                                >= MIN_BIN_COUNT
                            ):
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
    # -----------------------------------------------------------------------------

    logger.info(
        "Deteremining vertices which are not isolated and not in components without any labels"
    )

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

            while length != len(component):

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

    logger.info("Number of non-isolated contigs: " + str(len(non_isolated)))

    # Run label propagation
    # -----------------------

    data = []

    for contig in range(node_count):

        # Consider vertices that are not isolated

        if contig in non_isolated:
            line = []
            line.append(contig)

            assigned = False

            for i in range(n_bins):
                if contig in bins[i]:
                    line.append(i + 1)
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
        logger.error(
            "Initial binning result consists of contigs belonging to multiple bins. Please make sure that each contig in the initial binning result belongs to only one bin."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    # Label propagation

    lp = LabelProp()

    lp.load_data_from_mem(data)

    logger.info(
        "Starting label propagation with eps="
        + str(diff_threshold)
        + " and max_iteration="
        + str(max_iteration)
    )

    ans = lp.run(diff_threshold, max_iteration, show_log=True, clean_result=False)

    logger.info("Obtaining Label Propagation result")

    for l in ans:
        for i in range(n_bins):
            if l[1] == i + 1 and l[0] not in bins[i]:
                bins[i].append(l[0])

    # Remove labels of ambiguous vertices
    # -------------------------------------

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
    logger.info("Elapsed time: " + str(elapsed_time) + " seconds")

    # Write result to output file
    # -----------------------------

    logger.info("Writing the Final Binning result to file")

    output_bins = []

    output_bins_path = output_path + prefix + "bins/"
    output_file = output_path + prefix + "graphbin_output.csv"

    if not os.path.isdir(output_bins_path):
        subprocess.run("mkdir -p " + output_bins_path, shell=True)

    bin_files = {}

    for bin_name in set(final_bins.values()):
        bin_files[bin_name] = open(
            output_bins_path + prefix + "bin_" + bin_name + ".fasta", "w+"
        )

    for label, seq in MinimalFastaParser(contigs_file):

        contig_num = contig_names_rev[label]

        if contig_num in final_bins:
            bin_files[final_bins[contig_num]].write(f">{label}\n{seq}\n")

    # Close output files
    for c in set(final_bins.values()):
        bin_files[c].close()

    for b in range(len(bins)):

        # with open(output_bins_path + "bin_" + str(b+1) + "_ids.txt", "w") as bin_file:
        for contig in bins[b]:
            line = []
            line.append(str(contig_names[contig]))
            line.append(bins_list[b])
            output_bins.append(line)
        #         bin_file.write(str(contig_names[contig])+"\n")

        # subprocess.run("awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' " + output_bins_path + "bin_" + str(b+1) + "_ids.txt " + contigs_file + " > " + output_bins_path + "bin_" +str(b+1) +"_seqs.fasta", shell=True)

    with open(output_file, mode="w") as out_file:
        output_writer = csv.writer(
            out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
        )
        for row in output_bins:
            output_writer.writerow(row)

    logger.info("Final binning results can be found in " + str(output_bins_path))

    unbinned_contigs = []

    for i in range(node_count):
        if i in remove_labels or i not in non_isolated:
            line = []
            line.append(str(contig_names[i]))
            unbinned_contigs.append(line)

    if len(unbinned_contigs) != 0:
        unbinned_file = output_path + prefix + "graphbin_unbinned.csv"

        with open(unbinned_file, mode="w") as out_file:
            output_writer = csv.writer(
                out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
            )

            for row in unbinned_contigs:
                output_writer.writerow(row)

        logger.info("Unbinned contigs can be found at " + unbinned_file)

    # Exit program
    # --------------

    logger.info("Thank you for using GraphBin! Bye...!")

    logger.removeHandler(fileHandler)
    logger.removeHandler(consoleHeader)


def main():
    # Setup argument parser
    # -----------------------
    ap = PARSER
    args = ap.parse_args()
    run(args)


if __name__ == "__main__":
    main()
