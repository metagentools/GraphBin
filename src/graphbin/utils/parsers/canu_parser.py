#!/usr/bin/env python3

import csv
import logging
import os
import subprocess
import sys

from cogent3.parse.fasta import MinimalFastaParser
from igraph import *

from graphbin.utils.bidirectionalmap.bidirectionalmap import BidirectionalMap


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.6.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"


logger = logging.getLogger("GraphBin %s" % __version__)


def get_initial_binning_result(
    n_bins, bins_list, contig_bins_file, contigs_map_rev, delimiter
):
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
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that you have provided the correct assembler type and the correct path to the binning result file in the correct format."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    return bins


def parse_graph(assembly_graph_file):
    # Get the links from the .gfa file
    # -----------------------------------

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
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the assembly graph file is provided."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    contigs_map = my_map
    contigs_map_rev = my_map.inverse

    logger.info("Total number of contigs available: " + str(node_count))

    ## Construct the assembly graph
    # -------------------------------

    try:
        # Create the graph
        assembly_graph = Graph()

        # Create list of edges
        edge_list = []

        # Add vertices
        assembly_graph.add_vertices(node_count)

        # Name vertices
        for i in range(len(assembly_graph.vs)):
            assembly_graph.vs[i]["id"] = i
            assembly_graph.vs[i]["label"] = str(contigs_map[i])

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
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the assembly graph file is provided."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    logger.info("Total number of edges in the assembly graph: " + str(len(edge_list)))

    return assembly_graph, contigs_map, node_count


def write_output(
    output_path,
    prefix,
    final_bins,
    contigs_file,
    contigs_map_rev,
    bins,
    contigs_map,
    bins_list,
    delimiter,
    node_count,
    remove_labels,
    non_isolated,
):
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

    for label, seq in MinimalFastaParser(
        contigs_file, label_to_name=lambda x: x.split()[0]
    ):
        contig_num = contigs_map_rev[label]

        if contig_num in final_bins:
            bin_files[final_bins[contig_num]].write(f">{label}\n{seq}\n")

    # Close output files
    for c in set(final_bins.values()):
        bin_files[c].close()

    for b in range(len(bins)):
        for contig in bins[b]:
            line = []
            line.append(str(contigs_map[contig]))
            line.append(bins_list[b])
            output_bins.append(line)

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
            line.append(str(contigs_map[i]))
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
