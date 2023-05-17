#!/usr/bin/env python3

import csv
import logging
import os
import re
import subprocess
import sys

from cogent3.parse.fasta import MinimalFastaParser
from igraph import *

from graphbin.utils.bidirectionalmap.bidirectionalmap import BidirectionalMap


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.7.0"
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
                start = "contig-"
                end = ""
                contig_num = contigs_map_rev[
                    int(re.search("%s(.*)%s" % (start, end), row[0]).group(1))
                ]

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
    links = []

    contig_names = BidirectionalMap()

    my_map = BidirectionalMap()

    node_count = 0

    try:
        # Get contig connections from .asqg file
        with open(assembly_graph_file) as file:
            for line in file.readlines():
                line = line.strip()

                # Count the number of contigs
                if line.startswith("VT"):
                    start = "contig-"
                    end = ""
                    contig_name = str(line.split()[1])
                    contig_num = int(
                        re.search("%s(.*)%s" % (start, end), contig_name).group(1)
                    )
                    my_map[node_count] = contig_num
                    contig_names[node_count] = contig_name.strip()
                    node_count += 1

                # Identify lines with link information
                elif line.startswith("ED"):
                    link = []
                    strings = line.split("\t")[1].split()
                    link.append(int(strings[0][7:]))
                    link.append(int(strings[1][7:]))
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

    contig_names_rev = contig_names.inverse

    logger.info(f"Total number of contigs available: {node_count}")

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

    logger.info(f"Total number of edges in the assembly graph: {len(edge_list)}")

    return assembly_graph, contigs_map, contig_names, node_count


def write_output(
    output_path,
    prefix,
    final_bins,
    contigs_file,
    contig_names_rev,
    bins,
    contig_names,
    bins_list,
    delimiter,
    node_count,
    remove_labels,
    non_isolated,
    contig_descriptions,
):
    logger.info("Writing the Final Binning result to file")

    output_bins = []

    output_bins_path = f"{output_path}{prefix}bins/"
    output_file = f"{output_path}{prefix}graphbin_output.csv"

    if not os.path.isdir(output_bins_path):
        subprocess.run(f"mkdir -p {output_bins_path}", shell=True)

    bin_files = {}

    for bin_name in set(final_bins.values()):
        bin_files[bin_name] = open(
            f"{output_bins_path}{prefix}bin_{bin_name}.fasta", "w+"
        )

    for label, seq in MinimalFastaParser(
        contigs_file, label_to_name=lambda x: x.split()[0]
    ):
        contig_num = contig_names_rev[label]

        if contig_num in final_bins:
            bin_files[final_bins[contig_num]].write(f">{label}\n{seq}\n")

    # Close output files
    for c in set(final_bins.values()):
        bin_files[c].close()

    for b in range(len(bins)):
        for contig in bins[b]:
            line = []
            line.append(contig_descriptions[contig_names[contig]])
            line.append(bins_list[b])
            output_bins.append(line)

    with open(output_file, mode="w") as out_file:
        output_writer = csv.writer(
            out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
        )
        for row in output_bins:
            output_writer.writerow(row)

    logger.info(f"Final binning results can be found in {output_bins_path}")

    unbinned_contigs = []

    for i in range(node_count):
        if i in remove_labels or i not in non_isolated:
            line = []
            line.append(contig_names[i])
            unbinned_contigs.append(line)

    if len(unbinned_contigs) != 0:
        unbinned_file = f"{output_path}{prefix}graphbin_unbinned.csv"

        with open(unbinned_file, mode="w") as out_file:
            output_writer = csv.writer(
                out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
            )

            for row in unbinned_contigs:
                output_writer.writerow(row)

        logger.info(f"Unbinned contigs can be found at {unbinned_file}")


def get_contig_descriptions(contigs_file):
    contig_descriptions = {}

    for label, _ in MinimalFastaParser(contigs_file):
        name = label.split()[0]
        contig_descriptions[name] = label

    return contig_descriptions
