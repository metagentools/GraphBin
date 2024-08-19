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
__version__ = "1.7.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


logger = logging.getLogger(f"GraphBin {__version__}")


def get_initial_binning_result(
    n_bins,
    bins_list,
    contig_bins_file,
    contigs_map_rev,
    graph_to_contig_map_rev,
    delimiter,
):
    logger.info("Obtaining the initial binning result")

    bins = [[] for x in range(n_bins)]

    try:
        with open(contig_bins_file) as contig_bins:
            readCSV = csv.reader(contig_bins, delimiter=delimiter)
            for row in readCSV:
                contig_num = contigs_map_rev[int(graph_to_contig_map_rev[row[0]])]

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


def parse_graph(assembly_graph_file, original_contigs):
    node_count = 0

    graph_contigs = {}

    links = []

    my_map = BidirectionalMap()

    try:
        # Get links from .gfa file
        with open(assembly_graph_file) as file:
            for line in file.readlines():
                line = line.strip()

                # Identify lines with link information
                if line.startswith("L"):
                    link = []

                    strings = line.split("\t")

                    start_1 = "NODE_"
                    end_1 = "_length"

                    link1 = int(
                        re.search("%s(.*)%s" % (start_1, end_1), strings[1]).group(1)
                    )

                    start_2 = "NODE_"
                    end_2 = "_length"

                    link2 = int(
                        re.search("%s(.*)%s" % (start_2, end_2), strings[3]).group(1)
                    )

                    link.append(link1)
                    link.append(link2)
                    links.append(link)

                elif line.startswith("S"):
                    strings = line.split()

                    start = "NODE_"
                    end = "_length"

                    contig_num = int(
                        re.search("%s(.*)%s" % (start, end), strings[1]).group(1)
                    )

                    my_map[node_count] = int(contig_num)

                    graph_contigs[contig_num] = strings[2]

                    node_count += 1

        logger.info(f"Total number of contigs available: {node_count}")

        contigs_map = my_map
        contigs_map_rev = my_map.inverse

        # Create graph
        assembly_graph = Graph()

        # Add vertices
        assembly_graph.add_vertices(node_count)

        # Create list of edges
        edge_list = []

        for i in range(node_count):
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

    # Map original contig IDs to contig IDS of assembly graph
    # --------------------------------------------------------

    graph_to_contig_map = BidirectionalMap()

    for (n, m), (n2, m2) in zip(graph_contigs.items(), original_contigs.items()):
        if m == m2:
            graph_to_contig_map[n] = n2

    return assembly_graph, graph_to_contig_map, contigs_map, node_count


def write_output(
    output_path,
    prefix,
    final_bins,
    contigs_file,
    graph_to_contig_map,
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
    graph_to_contig_map_rev = graph_to_contig_map.inverse
    contigs_map_rev = contigs_map.inverse

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
        contig_num = contigs_map_rev[graph_to_contig_map_rev[label]]

        if contig_num in final_bins:
            bin_files[final_bins[contig_num]].write(f">{label}\n{seq}\n")

    # Close output files
    for c in set(final_bins.values()):
        bin_files[c].close()

    for b in range(len(bins)):
        for contig in bins[b]:
            line = []
            line.append(graph_to_contig_map[contigs_map[contig]])
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
            line.append(graph_to_contig_map[contigs_map[i]])
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


def get_contig_descriptors(contigs_file):
    original_contigs = {}
    contig_descriptions = {}

    for label, seq in MinimalFastaParser(contigs_file):
        name = label.split()[0]
        original_contigs[name] = seq
        contig_descriptions[name] = label

    return original_contigs
