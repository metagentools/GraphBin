#!/usr/bin/env python3

import csv
import logging
import os
import re
import subprocess
import sys

from collections import defaultdict

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
                start = "NODE_"
                end = "_length_"
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


def parse_graph(assembly_graph_file, contig_paths):
    paths = {}
    segment_contigs = {}
    node_count = 0

    contig_names = BidirectionalMap()

    my_map = BidirectionalMap()

    current_contig_num = ""

    try:
        with open(contig_paths) as file:
            name = file.readline()
            path = file.readline()

            while name != "" and path != "":
                while ";" in path:
                    path = path[:-2] + "," + file.readline()

                start = "NODE_"
                end = "_length_"
                contig_num = str(
                    int(re.search("%s(.*)%s" % (start, end), name).group(1))
                )

                segments = path.rstrip().split(",")

                if current_contig_num != contig_num:
                    my_map[node_count] = int(contig_num)
                    current_contig_num = contig_num
                    contig_names[node_count] = name.strip()
                    node_count += 1

                if contig_num not in paths:
                    paths[contig_num] = [segments[0], segments[-1]]

                for segment in segments:
                    if segment not in segment_contigs:
                        segment_contigs[segment] = set([contig_num])
                    else:
                        segment_contigs[segment].add(contig_num)

                name = file.readline()
                path = file.readline()

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the contig paths file is provided."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    contigs_map = my_map
    contigs_map_rev = my_map.inverse

    logger.info("Total number of contigs available: " + str(node_count))

    links = []
    links_map = defaultdict(set)

    ## Construct the assembly graph
    # -------------------------------

    try:
        # Get links from assembly_graph_with_scaffolds.gfa
        with open(assembly_graph_file) as file:
            for line in file.readlines():
                line = line.strip()

                # Identify lines with link information
                if line.startswith("L"):
                    strings = line.split("\t")
                    f1, f2 = strings[1] + strings[2], strings[3] + strings[4]
                    links_map[f1].add(f2)
                    links_map[f2].add(f1)
                    links.append(
                        strings[1] + strings[2] + " " + strings[3] + strings[4]
                    )

        # Create graph
        assembly_graph = Graph()

        # Add vertices
        assembly_graph.add_vertices(node_count)

        # Create list of edges
        edge_list = []

        # Name vertices
        for i in range(node_count):
            assembly_graph.vs[i]["id"] = i
            assembly_graph.vs[i]["label"] = str(contigs_map[i])

        for i in range(len(paths)):
            segments = paths[str(contigs_map[i])]

            start = segments[0]
            start_rev = ""

            if start.endswith("+"):
                start_rev = start[:-1] + "-"
            else:
                start_rev = start[:-1] + "+"

            end = segments[1]
            end_rev = ""

            if end.endswith("+"):
                end_rev = end[:-1] + "-"
            else:
                end_rev = end[:-1] + "+"

            new_links = []

            if start in links_map:
                new_links.extend(list(links_map[start]))
            if start_rev in links_map:
                new_links.extend(list(links_map[start_rev]))
            if end in links_map:
                new_links.extend(list(links_map[end]))
            if end_rev in links_map:
                new_links.extend(list(links_map[end_rev]))

            for new_link in new_links:
                if new_link in segment_contigs:
                    for contig in segment_contigs[new_link]:
                        if i != contigs_map_rev[int(contig)]:
                            # Add edge to list of edges
                            edge_list.append((i, contigs_map_rev[int(contig)]))

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

    for label, seq in MinimalFastaParser(contigs_file):
        contig_num = contig_names_rev[label]

        if contig_num in final_bins:
            bin_files[final_bins[contig_num]].write(f">{label}\n{seq}\n")

    # Close output files
    for c in set(final_bins.values()):
        bin_files[c].close()

    for b in range(len(bins)):
        for contig in bins[b]:
            line = []
            line.append(contig_names[contig])
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
            line.append(contig_names[i])
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
