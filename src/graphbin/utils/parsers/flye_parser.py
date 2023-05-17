#!/usr/bin/env python3

import csv
import logging
import os
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
__version__ = "1.7.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"


logger = logging.getLogger("GraphBin %s" % __version__)


def get_initial_binning_result(
    n_bins, bins_list, contig_bins_file, contig_names_rev, delimiter
):
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

    return bins


def parse_graph(assembly_graph_file, contig_paths):
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

    logger.info(f"Total number of contigs available: {node_count}")

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

    logger.info(f"Total number of edges in the assembly graph: {len(edge_list)}")

    return assembly_graph, contig_names, node_count


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

    output_bins_path = f"{output_path}{prefix}bins/"
    output_file = f"{output_path}{prefix}graphbin_output.csv"

    if not os.path.isdir(output_bins_path):
        subprocess.run(f"mkdir -p {output_bins_path}", shell=True)

    bin_files = {}

    for bin_name in set(final_bins.values()):
        bin_files[bin_name] = open(
            f"{output_bins_path}{prefix}bin_{bin_name}.fasta", "w+"
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

    logger.info(f"Final binning results can be found in {output_bins_path}")

    unbinned_contigs = []

    for i in range(node_count):
        if i in remove_labels or i not in non_isolated:
            line = []
            line.append(str(contig_names[i]))
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
