#!/usr/bin/env python3

"""graphbin_Canu.py: Refined binning of metagenomic contigs using Canu assembly graphs.

GraphBin is a metagenomic contig binning tool that makes use of the contig 
connectivity information from the assembly graph to bin contigs. It utilizes 
the binning result of an existing binning tool and a label propagation algorithm 
to correct mis-binned contigs and predict the labels of contigs which are 
discarded due to short length.

graphbin_Canu.py makes use of the assembly graphs produced by Canu long read assembler.
"""

import logging
import time

from graphbin.graphbin_Func import graphbin_main
from graphbin.parsers import get_initial_bin_count
from graphbin.parsers.canu_parser import (
    get_initial_binning_result,
    parse_graph,
    write_output,
)


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.7.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"

# create logger
logger = logging.getLogger(f"GraphBin {__version__}")


def run(args):
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

    logger.info(
        "Welcome to GraphBin: Refined Binning of Metagenomic Contigs using Assembly Graphs."
    )
    logger.info(
        "This version of GraphBin makes use of the assembly graph produced by Canu which is a long reads assembler based on the OLC approach."
    )

    logger.info(f"Assembly graph file: {assembly_graph_file}")
    logger.info(f"Existing binning output file: {contig_bins_file}")
    logger.info(f"Final binning output file: {output_path}")
    logger.info(f"Maximum number of iterations: {max_iteration}")
    logger.info(f"Difference threshold: {diff_threshold}")

    logger.info("GraphBin started")

    # Get the number of bins from the initial binning result
    # --------------------------------------------------------

    n_bins, bins_list = get_initial_bin_count(contig_bins_file, delimiter)

    # Get assembly graph
    # --------------------

    assembly_graph, contigs_map, node_count = parse_graph(assembly_graph_file)

    # Get initial binning result
    # ----------------------------

    bins = get_initial_binning_result(
        n_bins, bins_list, contig_bins_file, contigs_map.inverse, delimiter
    )

    # Run GraphBin logic
    # -------------------------------------

    final_bins, remove_labels, non_isolated = graphbin_main(
        n_bins,
        bins,
        bins_list,
        assembly_graph,
        node_count,
        diff_threshold,
        max_iteration,
    )

    elapsed_time = time.time() - start_time

    # Print elapsed time for the process
    logger.info(f"Elapsed time: {elapsed_time} seconds")

    # Write result to output file
    # -----------------------------

    write_output(
        output_path,
        prefix,
        final_bins,
        contigs_file,
        contigs_map.inverse,
        bins,
        contigs_map,
        bins_list,
        delimiter,
        node_count,
        remove_labels,
        non_isolated,
    )


def main(args):
    run(args)


if __name__ == "__main__":
    main()
