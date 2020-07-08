#!/usr/bin/env python3

import argparse

PARSER = argparse.ArgumentParser(
    description="""GraphBin Help. GraphBin is a metagenomic contig binning tool
    that makes use of the contig connectivity information from the assembly graph to bin contigs. It utilizes the
    binning result of an existing binning tool and a label propagation algorithm to correct mis-binned contigs
    and predict the labels of contigs which are discarded due to short length.""",
    prog="graphbin",
)

PARSER.add_argument("--version", default=False, action="store_true")

PARSER.add_argument("--graph", type=str, help="path to the assembly graph file")

PARSER.add_argument(
    "--binned",
    type=str,
    help="path to the .csv file with the initial binning output from an existing tool",
)

PARSER.add_argument("--output", type=str, help="path to the output folder")

PARSER.add_argument("--prefix", type=str, default="", help="prefix for the output file")

PARSER.add_argument(
    "--max_iteration",
    type=int,
    default=100,
    help="maximum number of iterations for label propagation algorithm. [default: 100]",
)

PARSER.add_argument(
    "--diff_threshold",
    type=float,
    default=0.1,
    help="difference threshold for label propagation algorithm. [default: 0.1]",
)
