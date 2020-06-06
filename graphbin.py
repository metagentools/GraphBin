#!/usr/bin/python3

"""graphbin.py: Refined binning of metagenomic contigs using assembly graphs."""

import argparse
import os
import sys
import subprocess

__author__ = "Vijini Mallawaarachchi, Anuradha Wickramarachchi, and Yu Lin"
__copyright__ = "Copyright 2020, GraphBin Project"
__credits__ = "Benjamin Kaehler and Gavin Huttley"
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Stable Release"

parser = argparse.ArgumentParser(description="""GraphBin Help. GraphBin is a metagenomic contig binning tool 
that makes use of the contig connectivity information from the assembly graph to bin contigs. It utilizes the 
binning result of an existing binning tool and a label propagation algorithm to correct mis-binned contigs 
and predict the labels of contigs which are discarded due to short length.""")

parser.add_argument("--assembler", 
                    required=True, 
                    type=str,
                    help="name of the assembler used (SPAdes, SGA or MEGAHIT). GraphBin supports Flye, Canu and Miniasm long-read assemblies as well.")

parser.add_argument("--graph", 
                    required=True,
                    type=str,
                    help="path to the assembly graph file")

parser.add_argument("--paths", 
                    required=False,
                    type=str,
                    help="path to the contigs.paths file")

parser.add_argument("--binned", 
                    required=True,
                    type=str,
                    help="path to the .csv file with the initial binning output from an existing tool")

parser.add_argument("--output", 
                    required=True,
                    type=str,
                    help="path to the output folder")

parser.add_argument("--prefix", 
                    required=False,
                    type=str,
                    default='',
                    help="prefix for the output file")

parser.add_argument("--max_iteration", 
                    required=False, 
                    type=int, 
                    default=100, 
                    help="maximum number of iterations for label propagation algorithm. [default: 100]")

parser.add_argument("--diff_threshold", 
                    required=False, 
                    type=float, 
                    default=0.1, 
                    help="difference threshold for label propagation algorithm. [default: 0.1]")

args = vars(parser.parse_args())


assembler = args["assembler"]
assembly_graph_file = args["graph"]
contig_paths = args["paths"]
contig_bins_file = args["binned"]
output_path = args["output"]
prefix = args["prefix"]
max_iteration = args["max_iteration"]
diff_threshold = args["diff_threshold"]


# Validation of inputs
#---------------------------------------------------

# Check assembler type
if not (assembler.lower() == "spades" or assembler.lower() == "sga" or assembler.lower() == "megahit" or assembler.lower() == "flye" or assembler.lower() == "canu" or assembler.lower() == "miniasm"):
    print("\nPlease make sure to provide the correct assembler type (SPAdes, SGA or MEGAHIT). GraphBin supports Flye, Canu and Miniasm long-read assemblies as well.")
    print("Exiting GraphBin...\nBye...!\n")
    sys.exit(1)

# Check assembly graph file
if not os.path.isfile(assembly_graph_file):
    print("\nFailed to open the assembly graph file.")
    print("Exiting GraphBin...\nBye...!\n")
    sys.exit(1)

# Check if paths files is provided when the assembler type is SPAdes
if assembler.lower() == "spades" and contig_paths is None:
    print("\nPlease make sure to provide the path to the contigs.paths file.")
    print("Exiting GraphBin...\nBye...!\n")
    sys.exit(1)

# Check contigs.paths file for SPAdes
if assembler.lower() == "spades" and not os.path.isfile(contig_paths):
    print("\nFailed to open the contigs.paths file.")
    print("Exiting GraphBin...\nBye...!\n")
    sys.exit(1)

# Check the file with the initial binning output
if not os.path.isfile(contig_bins_file):
    print("\nFailed to open the file with the initial binning output.")
    print("Exiting GraphBin...\nBye...!\n")
    sys.exit(1)

# Handle for missing trailing forwardslash in output folder path
if output_path[-1:] != "/":
    output_path = output_path + "/"

# Create output folder if it does not exist
if not os.path.isdir(output_path):
    subprocess.run("mkdir -p "+output_path, shell=True)

# Validate prefix
if args["prefix"] != '':
    if args["prefix"].endswith("_"):
        prefix = args["prefix"]
    else:
        prefix = args["prefix"]+"_"
else:
    prefix = ''

# Validate max_iteration
if args["max_iteration"] <= 0:
    print("\nPlease enter a valid number for max_iteration")
    print("Exiting GraphBin...\nBye...!\n")
    sys.exit(1)

# Validate diff_threshold
if args["diff_threshold"] < 0:
    print("\nPlease enter a valid number for diff_threshold")
    print("Exiting GraphBin...\nBye...!\n")
    sys.exit(1)


# Run GraphBin
#---------------------------------------------------
if assembler.lower() == "spades":
    cmdGraphBin = """python "{0}src/graphbin_SPAdes.py" --graph "{1}" --paths "{2}" --binned "{3}" --output "{4}" --prefix "{5}" --max_iteration "{6}" --diff_threshold "{7}" """.format(
        os.path.dirname(__file__), 
        assembly_graph_file, 
        contig_paths, 
        contig_bins_file, 
        output_path,
        prefix,
        max_iteration,
        diff_threshold)

elif assembler.lower() == "sga":
    cmdGraphBin = """python "{0}src/graphbin_SGA.py" --graph "{1}" --binned "{2}" --output "{3}" --prefix "{4}" --max_iteration "{5}" --diff_threshold "{6}" """.format(
        os.path.dirname(__file__), 
        assembly_graph_file,
        contig_bins_file, 
        output_path,
        prefix,
        max_iteration,
        diff_threshold)

elif assembler.lower() == "megahit":
    cmdGraphBin = """python "{0}src/graphbin_MEGAHIT.py" --graph "{1}" --binned "{2}" --output "{3}" --prefix "{4}" --max_iteration "{5}" --diff_threshold "{6}" """.format(
        os.path.dirname(__file__), 
        assembly_graph_file,
        contig_bins_file, 
        output_path,
        prefix,
        max_iteration,
        diff_threshold)

elif assembler.lower() == "flye":
    cmdGraphBin = """python "{0}src/graphbin_Flye.py" --graph "{1}" --binned "{2}" --output "{3}" --prefix "{4}" --max_iteration "{5}" --diff_threshold "{6}" """.format(
        os.path.dirname(__file__), 
        assembly_graph_file,
        contig_bins_file, 
        output_path,
        prefix,
        max_iteration,
        diff_threshold)

elif assembler.lower() == "canu":
    cmdGraphBin = """python "{0}src/graphbin_Canu.py" --graph "{1}" --binned "{2}" --output "{3}" --prefix "{4}" --max_iteration "{5}" --diff_threshold "{6}" """.format(
        os.path.dirname(__file__), 
        assembly_graph_file,
        contig_bins_file, 
        output_path,
        prefix,
        max_iteration,
        diff_threshold)

elif assembler.lower() == "miniasm":
    cmdGraphBin = """python "{0}src/graphbin_Miniasm.py" --graph "{1}" --binned "{2}" --output "{3}" --prefix "{4}" --max_iteration "{5}" --diff_threshold "{6}" """.format(
        os.path.dirname(__file__), 
        assembly_graph_file,
        contig_bins_file, 
        output_path,
        prefix,
        max_iteration,
        diff_threshold)

os.system(cmdGraphBin)
