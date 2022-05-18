#!/usr/bin/python3

"""miniasm_gfa2fasta.py: Obtain the sequences corresponding to edges in the Miniasm assembly graphs in FASTA format.

The assembly graph file (assembly_graph.gfa) should be provided as inputs.

"""

import sys
import os
import argparse
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "GPL"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"


# Sample command
# -------------------------------------------------------------------
# python miniasm_gfa2fasta.py  --graph /path/to/folder_with_binning_result
#                                   --output /path/to/output_folder
# -------------------------------------------------------------------


# Setup argument parser
#-----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--graph", required=True, help="path to the assembly graph file")
ap.add_argument("--assembler", required=True, type=str, default='flye', help="type of the assembler (Flye or Miniasm)")
ap.add_argument("--output", required=True, type=str, help="path to the output folder")
ap.add_argument("--prefix", required=False, type=str, default='', help="prefix for the output file")

args = vars(ap.parse_args())

assembler = args["assembler"]
assembler_name = ""
assembly_graph_file = args["graph"]
output_path = args["output"]
prefix = ""

# Check assembly graph file
if not os.path.isfile(assembly_graph_file):
    print("\nFailed to open the assembly graph file.")
    print("Exiting miniasm_gfa2fasta.py...\nBye...!\n")
    sys.exit(1)

# Check assembler type
if assembler.lower() == "flye":
    assembler_name = "Flye"
elif assembler.lower() == "miniasm":
    assembler_name = "Miniasm"


# Check if output folder exists
#---------------------------------------------------

# Handle for missing trailing forwardslash in output folder path
if output_path[-1:] != "/":
    output_path = output_path + "/"

# Create output folder if it does not exist
if not os.path.isdir(output_path):
    subprocess.run("mkdir -p "+output_path, shell=True)

# Validate prefix
#---------------------------------------------------
try:

    if args["prefix"] != '':
        if args["prefix"].endswith("_"):
            prefix = args["prefix"]
        else:
            prefix = args["prefix"]+"_"
    else:
        prefix = ""

except:
    print("\nPlease enter a valid string for prefix")
    print("Exiting miniasm_gfa2fasta.py...\n")
    sys.exit(1)


# Get the sequences corresponding to edges of the graph.
#---------------------------------------------------

print("\nObtaining edge sequences")

sequenceset = []

with open(assembly_graph_file) as file:
    
    for line in file.readlines():
        line = line.strip()

        if line.startswith("S"):
            
            strings = line.split("\t")

            print(strings)
            
            record = SeqRecord(Seq(re.sub('[^GATC]',"",str(strings[2]).upper())), id=str(strings[1]), name=str(strings[1]), description="")
            
            sequenceset.append(record)

print("\nWriting edge sequences to FASTA file")

if assembler.lower() == "flye":
    final_file = "edges.fasta"
elif assembler.lower() == "miniasm":
    final_file = "unitigs.fasta"

with open(output_path + prefix + final_file, "w") as output_handle:
    SeqIO.write(sequenceset, output_handle, "fasta")

print("\nThe FASTA file with", assembler_name, "sequences can be found at", output_handle.name)


# Exit program
#--------------

print("\nThank you for using miniasm_gfa2fasta for GraphBin!\n")