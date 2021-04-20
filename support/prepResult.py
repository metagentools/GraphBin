#!/usr/bin/python3

"""prepResult.py: Format the initial binning result from an existing binning tool.

Format the initial binning result from an existing binning tool in the .csv format
with contig ID and bin ID. Contigs are numbered starting from 0 and bins are 
numbered starting from 1.

"""

import sys
import os
import csv
import argparse
import re
import subprocess

from Bio import SeqIO

__author__ = "Vijini Mallawaarachchi, Anuradha Wickramarachchi, and Yu Lin"
__copyright__ = "Copyright 2019, GraphBin Project"
__credits__ = "Benjamin Kaehler and Gavin Huttley"
__license__ = "GPL"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"


# Sample command
# -------------------------------------------------------------------
# python prepResult.py     --binned /path/to/folder_with_binning_result
#                          --output /path/to/output_folder
# -------------------------------------------------------------------


# Setup argument parser
#-----------------------

ap = argparse.ArgumentParser()

ap.add_argument("--binned", required=True, type=str, help="path to the folder containing the initial binning result from an existing tool")
ap.add_argument("--output", required=True, type=str, help="path to the output folder")
ap.add_argument("--prefix", required=False, type=str, default='', help="prefix for the output file")

args = vars(ap.parse_args())

contig_bins_folder = args["binned"]
output_path = args["output"]
prefix = ""


# Check if folder to initial binning result exists
#---------------------------------------------------

# Handle for missing trailing forwardslash in folder path of binning result
if contig_bins_folder[-1:] != "/":
    contig_bins_folder = contig_bins_folder + "/"

# Throw an error if folder path of binning result does not exist.
if not os.path.isdir(contig_bins_folder):
    print("\nPlease enter a valid path to the folder containing the initial binning result.")
    print("\nExiting prepResult.py...\nBye...!\n")
    sys.exit(1)

# Get list of files in the folder path of binning result.
files = os.listdir(contig_bins_folder)


# Check if folder path of binning result is empty.
#---------------------------------------------------
if len(files) == 0:
    print("\nFolder containing the initial binning result is empty. Please enter a valid path to the folder containing the initial binning result.")
    print("\nExiting prepResult.py...\nBye...!\n")
    sys.exit(1)


# Check if binning result folder contains fasta files.
#---------------------------------------------------
isFasta = False
for myfile in files:
    if myfile.lower().endswith(('.fasta', '.fa', '.fna')):
        isFasta = True
    
if not isFasta:
    print("\nMake sure the folder containing the initial binning result contains fasta files (.fasta, .fa or .fna).")
    print("\nExiting prepResult.py...\nBye...!\n")
    sys.exit(1)


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
    print("Exiting prepResult.py...\n")
    sys.exit(1)


# Format binning results.
#---------------------------------------------------

print("\nFormatting initial binning results")

i = 1
contig_bins = []
bin_ids = []

for bin_file in files:

    if bin_file.lower().endswith(('.fasta', '.fa', '.fna')):

        bin_line = []
        bin_line.append(str(bin_file))
        bin_line.append(str(i))
        bin_ids.append(bin_line)

        for index, record in enumerate(SeqIO.parse(contig_bins_folder+bin_file, "fasta")):
            contig_name = str(record.id)

            line = []
            line.append(str(contig_name))

            except:
                print("\nContig naming does not match with the assembler type provided. Please make sure to provide the correct assembler type.")
                print("\nExiting prepResult.py...\nBye...!\n")
                sys.exit(1)

            line.append(str(i))
            contig_bins.append(line)

        i = i + 1


# Write binning results to output file.
#---------------------------------------------------

print("\nWriting initial binning results to output file")

with open(output_path + prefix + 'initial_contig_bins.csv', mode='w') as contig_bins_file:
    contig_writer = csv.writer(contig_bins_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in contig_bins:
        contig_writer.writerow(row)

with open(output_path + prefix + 'bin_ids.csv', mode='w') as bin_ids_file:
    bin_id_writer = csv.writer(bin_ids_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    for row in bin_ids:
        bin_id_writer.writerow(row)

print("\nFormatted initial binning results can be found at", contig_bins_file.name)
print("\nBin IDs and corresponding names of fasta files can be found at", bin_ids_file.name)


# Exit program
#--------------

print("\nThank you for using prepResult for GraphBin!\n")