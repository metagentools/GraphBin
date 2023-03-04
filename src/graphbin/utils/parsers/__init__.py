#!/usr/bin/env python3

import csv
import logging
import re
import sys


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.6.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"


logger = logging.getLogger("GraphBin %s" % __version__)


def get_initial_bin_count(contig_bins_file, delimiter):
    n_bins = 0

    try:
        all_bins_list = []

        with open(contig_bins_file) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=delimiter)
            for row in readCSV:
                all_bins_list.append(row[1])

        bins_list = list(set(all_bins_list))
        bins_list.sort()

        n_bins = len(bins_list)
        logger.info(
            "Number of bins available in the initial binning result: " + str(n_bins)
        )

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the initial binning result file is provided and it is having the correct format."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    return n_bins, bins_list
