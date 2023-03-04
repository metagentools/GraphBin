#!/usr/bin/env python3

import logging
import sys

from graphbin.utils.labelpropagation.labelprop import LabelProp

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.6.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Production"

logger = logging.getLogger("GraphBin %s" % __version__)

MIN_BIN_COUNT = 10


def getClosestLabelledVertices(graph, node, binned_contigs):
    # Remove labels of ambiguous vertices
    # -------------------------------------

    queu_l = [graph.neighbors(node, mode="ALL")]
    visited_l = [node]
    labelled = []

    while len(queu_l) > 0:
        active_level = queu_l.pop(0)
        is_finish = False
        visited_l += active_level

        for n in active_level:
            if n in binned_contigs:
                is_finish = True
                labelled.append(n)
        if is_finish:
            return labelled
        else:
            temp = []
            for n in active_level:
                temp += graph.neighbors(n, mode="ALL")
                temp = list(set(temp))
            temp2 = []

            for n in temp:
                if n not in visited_l:
                    temp2.append(n)
            if len(temp2) > 0:
                queu_l.append(temp2)
    return labelled


def graphbin_main(
    n_bins, bins, bins_list, assembly_graph, node_count, diff_threshold, max_iteration
):
    logger.info("Determining ambiguous vertices")

    remove_by_bin = {}

    remove_labels = []

    neighbours_have_same_label_list = []

    for b in range(n_bins):
        for i in bins[b]:
            my_bin = b

            # Get set of closest labelled vertices with distance = 1
            closest_neighbours = assembly_graph.neighbors(i, mode="all")

            # Determine whether all the closest labelled vertices have the same label as its own
            neighbours_have_same_label = True

            neighbours_binned = False

            for neighbour in closest_neighbours:
                for k in range(n_bins):
                    if neighbour in bins[k]:
                        neighbours_binned = True
                        if k != my_bin:
                            neighbours_have_same_label = False
                            break

            if not neighbours_have_same_label:
                if my_bin in remove_by_bin:
                    if len(bins[my_bin]) - len(remove_by_bin[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.append(i)
                        remove_by_bin[my_bin].append(i)
                else:
                    if len(bins[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.append(i)
                        remove_by_bin[my_bin] = [i]

            elif neighbours_binned:
                neighbours_have_same_label_list.append(i)

    for i in remove_labels:
        for n in range(n_bins):
            if i in bins[n]:
                bins[n].remove(i)

    # Further remove labels of ambiguous vertices
    binned_contigs = []

    for n in range(n_bins):
        binned_contigs = sorted(binned_contigs + bins[n])

    for b in range(n_bins):
        for i in bins[b]:
            if i not in neighbours_have_same_label_list:
                my_bin = b

                # Get set of closest labelled vertices
                closest_neighbours = getClosestLabelledVertices(
                    assembly_graph, i, binned_contigs
                )

                if len(closest_neighbours) > 0:
                    # Determine whether all the closest labelled vertices have the same label as its own
                    neighbours_have_same_label = True

                    for neighbour in closest_neighbours:
                        for k in range(n_bins):
                            if neighbour in bins[k]:
                                if k != my_bin:
                                    neighbours_have_same_label = False
                                    break

                    if not neighbours_have_same_label and i not in remove_labels:
                        if my_bin in remove_by_bin:
                            if (
                                len(bins[my_bin]) - len(remove_by_bin[my_bin])
                                >= MIN_BIN_COUNT
                            ):
                                remove_labels.append(i)
                                remove_by_bin[my_bin].append(i)
                        else:
                            if len(bins[my_bin]) >= MIN_BIN_COUNT:
                                remove_labels.append(i)
                                remove_by_bin[my_bin] = [i]

    logger.info("Removing labels of ambiguous vertices")

    # Remove labels of ambiguous vertices
    for i in remove_labels:
        for n in range(n_bins):
            if i in bins[n]:
                bins[n].remove(i)

    logger.info("Obtaining the refined binning result")

    # Get vertices which are not isolated and not in components without any labels
    # -----------------------------------------------------------------------------

    logger.info(
        "Deteremining vertices which are not isolated and not in components without any labels"
    )

    non_isolated = []

    for i in range(node_count):
        if i not in non_isolated and i in binned_contigs:
            component = []
            component.append(i)
            length = len(component)
            neighbours = assembly_graph.neighbors(i, mode="all")

            for neighbor in neighbours:
                if neighbor not in component:
                    component.append(neighbor)

            component = list(set(component))

            while length != len(component):
                length = len(component)

                for j in component:
                    neighbours = assembly_graph.neighbors(j, mode="all")

                    for neighbor in neighbours:
                        if neighbor not in component:
                            component.append(neighbor)

            labelled = False
            for j in component:
                if j in binned_contigs:
                    labelled = True
                    break

            if labelled:
                for j in component:
                    if j not in non_isolated:
                        non_isolated.append(j)

    logger.info("Number of non-isolated contigs: " + str(len(non_isolated)))

    # Run label propagation
    # -----------------------

    data = []

    for contig in range(node_count):
        # Consider vertices that are not isolated

        if contig in non_isolated:
            line = []
            line.append(contig)

            assigned = False

            for i in range(n_bins):
                if contig in bins[i]:
                    line.append(i + 1)
                    assigned = True

            if not assigned:
                line.append(0)

            neighbours = assembly_graph.neighbors(contig, mode="all")

            neighs = []

            for neighbour in neighbours:
                n = []
                n.append(neighbour)
                n.append(1.0)
                neighs.append(n)

            line.append(neighs)

            data.append(line)

    # Check if initial binning result consists of contigs belonging to multiple bins

    multiple_bins = False

    for item in data:
        if type(item[1]) is int and type(item[2]) is int:
            multiple_bins = True
            break

    if multiple_bins:
        logger.error(
            "Initial binning result consists of contigs belonging to multiple bins. Please make sure that each contig in the initial binning result belongs to only one bin."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    # Label propagation

    lp = LabelProp()

    lp.load_data_from_mem(data)

    logger.info(
        "Starting label propagation with eps="
        + str(diff_threshold)
        + " and max_iteration="
        + str(max_iteration)
    )

    ans = lp.run(diff_threshold, max_iteration, show_log=True, clean_result=False)

    logger.info("Obtaining Label Propagation result")

    for l in ans:
        for i in range(n_bins):
            if l[1] == i + 1 and l[0] not in bins[i]:
                bins[i].append(l[0])

    # Remove labels of ambiguous vertices
    # -------------------------------------

    logger.info("Determining ambiguous vertices")

    remove_by_bin = {}

    remove_labels = []

    for b in range(n_bins):
        for i in bins[b]:
            my_bin = b

            closest_neighbours = assembly_graph.neighbors(i, mode="all")

            # Determine whether all the closest labelled vertices have the same label as its own
            neighbours_have_same_label = True

            for neighbour in closest_neighbours:
                for k in range(n_bins):
                    if neighbour in bins[k]:
                        if k != my_bin:
                            neighbours_have_same_label = False
                            break

            if not neighbours_have_same_label:
                if my_bin in remove_by_bin:
                    if len(bins[my_bin]) - len(remove_by_bin[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.append(i)
                        remove_by_bin[my_bin].append(i)
                else:
                    if len(bins[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.append(i)
                        remove_by_bin[my_bin] = [i]

    logger.info("Removing labels of ambiguous vertices")

    # Remove labels of ambiguous vertices
    for i in remove_labels:
        for n in range(n_bins):
            if i in bins[n]:
                bins[n].remove(i)

    logger.info("Obtaining the Final Refined Binning result")

    final_bins = {}

    for i in range(n_bins):
        for contig in bins[i]:
            final_bins[contig] = bins_list[i]

    return final_bins, remove_labels, non_isolated
