#!/usr/bin/env python3


def getClosestLabelledVertices(graph, node, binned_contigs):
    # Remove labels of ambiguous vertices
    #-------------------------------------

    queu_l = [graph.neighbors(node, mode='ALL')]
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
                temp += graph.neighbors(n, mode='ALL')
                temp = list(set(temp))
            temp2 = []

            for n in temp:
                if n not in visited_l:
                    temp2.append(n)
            if len(temp2) > 0:
                queu_l.append(temp2)
    return labelled
