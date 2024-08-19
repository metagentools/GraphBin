import subprocess

from pathlib import Path

import pytest

from graphbin.bidirectionalmap.bidirectionalmap import *


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "1.7.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


def test_bidirectional_map_insert():
    try:
        my_map = BidirectionalMap()

        my_map[0] = 56
        my_map[7] = 56

    except BidirectionalError:
        print("BidirectionalError")


def test_bidirectional_map_delete():
    my_map = BidirectionalMap()

    my_map[0] = 56
    my_map[2] = 57
    my_map._del_item(2)
