import subprocess
from pathlib import Path

import pytest

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "1.6.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "vijini.mallawaarachchi@anu.edu.au"
__status__ = "Development"


TEST_ROOTDIR = Path(__file__).parent
EXEC_ROOTDIR = Path(__file__).parent.parent


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


def exec_wrong_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode == 1:
        return out.decode("utf8") if out is not None else None


def test_graphbin_assembler(tmp_dir):
    """test graphbin on wrong assembler"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spa --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_wrong_command(cmd)

def test_graphbin_graph(tmp_dir):
    """test graphbin on wrong graph file"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffold.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_wrong_command(cmd)

def test_graphbin_contigs(tmp_dir):
    """test graphbin on wrong contigs file"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contig.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_wrong_command(cmd)

def test_graphbin_paths(tmp_dir):
    """test graphbin on wrong paths file"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contig.paths"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_wrong_command(cmd)

def test_graphbin_binned(tmp_dir):
    """test graphbin on wrong binned result file"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_wrong_command(cmd)

def test_graphbin_spades_path(tmp_dir):
    """test graphbin spades without paths"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_wrong_command(cmd)

def test_graphbin_no_contigs(tmp_dir):
    """test graphbin with no contigs"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_wrong_command(cmd)

def test_graphbin_prefix(tmp_dir):
    """test graphbin with prefix"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir} --prefix test"
    exec_wrong_command(cmd)

def test_graphbin_delimiter(tmp_dir):
    """test graphbin with wrong delimiter"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    delimiter = "."
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir} --delimiter {delimiter}"
    exec_wrong_command(cmd)

def test_graphbin_max_iteration(tmp_dir):
    """test graphbin with wrong max_iteration"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    max_iteration = -10
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir} --max_iteration {max_iteration}"
    exec_wrong_command(cmd)

def test_graphbin_diff_threshold(tmp_dir):
    """test graphbin with wrong diff_threshold"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    diff_threshold = -10
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir} --max_iteration {diff_threshold}"
    exec_wrong_command(cmd)