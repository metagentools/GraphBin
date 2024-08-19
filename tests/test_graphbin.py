import subprocess

from pathlib import Path

import pytest


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "1.7.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
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


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0

    Parameters
    ----------

    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_graphbin_version():
    """test graphbin version"""
    cmd = "graphbin --version"
    exec_command(cmd)


def test_graphbin_on_spades_dataset(tmp_dir):
    """test graphbin on spades assembly"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_sga_dataset(tmp_dir):
    """test graphbin on sga assembly"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_SGA"
    graph = dir_name / "default-graph.asqg"
    contigs = dir_name / "default-contigs.fa"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler sga --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_megahit_dataset(tmp_dir):
    """test graphbin on megahit assembly"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_MEGAHIT"
    graph = dir_name / "final.gfa"
    contigs = dir_name / "final.contigs.fa"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler megahit --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_flye_dataset(tmp_dir):
    """test graphbin on flye assembly"""
    dir_name = TEST_ROOTDIR / "data" / "1Y3B_Flye"
    graph = dir_name / "assembly_graph.gfa"
    contigs = dir_name / "assembly.fasta"
    paths = dir_name / "assembly_info.txt"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler flye --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_canu_dataset(tmp_dir):
    """test graphbin on canu assembly"""
    dir_name = TEST_ROOTDIR / "data" / "1Y3B_Canu"
    graph = dir_name / "1y3b.contigs.gfa"
    contigs = dir_name / "1y3b.contigs.fasta"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler canu --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_miniasm_dataset(tmp_dir):
    """test graphbin on miniasm assembly"""
    dir_name = TEST_ROOTDIR / "data" / "1Y3B_Miniasm"
    graph = dir_name / "reads.gfa"
    contigs = dir_name / "unitigs.fasta"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler miniasm --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)
