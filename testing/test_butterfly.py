import os
import pytest
from pybustools.pybustools import _Bus
from pybustools import busio
from pybustools import butterfly


def test_make_ec_histograms(tmp_path, ec_matrix_file, transcript_file):
    records = [
        # Gene 1
        busio.Bus_record('ATAT', 'AAA', 1, 2, 1),
        busio.Bus_record('ATAT', 'GGG', 1, 1, 1),
        # Gene 2
        busio.Bus_record('TAGA', 'TAT', 2, 1, 13),
        busio.Bus_record('TAGA', 'TTT', 2, 1, 13),

        # Gene 3: multimapped, should be discarded
        busio.Bus_record('TTAT', 'AAA', 3, 1, 12),
        busio.Bus_record('TTAT', 'AAA', 4, 1, 12),
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)

    # gen = iterate_cells_of_busfile(fname)

    B = _Bus(fname, ec_name=ec_matrix_file, transcript_name=transcript_file)
    h = butterfly.make_ec_histograms(B, collapse_EC=False)
    assert h == {
        1: {1:1, 2:1},
        2: {1:2},
    }


def test_binomial_downsample():
    CUhist = {
        1: 10,  # 10 molecules with single read
        }
    assert butterfly.binomial_downsample(CUhist, 0.5) == {0:5, 1:5}

    CUhist = {
        2: 16,  # 15 molecules with two reads
        }
    assert butterfly.binomial_downsample(CUhist, 0.5) == {0:4, 1:8, 2:4}


def test_binomial_downsample_factors():
    # a simple amplificifaation that is strictly 1:1
    old_hist = {'gene1': {1:10}}
    new_hist = {'gene1': butterfly.binomial_downsample(old_hist['gene1'], 0.5)}
    assert butterfly.binomial_downsample_factors(old_hist, new_hist, CPM='reads') == {'gene1': 1}

    # a 2:1 amplificifaation
    old_hist = {'gene1': {2:10}}
    new_hist = {'gene1': butterfly.binomial_downsample(old_hist['gene1'], 0.5)}
    assert butterfly.binomial_downsample_factors(old_hist, new_hist, CPM='reads') == {'gene1': 1.5}


@pytest.fixture
def ec_matrix_file():
    "creates an ec_file with 10 entries"
    import tempfile
    fname = tempfile.mktemp()

    ECS = [
        [1], # EC1 maps to transcript 1
        [2], # EC2 maps to transcript 2
        [3], #3
        [4], #4
        [5], #5
        [6], #6
        [7], #7
        [8], #8
        [1,2], # EC9 maps to transcripts 1,2
        [1,3], # E10 maps to transcripts 1,3
    ]

    with open(fname, 'w') as fh:
        for i in range(10):
            fh.write(f'{i} 1,2,3,4\n')
    yield fname
    os.remove(fname)


@pytest.fixture
def transcript_file():
    "creates an transcript_file with 10 entries"
    import tempfile
    fname = tempfile.mktemp()
    with open(fname, 'w') as fh:
        for i in range(10):
            fh.write(f'ENST00000000000{i}.1\n')
    yield fname
    os.remove(fname)
