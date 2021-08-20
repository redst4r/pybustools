import os
import pytest
from pybustools.pybustools import _Bus
from pybustools import busio
from pybustools import butterfly
from pybustools.butterfly import CUHistogram, saturation_curve
import numpy as np


def test_CUHistogram_equals():
    b1 = CUHistogram({1:1, 2:1})
    b2 = CUHistogram({1:1, 2:1})
    b3 = CUHistogram({1:1, 2:2})

    assert b1 == b2 and b1 != b3 and b2 != b3

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
        1: butterfly.CUHistogram({1:1, 2:1}),
        2: butterfly.CUHistogram({1:2}),
    }

def _assert_dict_eq(d1, d2):
    assert d1.keys() == d2.keys()
    for k in d1.keys():
        np.testing.assert_almost_equal(d1[k], d2[k])


def test_binomial_downsample():
    CUhist = CUHistogram({
        1: 10,  # 10 molecules with single read
        })
    assert butterfly.binomial_downsample(CUhist, 0.5) == CUHistogram({0:5, 1:5})

    # =======================================================================
    CUhist = CUHistogram({
        2: 16,  # 16 molecules with two reads
        })
    # assert butterfly.binomial_downsample(CUhist, 0.5) ==
    _assert_dict_eq(butterfly.binomial_downsample(CUhist, 0.5).histogram,
                    CUHistogram({0:4, 1:8, 2:4}).histogram)

    # =======================================================================
    CUhist = CUHistogram({
        1: 4,
        2: 16,  # 15 molecules with two reads
        4: 20,  # 20 molecules with four reads
    })

    _assert_dict_eq(butterfly.binomial_downsample(CUhist, 0.5).histogram,
                    CUHistogram({0: 7.25, 1: 15, 2: 11.5, 3: 5, 4: 1.25}).histogram)


def test_binomial_downsample_all_genes():
    CUhist1 = CUHistogram({
        1: 10,  # 10 molecules with single read
        })
    CUhist2 = CUHistogram({
        2: 16,  # 16 molecules with two reads
        })
    d = {
        'gene1': CUhist1,
        'gene2': CUhist2
    }
    R = butterfly.binomial_downsample_all_genes(d, 0.5)

    R_expected = {g: butterfly.binomial_downsample(d[g], 0.5) for g in d.keys()}

    assert R == R_expected

def test_binomial_downsample_factors():
    # a simple amplificifaation that is strictly 1:1
    old_hist = {'gene1': CUHistogram({1:10})}
    new_hist = {'gene1': butterfly.binomial_downsample(old_hist['gene1'], 0.5)}
    assert butterfly.binomial_downsample_factors(old_hist, new_hist, CPM='reads') == {'gene1': 1}

    # a 2:1 amplificifaation
    old_hist = {'gene1': CUHistogram({2:10})}
    new_hist = {'gene1': butterfly.binomial_downsample(old_hist['gene1'], 0.5)}
    assert butterfly.binomial_downsample_factors(old_hist, new_hist, CPM='reads') == {'gene1': 1.5}


    # umi norm
    old_hist = {'gene1': CUHistogram({1:10})}
    new_hist = {'gene1': butterfly.binomial_downsample(old_hist['gene1'], 0.5)}
    assert butterfly.binomial_downsample_factors(old_hist, new_hist, CPM='umis') == {'gene1': 1}

    old_hist = {'gene1': CUHistogram({2:10})}
    new_hist = {'gene1': butterfly.binomial_downsample(old_hist['gene1'], 0.5)}
    assert butterfly.binomial_downsample_factors(old_hist, new_hist, CPM='umis') == {'gene1': 1.}


def test_aggregate():

    c1 = CUHistogram({1:10})
    c2 = CUHistogram({1:5, 10: 1})

    assert butterfly.aggregate_CUs({'g1': c1, 'g2':c2})  == CUHistogram({1:15, 10: 1})


def test_saturation():
    c1 = CUHistogram({1:10, 2:10})
    df = saturation_curve(c1, bins=20)

    c1 = CUHistogram({1:0, 2:0})
    df = saturation_curve(c1, bins=20)

    # c1 = CUHistogram({})
    # df = saturation_curve(c1, bins=20)


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
