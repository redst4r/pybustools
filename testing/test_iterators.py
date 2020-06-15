import pytest
from pybustools import busio, pybustools
import pathlib
import os


def test_iterate_cells(tmp_path):
    records = [
        busio.Bus_record('ATAT', 'AAA', 10, 20, 1),
        busio.Bus_record('ATAT', 'GGG', 10, 20, 1),
        busio.Bus_record('TAGA', 'TAT', 14, 250, 13),
        busio.Bus_record('TTAT', 'AAA', 13, 206, 12),
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)

    gen = pybustools.iterate_cells_of_busfile(fname)

    # the first record must have two uMIs
    cb1, list1 = next(gen)
    assert cb1 == 'ATAT' and len(list1) == 2

    cb2, list2 = next(gen)
    assert cb2 == 'TAGA' and len(list2) == 1

    cb3, list3 = next(gen)
    assert cb3 == 'TTAT' and len(list2) == 1


def test_iterate_cb_umi(tmp_path):
    records = [
        busio.Bus_record('ATAT', 'AAA', 10, 20, 1),
        busio.Bus_record('ATAT', 'GGG', 10, 20, 1),
        busio.Bus_record('ATAT', 'GGG', 11, 20, 1),
        busio.Bus_record('TAGA', 'TAT', 14, 250, 13),
        busio.Bus_record('TTAT', 'AAA', 13, 206, 12),
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)

    gen = pybustools.iterate_CB_UMI_of_busfile(fname)

    # the first record must have one entry
    cb1, list1 = next(gen)
    assert cb1 == ('ATAT', 'AAA') and len(list1) == 1

    cb2, list2 = next(gen)
    assert cb2 == ('ATAT', 'GGG') and len(list2) == 2

    cb3, list3 = next(gen)
    assert cb3 == ('TAGA', 'TAT') and len(list3) == 1

    cb4, list4 = next(gen)
    assert cb4 == ('TTAT', 'AAA') and len(list4) == 1
