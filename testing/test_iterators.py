from pybustools import busio, pybustools
import pytest


def test_return_busrecord(tmp_path):
    """
    make sure the functions return the namedtuple, not just the tuple
    """
    records = [
        busio.Bus_record('ATAT', 'AAA', 10, 20, 1),
        busio.Bus_record('ATAT', 'GGG', 10, 20, 1),
        busio.Bus_record('TAGA', 'TAT', 14, 250, 13),
        busio.Bus_record('TTAT', 'AAA', 13, 206, 12),
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)

    for cb, record_list in pybustools.iterate_cells_of_busfile(fname):
        for r in record_list:
            assert isinstance(r, busio.Bus_record)

    for cbumi, record_list in pybustools.iterate_CB_UMI_of_busfile(fname):
        for r in record_list:
            assert isinstance(r, busio.Bus_record)


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
    assert list1 == records[:2]

    cb2, list2 = next(gen)
    assert cb2 == 'TAGA' and len(list2) == 1
    assert list2 == [records[2]]

    cb3, list3 = next(gen)
    assert cb3 == 'TTAT' and len(list3) == 1
    assert list3 == [records[3]]


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


def test_iterate_cells_raise_unsorted(tmp_path):
    """
    iterate_cells must raise an error of the busfile is unsorted (in terms of CBs)
    """
    records = [
        busio.Bus_record('TTAT', 'AAA', 13, 206, 12),
        busio.Bus_record('ATAT', 'AAA', 10, 20, 1),
        busio.Bus_record('TAGA', 'TAT', 14, 250, 13),
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)

    with pytest.raises(ValueError):
        gen = pybustools.iterate_cells_of_busfile(fname)
        list(gen)


def test_iterate_cells_UMI_raise_unsorted(tmp_path):
    """
    iterate_CB_UMI_of_busfile must raise an error of the busfile is unsorted (in terms of CBs)
    """
    records = [
        busio.Bus_record('ATAT', 'TAT', 14, 250, 13),
        busio.Bus_record('ATAT', 'AAA', 10, 20, 1),
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)

    with pytest.raises(ValueError):
        gen = pybustools.iterate_CB_UMI_of_busfile(fname)
        list(gen)

    """
    also check that it raises when the CB is unsorted
    """
    records = [  # impotant: the UMI should be the same,
        busio.Bus_record('TTAT', 'TAT', 14, 250, 13),
        busio.Bus_record('ATAT', 'TAT', 10, 20, 1),
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)
    with pytest.raises(ValueError):
        gen = pybustools.iterate_CB_UMI_of_busfile(fname)
        list(gen)


def test_records_to_genes():

    ec2gene = {
        0: ['A','B'],
        1: ['B'],
        2: ['A'],
        3: ['C'],
    }

    records = [
        busio.Bus_record('TTAT', 'TAT', 0, 250, 13),  # maps to A,B
        busio.Bus_record('TTAT', 'TAT', 1, 250, 13),  # maps to B
    ]
    assert pybustools.records_to_gene(records, ec2gene) == list('B')

    # incompatible records
    records = [
        busio.Bus_record('TTAT', 'TAT', 0, 250, 13),  # maps to A,B
        busio.Bus_record('TTAT', 'TAT', 3, 250, 13),  # maps to C
    ]
    assert pybustools.records_to_gene(records, ec2gene) == list()
