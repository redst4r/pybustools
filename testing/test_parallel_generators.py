from pybustools.parallel_generators import ParallelCellGenerator, ParallelCellUMIGenerator, random_buslist
from pybustools import busio
from pybustools.pybustools import iterate_bus_cells_multiple, iterate_bus_cells_umi_multiple


def test_ParallelCellGenerator():
    """
    assert that it returns the same iterator as the serial version
    """
    records1 = random_buslist(1000, 4, 4, 20)
    records2 = random_buslist(1000, 4, 4, 20)
    records3 = random_buslist(1000, 4, 4, 20)

    import tempfile
    fname1 = tempfile.mktemp()
    fname2 = tempfile.mktemp()
    fname3 = tempfile.mktemp()

    busio.write_busfile(fname1, records1, cb_length=4, umi_length=4)
    busio.write_busfile(fname2, records2, cb_length=4, umi_length=4)
    busio.write_busfile(fname3, records3, cb_length=4, umi_length=4)

    pgen = ParallelCellGenerator({'sample1': fname1, 'sample2': fname2, 'sample3': fname3}, decode_seq=True, queue_size=10)
    pgen.start_queues()
    parallel_results = {cb: info for cb, info in pgen.iterate()}

    serial_results = {cb: info for cb, info in iterate_bus_cells_multiple(['sample1', 'sample2', 'sample3'], [fname1, fname2, fname3])}

    assert parallel_results == serial_results


def test_ParallelCellUMIGenerator():
    """
    assert that it returns the same iterator as the serial version
    """
    records1 = random_buslist(1000, 4, 4, 20)
    records2 = random_buslist(1000, 4, 4, 20)
    records3 = random_buslist(1000, 4, 4, 20)

    import tempfile
    fname1 = tempfile.mktemp()
    fname2 = tempfile.mktemp()
    fname3 = tempfile.mktemp()

    busio.write_busfile(fname1, records1, cb_length=4, umi_length=4)
    busio.write_busfile(fname2, records2, cb_length=4, umi_length=4)
    busio.write_busfile(fname3, records3, cb_length=4, umi_length=4)

    pgen = ParallelCellUMIGenerator({'sample1': fname1, 'sample2': fname2, 'sample3': fname3}, decode_seq=True, queue_size=10)
    pgen.start_queues()
    parallel_results = {cbumi: info for cbumi, info in pgen.iterate()}

    serial_results = {cbumi: info for cbumi, info in iterate_bus_cells_umi_multiple(['sample1', 'sample2', 'sample3'], [fname1, fname2, fname3])}

    assert parallel_results == serial_results
