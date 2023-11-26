from pybustools import busio, subsampling

def test_subsampling(tmp_path):

    # creating a total of 3 UMIs, nbut 10 counts
    records = [
        # CB UMI EC COUNT FLAG
        busio.Bus_record(0, 0, 1, 3, 1),
        busio.Bus_record(1, 0, 2, 3, 12),
        busio.Bus_record(2, 0, 3, 4, 13)
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=1, umi_length=2)

    fname_out = tmp_path / 'sub.bus'

    subsampling.subsample_busfile(fname, fname_out, fraction=0.5)

    # check the number of reads in the subsampled version
    nreads, nmol = subsampling.get_number_of_reads_and_molecules(fname_out)
    assert nreads == 5
