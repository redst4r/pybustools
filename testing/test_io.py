import pytest
from pybustools import busio, busio_old
import pathlib
import os


def test_encode():
    assert busio._encode_ACGT_to_int('A') == 0
    assert busio._encode_ACGT_to_int('AA') == 0

    assert busio._encode_ACGT_to_int('GCCA') == 148

    with pytest.raises(ValueError):
        busio._encode_ACGT_to_int('B')


def test_decode_int_to_ACGT():

    #  base order
    assert busio._decode_int_to_ACGT(0, seq_len=1) == 'A'
    assert busio._decode_int_to_ACGT(1, seq_len=1) == 'C'
    assert busio._decode_int_to_ACGT(2, seq_len=1) == 'G'
    assert busio._decode_int_to_ACGT(3, seq_len=1) == 'T'

    # padding leading A's
    assert busio._decode_int_to_ACGT(0, seq_len=3) == 'AAA'
    assert busio._decode_int_to_ACGT(1, seq_len=3) == 'AAC'
    assert busio._decode_int_to_ACGT(2, seq_len=3) == 'AAG'
    assert busio._decode_int_to_ACGT(3, seq_len=3) == 'AAT'

    assert busio._decode_int_to_ACGT(4, seq_len=2) == 'CA'  # note that AA is still 0
    assert busio._decode_int_to_ACGT(5, seq_len=2) == 'CC'
    assert busio._decode_int_to_ACGT(6, seq_len=2) == 'CG'
    assert busio._decode_int_to_ACGT(7, seq_len=2) == 'CT'

    assert busio._decode_int_to_ACGT(148, seq_len=4) == 'GCCA'

    # make sure to raise an error when the decoded string is actually longer
    # then requested (since its probably a bug in the code calling _decode_int_to_ACGT)
    with pytest.raises(AssertionError):
        busio._decode_int_to_ACGT(148, seq_len=2)
    with pytest.raises(AssertionError):
        busio._decode_int_to_ACGT(-1, seq_len=1)


def test_encode_decode_ordering():
    """
    ensure that the alphabetical ordering of sequences is that same as he ordering
    of their encodings
    """
    N = 8
    decoded = [busio._decode_int_to_ACGT(encoding, N) for encoding in range(4**N)]
    # since the encoding is sorted, so should be the decoding
    assert decoded == sorted(decoded)


def test_encode_decode():
    i = 4723452
    assert busio._encode_ACGT_to_int(busio._decode_int_to_ACGT(i, seq_len=13)) == i

    seq = 'AAAATTTTGGGGCCCC'
    assert busio._decode_int_to_ACGT(busio._encode_ACGT_to_int(seq), len(seq)) == seq


def test_get_header(tmp_path):

    records = [
        busio.Bus_record(0, 0, 10, 20, 1),
        busio.Bus_record(1, 0, 13, 206, 12),
        busio.Bus_record(2, 0, 14, 250, 13)
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=12, umi_length=5)

    _, cb, umi, _ = busio.get_header_info(fname)
    assert cb == 12 and umi == 5


def test_read_write(tmp_path):

    records = [
        busio.Bus_record(0, 0, 10, 20, 1),
        busio.Bus_record(1, 0, 13, 206, 12),
        busio.Bus_record(2, 0, 14, 250, 13)
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=1, umi_length=2)

    # check that file got created
    assert pathlib.Path(fname).exists()

    # check that read/write are inverses of each other
    # buffersize is on purpose smaller then len(records)
    new_records = list(busio.read_binary_bus(fname, decode_seq=False, buffersize=2))
    assert new_records == records

    # buffersize larger
    new_records = list(busio.read_binary_bus(fname, decode_seq=False, buffersize=20))
    assert new_records == records

    # check the decode_Seq works:
    record = next(busio.read_binary_bus(fname, decode_seq=False, buffersize=2))
    assert isinstance(record.CB, int) and isinstance(record.UMI, int)
    record = next(busio.read_binary_bus(fname, decode_seq=True, buffersize=2))
    assert isinstance(record.CB, str) and isinstance(record.UMI, str)


def test_read_write_str(tmp_path):
    """
    write records with strings instead of ints for CB/UMI
    """
    records = [
        busio.Bus_record('ATAT', 'GGG', 10, 20, 1),
        busio.Bus_record('TTAT', 'AAA', 13, 206, 12),
        busio.Bus_record('TAGA', 'TAT', 14, 250, 13)
    ]
    fname = tmp_path / 'some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)
    # read and compare to originl
    new_records = list(busio.read_binary_bus(fname, decode_seq=True, buffersize=2))
    assert new_records == records


def test_write_check_cb_umi_length(tmp_path):
    """
    make sure an exception is thrown if the records dont match the given UMI/CB length
    """
    records = [
        busio.Bus_record('ATAT', 'GGG', 10, 20, 1),
        busio.Bus_record('TTA', 'AAA', 13, 206, 12),
    ]
    fname = tmp_path / 'some.bus'
    with pytest.raises(AssertionError):
        busio.write_busfile(fname, records, cb_length=4, umi_length=3)


def test_batch():
    L = list(range(10))
    B = busio_old.batch(L, n=4)

    result = []
    for b in B:
        result.extend(b)
    assert result == L


@pytest.fixture
def ec_matrix_file():
    "creates an ec_file with 10 entries"
    import tempfile
    fname = tempfile.mktemp()
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


def test_read_matrix_ec(ec_matrix_file):
    ec_dict = busio.read_matrix_ec(ec_matrix_file)
    assert len(ec_dict) == 10
    assert ec_dict[0] == [1, 2, 3, 4]


def test_read_transcripts(transcript_file):
    t_dict = busio.read_transcripts(transcript_file)
    assert len(t_dict) == 10
    assert t_dict[0] == 'ENST000000000000.1'
