"""
storage for some older/slower busio functions
in particular, the old `read_binary_bus` used a custom function to load a chunk
of the binary and iterate over that chunk.
this got replaced by struct.iter_unpack, which is about 2x faster
"""
from functools import partial
import struct
from pybustools.busio import _decode_int_to_ACGT, Bus_record


def batch(iterable, n):
    """
    breaking an iterable into chunks of size n
    """
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def _iterate_busrecords_in_chunks(fh, n_chunks):
    """
    instead of loading one bus-record in an IO operation,
    we load a bunch of them and iterate over the individaul reacords later
    """
    BUS_ENTRY_SIZE = 32  # each entry is 32 bytes!!

    """
    the outer loop pulls 32bytes * N bytes out out the file, corresponding to N
    different BUS-records
    the inner loop goes over the individual single BUS-records
    """

    """
    we could do a
    while True:
        bus = fh.read(32)
        if bus == b'':
            break

    but the iter(..., stop_element) is alot nicer
    """
    for bus_entry_chunk in iter(partial(fh.read, BUS_ENTRY_SIZE * n_chunks), b''):
        for bus_entry in batch(bus_entry_chunk, BUS_ENTRY_SIZE):
            yield bus_entry


def read_binary_bus(fname, decode_seq:bool=True, buffersize=1000):
    """
    iterating over a binary busfile,
    yielding (CB,UMI,EC,Counts,Flag)

    decode_seq: CB/UMI are encoded as ints. decode them while going through the
    file (takes considerable time, so if the actual sequence is not of interest,
    dealing with the ints is alot faster)

    buffersize: how many busrecords to load in a single IO operation. ~10 seems to be a
    sweeyt spot, having decent speedups. increasing it more doesnt do much
    """

    with open(fname, 'rb') as fh:
        # read the header
        # Magic (4bytes)
        # version int
        # CB len
        # umi len
        # freetext header len
        header = fh.read(20)  # header is 20bytes
        magic, version, cb_len, umi_len, tlen = struct.unpack('4sIIII', header)
        assert magic == b'BUS\x00', "MAGIC doesnt match, wrong filetype??"

        # read the free header
        free_header = struct.unpack(f'{tlen}s', fh.read(tlen))

        print(f'Bustools {version}, CB length {cb_len}, UMI length {umi_len}')
        print(f'Free header  {free_header}')

        # read the bus entries and decode these bytes
        unpack_str = 'QQiIII'
        for bus_entry in _iterate_busrecords_in_chunks(fh, n_chunks=buffersize):
            cb, umi, ec, count, flags, pad = struct.unpack(unpack_str, bus_entry)
            assert pad == 0
            if decode_seq:
                cb = _decode_int_to_ACGT(cb, cb_len)
                umi = _decode_int_to_ACGT(umi, umi_len)

            record = Bus_record(CB=cb, UMI=umi, EC=ec, COUNT=count, FLAG=flags)
            yield record


def read_text_bus(fname):
    """
    iterating over a plaintext busfile,
    yielding CB/UMI/EC/Counts/Flag
    """
    with open(fname, 'r') as fh:
        for line in fh:
            cb, umi, ec, count = line.split()
            ec = int(ec)
            count = int(count)
            flag = 0  # TODO Flag fixed!!
            record = Bus_record(CB=cb, UMI=umi, EC=ec, COUNT=count, FLAG=flag)
            yield record
