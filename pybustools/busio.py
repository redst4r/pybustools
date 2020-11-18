import struct
import gmpy2
import collections
from functools import partial
Bus_record = collections.namedtuple('BUSRecord', 'CB UMI EC COUNT FLAG')


def _encode_ACGT_to_int(seq):
    intstr = seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3')
    the_int = int(intstr, base=4)
    return the_int


def _decode_int_to_ACGT(the_int, seq_len):
    """
    bustools encodes sequences as integers, this function decodes the ints
    into their sequence
    """

    assert the_int >= 0
    # the follwing line is a major performance HOG:
    # seq_decode = np.base_repr(the_int, 4) # seq_decode is a str: '10012031...'
    seq_decode = gmpy2.digits(the_int, 4)

    # the coding does not recognize leading 0 (adenines)
    # hence we pad 0 until we have the correct number of bases
    seq_decode_pad = seq_decode.zfill(seq_len)

    assert len(seq_decode_pad) == seq_len, f'{len(seq_decode_pad)} vs {seq_len}'

    seq_str = seq_decode_pad.replace('0', 'A').replace('1', 'C').replace('2', 'G').replace('3', 'T')
    return seq_str


def get_header_info(fname):
    """
    retrieve the bus header info (version, CB/UMI length, and the freeform header)
    """
    with open(fname, 'rb') as fh:
        header = fh.read(20)  # header is 20bytes
        version, cb_len, umi_len, tlen = parse_header_info(header)

    return version, cb_len, umi_len, tlen


def parse_header_info(header_bytes):
    # read the header
    # Magic (4bytes)
    # version int
    # CB len
    # umi len
    # freetext header len
    assert len(header_bytes) == 20
    magic, version, cb_len, umi_len, tlen = struct.unpack('4sIIII', header_bytes)
    assert magic == b'BUS\x00', "MAGIC doesnt match, wrong filetype??"

    return version, cb_len, umi_len, tlen


def read_binary_bus(fname, decode_seq:bool=True, buffersize=1000):
    """
    iterating over a binary busfile, yielding (CB,UMI,EC,Counts,Flag)

    decode_seq: CB/UMI are encoded as ints. decode them while going through the
    file (takes considerable time, so if the actual sequence is not of interest,
    dealing with the ints is alot faster)

    buffersize: how many busrecords to load in a single IO operation. ~1000 seems to be a
    sweeyt spot, having decent speedups. increasing it more doesnt do much
    """

    with open(fname, 'rb') as fh:
        header = fh.read(20)  # header is 20bytes

        version, cb_len, umi_len, tlen = parse_header_info(header)
        # read the free header
        free_header = struct.unpack(f'{tlen}s', fh.read(tlen))

        print(f'Bustools {version}, CB length {cb_len}, UMI length {umi_len}')
        print(f'Free header  {free_header}')

        BUS_ENTRY_SIZE = 32  # each entry is 32 bytes!!
        unpack_str = 'QQiIII'
        for bus_chunk in iter(partial(fh.read, BUS_ENTRY_SIZE * buffersize), b''):
            # iterate over each single entry in the loaded chunk.
            # mote that struct.iter_unpack neatly turns the chunk into an
            # iterator
            for cb, umi, ec, count, flags, pad in struct.iter_unpack(unpack_str, bus_chunk):
                assert pad == 0
                if decode_seq:
                    cb = _decode_int_to_ACGT(cb, cb_len)
                    umi = _decode_int_to_ACGT(umi, umi_len)

                yield Bus_record(CB=cb, UMI=umi, EC=ec, COUNT=count, FLAG=flags)


def read_matrix_ec(fname):
    """
    parsing the Equivalence Class file, which maps EC->list(int_transcript_ids)
    """
    D = {}
    with open(fname, 'r') as fh:
        for line in fh:
            ec, transcript_list = line.split()
            ec = int(ec)
            transcripts = [int(_) for _ in transcript_list.split(',')]
            D[ec] = transcripts
    return D


def read_transcripts(fname):
    """
    parsing the transcripts file,
    mapping int_transcript_ids -> ensemble_transcript_id

    int_transcript_ids is line number (0-based!)
    """
    D = {}
    with open(fname, 'r') as fh:
        for i, line in enumerate(fh):
            D[i] = line.strip()
    return D


def write_busfile(fname, bus_records:list, cb_length, umi_length):
    """
    write bus_records to disk.
    automatically checks if the CB/UMIs are given as ints or strings and encodes
    them if needed
    """
    with open(fname, 'wb') as fh:
        fre_header = 'my_bus'
        tlen = len(fre_header)
        header_bytes = struct.pack('4sIIII', b'BUS\x00', 1, cb_length, umi_length, tlen)
        free_header_bytes = struct.pack(f'{tlen}s', fre_header.encode())

        fh.write(header_bytes)
        fh.write(free_header_bytes)

        for record in bus_records:
            cb, umi, ec, count, flags = record
            # if the records have strings for CB/UMI, encode them into ints
            if isinstance(cb, str) and isinstance(umi, str):
                assert len(cb) == cb_length
                assert len(umi) == umi_length
                cb = _encode_ACGT_to_int(cb)
                umi = _encode_ACGT_to_int(umi)
            elif isinstance(cb, int) and isinstance(umi, int):
                pass
            else:
                raise ValueError('CB/UMI must be both str or both int')

            bytes_record = struct.pack('QQiIIxxxx', cb, umi, ec, count, flags)
            fh.write(bytes_record)


if __name__ == '__main__':
    import tqdm
    folder = '/run/media/michi/42506642-b470-4238-be14-bb0c303b3682/cruk/bustools_testing/E14B_mgi_bus_sorted'
    bin_bus = f'{folder}/output.corrected.sort.bus'
    I_bin = read_binary_bus(bin_bus)
    for _ in tqdm.tqdm(I_bin):
        pass
