import struct
from functools import partial
import gmpy2

def _decode_int_to_ACGT(the_int, seq_len):
    """
    bustools encodes sequences as integers, this function decodes the ints
    into their sequence
    """
    # the follwing line is a major performance HOG:
    # seq_decode = np.base_repr(the_int, 4) # seq_decode is a str: '10012031...'
    seq_decode = gmpy2.digits(the_int, 4)

    # the coding does not recognize leading 0 (adenines)
    # hence we pad 0 until we have the correct number of bases
    seq_decode_pad = seq_decode.zfill(seq_len)
    seq_str = seq_decode_pad.replace('0', 'A').replace('1', 'C').replace('2', 'G').replace('3', 'T')
    return seq_str


def read_binary_bus(fname):
    """
    iterating over a binary busfile,
    yielding (CB,UMI,EC,Counts,Flag)
    """

    with open(fname, 'rb') as fh:
        # read the header
        # Magic (4bytes)
        # version int
        # CB len
        # umi len
        # freetext header len
        header = fh.read(20)
        magic, version, cb_len, umi_len, tlen = struct.unpack('4sIIII', header)
        assert magic == b'BUS\x00', "MAGIC doesnt match, wrong filetype??"

        # read the free header
        free_header = struct.unpack(f'{tlen}s', fh.read(tlen))

        print(f'Bustools {version}, CB length {cb_len}, UMI length {umi_len}')
        print(f'Free header  {free_header}')

        # read the bus entries
        """
        we could do a
        while True:
            bus = fh.read(32)
            if bus == b'':
                break

        but the iter(..., stop_element) is alot nicer
        """
        for bus_entry in iter(partial(fh.read, 32), b''):
            cb, umi, ec, count, flags, pad = struct.unpack('QQiIII', bus_entry)
            assert pad == 0
            cb_str = _decode_int_to_ACGT(cb, cb_len)
            umi_str = _decode_int_to_ACGT(umi, umi_len)
            yield cb_str, umi_str, ec, count, flags


def read_text_bus(fname):
    """
    iterating over a plaintext busfile,
    yielding CB/UMI/EC/Counts/Flag
    """
    with open(fname, 'r') as fh:
        for line in fh:
            cb, umi, ec, count, flag = read_text_bus_entry(line)
            yield cb, umi, ec, count, flag


def read_text_bus_entry(line):
    """
    read a plaintext bus
    """
    cb, umi, ec, count = line.split()
    ec = int(ec)
    count = int(count)
    flag = 0  # TODO Flag fixed!!
    return cb, umi, ec, count, flag


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
            D[i] = line
    return D


if __name__ == '__main__':

    ## some speedtest
    folder = '/run/media/michi/42506642-b470-4238-be14-bb0c303b3682/cruk/bustools_testing/E14B_mgi_bus_sorted'
    text_bus = f'{folder}/output.corrected.sort.txt.bus'
    bin_bus = f'{folder}/output.corrected.sort.bus'

    I_plain = read_text_bus(text_bus)
    I_bin = read_binary_bus(bin_bus)

    import tqdm
    for _ in tqdm.tqdm(I_plain):
        pass

    for _ in tqdm.tqdm(I_bin):
        pass
