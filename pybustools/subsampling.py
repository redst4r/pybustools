import numpy as np
from pybustools.busio import read_binary_bus, Bus_record, write_busfile
import tqdm


def get_number_of_reads_and_molecules(fname):
    """
    similar to bustools inspect
    Gets the number of reads and number of bus-records (approx #UMI) in a busfile

    :param fname: filename of the busfile
    """
    # TODO: technically not correct: this count the number of reads (fine) and number of ENTRIES in the busfile.
    # ideally ach molecule would have a single entry.
    # However: a single molecule might map to two different EC classes. The molecule got fragmented into two places, mapping it to two diff locations
    total_reads = 0
    total_molecules = 0
    for record in tqdm.tqdm(read_binary_bus(fname, decode_seq=False), desc='counting reads'):
        total_reads += record.COUNT
        total_molecules += 1
    return total_reads, total_molecules


def subsample_busfile(fname_in, fname_out, fraction):
    """
    subsample the reads of an existing busfile by `fraction`, writing the Result
    into a new busfile. The major effect is that some entries will recieve
    0 reads and hence disappear from the file!

    :param fname_in: Filename of the input busfile
    :param fname_out: Filename of the resulting, subsampled busfile
    :param fraction: 0<fraction<1, the fraction of subsampling
    """
    assert 0 < fraction < 1, "fraction must be in [0,1]"

    # for this to work, we have to pass the input file twice:
    # 1. we  have to collect ALL the counts (for each record)
    #    this will then be jointly subsampled
    # 2. we iterate the inputfile again, just now we write out each reocrdin into a
    # different busfile with adjusted count
    huge_array = []
    I = read_binary_bus(fname_in, decode_seq=False)
    for record in tqdm.tqdm(I, desc='First pass'):
        huge_array.append(record.COUNT)

    huge_array = np.array(huge_array)
    n_total = np.sum(huge_array)
    n_target = int(n_total * fraction)
    print(f'Subsampling from {n_total} to {n_target}')
    x = _downsample_array(huge_array, target=n_target, replace=False, inplace=False)
    print(f'Downsampled reads: {x.sum()}')

    # create a generator for the bus-records
    def _helper_gen():
        I = read_binary_bus(fname_in, decode_seq=False)
        for i, record in tqdm.tqdm(enumerate(I), desc='Second pass'):
            if x[i] > 0:
                r = Bus_record(record.CB, record.UMI, record.EC, x[i], record.FLAG)
                yield r

    G = _helper_gen()
    write_busfile(fname_out, G, cb_length=16, umi_length=10)


def _downsample_array(
    col: np.ndarray,
    target: int,
    # random_state: AnyRandom = 0,
    replace: bool = True,
    inplace: bool = False,
):
    """\
    From scanpy!
    Evenly reduce counts in cell to target amount.
    This is an internal function and has some restrictions:
    * total counts in cell must be less than target
    """
    # np.random.seed(random_state)
    cumcounts = col.cumsum()
    if inplace:
        col[:] = 0
    else:
        col = np.zeros_like(col)
    total = np.int_(cumcounts[-1])
    sample = np.random.choice(total, target, replace=replace)
    sample.sort()
    geneptr = 0
    for count in sample:
        while count >= cumcounts[geneptr]:
            geneptr += 1
        col[geneptr] += 1
    return col
