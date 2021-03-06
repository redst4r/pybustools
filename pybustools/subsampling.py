import numpy as np
from pybustools.busio import read_binary_bus, Bus_record, write_busfile, get_header_info
from pybustools.pybustools import iterate_CB_UMI_of_busfile
import tqdm
import os
import shutil
import numba
import gc
import collections
import time
import pandas as pd


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


def subsample_bus_unseens_species(busfile, fractions):
    """
    subsamples the given bus file at different depth, and records the number of reads/umis obtained
    returns a DataFrame with the nReads and nUMIs as well as the UMI prevalences (duplicity of an umi -> #observations of that duplicity)
    """
    gc.collect()

    assert all([0 <= _ <= 1 for _ in fractions]), "Fractions must be in [0,1]"
    reads = []

    # lets group entires by CB/UMI instead of CB/UMI/EC (as done via B.iterate_bus())
    # and sum the counts over ECs
    for (cb, umi), record_list in tqdm.tqdm(iterate_CB_UMI_of_busfile(busfile, decode_seq=False)):
        s = sum([r.COUNT for r in record_list])
        reads.append(s)

    reads = np.array(reads)
    prevalences = collections.Counter(reads)

    df = []
    total = reads.sum()
    for percent in fractions:
        print(percent)
        target = int(percent * total)
        y = _downsample_array(reads, target=target, random_state=int(time.time()*1000), replace=False)

        numis = (y > 0).sum()
        df.append({
            'nUMIs': numis,
            'percent': percent,
            'nReads':  target,
            'nReads2':  y.sum()
        })
    df = pd.DataFrame(df)
    return df, prevalences


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
    x = _downsample_array(huge_array, target=n_target, random_state=int(time.time()*1000), replace=False, inplace=False)
    print(f'Downsampled reads: {x.sum()}')

    # create a generator for the bus-records
    def _helper_gen():
        I = read_binary_bus(fname_in, decode_seq=False)
        for i, record in tqdm.tqdm(enumerate(I), desc='Second pass'):
            if x[i] > 0:
                r = Bus_record(record.CB, record.UMI, record.EC, x[i], record.FLAG)
                yield r

    G = _helper_gen()
    # we need to write the correct header
    _, cb_len, umi_len, _ = get_header_info(fname_in)
    write_busfile(fname_out, G, cb_length=cb_len, umi_length=umi_len)


@numba.njit(cache=True)
def _downsample_array(
    col: np.ndarray,
    target: int,
    random_state=0,
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


def subsample_kallisto_bus(busdir, outdir, fraction):
    """
    not only subsamples the .bus file, but also sets up a folder that looks like a original kallisto quantification
    so that `bustools count` can be run on it

    :param busdir: directory containing the bustools quantification (if nextflow: .../kallisto/sort_bus/bus_output).
                   must contain `output.corrected.sort.bus`, `matrix.ec`, `transcripts.txt`

    :param outdir: must exist. puts the downsampled kallisto .bus and suppl.files here
    :param fraction: 0<fraction<1: fraction of reads to sample

    """
    assert os.path.exists(outdir)
    infile = f'{busdir}/output.corrected.sort.bus'
    outfile = f'{outdir}/output.corrected.sort.bus'

    subsample_busfile(infile, outfile, fraction=fraction)

    # copy the other parts of the kallisto output, they remain unchanged
    ec = f'{busdir}/matrix.ec'
    transcripts = f'{busdir}/transcripts.txt'

    shutil.copy(ec, outdir)
    shutil.copy(transcripts, outdir)
