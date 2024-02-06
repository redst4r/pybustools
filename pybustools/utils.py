from pybustools.pybus import iterate_CB_UMI_of_busfile, iterate_cells_of_busfile
import tqdm
import os
import tempfile
import subprocess
import sys
import h5py
from pybustools.busio import Bus_record, write_busfile, _encode_ACGT_to_int
import pandas as pd
from multiprocessing import Pool
import toolz

def map_ec2gene(bus, t2gfile):
    """
    create a (busfile specific) dict to resolve ECs to gene
    """
    t2g = read_t2g(t2gfile).set_index('transcript_id')['gene_symbol'].to_dict()
    ec2tr = {EC: [bus.transcript_dict[_] for _ in bus.ec_dict[EC]] for EC in bus.ec_dict.keys()}  # EC -> transcript
    ec2gene = {EC: list(set(t2g[t] if t in t2g else t for t in ec2tr[EC])) for EC in bus.ec_dict.keys()}
    return ec2gene

def read_t2g(t2g_file):
    """
    reading the kallisto trnascript2gene file
    :returns: a dataframe with transcript id, gene id and gene symbol
    """
    df = pd.read_csv(t2g_file, sep='\t', header=None, names=['transcript_id', 'ensembl_id', 'gene_symbol'])
#     t2g_dict = df.set_index('transcript_id')['gene_symbol'].to_dict()
    # t2g_dict = df.set_index('transcript_id')['ensembl_id'].to_dict()
    return df


def count_cb_umi_pairs(busfile):
    counter = 0
    for _ in tqdm.tqdm(iterate_CB_UMI_of_busfile(busfile, decode_seq=False)):
        counter += 1
    return counter


def count_cb_umi_pairs2(busfile):
    counter = 0
    for cb, info in tqdm.tqdm(iterate_cells_of_busfile(busfile, decode_seq=False)):
        s = set()
        for umi, ec, count, flag in info:
            s.add(umi)
        counter += len(s)
    return counter


def bus_to_countmatrix(busfile, ec_file, transcript_file, outfolder):
    """
    sort a busfile and call `bustools count` on it
    """
    TMPDIR = '/home/mstrasse/mountSSD/tmp'
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)

    # create a sorted version
    sorted_name = tempfile.mkstemp('.bus', 'bussort', TMPDIR)[1]
    print('sorting')
    ret = subprocess.run(["bustools", "sort", '-o', sorted_name, busfile])
    if ret.returncode != 0:
        print("Child was terminated by signal", ret, file=sys.stderr)
        raise ValueError()

    print('counting')

    t2g = '/home/mstrasse/resources/transcripts_to_genes.txt'
#     cmd = f'mkdir {outfolder} && bustools count -o {outfolder}/gene -g {t2g} -e {ec_file} -t {transcript_file} --genecounts {sorted_name}'
#     os.system(cmd)

    ret = subprocess.run(["bustools", "count", '-o', f'{outfolder}/gene', '-g', t2g, '-e', ec_file, '-t', transcript_file, '--genecounts', sorted_name])

    if ret.returncode != 0:
        print("Child was terminated by signal", ret, file=sys.stderr)
        raise ValueError()

    os.unlink(sorted_name)

    """
    ret = subprocess.run(["gsutil", "cp", '-r', gs_url, target_folder])
    if ret.returncode != 0:
    print("Child was terminated by signal", ret, file=sys.stderr)
    """


def h5_to_bus(h5filename, busfile_output, TMPDIR=None):
    """
    turns a 10x molecule_info.h5 into a "fake" bus file:
    Instead of EC, we just report the actual gene cellranger mapped the read.

    :param h5filename: path to the h5file
    :param busfile_output: path to the busfile to be created
    :param TMPDIR: optional, path to store temporary files
    """
    fh = h5py.File(h5filename, mode='r')
    CB_list = [_.decode() for _ in fh['/barcodes'][:]]
    # gene_list = [_.decode() for _ in fh['/features/name'][:]]
    n_entries = fh['/barcode_idx'].shape[0]

    cbs_idx = fh['/barcode_idx'][:]
    gene_idx = fh['/feature_idx'][:]
    counts = fh['/count'][:]
    umis = fh['/umi'][:]

    def _gen():
        for i in tqdm.trange(n_entries):
            cb = CB_list[cbs_idx[i]]
            cb = _encode_ACGT_to_int(cb)
            b = Bus_record(cb,
                           int(umis[i]),
                           gene_idx[i],
                           counts[i],
                           0)
            yield b

    unsorted_name = tempfile.mkstemp('.bus', 'unsorted_', TMPDIR)[1]

    write_busfile(unsorted_name, _gen(), cb_length=16, umi_length=12)

    fh.close()


    print('sorting')
    ret = subprocess.run(["bustools", "sort", '-o', busfile_output, unsorted_name])
    if ret.returncode != 0:
        print("Child was terminated by signal", ret, file=sys.stderr)
        raise ValueError()

    # note that the tmp file wont be deleted if an exception happens above!
    os.unlink(unsorted_name)


    # in addition create matrix.ec and transcripts.txt
    # Matrix.ec is a two column text file with
    # EC_id  --> [transcript_ID]
    #
    # nd transcripts.txt is just transcript names (ENST....), one per line (linenumber == transcript_ID)
    #
    # TODO: DO THIS, currently not even enccessary, as we only use this to determine non-unique CB/UMIs across multiple bus files


def _dummy(func, k, v):
    """
    helper for the parallel valmap, wrapping the val-function so that it also returns the key

    d = {'a':1,'b':2}
    def f(x):
        return x+1
    parallel_valmap(f, d, cores=8)
    """
    return k, func(v)


def parallel_valmap(func, dictionary, cores=8):
    """
    just as toolz.valmap, but do it in paralell across entries
    """
    kv_pairs = [(k,v) for k,v in dictionary.items()]
    s = toolz.partial(_dummy, func) 

    with Pool(cores) as p:
        res = p.starmap(s, kv_pairs)

    res2 = {k: v for k,v in res}
    return res2

