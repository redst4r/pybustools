from pybustools.pybustools import iterate_CB_UMI_of_busfile, iterate_cells_of_busfile
import tqdm
import os
import tempfile
import subprocess
import sys


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


import tqdm
import h5py
from pybustools.busio import Bus_record, write_busfile, _encode_ACGT_to_int
def h5_to_bus(h5filename, busfile_output):
    """
    turns a 10x molecule_info.h5 into a "fake" bus file:
    Instead of EC, we just report the actual gene cellranger mapped the read.
    """
    fh = h5py.File(h5filename)
    CB_list = [_.decode() for _ in fh['/barcodes'][:]]
    gene_list = [_.decode() for _ in fh['/features/name'][:]]
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

    TMPDIR = '/home/mstrasse/mountSSD/tmp'
    unsorted_name = tempfile.mkstemp('.bus', 'unsorted_', TMPDIR)[1]

    write_busfile(unsorted_name, _gen(), cb_length=16, umi_length=12)

    fh.close()
    
    
    print('sorting')
    ret = subprocess.run(["bustools", "sort", '-o', busfile_output, unsorted_name])
    if ret.returncode != 0:
        print("Child was terminated by signal", ret, file=sys.stderr)
        raise ValueError()
