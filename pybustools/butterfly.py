from scipy.stats import binom
import collections
import tqdm
import numpy as np
from pybustools.pybustools import iterate_CB_UMI_of_busfile
import pandas as pd

def make_ec2gene_dict(bus, t2gfile):
    """
    a dict mapping EC to a list of genes

    :params bus: a pybustools.Bus object
    :params t2gfile: filename of the transcript-to-gene mapping
    """
    t2g_df = pd.read_csv(t2gfile,
                         sep='\t', header=None,
                         names=['transcript_id', 'gene_id', 'symbol'])
    t2g = {row['transcript_id']: row['symbol'] for _, row in t2g_df.iterrows()}

    def genes_for_ec(ec_id):
        "get all genes for a particular EC"
        transcripts = [bus.transcript_dict[_] for _ in bus.ec_dict[ec_id]]
        genes = {t2g[trans] for trans in transcripts if trans in t2g}
        return genes

    ec2gene = {ec :genes_for_ec(ec) for ec in tqdm.tqdm(bus.ec_dict.keys())}
    return ec2gene


def make_ec_histograms(bus):
    """
    for all ECs in the busfile create a CU-histogram (a histogram of the
    amplification, ie. number of reads per UMI)

    :params bus: a pybustools.Bus object
    :returns: a dictionary of CU histograms (one item per EC)
    """
    ec_hists = {}
    n_multimapped = 0

    for (cb, umi), record_list in tqdm.tqdm(iterate_CB_UMI_of_busfile(bus.bus_file)):
        if len(record_list) == 1:
            r = record_list[0]
            if r.EC in ec_hists:
                ec_hists[r.EC][r.COUNT] += 1
            else:
                ec_hists[r.EC] = collections.defaultdict(int)
                ec_hists[r.EC][r.COUNT] += 1
        else:
            n_multimapped += 1
    return ec_hists


def binomial_downsample(CU_histogram, fraction):
    """
    given a histogram of counts per UMI for a particular gene
    calculate the histogram after downsampling the reads to `fraction`

    :params CU_histogram: a histogram (dictionary of amplification -> frequecny)
    :params fraction: fraction of reads to downsample to

    :returns: another CU histrogram, downsampled accordingly
    """
    jmax = np.max(list(CU_histogram.keys()))

    hnew = np.zeros(jmax+1)
    j = np.arange(jmax+1)

    for i, hi in CU_histogram.items():
        _t = binom.pmf(k=j, n=i, p=fraction) * hi
        hnew = hnew + _t
    return dict(zip(j, hnew))


def binomial_downsample_factors(old_CU_hist, new_CU_hist, CPM='reads'):
    """
    getting the scale factors for each gene

    :params old_CU_hist: a histogram (dictionary of amplification -> frequecny)
    before downsampling
    :params new_CU_hist: a histogram (dictionary of amplification -> frequecny)
    after downsampling

    :returns: a dict which has the scaling factor for each gene/ec
    """

    assert CPM in ['reads', 'umis'], "CPM must be either 'reads' or 'umis'"
    ecs = []
    counts_before = []
    counts_after = []

    # TODO this could be done more elegantly with toolz.valmap
    for ec in tqdm.tqdm(old_CU_hist.keys()):
        old_CU = old_CU_hist[ec]
        new_CU = new_CU_hist[ec]
        ecs.append(ec)
        # how many UMIs survived the downsampling
        counts_before.append(sum([freq for nreads, freq in old_CU.items() if nreads>0]))
        counts_after.append(sum([freq for nreads, freq in new_CU.items() if nreads>0]))

    counts_before = np.array(counts_before)
    counts_after  = np.array(counts_after)

    def nreads_from_CU(CU):
        "how many reads in total are in the CU histogram"
        return sum([freq*reads for reads, freq in CU.items()])

    def numi_from_CU(CU):
        "how many UMIs  in total are in the CU histogram"
        return sum([freq for reads, freq in CU.items() if reads>0])

    n_reads_before = sum(toolz.valmap(nreads_from_CU, old_CU_hist).values())
    n_reads_after  = sum(toolz.valmap(nreads_from_CU, new_CU_hist).values())

    n_umi_before = sum(toolz.valmap(numi_from_CU, old_CU_hist).values())
    n_umi_after = sum(toolz.valmap(numi_from_CU, new_CU_hist).values())

    assert counts_before.sum() == n_umi_before
    assert int(counts_after.sum()) == int(n_umi_after)

    print(f'UMIs: {n_umi_before} - {n_umi_after} ({100*n_umi_after/n_umi_before} %)')
    print(f'Reads: {n_reads_before} - {n_reads_after} ({100*n_reads_after/n_reads_before} %)')

    if CPM == 'reads':
        before_cpm = n_reads_before
        after_cpm = n_reads_after
    elif CPM == 'umis':
        before_cpm = n_umi_before
        after_cpm = n_umi_after
    else:
        raise ValueError(f'unknown CPM: {CPM}')

    CPM_before = counts_before / before_cpm
    CPM_after = counts_after / after_cpm

    factor = CPM_after / CPM_before

    return dict(zip(ecs, factor))


"""
ec2g = make_ec2gene_dict(bus, t2gfile)
h = make_ec_histograms(bus)

import toolz
h_new = toolz.valmap(lambda CUhist: binomial_downsample(CUhist, 0.50), h)

fg = binomial_downsample_factors(h2, h2_new)
"""
