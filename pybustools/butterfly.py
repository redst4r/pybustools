import collections
import pandas as pd
import toolz
from scipy.stats import binom
import tqdm
import numpy as np
import matplotlib.pyplot as plt
from pybustools.pybustools import iterate_CB_UMI_of_busfile, Bus
from pybustools.busio import Bus_record
# t2gfile='/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt'

class CUHistogram():
    """
    Histogram of the Copies per UMI (CU) with some convenience functions

    dictionary of amplification -> frequecny
    """
    def __init__(self, CU_dict):
        self.histogram = CU_dict  # copies per UMI histogram

    def get_nreads(self):
        """
        return the number of reads in the histogram
        """
        return sum([copies*freq for copies, freq in self.histogram.items()])

    def get_nUMI(self):
        """
        return the number of UMIs in the histogram
        """
        return sum([freq for copies, freq in self.histogram.items() if copies > 0])

    def toVectors(self):
        """
        return copies and frequencies as vectors
        """
        copies, freq = zip(*self.histogram.items())
        return np.array(copies), np.array(freq)

    def FSCM(self):
        """
        return fraction of single copy molecules
        """
        n1 = self.histogram[1]
        nTotal = self.get_nUMI()
        return n1/nTotal

    def __eq__(self, other):
        """Overrides the default implementation"""
        if isinstance(other, CUHistogram):
            return self.histogram == other.histogram
        return False

def _collapsed_gene_busiterator(bus_file: str, ec2gene_dict, verbose=False):
    """
    iterating over a busfile, grouping CB/UMI and gene, i.e. all busrecords
    matchibg the same CB/UMI and gene (EC point to the same gene)
    """
    discarded = 0
    multimapped = 0
    total = 0
    for (cb, umi), record_list in iterate_CB_UMI_of_busfile(bus_file):
        total += 1
        # whats the ECs and are they compatible with a single gene
        # ecs can map to multiple genes
        gene_sets = [set(ec2gene_dict[r.EC]) for r in record_list]
        common_gene = set.intersection(*gene_sets)
        if len(common_gene) == 0:
            # no single gene that can explain it
            discarded += 1
        elif len(common_gene) == 1:
            counts = sum([r.COUNT for r in record_list])
            common_gene = common_gene.pop()  # get the single element of the set
            yield Bus_record(cb, umi, common_gene, counts, record_list[0].FLAG)
        else:
            # more than one gene
            multimapped += 1
#             print(common_gene)
        if total % 5_000_000 == 0 and verbose:
            print(f'Total CB/UMI: {total}\nMultimapped: {multimapped} ({100*multimapped/total:.2f}%)\nDiscarded:{discarded} ({100*discarded/total:.2f}%)')

    print(f'Total CB/UMI: {total}\nMultimapped: {multimapped} ({100*multimapped/total:.2f}%)\nDiscarded:{discarded} ({100*discarded/total:.2f}%)')


def _make_transcript2gene_df(t2gfile: str):
    """
    loading a transcript->gene file, turning it into a proper dataframe
    """
    
    t2g_df = pd.read_csv(t2gfile,
                         sep='\t', header=None,
                         names=['transcript_id', 'gene_id', 'symbol'])
    return t2g_df


def _make_ec2gene_dict(bus: Bus, t2gfile):
    """
    a dict mapping EC to a list of genes

    :params bus: a pybustools.Bus object
    :params t2gfile: filename of the transcript-to-gene mapping
    """
    t2g_df = _make_transcript2gene_df(t2gfile)
    t2g = {row['transcript_id']: row['symbol'] for _, row in t2g_df.iterrows()}

    def genes_for_ec(ec_id):
        "get all genes for a particular EC"
        transcripts = bus.resolve_EC_to_transcripts(ec_id)
        genes = {t2g[trans] for trans in transcripts if trans in t2g}
        return genes

    ec2gene = {ec :genes_for_ec(ec) for ec in tqdm.tqdm(bus.ec_dict.keys(), desc='translating EC to gene')}
    return ec2gene


def make_ec_histograms(bus: Bus, collapse_EC=False, t2gfile=None):
    """
    for all ECs in the busfile create a CU-histogram (a histogram of the
    amplification, ie. number of reads per UMI)

    :params bus: a pybustools.Bus object
    :returns: a dictionary of CU histograms (one item per EC)
    """
    ec_hists = {}

    if collapse_EC:
        assert t2gfile is not None
        ec2g = _make_ec2gene_dict(bus, t2gfile)
        I = _collapsed_gene_busiterator(bus.bus_file, ec2g)
    else:
        I = iterate_CB_UMI_of_busfile(bus.bus_file)
        # filtering records that map to more than one EC
        # this is not totally correct (the two ECs might overlap in a single gene)
        I = (record_list[0] for (cb, umi), record_list in I if len(record_list) == 1)

    for r in tqdm.tqdm(I):
        if not r.EC in ec_hists:
            ec_hists[r.EC] = collections.defaultdict(int)

        ec_hists[r.EC][r.COUNT] += 1
    ec_hists = toolz.valmap(CUHistogram, ec_hists)  # turn into class
    return ec_hists

def make_ec_histogram_across_genes(bus: Bus, collapse_EC=False, t2gfile=None):
    """
    create one big CU-histogram (pooled over genes). In contrast make_ec_histograms() does this gene by gene
    This method here has a much smaller memory footprint and for things like unseen species it is enough 
    (we'd aggregate the gene wise CUs anyway)!
    
    Note that theres still a collapse_EC argument: This determines how we count: Suppose there's two entries of CB/UMI/EC1 and CB/UMI/EC2.
    If Gene[EC1] == Gene[EC2], do we count this a single molecule with two copies (collapse_EC=True) or two molecules with a single copy (collapse_EC=False)
    Actually, the collapse_EC=False would be dropped as it is multimapped!!
    
    :params bus: a pybustools.Bus object
    :returns: a dictionary of CU histograms (one item per EC)
    """
    ec_hist = collections.defaultdict(int)

    if collapse_EC:
        assert t2gfile is not None
        ec2g = _make_ec2gene_dict(bus, t2gfile)
        I = _collapsed_gene_busiterator(bus.bus_file, ec2g)
    else:
        I = iterate_CB_UMI_of_busfile(bus.bus_file)
        # filtering records that map to more than one EC
        # this is not totally correct (the two ECs might overlap in a single gene)
        I = (record_list[0] for (cb, umi), record_list in I if len(record_list) == 1)

    for r in tqdm.tqdm(I):
        ec_hist[r.COUNT] += 1
        
    ec_hist = CUHistogram(ec_hist)  # turn into class
    return ec_hist




def binomial_downsample_all_genes(CU_histogram_dict, fraction):
    """
    instead of doing it gene by gene [as binomial_downsample()] do it for
    all genes at once. should be faster!
    """

    jmax_per_gene = toolz.valmap(lambda CU_hist: np.max(list(CU_hist.histogram.keys())), CU_histogram_dict)
    jmax = np.max(list(jmax_per_gene.values()))

    j = np.arange(jmax+1)

    # precalculate this binomial distribution for the various values of n
    binomial_vector = {i: binom.pmf(k=j, n=i, p=fraction) for i in range(jmax+1)}

    downsampled_dict = {}  # TODO toolz.valmap instead
    for gene, CU_histogram in CU_histogram_dict.items():
        hnew = np.zeros(jmax+1)
        for i, hi in CU_histogram.histogram.items():
            _t2 = binomial_vector[i] * hi
            hnew = hnew + _t2
        res = dict(zip(j, hnew))

        # sometimes values are 0, these are useless so filter
        res = toolz.valfilter(lambda v: v != 0, res)
        downsampled_dict[gene] = CUHistogram(res)
    return downsampled_dict

def binomial_downsample(CU: CUHistogram, fraction):
    """
    given a histogram of counts per UMI for a particular gene
    calculate the histogram after downsampling the reads to `fraction`

    :params CU: a histogram (dictionary of amplification -> frequecny)
    :params fraction: fraction of reads to downsample to

    :returns: another CU histrogram, downsampled accordingly
    """
    jmax = np.max(list(CU.histogram.keys()))

    hnew = np.zeros(jmax+1)
    j = np.arange(jmax+1)

    for i, hi in CU.histogram.items():
        _t = binom.pmf(k=j, n=i, p=fraction) * hi
        hnew = hnew + _t
    return CUHistogram(dict(zip(j, hnew)))


def binomial_downsample_factors(old_CU_hist, new_CU_hist, CPM='reads'):
    """
    getting the scale factors for each gene

    :params old_CU_hist: a histogram (dictionary of amplification -> frequecny)
    before downsampling
    :params new_CU_hist: a histogram (dictionary of amplification -> frequecny)
    after downsampling

    :returns: a dict which has the scaling factor for each gene/ec
    """

    assert CPM in ['reads', 'umis', 'raw'], "CPM must be either 'reads' or 'umis' or 'raw'"

    # per gene/EC UMI counts
    counts_before = toolz.valmap(lambda CU: CU.get_nUMI(), old_CU_hist)
    counts_after  = toolz.valmap(lambda CU: CU.get_nUMI(), new_CU_hist)
    ecs = list(sorted(old_CU_hist.keys()))
    counts_before = np.array([counts_before[ec] for ec in ecs])
    counts_after  = np.array([counts_after[ec]  for ec in ecs])

    # reads across all genes
    n_reads_before = sum(toolz.valmap(lambda CU: CU.get_nreads(), old_CU_hist).values())
    n_reads_after  = sum(toolz.valmap(lambda CU: CU.get_nreads(), new_CU_hist).values())

    n_umi_before = counts_before.sum()
    n_umi_after = counts_after.sum()

    print(f'UMIs: {n_umi_before} - {n_umi_after} ({100*n_umi_after/n_umi_before:.3f} %)')
    print(f'Reads: {n_reads_before} - {n_reads_after} ({100*n_reads_after/n_reads_before:.3f} %)')

    if CPM == 'reads':
        before_cpm = n_reads_before
        after_cpm = n_reads_after
    elif CPM == 'umis':
        before_cpm = n_umi_before
        after_cpm = n_umi_after
    elif CPM == 'raw':
        before_cpm = 1
        after_cpm = 1
    else:
        raise ValueError(f'unknown CPM: {CPM}')

    CPM_before = counts_before / before_cpm
    CPM_after = counts_after / after_cpm

    factor = CPM_after / CPM_before

    return dict(zip(ecs, factor))


def aggregate_CUs(CU_dict):
    """
    aggregates (adds up) all the CU_dicts (each dict for a different gene) into a
    single CU dict (agnostic of the gene). Useful for the sequencing depth saturation plots etc
    """
    histograms = (h.histogram for h in CU_dict.values())  # generator to save mem
    hall = toolz.merge_with(sum, histograms)
    return CUHistogram(hall)


from scipy.spatial.distance import cdist
def compare_histograms_OT(h1, h2):
    import ot
    """
    optimal transport distance on the CU histograms. Better than KL since the
    domains might not overlap!
    """
    x_min = np.minimum(min(h1.keys()), min(h2.keys()))
    x_max = np.maximum(max(h1.keys()), max(h2.keys()))
    xrange = np.arange(x_min, x_max)
    M = cdist(xrange.reshape(-1,1), xrange.reshape(-1,1))

    a = np.array([h1[x] for x in xrange])
    b = np.array([h2[x] for x in xrange])
    a = a/a.sum()
    b = b/b.sum()
    return ot.emd2(a, b, M)


def plot_CU(CU, norm=True):
    q = pd.Series(CU.histogram).sort_index()

    values = q.values
    if norm:
        values = values/values.sum()
    ix = values!=0
    plt.scatter(q.index[ix], values[ix])
    plt.xlabel('Reads per UMI')
    plt.ylabel('Fraction')


def saturation_curve(CU_aggr, bins=20):
    down_percentages = np.linspace(0.01, 1, bins)
    df_down2 = []
    for f in tqdm.tqdm(down_percentages):
        hdown = binomial_downsample(CU=CU_aggr, fraction=f)
        n_reads = hdown.get_nreads()# sum([k*v for k,v in hdown.items()])
        n_umi = hdown.get_nUMI()# sum([v for k,v in hdown.items() if k > 0])
        df_down2.append({
            'n_reads': n_reads,
            'n_umi': n_umi,
            'f': f,
        })
    df_down2 = pd.DataFrame(df_down2)
    return df_down2


import multiprocessing as mp
def _parralel_helper_saturation(gene, f, CU):
    """
    given a CU, downsample it to fraction f
    gene is return to emulate a parallel dict.valmap
    """
    return gene, binomial_downsample(CU, fraction=f)


def saturation_curve_per_gene(CU_dict, fractions=None, cores=1):
    if fractions is None:
        down_fractions = np.linspace(0.01,1,10)
    else:
        down_fractions = fractions
    df_down = []
    for f in tqdm.tqdm(down_fractions):
    #     hdown = toolz.valmap(lambda CU: binomial_downsample(CU, fraction=f), h)

        with mp.Pool(cores) as pool:
            inputs = [(gene, f, CU) for gene,CU in CU_dict.items()]
            hdown = pool.starmap(_parralel_helper_saturation, inputs)
            hdown = dict(hdown)

        n_reads = sum(toolz.valmap(lambda CU: sum([k*v for k,v in CU.items()]), hdown).values())
        n_umi = sum(toolz.valmap(lambda CU: sum([v for k,v in CU.items() if k>0]), hdown).values())
        df_down.append({
            'n_reads': n_reads,
            'n_umi': n_umi,
            'f': f,
        })
    df_down = pd.DataFrame(df_down)
    return df_down

"""
ec2g = _make_ec2gene_dict(bus, t2gfile)
h = make_ec_histograms(bus)

import toolz
h_new = toolz.valmap(lambda CUhist: binomial_downsample(CUhist, 0.50), h)

fg = binomial_downsample_factors(h2, h2_new)
"""
