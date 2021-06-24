import collections
import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pybustools.pybustools import Bus


"""
a few functions to measure the 10x library amplification in kallisto-bus files
"""


def create_per_EC_amp(bus: Bus):
    """
    quantify the cDNA amplification of each molecule, grouped by EC
    (equivalence class, related to gene). The resulting dataframe will have one
    row per EC/gene and this ECs's amplification characteristics
    """
    ec_counts = collections.defaultdict(list)
    for record in tqdm.tqdm(bus.iterate_bus()):
        ec_counts[record.EC].append(record.COUNT)

    df = [
        {'EC': ec,
         'nUMI': len(ec_counts[ec]),
         'mean_amp': np.mean(ec_counts[ec]),
         'std_amp': np.std(ec_counts[ec]),
         'total_counts': np.sum(ec_counts[ec]),
         'max_amp': np.max(ec_counts[ec]),
         'percent_single_read': np.mean(np.array(ec_counts[ec]) == 1)
         } for ec in ec_counts.keys()]
    df = pd.DataFrame(df)

    # add some more info to the EC
    t2g_df = pd.read_csv('/home/mstrasse/resources/transcripts_to_genes.txt',
                         sep='\t', header=None,
                         names=['transcript_id', 'gene_id', 'symbol'])
    t2g = {row['transcript_id']: row['symbol'] for _, row in t2g_df.iterrows()}

    def genes_for_ec(ec_id):
        transcripts = [bus.transcript_dict[_] for _ in bus.ec_dict[ec_id]]
        genes = {t2g[trans] for trans in transcripts if trans in t2g}
        return genes

    df['nGenes'] = df.EC.apply(lambda x: len(genes_for_ec(x)))

    # if the EC is a unique gene, annotate taht too
    def _h(x):
        genes = genes_for_ec(x)
        return list(genes)[0] if len(genes) == 1 else 'multiple'
    df['the_gene'] = df.EC.apply(_h)
    return df


def create_per_cell_amp(bus: Bus):
    """
    quantify the cDNA amplification of each molecule, grouped by cells.
    the resulting dataframe will have a row per cell and this cell's
    amplification characteristics
    """
    cb_counts = collections.defaultdict(list)
    for record in tqdm.tqdm(bus.iterate_bus()):
        cb_counts[record.CB].append(record.COUNT)
    df_cb = [
        {'cb': cb,
         'nUMIs': len(cb_counts[cb]),
         'mean_amp': np.mean(cb_counts[cb]),
         'std_amp': np.std(cb_counts[cb]),
         'total_counts': np.sum(cb_counts[cb]),
         'max_amp': np.max(cb_counts[cb]),
         'percent_single_read': np.mean(np.array(cb_counts[cb]) == 1)
         } for cb in cb_counts.keys()]
    df_cb = pd.DataFrame(df_cb)
    return df_cb


def plot_amplification(df_amp, label=None):

    plt.scatter(df_amp['nUMIs'], df_amp['mean_amp'], s=1, alpha=0.5, label=label)

    plt.xscale('log')
#     plt.xlim(10,60000)
#     plt.ylim(0.99,2)
#     plt.title(sname)
    plt.xlabel('#UMIs')
    plt.ylabel('#reads/#UMIs')
    plt.grid()
