import collections
import scanpy as sc
import pandas as pd
import itertools
from scipy.sparse import csr_matrix
from pybustools.pybus import Bus, records_to_gene

import tqdm
from sctools.kallisto import annotate_gene_symbols
from pybustools.utils import read_t2g
from pathlib import Path

def bustools_count_call(bus: Bus, outdir: str, t2g_file: str, optional_args:str=''):
    """
    creates the call to turn a busfile into a count matrix via bustools count
    """
    assert Path(bus.ec_file).exists() and Path(bus.transcript_file).exists() and Path(bus.bus_file).exists()
    return f'bustools count -o {outdir} -e {bus.ec_file} -g {t2g_file} -t {bus.transcript_file} --genecounts {optional_args} {bus.bus_file}'


def _list_of_expression_vectors_to_matrix(expressionvectors, all_genes):
    """
    turns a list of expression vectors (dicts of gene->abundance) into a
    sparse matrix. Order of the columns is determined by `all_genes`
    """
    ii, jj, vv = [], [], []
    gene_to_ix = {g: i for i, g in enumerate(all_genes)}
    for i, ev in enumerate(expressionvectors):
        for gene, exp in ev.items():
            ii.append(i)
            jj.append(gene_to_ix[gene])
            vv.append(exp)

    X = csr_matrix((vv, (ii, jj)), shape=(len(expressionvectors), len(all_genes)))
    return X


def _group_records_by_cb_umi(record_list):
    D = collections.defaultdict(list)

    for r in record_list:
        D[(r.CB, r.UMI)].append(r)

    for (cb, umi), records in D.items():
        yield (cb, umi), records


def _records2genevector_better(records, ec2gene_dict):
    """
    turns a set of bus-records (from one cell) into a dict: gene->abundance
    multimapped records (EC mapping to more than one gene) are discarded

    TODO this doesnt deal with the same UMI occuring with different ECs
    """
    expr_vector = collections.defaultdict(int)
    n_multimapped = 0
    n_unmappable = 0

    for (_, _), rlist in _group_records_by_cb_umi(records):
        genes = records_to_gene(rlist, ec2gene_dict)
        if len(genes) > 1:
            # multimapped
            n_multimapped += 1
        elif len(genes) == 1:
            genes = list(genes)[0]
            expr_vector[genes] += 1
        else:
            #unmappable
            n_unmappable += 1
    return expr_vector, n_multimapped, n_unmappable



def _records2genevector(records, ec2gene_dict):
    """
    turns a set of bus-records (from one cell) into a dict: gene->abundance
    multimapped records (EC mapping to more than one gene) are discarded

    TODO this doesnt deal with the same UMI occuring with different ECs
    """
    expr_vector = collections.defaultdict(int)
    n_multimapped = 0
    for r in records:
        genes = ec2gene_dict[r.EC]
        if len(genes) > 1:
            pass  # multimapped
            n_multimapped += 1
        else:
            genes = list(genes)[0]
            expr_vector[genes] += 1
    return expr_vector, n_multimapped


def _kallisto_count(cell_iterator, ec_dict, transcript_dict, t2g_file):
    """
    python version of the kallisto count command. Turns a busfile into a adata object
    """
    t2g_df = read_t2g(t2g_file)  # transcripts -> gene name
    t2g = t2g_df.set_index('transcript_id')['ensembl_id'].to_dict()

    ec2tr = {EC: [transcript_dict[_] for _ in ec_dict[EC]] for EC in ec_dict.keys()}  # EC -> transcript
    ec2gene = {EC: set(t2g[t] if t in t2g else t for t in ec2tr[EC]) for EC in ec_dict.keys()}

    all_genes = set(t2g.values()) | set(itertools.chain.from_iterable(ec2gene.values()))
    all_genes = sorted(list(all_genes))

    expressionvectors = []
    cbs = []
    n_multimapped = 0
    n_total = 0
    for cb, recordlist in tqdm.tqdm(cell_iterator, desc='Iterating cells'):
        expr_vector, n_multi = _records2genevector(recordlist, ec2gene)
        n_multimapped += n_multi
        n_total += len(recordlist)
        expressionvectors.append(expr_vector)
        cbs.append(cb)

#     # build the count matrix and adata
    X = _list_of_expression_vectors_to_matrix(expressionvectors, all_genes)
    adata = sc.AnnData(
        X,
        obs=pd.DataFrame(cbs, columns=['CB']).set_index('CB'),
        var=pd.DataFrame(all_genes, columns=['var_index']).set_index('var_index'))

    print('multimapped', n_multimapped)
    print('total', n_total)

    adata = annotate_gene_symbols(adata)
    return adata


def kallisto_count(bus: Bus, t2g_file):
    """
    python version of the kallisto count command. Turns a busfile into a adata object
    """
    cell_iterator = bus.iterate_cells()
    return _kallisto_count(cell_iterator, bus.ec_dict, bus.transcript_dict, t2g_file)
