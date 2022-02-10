import collections
import scanpy as sc
import pandas as pd
import itertools
from scipy.sparse import csr_matrix
from pybustools.pybustools import Bus
import tqdm
from sctools.kallisto import annotate_gene_symbols


def read_t2g(t2g_file):
    """
    reading the kallisto trnascript2gene file
    :returns: a dictionary, resolving each transcript to a ensembl_gene_id
    """
    df = pd.read_csv(t2g_file, sep='\t', header=None, names=['transcript_id', 'ensembl_id', 'gene_symbol'])
#     t2g_dict = df.set_index('transcript_id')['gene_symbol'].to_dict()
    t2g_dict = df.set_index('transcript_id')['ensembl_id'].to_dict()
    return t2g_dict


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


def _records2genevector(records, ec2gene_dict):
    """
    turns a set of bus-records (from one cell) into a dict: gene->abundance
    multimapped records (EC mapping to more than one gene) are discarded
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


def kallisto_count(bus: Bus, t2g_file):
    """
    python version of the kallisto count command. Turns a busfile into a adata object
    """

    # establish the mapping from EC->set(gene_ids):
    # -------------------------------------------
    t2g = read_t2g(t2g_file)  # transcripts -> gene name
    ec2tr = {EC: [bus.transcript_dict[_] for _ in bus.ec_dict[EC]] for EC in bus.ec_dict.keys()}  # EC -> transcript
    ec2gene = {EC: set(t2g[t] if t in t2g else t for t in ec2tr[EC]) for EC in bus.ec_dict.keys()}

    all_genes = set(t2g.values()) | set(itertools.chain.from_iterable(ec2gene.values()))
    all_genes = sorted(list(all_genes))

    expressionvectors = []
    cbs = []
    n_multimapped = 0
    n_total = 0
    for cb, recordlist in tqdm.tqdm(bus.iterate_cells()):
        expr_vector, n_multi = _records2genevector(recordlist, ec2gene)
        n_multimapped += n_multi
        n_total += len(recordlist)
        expressionvectors.append(expr_vector)
        cbs.append(cb)

#     # build the count matrix and adata
    X = _list_of_expression_vectors_to_matrix(expressionvectors, all_genes)
    adata = sc.AnnData(
        X,
        bs=pd.DataFrame(cbs, columns=['CB']).set_index('CB'),
        var=pd.DataFrame(all_genes, columns=['var_index']).set_index('var_index'))

    print('multimapped', n_multimapped)
    print('total', n_total)

    adata = annotate_gene_symbols(adata)
    return adata
