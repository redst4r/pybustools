import collections
import tqdm
from scipy import sparse
from pybustools.pybus import Bus
import pandas as pd
import scanpy as sc
import numpy as np


def busfolder_to_ec_transcript_compatibility(bus: Bus):
    """
    create a EC-transcript compatibility matrix
    transcripts are ordered as per bus.transcript_dict
    """
    rows = []
    cols = []
    data = []
    transcripts = [bus.transcript_dict[_] for _ in range(len(bus.transcript_dict))]

    for ec, tid in tqdm.tqdm(bus.ec_dict.items()):
        for t in tid:
            rows.append(ec)
            cols.append(t)
            data.append(1)

    nECs = len(bus.ec_dict)
    nTranscripts = len(transcripts)

    # note: sometimes a transcript has NO EC! WEIRD!
    # but for consistency, we include a empty column for those (via the explicit `shape`)
    compatibility_matrix = sparse.coo_matrix((data, (rows, cols)), shape=(nECs,nTranscripts)).tocsr()
    return compatibility_matrix, transcripts

def t2g_into_matrix(t2g):
    """
    create a EC-transcript compatibility matrix
    transcripts are ordered as per bus.transcript_dict

    actually returns the matrix as an annotated dataframe (AnnData)
    to make it easy to resort rows and columns
    """
    rows = []
    cols = []
    data = []
    transcripts = sorted(set(list(t2g.keys())))
    genes = sorted(list(set(t2g.values())))

    print(len(transcripts))
    print(len(genes))

    transcripts2index = {t: i for i,t in enumerate(transcripts)}
    gene2index = {g: i for i,g in enumerate(genes)}

    for tname, gname in tqdm.tqdm(t2g.items()):
        rows.append(transcripts2index[tname])
        cols.append(gene2index[gname])
        data.append(1)

    # note: sometimes a transcript has NO EC! WEIRD!
    # but for consistency, we include a empty column for those (via the explicit `shape`)
    compatibility_matrix = sparse.coo_matrix((data, (rows, cols)), shape=(len(transcripts), len(genes))).tocsr()

    adata = sc.AnnData(compatibility_matrix, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=transcripts))

    return adata

def load_t2g(tx2gene_file: str):
    """
    turns the transcript_to_gene file into a dict
    ENSTxxxx -> ENSGyyyy
    Note: this is a 1:1 mapping (each transcript has exactly one gene)
    """
    t2g = {}
    df_t2g = pd.read_csv(tx2gene_file, sep='\t', header=None)
    for _, row in df_t2g.iterrows():
        assert row[0] not in t2g
        t2g[row[0]] = row[1]
    return t2g

def busfolder_to_ec_gene_compatibility(bus, t2g: dict):
    """
    returns a EC-> gene compatibility matrix.
    rows are ECs, cols are genes (order according to all_genes
    """
    all_genes = set()
    ec2gene = collections.defaultdict(lambda: set())
    for ec in tqdm.tqdm(bus.ec_dict.keys()):
        transcripts = bus.resolve_EC_to_transcripts(ec)
        for t in transcripts:
            g = t2g[t]
            ec2gene[ec].add(g)
            all_genes.add(g)

    all_genes = sorted(list(all_genes))

    rows = []
    cols = []
    data = []
    for ec, genes in tqdm.tqdm(ec2gene.items()):
        for g in genes:
            gid = all_genes.index(g)
            rows.append(ec)
            cols.append(gid)
            data.append(1)

    compatibility_matrix = sparse.coo_matrix((data, (rows, cols)) ).tocsr()
    return compatibility_matrix, all_genes


def do_em(ec_vector, y, n_iter=50):
    """
    just a wrapper around the EM which accounts for effective transcript length.
    we subsitute each transcript length as `1`, so no normlaization of transcript length!
    """
    eff_len = np.ones(y.shape[1])
    return do_em_efflen(ec_vector, y, eff_len, n_iter)

def alpha_to_rho(alpha, eff_lengths):
    """
    transform from probabilty of observing a read from a feature (alpha) to "relative abundance of that feature" (rho)
    """
    nom = alpha / eff_lengths
    rho = nom / nom.sum()
    return rho

def rho_to_alpha(rho, eff_lengths):
    """
    transform from "relative abundance of that feature" (rho) to probabilty of observing a read from a feature (alpha)
    """
    nom = rho * eff_lengths
    alpha = nom / nom.sum()
    return alpha

def loglike_fn_efflen(reads, y, rho, eff_len):
    """
    # Parameters:
    :param reads: vector of length == n_ECs; how many reads are consistent with the respective EC
    :param y: compatibility matrix (n_EC x nTranscripts or nGenes): 1 if EC_i is mappable to target_j
    :param rho: relative abundance of each transcript/gene
    :param eff_len: the ffective length of each transcript

    TODO: why does the length not enter here??!!

    Essentially calculates
    log[ (sum_k y_ik rho_k)**n_obs]
    """
    rho_norm = rho/ np.sum(rho)

    # the innter term: sum over all transcripts of an EC, weighted by the estimated abundance
    # this should be nEC x 1
    M = y @ sparse.diags(rho_norm).sum(1).A1

    # something silly: if rho is set to zero, logM will return NaN
    # usually what happened is that some transcript is not supported by a single read
    # hence the EM sets it to zero
    ix_discard = np.logical_and(reads==0, M==0)
    # now weight by each datapoint
    return np.sum(reads[~ix_discard] * np.log(M[~ix_discard]))

def do_em_efflen(ec_vector, y, eff_len, n_iter=50):
    """
    EM algorithm for abundace estimation of trasncripts from EC classes

    :param ec_vector: for each EC, how many counts are observed [np.ndarray]
    :param y: binary compatibility matrix (#EC x #transcripts)
    :param eff_len: effective length of each transcript [np.ndarray]
    :param n_iter: number of EM iterations
    """
    nEC, nGenes = y.shape

    assert len(ec_vector) == nEC
    assert len(eff_len) == nGenes
    loglike = []
    rho = np.ones(nGenes)/nGenes
    updates = []

    # some elements of ec_vector might be zero: Those will never influence the calculations
    # and can be removed (actually they screw things up numerically
    ix_zero = ec_vector == 0
    y = y[~ix_zero, :]
    ec_vector = ec_vector[~ix_zero]

    # print(y.shape, ec_vector.shape)
    for i in tqdm.trange(n_iter):

        # only keep every 50th update, it'll get to much otherwise
        if i % 50 == 0:
            updates.append(rho)
        loglike.append(loglike_fn_efflen(ec_vector, y, rho, eff_len))

        # how to split reads across transcripts with the current transcript frequencies
        splitting_matrix = y @ sparse.diags(rho / eff_len)
        splitting_matrix = splitting_matrix / splitting_matrix.sum(1)
        # print(splitting_matrix)

        # split the reads, caluclate new freuqneices
        rho_new = ec_vector @ splitting_matrix
        rho_new = rho_new/rho_new.sum()

        abs_change = np.abs(rho_new - rho).sum()
        # if i %100 == 0:
        #     print("Absolute change", abs_change)

        rho = rho_new

    # the final result
    updates.append(rho)
    updates = np.stack(updates)

    print("Absolute change", abs_change)
    return updates, loglike

