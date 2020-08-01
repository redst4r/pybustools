"""
parsimonous UMI graph ala salmon/alevin
"""
import collections
import itertools
import numpy as np
import networkx as nx
import Levenshtein
import pybktree
from pybustools.pybustools import Bus
from pybustools.busio import Bus_record, write_busfile
from pybustools.pybustools import iterate_cells_of_busfile_new
import tqdm
from disjoint_set import DisjointSet

"""
ds = DisjointSet()
for r in cb_records:
    the_set = ec_dict[r.EC]
    ds.find(the_set[0])
    for s in the_set[1:]:
        ds.union(the_set[0], s)

# group the records based on the connected component they fall into
grouping = collections.defaultdict(list)
for r in cb_records:
    the_set = ec_dict[r.EC]
    roots = [ds.find(_) for _ in the_set]
    assert np.all(np.array(roots) == roots[0])
    grouping[roots[0]].append(r)

# within each group, theres potential UMI that are duplicated. Their ECs might still not intersect directly!!
"""


class DisjointSets():

    def __init__(self):
        self.disjoint_sets = {}

    def add_set(self, aset, name):

        # look for disjoint sets that share a member with the current set
        candidate_setix = []
        for n, s in self.disjoint_sets.items():
            # if len(s & aset) > 0: # there;s a shared element
            if not s.isdisjoint(aset): # there;s a shared element
                candidate_setix.append(n)

        # now, anything thats in the candidate set (+aset iself)
        # for a new aggregated disjoint set!
        newset = aset
        # new_name = name
        new_name = (name,)  # turn into tuple
        for n in candidate_setix:
            newset = newset | self.disjoint_sets[n]
            # new_name = new_name + '_' + n
            new_name = new_name +  n
            del self.disjoint_sets[n]

        self.disjoint_sets[new_name] = newset

    def n_disjoint_sets(self):
        return len(self.disjoint_sets)


def hamming_distance(first, second):
    ''' returns the edit distance/hamming distances between
    its two arguements '''

    # dist = sum([not a == b for a, b in zip(first, second)])
    # return dist
    return Levenshtein.hamming(first, second)


def set_intersect(set_list):
    l = set_list[0]
    for s in set_list[1:]:
        l = l & s
    return l


def create_PUG_umi_based(cb_records:list, ec_dict):

    umi_dict = collections.defaultdict(list)
    for record in cb_records:
        umi_dict[record.UMI].append(record)

    # a BKTree of all UMIs in that cell
    tree = pybktree.BKTree(hamming_distance, list(umi_dict.keys()))

    nodes = set()
    edges = []
    for record in cb_records:

        nodes.add(record)
        # any sequence neighbours?
        for distance, umi_neighbor in tree.find(record.UMI, 1):
            # this particular UMI might have multiple records:
            for neighbor_record in umi_dict[umi_neighbor]:
                if record == neighbor_record:
                    continue  # due to d==0 this can be the same record
                # check EC overlap
                T1 = set(ec_dict[record.EC])
                T2 = set(ec_dict[neighbor_record.EC])
                e1 = (record, neighbor_record)
                e2 = (neighbor_record, record)
                if len(T1 & T2) > 0:
                    if record.COUNT > 2 * neighbor_record.COUNT - 1:
                        edges.append(e1)
                    elif neighbor_record.COUNT > 2 * record.COUNT - 1:
                        edges.append(e2)
                    else:
                        edges.append(e1)
                        edges.append(e2)
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return G


def create_PUG2(cb_records:list, ec_dict):
    """
    first partitions the records according to their ECs (and potnetial transcript overlap thereoff)
    then creates the PUG
    """
    ds = DisjointSet()
    for r in cb_records:
        the_set = ec_dict[r.EC]
        ds.find(the_set[0])
        for s in the_set[1:]:
            ds.union(the_set[0], s)

    # group the records based on the connected component they fall into
    grouping = collections.defaultdict(list)
    for r in cb_records:
        the_set = ec_dict[r.EC]
        roots = [ds.find(_) for _ in the_set]
        assert np.all(np.array(roots) == roots[0])
        grouping[roots[0]].append(r)

    nodes = set()
    edges = []
    for records in grouping.values():
        if len(records) == 1:
            nodes.add(records[0])
        else:
            for r1, r2 in itertools.combinations(records, 2):
                node1 = r1
                node2 = r2
                nodes.add(node1)
                nodes.add(node2)

                hd = hamming_distance(r1.UMI, r2.UMI)
                if hd > 1:
                    continue
                if set(ec_dict[r1.EC]).isdisjoint(ec_dict[r2.EC]):
                    continue

                if hd <= 1 and r1.COUNT > 2 * r2.COUNT - 1:
                    # print(r1, r2)
                    e = (node1, node2)
                    edges.append(e)
                elif hd <= 1 and r2.COUNT > 2 * r1.COUNT - 1:
                    e = (node2, node1)
                    edges.append(e)
                elif hd <= 1:
                    e1 = (node1, node2)
                    e2 = (node2, node1)
                    edges.append(e1)
                    edges.append(e2)
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return G


def create_PUG(cb_records:list, ec_dict):
    """
    take all the BUSrecords from a single cell and turn them into the PUG
    """
    raise ValueError("not correct: Two members of a Disjoint set might not overlap in EC!!")
    # group all these records according to if they have overlapping EC
    DS = DisjointSets()
    for r in cb_records:
        DS.add_set(set(ec_dict[r.EC]), name=r)

    nodes = set()
    edges = []
    for busrecord_tuple in DS.disjoint_sets.keys(): # all the different overlapping ecs

        # if its a single record as the member of the disjoint set, create
        # its own node, connect to nothing (it has no overlap with and other EC)
        if len(busrecord_tuple) == 1:
            nodes.add(busrecord_tuple[0])
        else:
            # with multiple member in the set, we have to check which ones are close in
            # UMI seqeuence space
            for r1, r2 in itertools.combinations(busrecord_tuple, 2):
                node1 = r1
                node2 = r2
                nodes.add(node1)
                nodes.add(node2)

                hd = hamming_distance(r1.UMI, r2.UMI)
                if hd <= 1 and r1.COUNT > 2 * r2.COUNT - 1:
                    # print(r1, r2)
                    e = (node1, node2)
                    edges.append(e)
                elif hd <= 1 and r2.COUNT > 2 * r1.COUNT - 1:
                    e = (node2, node1)
                    edges.append(e)
                elif hd <= 1:
                    e1 = (node1, node2)
                    e2 = (node2, node1)
                    edges.append(e1)
                    edges.append(e2)

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    return G


def PUG_to_records(G, inverse_ec):
    # now we have to find a set of EC that cover the entire CC
    # note that the records have been joined based on EC union, not intersect,
    # i.e. two elements in the same CC might not share a EC
    for nodes in nx.weakly_connected_components(G):
        if len(nodes) == 1:
            # just yield the single element of the set
            # yield list(nodes)[0]
            newr = list(nodes)[0]
        elif len(nodes) == 2:
            # these MUSt overlap in some transcrirpt
            n1, n2 = nodes
            overlapping_transcripts = frozenset(bus.ec_dict[n1.EC]) & set(bus.ec_dict[n2.EC])
            if overlapping_transcripts in inverse_ec:
                EC = inverse_ec[overlapping_transcripts]
                newr = Bus_record(n1.CB, n1.UMI, EC, n1.COUNT+n2.COUNT, n1.FLAG)
            else:
                EC = n1.EC if n1.COUNT > n2.COUNT else n2.EC
                UMI = n1.UMI if n1.COUNT > n2.COUNT else n2.UMI
                FLAG = n1.FLAG if n1.COUNT > n2.COUNT else n2.FLAG
                newr = Bus_record(n1.CB, UMI, EC, n1.COUNT+n2.COUNT, FLAG)

        else:
            overlapping_transcripts = frozenset(set_intersect([set(bus.ec_dict[n.EC]) for n in nodes]))
            if overlapping_transcripts in inverse_ec:
                EC = inverse_ec[overlapping_transcripts]
                umis = [_.UMI for _ in nodes]
                counts = [_.COUNT for _ in nodes]
                flags = [_.FLAG for _ in nodes]
                cbs = [_.CB for _ in nodes]
                newr = Bus_record(cbs[0], umis[0], EC, np.sum(counts), flags[0])
            else:  # theres either no overlapping EC, or just no EC correpsonding.
                # lets just take the EC what has the most counts
                L = list(nodes)
                counts = [_.COUNT for _ in L]
                ix = np.argmax(counts)
                the_node = L[ix]
                newr = Bus_record(the_node.CB, the_node.UMI, the_node.EC, np.sum(counts), the_node.FLAG)
        yield newr


def alevin_style_umi_correction(bus):
    """
    removing redundant bus records (same cell, close UMI and overlapping EC)
    from the bus
    """

    inverse_ec = {frozenset(transcript_list): ec for ec, transcript_list in bus.ec_dict.items()}

    for cb, record_list in tqdm.tqdm(iterate_cells_of_busfile_new(bus.bus_file, decode_seq=True)):

        # G = create_PUG_umi_based(record_list, bus.ec_dict)
        G = create_PUG2(record_list, bus.ec_dict)

        # iterate over the connected commponents and yield a set of records or each CC
        for r in PUG_to_records(G, inverse_ec):
            yield r


if __name__ == '__main__':

    bus = Bus('/home/mstrasse/mountSSD/kallisto_out/trimmed-kallisto/novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs/kallisto/sort_bus/bus_output')
    import tqdm
    for record in alevin_style_umi_correction(bus):
        pass

    for cb, record_list in iterate_cells_of_busfile_new(bus.bus_file, decode_seq=True):
        if len(record_list) > 1000: #  and len(record_list) < 1000:
            break


    G1 = create_PUG_umi_based(record_list, bus.ec_dict)
    G3 = create_PUG2(record_list, bus.ec_dict)
    # G2 = create_PUG(record_list, bus.ec_dict)
