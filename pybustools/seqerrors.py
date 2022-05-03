import tqdm
import collections
import toolz
from pybustools.busio import Bus_record, get_header_info, _encode_ACGT_to_int, _decode_int_to_ACGT
from pybustools.pybustools import iterate_CB_UMI_of_busfile
from pybustools.butterfly import CUHistogram
import numpy as np
import random


"""
some code to artificially introduce errors into the UMIs
and check their effect on the CU histograms (unseen species)
"""
def _group_busrecords_by_cb_umi_ec(record_list):
    new_list = {} #indexed by cb,umi,ec
    for r in record_list:
        CUE = (r.CB,r.UMI,r.EC)
        if CUE not in new_list: # create a new record
            new_list[CUE] = r
        else:
            # update the counts of the existing record
            c = new_list[CUE]
            new_record = Bus_record(c.CB, c.UMI, c.EC, c.COUNT+r.COUNT, c.FLAG)
            new_list[CUE] = new_record
    return list(new_list.values())

def errorneous_CB_UMI_iterator(bus_file, R1error, decode_seq):
 
    """ 
    some problems:
        - this yields as  (CB,UMI). Errors in the CB will lead to out of order CBs
        - CB/UMI errors are lumped together
    """
    version, cb_len, umi_len, tlen = get_header_info(bus_file)

    for (cb, umi), record_list in iterate_CB_UMI_of_busfile(bus_file, decode_seq):
        n_reads = sum(r.COUNT for r in record_list)

        n_errorneous_reads = np.random.binomial(n_reads, R1error)
        if n_errorneous_reads == 0:
            yield (cb, umi), record_list
        else:
            unfolded_records = []
            for r in record_list:
                # unfold the reads
                unfolded_records.extend([Bus_record(r.CB, r.UMI, r.EC, 1, r.FLAG)]*r.COUNT)

            # first yield the undamaged ones
            # this shuffles in place
            random.shuffle(unfolded_records)
            #  unfolded_records = random.choices(unfolded_records, k=len(unfolded_records))  # doesnt work, samples with REPLACEMENT
            #  unfolded_records = np.random.permutation(unfolded_records) # doesnt work, turns it into an array
            n_correct_reads = n_reads - n_errorneous_reads
            # groub by CB,UMI,EC
            correct_reads = _group_busrecords_by_cb_umi_ec(unfolded_records[:n_correct_reads])
            if n_correct_reads>0:
#                 print('yielding correct')
                yield (cb, umi), correct_reads
            
            # now the damaged ones
            # note that we just add an artifical base X. In theory, that CB/UMI could mutate into an already existing one, adding a count there
            # but the sequential way of the busfile-iterator cant hanlde that anyway. As long as the error rate is small this shouldt matter!
            # 
            # also, to introduce the error, we have to decode/encode the CB/UMI (if the iterator doesnt decode already)
            damaged_reads = []
            for r in unfolded_records[n_correct_reads:]: # all the ones with mistakes in R1
                if not decode_seq:
                    cb_tmp = _decode_int_to_ACGT(cb, cb_len)
                    umi_tmp = _decode_int_to_ACGT(umi, umi_len)
                    cb_umi = cb_tmp+umi_tmp
                else:
                    cb_umi = cb+umi

                cb_umi = _mutate_randomly(cb_umi)
                cb_new, umi_new = cb_umi[:cb_len], cb_umi[cb_len:]

                if not decode_seq: # reencode
                    cb_new = _encode_ACGT_to_int(cb_new)
                    umi_new = _encode_ACGT_to_int(umi_new)

                damaged_reads.append(Bus_record(cb_new, umi_new, r.EC, r.COUNT, r.FLAG))

            damaged_reads = _group_busrecords_by_cb_umi_ec(damaged_reads) # this should rarely group anything!

            for r in damaged_reads:
#                 print('yielding damaged')
                yield (r.CB, r.UMI), [r]

def _mutate_randomly(seq):
    error_position = np.random.choice(len(seq))
    current_base = seq[error_position]
    subs = set(['A','T','G','C']) - set(current_base) # any of the other 3 bases
    new_base = np.random.choice(list(subs), 1)[0]
    seq = seq[0:error_position] + new_base + seq[error_position+1: ] # seq[error_position] = 'X' doesnt work
    return seq

def make_ec_histograms_adding_umi_error(bus, t2gfile, R1error):
    """
    :param R1error: the probability of making an error in the R1 reads (this is not per base!)
    """
    ec_hists = {}
    I = errorneous_CB_UMI_iterator(bus.bus_file, R1error, bus.decode_seq)
    # filtering records that map to more than one EC
    I = (record_list[0] for (cb, umi), record_list in I if len(record_list) == 1)
    counter = 0
    for r in tqdm.tqdm(I):
        counter += 1
        if not r.EC in ec_hists:
            ec_hists[r.EC] = collections.defaultdict(int)

        ec_hists[r.EC][r.COUNT] += 1
    ec_hists = toolz.valmap(CUHistogram, ec_hists)

    return ec_hists

if False:
    from pybustools import busio, pybustools
    import pytest


    # some testing
    records = [
        busio.Bus_record('ATAT', 'AAA', 10, 5, 1),
        busio.Bus_record('ATAT', 'GGG', 10, 3, 1),
        busio.Bus_record('TAGA', 'TAT', 14, 2, 13),
        busio.Bus_record('TTAT', 'AAA', 13, 100, 12),
    ]
    fname = '/tmp/some.bus'
    busio.write_busfile(fname, records, cb_length=4, umi_length=3)
    n_umi = 0
    n_read = 0
    for (cb, umi), records in errorneous_CB_UMI_iterator(fname, R1error=0.5, decode_seq=True):
        n_umi += 1
        n_read += sum(r.COUNT for r in records)
    n_umi,n_read
