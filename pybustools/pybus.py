import pathlib
import toolz
from pybustools import busio

# TODO: multiple dispatch instead if inheritance?
class _Bus():
    """
    Main class to represent a kallisto-Bus file and its associated metadata (EC and transcript files)
    """

    def __init__(self,
                 bus_name,
                 ec_name,
                 transcript_name,
                 decode_seq=True):
        self.bus_file = bus_name
        self.ec_file = ec_name
        self.transcript_file = transcript_name
        self.decode_seq = decode_seq
        self.ec_dict = busio.read_matrix_ec(self.ec_file)
        self.transcript_dict = busio.read_transcripts(self.transcript_file)

    def iterate_bus(self):
        return busio.read_binary_bus(self.bus_file, self.decode_seq)

    def iterate_cells(self):
        return iterate_cells_of_busfile(self.bus_file, decode_seq=self.decode_seq)

    def resolve_EC_to_transcripts(self, EC):
        "which transcript IDs correspond to the Equivalence class"
        transcripts = [self.transcript_dict[_] for _ in self.ec_dict[EC]]
        return transcripts


class Bus(_Bus):
    """
    this is just a convenience wrapper around _Bus where we only specify
    the bustools folder and it infers all the filenames
    """

    def __init__(self, folder,
                 bus_name='output.corrected.sort.bus',
                 ec_name='matrix.ec',
                 transcript_name='transcripts.txt',
                 decode_seq=True):

        folder = pathlib.Path(folder)
        super().__init__(
            bus_name=str(folder / bus_name),
            ec_name=str(folder / ec_name),
            transcript_name=str(folder / transcript_name),
            decode_seq=decode_seq
        )


def records_to_gene(records: list, ec2g) -> list:
    """
    is there one or more genes that are consistent across all the records of this CB/UMI
    returns the set of genes that are shared by all the records

    :param records: list of BusRecords, usually from a single CB/UMI
    :param ec2g: a dictionary mapping each EC to a list of genes

    :returns: a list of genes that are compatible with the records (i.e. each record supports all those genes)
    """
    if len(records) == 1:
        genelist = ec2g[records[0].EC]
        assert isinstance(genelist, list), "ec2g didnt return a list"
        return genelist

    genes_per_record = []
    for r in records:
        genelist = ec2g[r.EC]
        assert isinstance(genelist, list), "ec2g didnt return a list"
        genes_per_record.append(set(genelist))

    # si there an intersection?
    intersect = set.intersection(*genes_per_record)
    return list(intersect)


def iterate_cells_of_busfile(fname, decode_seq=True):
    """
    runs over the !!!SORTED!!! busfile, collecting all entries for a single CB
    and yield it as `cb,info_list`

    this one returns a list of BusRecords
    """
    bus_iterator = busio.read_binary_bus(fname, decode_seq)

    # get the first entry to get started
    record = next(bus_iterator)
    current_cell = record.CB
    current_recordlist = [record]

    for record in bus_iterator:
        if record.CB > current_cell:
            # we're finished with one cells, yield it and start the next
            yield current_cell, current_recordlist

            # reset for the next cell
            # process results and reset
            current_cell = record.CB
            current_recordlist = [record]
        elif record.CB == current_cell:
            current_recordlist.append(record)
        else:
            raise ValueError(f'Bus file not sorted!! {record.CB} vs {current_cell}')

    # emitting the final cell
    yield current_cell, current_recordlist


def iterate_CB_UMI_of_busfile(fname, decode_seq=True):
    """
    iterates over CB/UMI entries, i.e. all entries with the same CB/UMI
    are emitted together.
    ideally, there'd only be one entry per CB/UMI, but sometimes thers doublets:
    This happens when two reads have the same CB/UMI but map to different genes (ECs)
    """
    bus_iterator = busio.read_binary_bus(fname, decode_seq)

    record = next(bus_iterator)
    current_cell = record.CB
    current_umi = record.UMI
    current_recordlist = [record]
    for record in bus_iterator:
        if record.CB > current_cell or (record.CB == current_cell and record.UMI > current_umi):

            yield (current_cell, current_umi), current_recordlist

            # reset for the next cell/UMI
            # process results and reset
            current_cell = record.CB
            current_umi = record.UMI
            current_recordlist = [record]
        elif record.CB == current_cell and record.UMI == current_umi:
            current_recordlist.append(record)
        else:
            raise ValueError(f'bsufile unsorted:  {record.CB}/{record.UMI}  vs {current_cell}/{current_umi}')

    yield (current_cell, current_umi), current_recordlist


def merge_iterators(dict_of_iterators):
    """
    given a dictionary  iname->Iterator
    where each iterator yields a key,value pair,
    this function aggregates items across iterators that have the same key
    {
        I1 : ['A': 1, 'B':2, 'C':3]
        I2 : ['A': 10, 'B':20, 'D':30]
    }

    will become:
    'A': {I1: 1, I2:10},
    'B': {I1: 2, I2:20}
    'C': {I1: 3}
    'D': {I2: 30}

    This assumes the ITERATORS ARE KEY SORTED!!!
    Keys also must be strings
    """
    def minimum_str(str_list):
        return toolz.reduce(lambda x, y: x if x < y else y, str_list)

    # first elements:
    elements = {}
    for n, it in dict_of_iterators.items():
        key, value = next(it)
        elements[n] = (key, value)
    current_keys = toolz.valmap(lambda keyvalue: keyvalue[0], elements)
    current_min = minimum_str(current_keys.values())  # the smallest CB/UMI tuple

    while len(dict_of_iterators) > 0:
        # whichever iterators have the minimum value:
        to_emit = []  # record the names of iterators that will emit an item
        for n, key in current_keys.items():
            if key == current_min:
                to_emit.append(n)

        # emit all candidates,
        emit_values = {}
        for candidate_name in to_emit:
            key, value = elements[candidate_name]
            emit_values[candidate_name] = value

        yield current_min, emit_values

        #  advance their iterators
        for candidate_name in to_emit:
            try:
                key, value = next(dict_of_iterators[candidate_name])
                elements[candidate_name] = (key, value)
            except StopIteration:
                del dict_of_iterators[candidate_name]
                del current_keys[candidate_name]
                del elements[candidate_name]

        # new minimum!
        if len(dict_of_iterators) > 0:
            current_keys = toolz.valmap(lambda key_value: key_value[0], elements)
            new_min = minimum_str(current_keys.values())
            assert new_min > current_min, "some iterator is not sorted!"
            current_min = new_min
        else:
            break


def iterate_bus_cells_umi_multiple(fname_dict: dict, decode_seq=True):
    """
    simultaniously iterating multiple busfiles, yielding groups of
    records with the same CB/UMI in the busfiles

    :param bus_dict: a dictionary of name->bus_filename. Each entry is one busfile to be iterated. Names can be arbitrary, but will be used in the yielded results
    :param decode_seq: if True, decodes the int represnetation of the CB/UMI into bases
    :return: yields (cb, umi) and a info_dict (name->list[records])
    """
    # a dict of all the busfile-iterators
    iterators = {n: iterate_CB_UMI_of_busfile(fname, decode_seq)
                 for n, fname in fname_dict.items()}

    for (cb, umi), info in merge_iterators(iterators):
        yield (cb, umi), info


def iterate_bus_cells_multiple(fname_dict: dict, decode_seq=True):
    """
    iterates over multiple busfiles, emitting an entries grouped by CB:
    If a CB is present in multiple bus-files, this will yield:
        cb, {n1: [umi, ...,
                  umi ,...]
             n2: [umi,....
                  umi,....]
        }
    """

    # a dict of all the busfile-iterators
    iterators = {n: iterate_cells_of_busfile(fname, decode_seq)
                 for n, fname in fname_dict.items()}

    for cb, info in merge_iterators(iterators):
        yield cb, info
