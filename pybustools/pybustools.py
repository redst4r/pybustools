from pybustools import busio
import pathlib
import toolz


class Bus():

    def __init__(self, folder,
                 bus_name='output.corrected.sort.bus',
                 ec_name='matrix.ec',
                 transcript_name='transcripts.txt',
                 decode_seq=True):
        self.folder = pathlib.Path(folder)
        self.bus_file = self.folder / bus_name
        self.ec_file = self.folder / ec_name
        self.transcript_file = self.folder / transcript_name
        self.decode_seq = decode_seq
        self.ec_dict = busio.read_matrix_ec(self.ec_file)
        self.transcript_dict = busio.read_transcripts(self.transcript_file)

    def iterate_bus(self):
        return busio.read_binary_bus(self.bus_file, self.decode_seq)

    def iterate_cells(self):
        return iterate_cells_of_busfile(self.bus_file, decode_seq=self.decode_seq)


def iterate_cells_of_busfile(fname, decode_seq=True):
    """
    runs over the !!!SORTED!!! busfile, collecting all entries for a single CB
    and yield it as `cb,info_list`
    """
    bus_iterator = busio.read_binary_bus(fname, decode_seq)

    # get the first entry to get started
    cb, umi, ec, count, flag = next(bus_iterator)
    current_cell = cb
    current_info = [(umi, ec, count, flag)]

    for cb, umi, ec, count, flag in bus_iterator:
        if cb > current_cell:
            # we're finished with one cells, yield it and start the next
            yield current_cell, current_info

            # reset for the next cell
            # process results and reset
            current_cell = cb
            current_info = [(umi, ec, count, flag)]
        elif cb == current_cell:
            current_info.append((umi, ec, count, flag))
        else:
            raise ValueError(f'Bus file not sorted!! {cb} vs {current_cell}')

    # emitting the final cell
    yield current_cell, current_info


def iterate_CB_UMI_of_busfile(fname, decode_seq=True):
    """
    iterates over CB/UMI entries, i.e. all entries with the same CB/UMI
    are emitted together.
    ideally, there'd only be one entry per CB/UMI, but sometimes thers doublets
    """
    bus_iterator = busio.read_binary_bus(fname, decode_seq)

    cb, umi, ec, count, flag = next(bus_iterator)
    current_cell = cb
    current_umi = umi
    current_info = [(ec, count, flag)]
    for cb, umi, ec, count, flag in bus_iterator:
        if cb > current_cell or (cb == current_cell and umi > current_umi):

            yield (current_cell, current_umi), current_info

            # reset for the next cell/UMI
            # process results and reset
            current_cell = cb
            current_umi = umi
            current_info = [(ec, count, flag)]
        elif cb == current_cell and umi == current_umi:
            current_info.append((ec, count, flag))
        else:
            raise ValueError(f'bsufile unsorted:  {cb}/{umi}  vs {current_cell}/{current_umi}')

    yield (current_cell, current_umi), current_info


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
            current_min = minimum_str(current_keys.values())
        else:
            break


def iterate_bus_cells_umi_multiple(names, fname_list, decode_seq=True):

    # a dict of all the busfile-iterators
    iterators = {n: iterate_CB_UMI_of_busfile(fname, decode_seq)
                 for n, fname in zip(names, fname_list)}

    for (cb, umi), info in merge_iterators(iterators):
        yield (cb, umi), info


def iterate_bus_cells_multiple(names, fname_list, decode_seq=True):
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
                 for n, fname in zip(names, fname_list)}

    for cb, info in merge_iterators(iterators):
        yield cb, info
