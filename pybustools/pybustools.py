from pybustools import busio, busio_old
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
        return iterate_cells_of_busfile(self.bus_file, is_binary=True, decode_seq=self.decode_seq)


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
        if cb != current_cell:
            # we're finished with one cells, yield it and start the next
            yield current_cell, current_info

            # reset for the next cell
            # process results and reset
            current_cell = cb
            current_info = [(umi, ec, count, flag)]
        else:
            current_info.append((umi, ec, count, flag))

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
        if cb != current_cell or umi != current_umi:
            yield (current_cell, current_umi), current_info

            # reset for the next cell/UMI
            # process results and reset
            current_cell = cb
            current_umi = umi
            current_info = [(ec, count, flag)]

        else:
            current_info.append((ec, count, flag))

    yield (current_cell, current_umi), current_info


# def iterate_bus_cells_pairs(fname1, fname2, is_binary=True, decode_seq=True):
#     """
#     busfiles must be sorted!!

#     iterate over the files until we have a matchign pair of cells
#     collect the info on both cells and advance to the next pair

#     the pairing is tricky, given two CBs cb1, cb2:
#     - if cb1 > cb2 : advance the 2nd iterator (since cb1 is in the future of that 2nd iterator)
#     - if cb2 > cb1 : advance the 1st iterator (since cb2 is in the future of that 1nd iterator)
#     """
#     I1 = iterate_cells_of_busfile(fname1, is_binary, decode_seq)
#     I2 = iterate_cells_of_busfile(fname2, is_binary, decode_seq)

#     try:
#         # get it started outside the loop
#         cb1, info1 = next(I1)
#         cb2, info2 = next(I2)
#         while True:
#             if cb1 == cb2:
#                 yield cb1, info1, info2
#                 # advancing both iterators
#                 cb1, info1 = next(I1)
#                 cb2, info2 = next(I2)
#             elif cb1 > cb2:
#                 # get the next cell in I2
#                 cb2, info2 = next(I2)
#             elif cb2 > cb1:
#                 cb1, info1 = next(I1)
#             else:
#                 raise ValueError('cant happen')
#     # the next() will throw an exception if one generator runs out
#     # thats the signal that we're done with the pairs
#     except StopIteration:
#         print('One iterator finished!')
#         return


# def iterate_bus_cells_pairs_basic(fname1, fname2, is_binary=True, decode_seq=True):
#     """
#     opposed to iterate_bus_cells_pairs, this yields an entry for each CB
#     (wether its in both or only a single file).
#     If the CB is not present in the other file, it'll yield `None`


#     WARNING: this skipps the last elements if one iterator finishes
#     """
#     I1 = iterate_cells_of_busfile(fname1, is_binary, decode_seq)
#     I2 = iterate_cells_of_busfile(fname2, is_binary, decode_seq)

#     try:
#         # get it started outside the loop
#         cb1, info1 = next(I1)
#         cb2, info2 = next(I2)
#         while True:
#             if cb1 == cb2:
#                 yield cb1, info1, info2
#                 # advancing both iterators
#                 cb1, info1 = next(I1)
#                 cb2, info2 = next(I2)
#             elif cb1 > cb2:
#                 yield cb2, None, info2
#                 # get the next cell in I2
#                 cb2, info2 = next(I2)
#             elif cb2 > cb1:
#                 yield cb1, info1, None
#                 cb1, info1 = next(I1)
#             else:
#                 raise ValueError('cant happen')
#     # the next() will throw an exception if one generator runs out
#     # thats the signal that we're done with the pairs
#     except StopIteration:
#         print('One iterator finished!')
#         return


def iterate_bus_cells_umi_multiple(names, fname_list, is_binary=True, decode_seq=True):

    def minimum_str(str_list):
        return toolz.reduce(lambda x, y: x if x < y else y, str_list)

    # a dict of all the busfile-iterators
    iterators = {n: iterate_CB_UMI_of_busfile(fname, is_binary, decode_seq)
                 for n, fname in zip(names, fname_list)}

    # the CB/UMI iterator will yield tuples of CB,UMI
    # luckily python does the righth thing for tuple comparison, i.e. > < work

    # first elements:
    elements = {}
    for n, it in iterators.items():
        cb_umi, info = next(it)
        elements[n] = (cb_umi, info)
    current_cbs = toolz.valmap(lambda cb_and_info: cb_and_info[0], elements)
    current_min = minimum_str(current_cbs.values())  # the smallest CB/UMI tuple

    while len(iterators) > 0:
        # whichever iterators have the minimum value:
        to_emit = []  # record the names of iterators that will emit an item
        for n, cb_umi in current_cbs.items():
            if cb_umi == current_min:
                to_emit.append(n)

        # emit all candidates,
        emit_infos = {}
        for candidate_name in to_emit:
            cb_umi, info = elements[candidate_name]
            emit_infos[candidate_name] = info

        yield current_min, emit_infos

        #  advance their iterators
        for candidate_name in to_emit:
            try:
                cb_umi, info = next(iterators[candidate_name])
                elements[candidate_name] = (cb_umi, info)
            except StopIteration:
                del iterators[candidate_name]
                del current_cbs[candidate_name]
                del elements[candidate_name]

        # new minimum!
        if len(iterators) > 0:
            current_cbs = toolz.valmap(lambda cb_and_info: cb_and_info[0], elements)
            current_min = minimum_str(current_cbs.values())
        else:
            break


def iterate_bus_cells_multiple(names, fname_list, is_binary=True, decode_seq=True):
    """
    iterates over multiple busfiles, emitting an entries grouped by CB:
    If a CB is present in multiple bus-files, this will yield:
        cb, {n1: [umi, ...,
                  umi ,...]
             n2: [umi,....
                  umi,....]
        }
    """

    def minimum_str(str_list):
        return toolz.reduce(lambda x, y: x if x < y else y, str_list)

    # a dict of all the busfile-iterators
    iterators = {n: iterate_cells_of_busfile(fname, is_binary, decode_seq)
                 for n, fname in zip(names, fname_list)}

    # first elements:
    elements = {}
    for n, it in iterators.items():
        cb, info = next(it)
        elements[n] = (cb, info)
    current_cbs = toolz.valmap(lambda cb_and_info: cb_and_info[0], elements)
    current_min = minimum_str(current_cbs.values())

    while len(iterators) > 0:
        # whichever iterators have the minimum value:
        to_emit = []  # record the names of iterators that will emit an item
        for n, cb in current_cbs.items():
            if cb == current_min:
                to_emit.append(n)

        # emit all candidates,
        emit_infos = {}
        for candidate_name in to_emit:
            cb, info = elements[candidate_name]
            emit_infos[candidate_name] = info

        yield current_min, emit_infos
        # print(current_min, emit_infos)

        #  advance their iterators
        for candidate_name in to_emit:
            try:
                cb, info = next(iterators[candidate_name])
                elements[candidate_name] = (cb, info)
            except StopIteration:
                del iterators[candidate_name]
                del current_cbs[candidate_name]
                del elements[candidate_name]

        # new minimum!
        if len(iterators) > 0:
            current_cbs = toolz.valmap(lambda cb_and_info: cb_and_info[0], elements)
            current_min = minimum_str(current_cbs.values())
        else:
            break


if __name__ == '__main__':
    pass
