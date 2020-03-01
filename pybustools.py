from busio import read_binary_bus, read_text_bus

def iterate_cells_of_busfile(fname, is_binary=True):
    """
    runs over the !!!SORTED!!! busfile, collecting all entries for a single CB
    and yield it as `cb,info_list`
    """
    if is_binary:
        bus_iterator = read_binary_bus(fname)
    else:
        bus_iterator = read_text_bus(fname)

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
    return


# def iterate_cells_of_busfile_old(fname):
#     """
#     runs over the !!!SORTED!!! busfile, collecting all entries for a single CB
#     and yield it as `cb,info_list`
#     """
#     with open(fname, 'r') as fh:
#         # first line to get going
#         line = fh.readline()
#         cb, umi, ec, count = read_text_bus_entry(line)
#
#         current_cell = cb
#         current_info = [(umi, ec, count)]
#
#         for line in fh:
#             cb, umi, ec, count = read_text_bus_entry(line)
#
#             if cb != current_cell:
#
#                 # we're finished with one cells, yield it and start the next
#                 yield current_cell, current_info
#
#                 # reset for the next cell
#                 # process results and reset
#                 current_cell = cb
#                 current_info = [(umi, ec, count)]
#
#             else:
#                 current_info.append((umi, ec, count))
#
#         # emitting the final cell
#         yield current_cell, current_info


def iterate_bus_cells_pairs(fname1, fname2, is_binary=True):
    """
    bustfiles must be sorted!!

    iterate over the files until we have a matchign pair of cells
    collect the info on both cells and advance to the next pair

    the pairing is tricky, given two CBs cb1, cb2:
    - if cb1 > cb2 : advance the 2nd iterator (since cb1 is in the future of that 2nd iterator)
    - if cb2 > cb1 : advance the 1st iterator (since cb2 is in the future of that 1nd iterator)
    """
    I1 = iterate_cells_of_busfile(fname1, is_binary)
    I2 = iterate_cells_of_busfile(fname2, is_binary)

    try:
        # get it started outside the loop
        cb1, info1 = next(I1)
        cb2, info2 = next(I2)
        while True:
            if cb1 == cb2:
                yield cb1, info1, info2
                # advancing both iterators
                cb1, info1 = next(I1)
                cb2, info2 = next(I2)
            elif cb1 > cb2:
                # get the next cell in I2
                cb2, info2 = next(I2)
            elif cb2 > cb1:
                cb1, info1 = next(I1)
            else:
                raise ValueError('cant happen')
    # the next() will throw an exception if one generator runs out
    # thats the signal that we're done with the pairs
    except StopIteration:
        print('One iterator finished!')
        return
