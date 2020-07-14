from pybustools.pybustools import iterate_cells_of_busfile, iterate_CB_UMI_of_busfile
from pybustools import busio
import multiprocessing as mp
# import time
import random
import toolz
"""
for each file, we use a separate process to load bus entries from the file
and group the entries.

The processes put the resulting items (same as the generator items)
onto a queue (bascially replaceing the generator)
Another process consumes the items on the queues

see https://stackoverflow.com/questions/43078980/python-multiprocessing-with-generator/43299825

"""

TERMINATOR = None  # queue item signalling done


class ParallelGenerator(object):
    """docstring for ParallelCellGenerator."""

    def __init__(self, busfile_dict, decode_seq, queue_size):
        super(ParallelGenerator, self).__init__()
        self.busfile_dict = busfile_dict
        self.queues = {name: mp.Queue(queue_size) for name in busfile_dict.keys()}
        self.tasks = {}
        self.decode_seq = decode_seq

    def start_queues(self):
        """
        starts the processes that fill up the queues
        """
        raise NotImplementedError('please derived from this class and implement')

    def print_queue_load(self):
        """
        print the number of items per queue currently loaded
        """
        for name, q in self.queues.items():
            print(f'Queue {name}: {q.qsize()} items')

    def iterate(self):
        """
        pulls items from the queues and emits them such that entries from the
        same group (e.g. same CB/UMI combination) are emitted together
        """
        # the CB/UMI iterator will yield tuples of CB,UMI
        # luckily python does the righth thing for tuple comparison,
        # i.e. > < work
        def minimum_str(str_list):
            return toolz.reduce(lambda x, y: x if x < y else y, str_list)

        # first elements:
        elements = {}
        for n, queue in self.queues.items():
            bus_group_tuple, info = _helper_queue_iterator(queue)
            elements[n] = (bus_group_tuple, info)
        current_cbs = toolz.valmap(lambda group_and_info: group_and_info[0], elements)
        current_min = minimum_str(current_cbs.values())  # the smallest CB/UMI tuple

        while len(self.queues) > 0:
            # whichever iterators have the minimum value:
            to_emit = []  # record the names of iterators that will emit an item
            for n, bus_group_tuple in current_cbs.items():
                if bus_group_tuple == current_min:
                    to_emit.append(n)

            # emit all candidates,
            emit_infos = {}
            for candidate_name in to_emit:
                bus_group_tuple, info = elements[candidate_name]
                emit_infos[candidate_name] = info

            yield current_min, emit_infos

            #  advance their iterators
            for candidate_name in to_emit:
                try:
                    # bus_group_tuple, info = next(iterators[candidate_name])
                    bus_group_tuple, info = _helper_queue_iterator(self.queues[candidate_name])
                    elements[candidate_name] = (bus_group_tuple, info)
                except StopIteration:
                    # cleaning up finished queues and jobs
                    self.cleanup(candidate_name)
                    del self.queues[candidate_name]
                    del self.tasks[candidate_name]

                    del current_cbs[candidate_name]
                    del elements[candidate_name]

            # new minimum!
            if len(self.queues) > 0:
                current_cbs = toolz.valmap(lambda group_and_info: group_and_info[0], elements)
                current_min = minimum_str(current_cbs.values())
            else:
                break

    def cleanup(self, samplename):
        "close the queue and task for a specific job"
        assert self.queues[samplename].qsize() == 0, "trying to close nonempty queue"
        self.queues[samplename].close()
        self.tasks[samplename].close()


class ParallelCellGenerator(ParallelGenerator):
    """docstring for ParallelCellGenerator."""

    def __init__(self, busfile_dict, decode_seq, queue_size):
        super(ParallelCellGenerator, self).__init__(busfile_dict, decode_seq, queue_size)

    def start_queues(self):
        # connect the queues to the files
        # this will start a process for each file, populating the queue
        for name, fname in self.busfile_dict.items():
            task = mp.Pool(1, initializer=cell_producer, initargs=(fname, self.queues[name], self.decode_seq))
            self.tasks[name] = task


class ParallelCellUMIGenerator(ParallelGenerator):
    """docstring for ParallelCellGenerator."""

    def __init__(self, busfile_dict, decode_seq, queue_size):
        super(ParallelCellUMIGenerator, self).__init__(busfile_dict, decode_seq, queue_size)

    def start_queues(self):
        # connect the queues to the files
        # this will start a process for each file, populating the queue
        for name, fname in self.busfile_dict.items():
            task = mp.Pool(1, initializer=cell_umi_producer, initargs=(fname, self.queues[name], self.decode_seq))
            self.tasks[name] = task


def cell_producer(fname, out_queue, decode_seq):
    """
    turns iterate_cells_of_busfile into a queue
    """
    for record in iterate_cells_of_busfile(fname, decode_seq):
        # print('Putting record:', fname, record[0])
        out_queue.put(record)
        # print('sleeping')
        # time.sleep(random.randrange(1,3))
    out_queue.put(TERMINATOR)


def cell_umi_producer(fname, out_queue, decode_seq):
    """
    turns iterate_CB_UMI_of_busfile into a queue
    """
    for record in iterate_CB_UMI_of_busfile(fname, decode_seq):
        # print('Putting record:', fname, record[0])
        out_queue.put(record)
        # print('sleeping')
        # time.sleep(random.randrange(1,3))
    out_queue.put(TERMINATOR)


def _helper_queue_iterator(queue):
    """
    iterate over the queue until we hit the stop item
    """
    result = queue.get()
    if result == TERMINATOR:
        raise StopIteration
    else:
        return result


def parallel_iterate_bus_cells_pairs_basic(fname1, fname2, decode_seq=True):
    """
    opposed to iterate_bus_cells_pairs, this yields an entry for each CB
    (wether its in both or only a single file).
    If the CB is not present in the other file, it'll yield `None`
    """
    raise ValueError('deprecated')
    QUEUE_LENGTH = 1
    queue1 = mp.Queue(QUEUE_LENGTH)
    queue2 = mp.Queue(QUEUE_LENGTH)

    task1 = mp.Pool(1, initializer=cell_producer, initargs=(fname1, queue1, True))
    task2 = mp.Pool(1, initializer=cell_producer, initargs=(fname2, queue2, True))

    try:
        # get it started outside the loop
        cb1, info1 = _helper_queue_iterator(queue1)
        cb2, info2 = _helper_queue_iterator(queue2)
        while True:
            print(f'Q1 {queue1.qsize()}, Q2 {queue1.qsize()}')
            if cb1 == cb2:
                yield cb1, info1, info2
                # advancing both iterators
                cb1, info1 = _helper_queue_iterator(queue1)
                cb2, info2 = _helper_queue_iterator(queue2)
            elif cb1 > cb2:
                yield cb2, None, info2
                # get the next cell in I2
                cb2, info2 = _helper_queue_iterator(queue2)

            elif cb2 > cb1:
                yield cb1, info1, None
                cb1, info1 = _helper_queue_iterator(queue1)

            else:
                raise ValueError('cant happen')
    # the next() will throw an exception if one generator runs out
    # thats the signal that we're done with the pairs
    except StopIteration:
        print('One iterator finished!')
        queue1.close()
        queue2.close()
        task1.close()
        task2.close()
        return


def random_busrecord(cb_length, umi_length, ngenes):
    d = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    cb = [random.randrange(0, 4) for _ in range(cb_length)]
    cb = ''.join([d[_] for _ in cb])
    umi = [random.randrange(0, 4) for _ in range(umi_length)]
    umi = ''.join([d[_] for _ in umi])
    gene = random.randrange(0, ngenes)
    counts = random.randint(0, 50)

    return busio.Bus_record(cb, umi, gene, counts, 1)


def random_buslist(n_records, cb_length, umi_length, ngenes):
    records = sorted([random_busrecord(cb_length, umi_length, ngenes)
                      for _ in range(n_records)])
    return records


if __name__ == '__main__':

    fname1 = '/tmp/some1.bus'
    fname2 = '/tmp/some2.bus'

    records1 = random_buslist(500, cb_length=4, umi_length=5, ngenes=10)
    records2 = random_buslist(500, cb_length=4, umi_length=5, ngenes=10)
    busio.write_busfile(fname1, records1, cb_length=4, umi_length=5)
    busio.write_busfile(fname2, records2, cb_length=4, umi_length=5)


    PCG = ParallelCellGenerator({'sample1': fname1, 'sample2': fname2}, decode_seq=True, queue_size=3)
    PCG.start_queues()

    results = {cb: info for cb, info in PCG.iterate()}

    from pybustools.pybustools import iterate_bus_cells_multiple
    results_serial = {cb: info for cb,info in iterate_bus_cells_multiple(['sample1', 'sample2'], [fname1, fname2])}
