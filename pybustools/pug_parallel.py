import multiprocessing as mp
import collections
import itertools
import os
from pybustools.pybus import Bus
from pybustools.busio import Bus_record, write_busfile
from pybustools.pybus import iterate_cells_of_busfile_new
import tqdm
from pybustools import busio
import tempfile

from pybustools.pug import create_PUG_umi_based, PUG_to_records, create_PUG2

TERMINATOR = None  # queue item signalling done


def cell_producer(bus, out_queue, n_workers):
    """
    turns iterate_cells_of_busfile into a queue
    """
    I = iterate_cells_of_busfile_new(bus.bus_file, decode_seq=True)
    # I = toolz.take(10000, I)
    for cb, record_list in tqdm.tqdm(I):
        # workaround since pickle doesnt like namedtuples
        record_list = [tuple(_) for _ in record_list]
        out_queue.put((cb, record_list))
        # print('sleeping')
        # time.sleep(random.randrange(1,3))

    # each worker needs to see/consume a terminator, otherwise they wont stop
    for i in range(n_workers):
        out_queue.put(TERMINATOR)


def _helper_queue_iterator(queue):
    """
    iterate over the queue until we hit the stop item
    """
    result = queue.get()
    if result == TERMINATOR:
        raise StopIteration
    else:
        cb, record_list = result
        record_list = [Bus_record(*_) for _ in record_list]

        return cb, record_list


def _gen(bus, queue):
    """
    picks elemets from the queue (CB, [list of records]) and turns them
    into a generator of new records (filtering UMI duplicates etc)
    """
    inverse_ec = {frozenset(transcript_list): ec for ec, transcript_list in bus.ec_dict.items()}

    try:
        while True:
            cb, record_list = _helper_queue_iterator(queue)
            # G = create_PUG_umi_based(record_list, bus.ec_dict)
            G = create_PUG2(record_list, bus.ec_dict)
            # iterate over the connected commponents and yield a set of records or each CC
            for r in PUG_to_records(G, inverse_ec):
                yield r
    except StopIteration:
        print('iterator finished!')
        return


def pug_writer(bus, queue, outfile):
    """
    wrapper for cont. taking an item from the queue and writing it to a busfile
    """
    generator = _gen(bus, queue)
    write_busfile(outfile, generator, cb_length=16, umi_length=12)


def in_parallel(bus, outfile, cores):
    """

    """
    # put al the intermediate results in here
    tmpfolder = tempfile.mkdtemp(prefix='pug_', dir='/tmp')

    QUEUE_LENGTH = 1000
    cb_queue = mp.Queue(QUEUE_LENGTH)  # has (CB, [records]) as elements

    worker_tasks = []

    # task1 = mp.Pool(1, initializer=cell_producer, initargs=(bus, cb_queue))
    task1 = mp.Process(target=cell_producer, args=(bus, cb_queue, cores))
    worker_tasks.append(task1)
    for i in range(cores):
        busfile = f'{tmpfolder}/{i}.bus'
        # t = mp.Pool(1, initializer=pug_writer, initargs=(bus, cb_queue, busfile))
        t = mp.Process(target=pug_writer, args=(bus, cb_queue, busfile))
        worker_tasks.append(t)

    for t in worker_tasks:
        t.start()

    for t in worker_tasks:
        t.join()

    # cleanup
    assert cb_queue.empty(), "queue not empty!!"
    cb_queue.close()

    for t in worker_tasks:
        t.close()

    # merge all busfiles into a big one!
    bus_iterators = []
    for i in range(cores):
        busfile = f'{tmpfolder}/{i}.bus'
        gen = busio.read_binary_bus(busfile)
        bus_iterators.append(gen)

    big_gen = itertools.chain.from_iterable(bus_iterators)
    unsorted_outfile = f'{tmpfolder}/unsorted.bus'
    write_busfile(unsorted_outfile, big_gen, cb_length=16, umi_length=12)

    # cleanup the parts
    for i in range(cores):
        busfile = f'{tmpfolder}/{i}.bus'
        os.remove(busfile)

    # sort the file
    import subprocess
    import sys
    ret = subprocess.run(["bustools", "sort", '-o', outfile, unsorted_outfile])
    if ret.returncode != 0:
        print("Child was terminated by signal", ret, file=sys.stderr)
        raise ValueError()

    # os.remove(busfile)

# bus = Bus('/home/mstrasse/mountSSD/kallisto_out/trimmed-kallisto/novaseq.190919_A00266_0278_BHFJY5DRXX_ParkFoulkesCRUK.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-3_0_0.outs/kallisto/sort_bus/bus_output')
in_parallel(bus, '/tmp/finally.bus', cores=6)
