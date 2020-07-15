if __name__ == '__main__':
    from pybustools.pybustools import Bus
    from pybustools.busio import read_binary_bus#, read_binary_bus2
    import collections
    import tqdm
    import toolz
    import time

    "gsutil -m cp -r gs://cruk-01-kallisto-nextflow/dnbseqg400.V300026370_88A.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs/ /tmp"

    f = '/run/media/michi/42506642-b470-4238-be14-bb0c303b3682/ISB_data/190524_10x_IRubin_HL60_groundZero/day_2_kallisto/kallisto/sort_bus/bus_output'
    f = '/tmp/dnbseqg400.V300026370_88A.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs/kallisto/sort_bus/bus_output'

    B1 = Bus(folder=f, decode_seq=True)
    B2 = Bus(folder=f, decode_seq=False)

    I = B1.iterate_bus()
    next(I)

    I = B2.iterate_bus()
    next(I)

    # speedtest
    %load_ext snakeviz
    %snakeviz list(toolz.take(1_000_000, B1.iterate_bus()))
    %snakeviz list(toolz.take(1_000_000, B2.iterate_bus()))


    # different buffesizes
    %snakeviz list(toolz.take(10_000_000, read_binary_bus(B1.bus_file,decode_seq=True, buffersize=1)))
    %snakeviz list(toolz.take(10_000_000, read_binary_bus(B1.bus_file,decode_seq=True, buffersize=10)))
    %snakeviz list(toolz.take(10_000_000, read_binary_bus(B1.bus_file,decode_seq=True, buffersize=100)))

    buffersize = [1,2,3,4,5,6,7,8,9,10,20,30,40.50,100,200,300,400,500, 1000, 10000, 100000, ]
    buffersize = [1,10,100, 1000, 10000, 100000, 1000000]
    import numpy as np
    buffersize = np.logspace(0, 6, 6).astype(int)

    def speed_test_buffer(buffersize):
        t1 = time.time()
        n_iter = 10_000_000
        gen = read_binary_bus(B1.bus_file, decode_seq=False, buffersize=buffersize)
        gen = toolz.take(n_iter, gen)
        for a in gen:
            pass
        t2 = time.time()
        return t2 - t1

    results = [{'buffer': b, 'time': speed_test_buffer(b)} for b in buffersize]

    import pandas as pd
    results = pd.DataFrame(results)

    %matplotlib
    import matplotlib.pyplot as plt
    plt.plot(results.buffer, results.time, 'x-')
    plt.xscale('log')
    plt.show()

    %load_ext snakeviz
    %%snakeviz
    for a in toolz.take(10_000_000, B2.iterate_bus()):
        pass


    """
    comparing the old read_binary_bus vs the new version using struct.iter_unpack
    """
    from pybustools.pybustools import Bus
    from pybustools import busio, busio_old

    import collections
    import tqdm
    import toolz
    import time

    "gsutil -m cp -r gs://cruk-01-kallisto-nextflow/dnbseqg400.V300026370_88A.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs/ /tmp"
    f = '/tmp/dnbseqg400.V300026370_88A.L05A_2-658952_cellranger_v3p0p1_refdata-cellranger-GRCh38-1_2_0.outs/kallisto/sort_bus/bus_output'

    B1 = Bus(folder=f, decode_seq=True)

    %load_ext snakeviz
    %%snakeviz
    for a in toolz.take(10_000_000, busio.read_binary_bus(B1.bus_file, decode_seq=False, buffersize=1000)):
        pass

    %%snakeviz
    for a in toolz.take(10_000_000, busio_old.read_binary_bus2(B1.bus_file, decode_seq=False, buffersize=1000)):
        pass

    import tqdm
    I1 = toolz.take(10_000_000, busio.read_binary_bus(B1.bus_file, decode_seq=False, buffersize=1000))
    I2 = toolz.take(10_000_000, busio_old.read_binary_bus(B1.bus_file, decode_seq=False, buffersize=1000))
    for a, b in tqdm.tqdm(zip(I1, I2)):
        assert a==b
