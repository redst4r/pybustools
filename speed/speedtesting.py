if __name__ == '__main__':
    from pybustools.pybus import Bus
    from pybustools.busio import read_binary_bus#, read_binary_bus2
    import collections
    import tqdm
    import toolz
    import time
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import plotnine as pn
    import tqdm

    f = '/tmp/bus_output'

    B1 = Bus(folder=f, decode_seq=True)
    B2 = Bus(folder=f, decode_seq=False)

    # speedtest
    %load_ext snakeviz
    %snakeviz list(toolz.take(1_000_000, B1.iterate_bus()))
    %snakeviz list(toolz.take(1_000_000, B2.iterate_bus()))


    %load_ext snakeviz
    %%snakeviz
    for a in toolz.take(10_000_000, B2.iterate_bus()):
        pass


    """
    comparing the old read_binary_bus vs the new version using struct.iter_unpack
    """
    from pybustools.pybus import Bus
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
