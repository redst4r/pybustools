from pybustools.pybus import Bus
from pybustools.busio import read_binary_bus#, read_binary_bus2
import collections
import toolz
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotnine as pn
import tqdm


f = '/tmp/bus_output'
B1 = Bus(folder=f, decode_seq=False)

# different buffesizes
buffersize = np.logspace(0, 7, 15).astype(int)


def speed_test_buffer(buffersize):
    t1 = time.time()
    n_iter = 10_000_000
    gen = read_binary_bus(B1.bus_file, decode_seq=False, buffersize=buffersize)
    gen = toolz.take(n_iter, gen)
    for a in gen:
        pass
    t2 = time.time()
    return t2 - t1


results = [
    {'buffer': b, 'time': speed_test_buffer(b)}
    for _ in tqdm.trange(20) for b in buffersize
]

results = pd.DataFrame(results)
results.to_csv('/home/michi/ms_python_packages/pybustools/speed/results.csv')
%matplotlib
plt.scatter(results.buffer, results.time)
#     plt.plot(results.buffer, results.time, 'x-')
plt.xscale('log')
plt.savefig('/home/michi/ms_python_packages/pybustools/speed/buffer.png', dpi=300)
plt.show()

p = pn.ggplot(results.query('buffer>=100'), pn.aes('factor(buffer)', 'time', color='factor(buffer)')) + pn.geom_boxplot()+ pn.geom_jitter() +pn.labs(title='Buffersize vs Time', x='Buffersize', y='Time (sec)')
p.save('/home/michi/ms_python_packages/pybustools/speed/buffer_vs_time.png', dpi=300)
