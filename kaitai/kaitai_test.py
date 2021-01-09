from pybustools import busio, busio_old
import pathlib
import os

records = [
    busio.Bus_record(0, 0, 10, 20, 1),
    busio.Bus_record(1, 0, 13, 206, 12),
    busio.Bus_record(2, 0, 14, 250, 13)
]
fname =  '/tmp/some.bus'
busio.write_busfile(fname, records, cb_length=12, umi_length=5)

# ~/kaitai-struct-compiler-0.9/bin/kaitai-struct-compiler -t python busformat.ksy
from busfile import Busfile

b = Busfile.from_file('/tmp/some.bus')
for r in b.records:
    print(r.cb, r.umi, r.count)
