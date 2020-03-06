# pybustools
[![Build Status](https://travis-ci.org/redst4r/pybustools.svg?branch=master)](https://travis-ci.org/redst4r/pybustools)

This is a python interface to the [bustools-format](https://github.com/BUStools/BUS-format) for single cell RNAseq usually produced using [kallisto-bustools](https://www.kallistobus.tools/).

## Installation
Note the required `gmpy2` library, which itself relies on some c-libraries. These can be installed (on ubuntu) via
`apt-get install libgmp3-dev libgmp-dev libmpfr-dev libmpc-dev`

```
git clone https://github.com/redst4r/pybustools
pip install -e pybustools
```

## Usage

```python

from pybustools.pybustools import Bus

B = Bus(folder='/path/to/bus/output', bus_name='sorted.bus')

# iterate the records one by one; it yields namedtuples 
for record in B.iterate_bus():
    record.CB, record.UMI, record.EC, record.COUNT, record.FLAG  # str, str, int int, int

# aggregate the output of entire cells (WARNING: bus-file must be SORTED)
for cb, umi_list  in B.iterate_cells():
    """
    umi_list contains all UMIs of that cell:
          CB                    UMI         EC  Count  FLAG
    ('AAACCCAAGAACCGCA', [('GCATTAAATCAG', 822, 1,     0),
                          ('GCATTAAATCCG', 822, 1,     0)])
    """

    for umi, ec, count, flag in umi_list:
        pass
```

## Some tricks for speed
The bus-format encodes sequences (CB/UMI) as ints (from the docs: ... the barcode `GCCA` corresponds to the bit code `10010100` or the integer `148`). 
Decoding these into sequences is currently a performance bottleneck (done via `gmpy2.digits()` atm, which is already a pretty fast C-implementation).
If you dont care about actual sequence, you can turn of the decoding step and iterating over the bus records will return ints for CB/UMI:

```python
B = Bus(folder='/path/to/bus/output', decode_seq=False)
for record in B.iterate_bus():
    print(record.CB, record.UMI)  # these will be python ints now instead of str 
``` 

## TODO
- Unit-tests
- parsing the binary format is slow; the bottleneck is `busio._decode_int_to_ACGT()` and `gmpy2.digits()` in there. This converts the int encoding into cDNA sequence
- legacy code for reading plaintext-bus files
