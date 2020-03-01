# pybustools

This is a python interface to the [bustools-format](https://github.com/BUStools/BUS-format) for single cell RNAseq usually produced using [kallisto-bustools](https://www.kallistobus.tools/).

## Installation

```
git clone https://github.com/redst4r/pybustools
pip install -e pybustools
```

## Usage

```python

from pybustools.pybustools import Bus

B = Bus(folder='/path/to/bus/output', bus_name='sorted.bus')

# iterate the records/moleculesone by one
for cb,umi,ec,count,flag in B.iterate_bus():
    pass

# aggregate the output of entire cells (WARNING: bus-file must be SORTED)
for cb, umi_list  in B.iterate_cells():
    """
    umi_list contains all UMIs of that cell:
    ('AAACCCAAGAACCGCA', [('GCATTAAATCAG', 822, 1, 0),
                          ('GCATTAAATCCG', 822, 1, 0)])
    """

    for umi, ec, count, flag in umi_list:
        pass
```
