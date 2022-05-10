import pathlib
import numpy as np
from pybustools.busio import read_binary_bus, Bus_record, write_busfile, get_header_info
from pybustools.pybustools import iterate_CB_UMI_of_busfile, iterate_cells_of_busfile
import tqdm
import shutil


def get_number_of_reads_and_molecules(fname):
    """
    similar to bustools inspect
    Gets the number of reads and number of bus-records (approx #UMI) in a busfile

    :param fname: filename of the busfile
    """
    # TODO: technically not correct: this count the number of reads (fine) and number of ENTRIES in the busfile.
    # ideally ach molecule would have a single entry.
    # However: a single molecule might map to two different EC classes.
    # The molecule got fragmented into two places, mapping it to two diff locations
    total_reads = 0
    total_molecules = 0
    for record in tqdm.tqdm(read_binary_bus(fname, decode_seq=False), desc='counting reads'):
        total_reads += record.COUNT
        total_molecules += 1
    return total_reads, total_molecules


def subsample_busfile(fname_in, fname_out, fraction):
    """
    subsample the reads of an existing busfile by `fraction`, writing the Result
    into a new busfile. The major effect is that some entries will recieve
    0 reads and hence disappear from the file!

    :param fname_in: Filename of the input busfile
    :param fname_out: Filename of the resulting, subsampled busfile
    :param fraction: 0<fraction<1, the fraction of subsampling
    """
    assert 1 == 0, "we can replace the entire thing with bustools count --downsample!!"
    assert 0 < fraction < 1, "fraction must be in [0,1]"

    G = tqdm.tqdm(subsample_busfile_cb_umi_generator(fname_in, fraction))
    # we need to write the correct header
    _, cb_len, umi_len, _ = get_header_info(fname_in)
    write_busfile(fname_out, G, cb_length=cb_len, umi_length=umi_len)


def subsample_busfile_cb_umi_generator(fname, fraction):
    """
    returns a generator across cb/umi but subsampled the reads (not the records!)
    essentially, eaech record is binomially subsampled. If it still has counts, its yieled
    """
    for (cb, umi), records in iterate_CB_UMI_of_busfile(fname, decode_seq=False):
        for r in records:
            c = np.random.binomial(r.COUNT, fraction)
            if c > 0:
                yield Bus_record(cb, umi, r.EC, c, r.FLAG)


def subsample_busfile_cb_generator(fname, fraction):
    """
    returns a generator across cb but subsampled the reads (not the records!)
    essentially, eaech record is binomially subsampled. If it still has counts, its yieled
    """
    for cb, records in iterate_cells_of_busfile(fname, decode_seq=False):
        subsampled_records = []
        for r in records:
            c = np.random.binomial(r.COUNT, fraction)
            if c > 0:
                subsampled_records.append(
                    Bus_record(cb, r.UMI, r.EC, c, r.FLAG)
                    )
        if len(subsampled_records) > 0:
            yield cb, subsampled_records


def subsample_kallisto_bus(busdir, outdir, fraction):
    """
    not only subsamples the .bus file, but also sets up a folder that looks like a original kallisto quantification
    so that `bustools count` can be run on it

    :param busdir: directory containing the bustools quantification (if nextflow: .../kallisto/sort_bus/bus_output).
                   must contain `output.corrected.sort.bus`, `matrix.ec`, `transcripts.txt`

    :param outdir: must exist. puts the downsampled kallisto .bus and suppl.files here
    :param fraction: 0<fraction<1: fraction of reads to sample

    """
    assert 1 == 0, "we can replace the entire thing with bustools count --downsample!!"
    assert pathlib.Path(outdir).exists()
    infile = f'{busdir}/output.corrected.sort.bus'
    outfile = f'{outdir}/output.corrected.sort.bus'

    subsample_busfile(infile, outfile, fraction=fraction)

    # copy the other parts of the kallisto output, they remain unchanged
    ec = f'{busdir}/matrix.ec'
    transcripts = f'{busdir}/transcripts.txt'

    shutil.copy(ec, outdir)
    shutil.copy(transcripts, outdir)


def subsample_bustools(bus, t2gfile, fraction, outdir_and_prefix):
    """
    creates the call to bustools count for subsampling a busfile. MUCH fast than doign it in python
    """
    p = pathlib.Path(outdir_and_prefix)
    assert p.parent.exists(), "outputfolder doenst exist. if not estining bustools will quitely fail"
    print(f"bustools count --genecounts -o {outdir_and_prefix} -d {fraction} -t {bus.transcript_file} -e {bus.ec_file} --genemap {t2gfile} {bus.bus_file}")
