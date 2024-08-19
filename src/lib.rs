mod python_module;
mod phantom;
mod counting;
mod butterfly;
//:w
//mod tests;


/// Get a `indicatif::ProgressBar` with unknown size.
fn get_spinner() -> indicatif::ProgressBar{
    let bar = indicatif::ProgressBar::new_spinner();
    bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {pos} {per_sec}")
            .unwrap()
            .progress_chars("##-"),
    );
    bar
}

/* some testing for maturin

from pybustools import pybustools
# uncomprssed
q1 = pybustools.make_ecs_across_cb(
    '/home/michi/bus_testing/bus_output/output.corrected.sort.bus', 
    '/home/michi/bus_testing/bus_output/matrix.ec', 
    '/home/michi/bus_testing/bus_output/transcripts.txt')

# comprssed
q2 = pybustools.make_ecs_across_cb(
    '/home/michi/bus_testing/bus_output/output.corrected.sort.busz', 
    '/home/michi/bus_testing/bus_output/matrix.ec', 
    '/home/michi/bus_testing/bus_output/transcripts.txt')

from pybustools import pybustools
df, multildict = pybustools.estimate_tgt(
    '/home/michi/bus_testing/bus_output/output.corrected.sort.busz', 
    '/home/michi/bus_testing/bus_output/matrix.ec',
    '/home/michi/bus_testing/bus_output/transcripts.txt',
    '/home/michi/bus_testing/transcripts_to_genes.txt'
)


from pybustools import pybustools
t2g_file = '/home/michi/mounts/TB4drive/kallisto_resources_v50/t2g.txt'

df, multildict, inconsistent_dict = pybustools.estimate_tgt(
    '/home/michi/mounts/TB4drive/ISB_data/240116_VH00715_135_AACGWHMHV/kallisto/HL60_S1/kallisto/sort_bus/bus_output/output.corrected.sort.busz', 
    '/home/michi/mounts/TB4drive/ISB_data/240116_VH00715_135_AACGWHMHV/kallisto/HL60_S1/kallisto/sort_bus/bus_output/matrix.ec',
    '/home/michi/mounts/TB4drive/ISB_data/240116_VH00715_135_AACGWHMHV/kallisto/HL60_S1/kallisto/sort_bus/bus_output/transcripts.txt',
    t2g_file
)
import toolz
toolz.keyfilter(lambda k: len(k)<2, multildict)  # better be empty!
import collections
collections.Counter(multildict).most_common(10)


from pybustools.butterfly import _make_ec2gene_dict
from pybustools import Bus
bus = Bus('/home/michi/mounts/TB4drive/ISB_data/240116_VH00715_135_AACGWHMHV/kallisto/HL60_S1/kallisto/sort_bus/bus_output')
_make_ec2gene_dict(bus,  t2g_file)
*/