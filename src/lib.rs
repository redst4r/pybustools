mod python_module;
// mod phantom;
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
*/