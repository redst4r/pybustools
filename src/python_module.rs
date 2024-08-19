use std::{collections::{HashMap, HashSet}, time::Instant};

use bustools_cli::correct::{build_correct_map, load_whitelist};
use pyo3::prelude::*;

/// Quck howto
/// ```python
/// import rustpybustools
/// rustpybustools.make_ecs('/home/michi/bus_testing/bus_output_shortest/', '/home/michi/bus_testing/transcripts_to_genes.txt', True)
/// 

/// A Python module implemented in Rust.
#[pymodule]
fn pybustools(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::make_ecs_cb, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::make_ecs_ec, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::make_ecs_across_cb, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::make_ecs, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::estimate_tgt, m)?)?;
    
    // m.add_function(wrap_pyfunction!(crate::phantom::phantom_fingerprint_cb, m)?)?;
    // m.add_function(wrap_pyfunction!(crate::phantom::phantom_fingerprint_cbumi, m)?)?;
    // m.add_function(wrap_pyfunction!(crate::phantom::rustphantom, m)?)?;
    // m.add_function(wrap_pyfunction!(crate::phantom::rustphantom_filter, m)?)?;

    m.add_function(wrap_pyfunction!(crate::counting::cbumi_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::cb_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::detect_cb_frequencies, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::kmer_counter_cb, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::count_reads_per_ec, m)?)?;
    
    m.add_function(wrap_pyfunction!(cb_correct, m)?)?;
    Ok(())
}

#[pyfunction]
/// Creates a error correction mapping of Cell barcodes (CBs) based on a whitelist
/// of CBs.
/// In particular it corrects every CB if its max 1-hamming away from one (*and only one*) whitelisted CB
/// (sometimes a seqeuence can be 1-hamming away from two whitelisted CBs )
pub (crate) fn cb_correct(_py: Python<'_>, cb_list: HashSet<String>, whitelist_fname: &str) -> PyResult<HashMap<String, String>> {

    println!("Loading whitelist");
    let whitelist = load_whitelist(whitelist_fname);
    println!("Done Loading whitelist");

    let t = Instant::now();
    let corrector = build_correct_map(&cb_list, &whitelist);

    println!("Build corrector: {} sec", t.elapsed().as_secs());
    // convert to strings
    let mut corrector_seqs = HashMap::new();
    for (cb_candidate, cb_true) in corrector {
        corrector_seqs.insert(
            bustools::utils::int_to_seq(cb_candidate, 16),
            bustools::utils::int_to_seq(cb_true, 16)
        );
    }

    Ok(corrector_seqs)
}
