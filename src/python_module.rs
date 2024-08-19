use pyo3::prelude::*;

/// Quck howto
/// ```python
/// import rustpybustools
/// rustpybustools.make_ecs('/home/michi/bus_testing/bus_output_shortest/', '/home/michi/bus_testing/transcripts_to_genes.txt', True)
/// 

/// A Python module implemented in Rust.
#[pymodule]
fn pybustools(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::make_ecs_cb, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::make_ecs_ec, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::make_ecs_across_cb, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::make_ecs, m)?)?;
    m.add_function(wrap_pyfunction!(crate::butterfly::estimate_tgt, m)?)?;
    
    m.add_function(wrap_pyfunction!(crate::phantom::phantom_fingerprint_cb, m)?)?;
    m.add_function(wrap_pyfunction!(crate::phantom::phantom_fingerprint_cbumi, m)?)?;
    m.add_function(wrap_pyfunction!(crate::phantom::rustphantom, m)?)?;
    m.add_function(wrap_pyfunction!(crate::phantom::rustphantom_filter, m)?)?;

    m.add_function(wrap_pyfunction!(crate::counting::cbumi_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::cb_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::detect_cb_frequencies, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::kmer_counter_cb, m)?)?;
    m.add_function(wrap_pyfunction!(crate::counting::count_reads_per_ec, m)?)?;
    
    Ok(())
}
