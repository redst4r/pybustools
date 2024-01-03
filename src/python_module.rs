use std::collections::HashMap;
use bustools::{io::BusFolder, consistent_genes::{MappingMode, InconsistentResolution}};
use bustools_cli::butterfly;
use pyo3::prelude::*;

/// Quck howto
/// ```python
/// import rustpybustools
/// rustpybustools.make_ecs('/home/michi/bus_testing/bus_output_shortest/', '/home/michi/bus_testing/transcripts_to_genes.txt', True)
/// 


/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn rustpybustools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;

    // m.add_class::<MyIterator>()?;

    #[pyfn(m)]
    #[pyo3(name = "testing")]
    fn testing() -> PyResult<()> {
        println!("Hello rust");
        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(name = "make_ecs")]
    fn make_ecs(foldername: &str, t2g_file: &str, collapse_ec:bool) -> PyResult<HashMap<usize, usize>> {
        let folder = BusFolder::new(foldername);
        
        let mode = if collapse_ec {
            println!("Making mapper");
            let mapper = folder.make_mapper(t2g_file);
            MappingMode::Gene(mapper, InconsistentResolution::IgnoreInconsistent)
        } else {
            MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
        };
        println!("Done Making mapper");


        let h = butterfly::make_ecs(&folder, mode);
        println!("umis: {}", h.get_numis());
        Ok(h.into())  // conversion into a hashmap/dict
    }

    Ok(())
}


// #[pyclass]
// struct MyIterator {
//     iter: Box<dyn Iterator<Item = usize> + Send>,
// }

// #[pymethods]
// impl MyIterator {

//     #[new]
//     fn new() -> Self {

//         let iter = vec![0_usize,1,2,3,4]
//             .iter()
//             .map(|x| *x);
//         MyIterator { iter: Box::new(iter) }
//     }

//     fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
//         slf
//     }
//     fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyObject> {
//         let x = slf.iter.next();
//         x.into()
//     }
// }
