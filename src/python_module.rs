use std::collections::HashMap;
use bustools::{consistent_genes::{MappingMode, InconsistentResolution}, io::{BusFolder, BusReader}, iterators::{CbUmiGroupIterator, CellGroupIterator}, merger::MultiIterator, utils::get_progressbar};
use bustools_cli::butterfly;
use pyo3::{prelude::*, types::{PyFrozenSet, PyDict}};

/// Quck howto
/// ```python
/// import rustpybustools
/// rustpybustools.make_ecs('/home/michi/bus_testing/bus_output_shortest/', '/home/michi/bus_testing/transcripts_to_genes.txt', True)
/// 


/// Formats the sum of two numbers as string.
// #[pyfunction]
// fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
//     Ok((a + b).to_string())
// }

/// A Python module implemented in Rust.
#[pymodule]
fn pybustools(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(cbumi_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(cb_overlap, m)?)?;

    // m.add_class::<MyIterator>()?;

    // #[pyfn(m)]
    // #[pyo3(name = "testing")]
    // fn testing() -> PyResult<()> {
    //     println!("Hello rust");
    //     Ok(())
    // }

    #[pyfn(m)]
    #[pyo3(name = "make_ecs")]
    // fn make_ecs(foldername: &str, t2g_file: &str, collapse_ec:bool) -> PyResult<HashMap<usize, usize>> {
    fn make_ecs(busfile: &str, matrix_file: &str, transcript_file: &str, t2g_file: &str, collapse_ec:bool) -> PyResult<HashMap<usize, usize>> {

        let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
        
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

#[pyfunction]
fn cbumi_overlap(py: Python<'_>, busfolders: HashMap<String, String>) -> PyResult<PyObject> {

    let overlap = detect_cbumi_overlap(busfolders);
    let pythondict = PyDict::new(py);

    for (k,v) in overlap {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyFrozenSet::new(py, k.iter())?;
        pythondict.set_item(fset, v)?;
    }
    Ok(pythondict.to_object(py))
}

#[pyfunction]
fn cb_overlap(py: Python<'_>, busfolders: HashMap<String, String>) -> PyResult<PyObject> {

    let overlap = detect_cb_overlap(busfolders);
    let pythondict = PyDict::new(py);

    for (k,v) in overlap {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyFrozenSet::new(py, k.iter())?;
        pythondict.set_item(fset, v)?;
    }
    Ok(pythondict.to_object(py))
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



pub fn detect_cbumi_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {

    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cbumi());
    }

    let multi_iter = MultiIterator::new(h);
    let mut counter: HashMap<Vec<String>, usize> = HashMap::new();

    let bar = get_spinner();
    let interval = 100_000;

    for (i, ((_cb, _umi), record_dict)) in multi_iter.enumerate() {
        let mut the_set: Vec<String> = record_dict.keys().cloned().collect();
        the_set.sort();
        let val = counter.entry(the_set).or_insert(0);
        *val += 1;

        if i % interval == 0 {
            bar.inc(interval as u64);
        }
    }
    counter
}

pub fn detect_cb_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {

    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cb());
    }

    let multi_iter = MultiIterator::new(h);

    let bar = get_spinner();
    let interval = 100_000;

    let mut counter: HashMap<Vec<String>, usize> = HashMap::new();

    for (i, (_cb, record_dict)) in multi_iter.enumerate() {
        let mut the_set: Vec<String> = record_dict.keys().cloned().collect();
        the_set.sort();
        let val = counter.entry(the_set).or_insert(0);
        *val += 1;

        if i % interval == 0 {
            bar.inc(interval as u64);
        }
    }
    counter
}


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