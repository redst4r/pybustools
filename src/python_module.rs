use std::collections::HashMap;
use bustools::{consistent_genes::{InconsistentResolution, MappingMode}, io::{BusFolder, BusReader, BusRecord}, iterators::{CbUmiGroupIterator, CellGroupIterator}, merger::MultiIterator, utils::get_progressbar};
use bustools_cli::butterfly::{self, CUHistogram};
use pyo3::{prelude::*, types::{PyDict, PyFrozenSet, PyTuple}};

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
    m.add_function(wrap_pyfunction!(make_ecs_cb, m)?)?;
    m.add_function(wrap_pyfunction!(phantom_fingerprint_cb, m)?)?;
    m.add_function(wrap_pyfunction!(phantom_fingerprint_cbumi, m)?)?;

    #[pyfn(m)]
    #[pyo3(name = "make_ecs")]
    // fn make_ecs(foldername: &str, t2g_file: &str, collapse_ec:bool) -> PyResult<HashMap<usize, usize>> {
    fn make_ecs(busfile: &str, matrix_file: &str, transcript_file: &str, t2g_file: &str, collapse_ec:bool) -> PyResult<HashMap<usize, usize>> {

        let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
        
        let mode = if collapse_ec {
            // println!("Making mapper");
            let mapper = folder.make_mapper(t2g_file);
            MappingMode::Gene(mapper, InconsistentResolution::IgnoreInconsistent)
        } else {
            MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
        };
        // println!("Done Making mapper");

        let h = butterfly::make_ecs(&folder, mode);
        // println!("umis: {}", h.get_numis());
        Ok(h.into())  // conversion into a hashmap/dict
    }

    Ok(())
}

/// counts frequencies of barcodes
#[pyfunction]
fn make_ecs_cb(
    busfile: &str, matrix_file: &str, transcript_file: &str, //t2g_file: &str, 
    // collapse_ec:bool,
) -> PyResult<(HashMap<usize, usize>, HashMap<usize, usize>)> {

    let folder = BusFolder::from_files(busfile, matrix_file, transcript_file);
    
    // let mode = if collapse_ec {
    //     // println!("Making mapper");
    //     let mapper = folder.make_mapper(t2g_file);
    //     MappingMode::Gene(mapper, InconsistentResolution::IgnoreInconsistent)
    // } else {
    //     MappingMode::EC(InconsistentResolution::IgnoreInconsistent)
    // };
    // println!("Done Making mapper");
    let mut h_umi: HashMap<usize, usize> = HashMap::new();
    let mut h_read: HashMap<usize, usize> = HashMap::new();

    for (_cb, records) in folder.get_iterator().groupby_cb() {
        let nreads: u32 = records.iter().map(|r| r.COUNT).sum();
        let v = h_read.entry(nreads as usize).or_insert(0);
        *v += 1;

        let numi = records.len();
        let v = h_umi.entry(numi).or_insert(0);
        *v += 1;
    }
    let cu_umi = CUHistogram::new(h_umi);
    let cu_reads = CUHistogram::new(h_read);

    Ok(
        (cu_umi.into(),cu_reads.into())  // conversion into a hashmap/dict
    ) 
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

// fn general_phantom_fingerprint<I,K>(h: HashMap<String, I>) -> (Vec<String>, HashMap<Vec<usize>, usize>, HashMap<Vec<usize>, usize>)

/// Given a set of iterators over busfiles (keyed/grouped by CB, or CB+UMI) count how many
/// keys overlap between busfiles
pub fn general_detect_overlap<I,K>(h: HashMap<String, I>) -> HashMap<Vec<String>, usize> 
where
    I: Iterator<Item = (K, Vec<BusRecord>)>,
    K: Ord + Eq + Debug + Copy,
{
    let multi_iter = MultiIterator::new(h);
    let mut counter: HashMap<Vec<String>, usize> = HashMap::new();

    let bar = get_spinner();
    let interval = 100_000;

    for (i, (_key, record_dict)) in multi_iter.enumerate() {
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
pub fn detect_cbumi_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {

    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cbumi());
    }
    general_detect_overlap(h)
}

pub fn detect_cb_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {

    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cb());
    }
    general_detect_overlap(h)
}

pub fn detect_cbumi_overlap_old(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {

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

pub fn detect_cb_overlap_old(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {

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


#[pyfunction]
fn phantom_fingerprint_cb(py: Python<'_>, busfolders: HashMap<String, String>) -> PyResult<(PyObject, PyObject, PyObject)> {

    let (samplenames, fp_umi, fp_read ) = _phantom_fingerprint_cb(busfolders);
    let pythondict_umi = PyDict::new(py);
    let pythondict_reads = PyDict::new(py);

    for (k,v) in fp_umi {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict_umi.set_item(fset, v)?;
    }
    for (k,v) in fp_read {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict_reads.set_item(fset, v)?;
    }
    Ok((samplenames.to_object(py), pythondict_umi.to_object(py), pythondict_reads.to_object(py)))
}

#[pyfunction]
fn phantom_fingerprint_cbumi(py: Python<'_>, busfolders: HashMap<String, String>) -> PyResult<(PyObject, PyObject, PyObject)> {

    let (samplenames, fp_umi, fp_read ) = _phantom_fingerprint_cbumi(busfolders);
    let pythondict_umi = PyDict::new(py);
    let pythondict_reads = PyDict::new(py);

    for (k,v) in fp_umi {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict_umi.set_item(fset, v)?;
    }
    for (k,v) in fp_read {
        // need to convert it to a frozenSet explicitly (regular sets cant be hashed in python)
        let fset = PyTuple::new(py, k.iter());
        pythondict_reads.set_item(fset, v)?;
    }
    Ok((samplenames.to_object(py), pythondict_umi.to_object(py), pythondict_reads.to_object(py)))
}


fn _phantom_fingerprint_cb(busfolders: HashMap<String, String>) -> (Vec<String>, HashMap<Vec<usize>, usize>, HashMap<Vec<usize>, usize>){
    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cb());
    }
    general_phantom_fingerprint(h)
}

fn _phantom_fingerprint_cbumi(busfolders: HashMap<String, String>) -> (Vec<String>, HashMap<Vec<usize>, usize>, HashMap<Vec<usize>, usize>){
    let mut h = HashMap::new();
    for (k,filename) in busfolders{
        h.insert(k, BusReader::new(&filename).groupby_cbumi());
    }
    general_phantom_fingerprint(h)
}

use std::fmt::Debug;

fn general_phantom_fingerprint<I,K>(h: HashMap<String, I>) -> (Vec<String>, HashMap<Vec<usize>, usize>, HashMap<Vec<usize>, usize>)
where
    I: Iterator<Item = (K, Vec<BusRecord>)>,
    K: Ord + Eq + Debug + Copy,
    
{
    let samplenames: Vec<_> = h.keys().cloned().collect();

    let multi_iter = MultiIterator::new(h);

    let bar = get_spinner();
    let interval = 100_000;

    let mut umi_hist : HashMap<Vec<usize>, usize> = HashMap::new();
    let mut read_hist : HashMap<Vec<usize>, usize> = HashMap::new();

    for (i, (_cb, record_dict)) in multi_iter.enumerate() {

        let mut  umi_fp = Vec::new();
        let mut read_fp = Vec::new();

        for s in samplenames.iter(){
            let (numi, nreads) =match record_dict.get(s) {
                Some(records) => {
                    let numi = records.len();
                    let nreads: u32 = records.iter().map(|r| r.COUNT).sum() ;
                    (numi, nreads as usize)
                },
                None => (0,0),
            };
            umi_fp.push(numi);
            read_fp.push(nreads);
        }
        let val = umi_hist.entry(umi_fp).or_insert(0);
        *val += 1;

        let val = read_hist.entry(read_fp).or_insert(0);
        *val += 1;

        if i % interval == 0 {
            bar.inc(interval as u64);
        }
    }
    (samplenames, umi_hist, read_hist)
}